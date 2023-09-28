""" Methods for FRB targeting """

from YSE_App.models import *

import pandas

from astroplan import Observer
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units

from YSE_App.chime import tags as chime_tags

from IPython import embed

def calc_airmasses(frb_fu, gd_frbs,
                    debug:bool=False):
    """ Calulate the minimum AM for a set of FRBs
    tied to a FRBFollowupResource

    Moved this code here to speed up development.  Could be
    a method on FRBFollowupResource someday..

    Args:
        frb_fu (FRBFollowupResource): _description_
        gd_frbs (QuerySet of FRBTransient): _description_
        debug (bool, optional): _description_. Defaults to False.

    Returns:
        np.ndarray: Minimum airmass values for the input FRBs
            during the observing period
    """
                    
    # Grab the telescope 
    telescope = frb_fu.instrument.telescope
    location = EarthLocation.from_geodetic(
        telescope.longitude*units.deg,telescope.latitude*units.deg,
        telescope.elevation*units.m)
    tel = Observer(location=location, timezone="UTC")

    nfrb = len(gd_frbs)

    # Coords
    ras = [frbt.ra for frbt in gd_frbs]
    decs = [frbt.dec for frbt in gd_frbs]
    frb_coords = SkyCoord(ra=ras, dec=decs, unit='deg')

    min_AM = np.array([1e9]*nfrb)

    # Loop on nights
    this_time = Time(frb_fu.valid_start)
    end_time = Time(frb_fu.valid_stop)
    
    # Loop on nights
    while(this_time < end_time):
        night_end = tel.twilight_morning_astronomical(this_time)
        this_obs_time = this_time.copy()


        # Loop on 30min intervals from 
        while(this_obs_time < min(end_time,night_end)):
            if debug:
                print(this_obs_time.datetime)

            # Calculate AM
            altaz = tel.altaz(this_obs_time, frb_coords)
            airmass = 1/np.cos((90.-altaz.alt.value)*np.pi/180.)
            # Below horizon
            airmass[altaz.alt.value < 0] = 1e9

            # Save minimum AM
            min_AM = np.minimum(min_AM, airmass)

            # Increment by 30min
            this_obs_time = this_obs_time + TimeDelta(1800, format='sec')

        # Add a day -- Added to night_end to avoid infinite while loop
        this_time = night_end + TimeDelta(1, format='jd')
        this_time = tel.twilight_evening_astronomical(this_time, which='previous')

    return min_AM

def target_table_from_frbs(frbs, mode:str):
    """ Generate a pandas table from a list or QuerySet of FRBs

    Args:
        frbs (QuerySet): List of FRBTransient objects
        mode (str): Observation mode (e.g. 'imaging')

    Returns:
        pandas.DataFrame: table of FRBTransient properties
    """


    # Loop in FRBTansients
    rows = []
    for itransient in frbs.all():
        rdict = {}
        rdict['TNS'] = itransient.name
        rdict['FRB_RA'] = itransient.ra
        rdict['FRB_Dec'] = itransient.dec
        rdict['FRB_DM'] = itransient.DM
        rdict['FRB_Survey'] = itransient.FRBSurveyString()
        rdict['FRB_tags'] = itransient.FRBTagsString()

        # Host candidates?
        if Path.objects.filter(transient=itransient).count() > 0:
            path_values, galaxies = itransient.get_Path_values()
            srt = np.argsort(path_values)[::-1]

            for kk in range(min(2,len(path_values))):
                galaxy = galaxies[srt[kk]]
                if kk == 0:
                    prefix = 'Pri'
                elif kk == 1:
                    prefix = 'Sec'

                # Take top two
                rdict[f'{prefix}_name'] = galaxy.name
                rdict[f'{prefix}_RA'] = galaxy.ra
                rdict[f'{prefix}_Dec'] = galaxy.dec
                rdict[f'{prefix}_POx'] = galaxy.P_Ox
                # Photometry
                ifilter, mag = galaxy.FilterMagString()
                rdict[f'{prefix}_mag'] = float(mag)
                rdict[f'{prefix}_filter'] = ifilter
        #
        rows.append(rdict)

    # Table
    df = pandas.DataFrame(rows)
    # Add mode
    if len(df) > 0:
        df['mode'] = mode

    # Return
    return df

def targetfrbs_for_fu(frb_fu):
    """ Provided an FRBFollowupResource, return the
    QuerySet of FRBTransient objects that are acceptible

    This may eventually be made a method on FRBFollowupResource

    Args:
        frb_fu (FRBFollowupResource): Follow-up resource

    Returns:
        QuerySet of FRBTransient:  All transients satisfying the criteria
    """
    # 
    # Check Surveys
    if frb_fu.frb_surveys == 'all':
        gd_frbs = FRBTransient.objects.all()
    else:
        gd_frbs = FRBTransient.objects.filter(frb_survey__in=FRBSurvey.objects.filter(
            name__in=frb_fu.frb_surveys.split(',')))

    # Status
    if frb_fu.frb_statuses:
        gd_frbs = gd_frbs.filter(status__in=TransientStatus.objects.filter(
            name__in=frb_fu.frb_statuses.split(',')))
    else:
        gd_frbs = gd_frbs.filter(
            status__in=TransientStatus.objects.filter(
                name__in=['NeedImage', 'NeedSpectrum']))

    # Tags? aka samples
    if frb_fu.frb_tags:
        gd_frbs = gd_frbs.filter(frb_tags__in=FRBTag.objects.filter(name__in=frb_fu.frb_tags.split(',')))
    
    return gd_frbs

def grab_targets_by_mode(frb_fu, frbs):
    """ Grab targets by observing mode

    Logic is applied as follows:

      -- Imaging: The FRB must not have a host
      -- Longslit: The FRB must have a host and the host must have a P_Ox > frb_fu.min_POx
      -- Mask: Not implemented yet

    Args:
        frb_fu (FRBFollowupResource): Follow-up resource
        frbs (QuerySet): List of FRBTransient objects

    Returns:
        dict: Dictionary of QuerySets for each observing mode
    """

    # Imaging
    if frb_fu.num_targ_img > 0:
        #imaging_frbs = frbs.filter(host__isnull=True)
        imaging_frbs = frbs.filter(
            status=TransientStatus.objects.get(name='Image'))
    else:
        imaging_frbs = FRBTransient.objects.none()


    # Longslit
    if frb_fu.num_targ_longslit > 0:
        longslit_frbs = frbs.filter(
            status=TransientStatus.objects.get(name='Spectrum'))
        if frb_fu.min_POx:
            gd_ids = []
            for frb in longslit_frbs:
                if frb.host.P_Ox is not None and frb.host.P_Ox > frb_fu.min_POx:
                    gd_ids.append(frb.id)
            longslit_frbs = longslit_frbs.filter(id__in=gd_ids)

    # Mask -- Not yet implemented
    mask_frbs = FRBTransient.objects.none()

    # Populate
    targets_by_mode = {}
    targets_by_mode['imaging'] = imaging_frbs
    targets_by_mode['longslit'] = longslit_frbs
    targets_by_mode['mask'] = mask_frbs

    # Return
    return targets_by_mode

def select_with_priority(frb_fu, frbs_by_mode:dict): 
    """ Select targets from the total set by priority

    Args:
        frb_fu (FRBFollowupResource): Follow-up resource
        frbs_by_mode (dict): Dictionary of QuerySets for each observing mode

    Returns:
        dict: final dictionary of selected QuerySets for each observing mode
    """

    # Final product
    selected_frbs = {}

    # Loop on all the modes
    for mode, num_targ in zip(['imaging', 'longslit', 'mask'], 
                              [frb_fu.num_targ_img, frb_fu.num_targ_longslit, 
                               frb_fu.num_targ_mask]):
        # Do we want any?
        if num_targ > 0:
            # Do we not have enough?
            if frb_fu.num_targ_img >= len(frbs_by_mode[mode]):
                selected_frbs[mode] = frbs_by_mode[mode]
            else:    
                # Assign priorities
                probs = assign_probs(frbs_by_mode[mode])

                # Init
                keep = []
                still_left = np.arange(len(frbs_by_mode[mode])).tolist()

                # Loop until we have enough
                while(len(keep) < num_targ):
                    need = num_targ - len(keep)
                    # Random numbers
                    rand = np.random.uniform(0.,1.,size=len(still_left))
                    left_probs = np.array(probs)[still_left]
                    gd = np.where(rand < left_probs)[0]
                    if len(gd) > need:
                        gd = gd[:need]
                    keep += list(np.array(still_left)[gd])
                    # Remove
                    for idx in gd:
                        still_left.remove(idx)
                # Finish
                gd_ids = [frbs_by_mode[mode].all()[int(ii)].id for ii in keep]
                selected_frbs[mode] = frbs_by_mode[mode].filter(id__in=gd_ids)
        else:
            # We don't want any
            selected_frbs[mode] = FRBTransient.objects.none()

    # Return
    return selected_frbs

        
def assign_probs(frbs):

    # Collect all possible samples
    all_samples = chime_tags.all_samples

    probs = []
    for frb in frbs.all():
        max_prob = 0.1 # Default
        tag_names = [itag.name for itag in frb.frb_tags.all()]

        # Query them all
        for sample in all_samples:
            if sample['name'] in tag_names:
                max_prob = max(max_prob, sample['prob'])
        # Save
        probs.append(max_prob)
    # Return
    return probs
