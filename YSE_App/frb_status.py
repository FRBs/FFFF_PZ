""" Code related to the Status of FRBs """

import numpy as np


from YSE_App import frb_tags
from YSE_App.chime import tags as chime_tags

from IPython import embed

# Create the FRB ones
all_status = [\
    'Unassigned', # Does not meet the criteria for FFFF FollowUp
    'BrightStar', # Localization is too close to a very bright star
    'TooDusty', # Sightline exceeds E(B-V) threshold
    'RunPublicPATH', # Needs to be run through PATH with public data
        # P_Ux is None
        # At least one frb_tag is in the list of run_public_path entries below
    'NeedImage', # Needs deeper imaging
        # P_Ux > maximum of all P_Ux_max for the frb_tags
        #  or r-mag of the top candidate is too faint for the public survey
        #   23.0 for Blanco/DECam, 21.0 for Pan-STARRS
        # No successful Image taken
        # No Image pending
    'NeedSpectrum', # Needs spectroscopy for redshift
        # P(O|x) of top 2 > P_Ox_min
        # No pending spectrum
        # No succesfully observed spectrum
    'RunDeepPATH', # Needs PATH run on deeper (typically private) imaging
    'ImagePending', # Pending deeper imaging with an FRBFollowUp
        # FRB appears in FRBFollowUpRequest with mode='image'
    'SpectrumPending', # Pending spectroscopy with an FRBFollowUp
        # FRB appears in FRBFollowUpRequest with mode='longslit','mask'
    'GoodSpectrum', # Observed with spectroscopy successfully
        # P(O|x) of top 2 > P_Ox_min
        # If Primary does not exceed min_POx, then the top two must have a spectrum
    'TooFaint', # Host is too faint for spectroscopy
        # r-magnitude (or equivalent; we use the PATH band) of the top host candidate
        #   is fainter than the maximum(mr_max) for the sample/surveys
    'AmbiguousHost',  # Host is considered too ambiguous for further follow-up
        #  For primary-only, it isn't satisfied and for top two, the top two P(O|x) are not satisfied
    'UnseenHost',  # Even with deep imaging, no compelling host was found
        # P(U|x) is set
        # deep PATH was run
        # P(U|x) > P_Ux_max
        # Note that this takes precedence over 'AmbiguousHost' 
    'Redshift', # Redshift measured
        # P(O|x) of top 2 > P_Ox_min
        # If Primary does not exceed min_POx, then the top two redshifts must be nearly the same
        # Else, take primary
]



# Add all of the chime
def set_status(frb):
    """ Set the status of an FRB transient 

    The frb is modified and saved

    Args:
        frb (FRBTransient): FRBTransient instance
    """
    log_message = ''

    # Hide here for circular imports
    from YSE_App.models import TransientStatus
    from YSE_App.models import Path
    from YSE_App.models import FRBFollowUpObservation
    from YSE_App.models import FRBFollowUpRequest

    # Check Criteria 
    criteria = frb_tags.chk_all_criteria(frb)
    PATH_run = False if frb.host is None else True

    # Is the top candidate too faint?
    r_too_faint = False  
    if frb.host is not None:
        POx_values, galaxies, _ = frb.get_Path_values()
        argsrt = np.argsort(POx_values)
        pri_gal = galaxies[argsrt[-1]]  # Primary galaxy
        # Check the top candidate magnitude
        rfilter = pri_gal.FilterMagString()
        if 'Blanco' in rfilter or 'DECam' in rfilter:
            if pri_gal.mag > 23.0:
                r_too_faint = True
        elif 'Pan-STARRS' in rfilter:
            if pri_gal.mag > 21.0:
                r_too_faint = True

    # Run in reverse order of completion

    # #########################################################
    # #########################################################
    # Bright star?
    # #########################################################
    if np.all(criteria['bright_star']):
        frb.status = TransientStatus.objects.get(name='BrightStar')
        frb.save()
        return

    # #########################################################
    # #########################################################
    # Too Dusty??
    # #########################################################
    if np.all(np.invert(criteria['bright_star']) & np.invert(criteria['EBV'])):
        frb.status = TransientStatus.objects.get(name='TooDusty')
        frb.save()
        return

    # #########################################################
    # Run Public PATH
    # #########################################################
    if np.any(np.invert(criteria['bright_star']) & criteria['EBV'] & \
        criteria['run_public_PATH']):
        if not PATH_run:
            frb.status = TransientStatus.objects.get(name='RunPublicPATH')
            frb.save()
            return
    else: # We have chosen not to proceed with this FRB; all items that follow require PATH
        frb.status = TransientStatus.objects.get(name='Unassigned')
        frb.save()
        return

    # Grab the sample that satisfy the criteria so far
    good = np.invert(criteria['bright_star']) & criteria['EBV'] & \
        criteria['run_public_PATH']
    good_idx = np.where(good)[0]

    # #########################################################
    # Need Image
    # #########################################################
    if (np.any(criteria['PUx'][good_idx]) or r_too_faint) and (
        not FRBFollowUpRequest.objects.filter(
            transient=frb,
            mode='image').exists()) and (
        not FRBFollowUpObservation.objects.filter(
            transient=frb,
            success=True,
            mode='image').exists()):

        frb.status = TransientStatus.objects.get(name='NeedImage')
        frb.save()
        return

    # #########################################################
    # Run deep PATH
    # #########################################################
    if FRBFollowUpObservation.objects.filter(
            transient=frb,
            success=True,
            mode='image').exists() and (not criteria['ran_deep_PATH'][0]):
        frb.status = TransientStatus.objects.get(name='RunDeepPATH')
        frb.save()
        return

    # #########################################################
    # Unseen host
    # #########################################################
    # In PATH table?
    if np.any(criteria['PUx'][good_idx] & criteria['ran_deep_PATH'][good_idx]):
        frb.status = TransientStatus.objects.get(
                        name='UnseenHost')
        frb.save()
        return

    # #########################################################
    # Ambiguous host
    # #########################################################
    if np.all(np.invert(criteria['POx'][good_idx])):
        frb.status = TransientStatus.objects.get(name='AmbiguousHost')
        frb.save()
        return


    # #########################################################
    # Redshift?
    # #########################################################
    if frb.host is not None and frb.host.redshift is not None and \
        np.any(criteria['POx'][good_idx]):  # This last query is superfluous but it is here for clarity 

        if np.any(criteria['z_done'][good_idx] & criteria['z_consistent'][good_idx]):
            # We have a redshift
            frb.status = TransientStatus.objects.get(name='Redshift')
            frb.save()
            return

        if np.any(criteria['z_done'][good_idx] & np.invert(criteria['z_consistent'][good_idx])):
            # Amibiguous host
            frb.status = TransientStatus.objects.get(name='AmbiguousHost')
            frb.save()
            return

        '''
        # Grab info
        path_values, galaxies, _ = frb.get_Path_values()
        argsrt = np.argsort(path_values)

        # Items to loop over
        if POx_satisfied_primary: 
            idxs = argsrt[-1:]  # Primary is the last one
        else:
            idxs = argsrt[-2:]  # Top 2

        # Loop on the info
        has_redshift = []
        source_ok = []
        for ss, idx in enumerate(idxs):
            # Grab the galaxy
            gal = galaxies[idx]
            # Check redshift
            if gal.redshift is None:
                has_redshift.append(False)
            else:
                has_redshift.append(True)
            # Check the redshift source
            tmp_ok = False
            for gd_source in good_z_sources:
                if gd_source in gal.redshift_source:
                    tmp_ok = True
            source_ok.append(tmp_ok)

        # Time to set the status
        if POx_satisfied_primary:
            log_message += "I AM A PRIMARY-"
            # Ok?
            if np.all(has_redshift) and np.all(source_ok):
                frb.status = TransientStatus.objects.get(name='Redshift')
                frb.save()
                return log_message
        else: # Top 2
            log_message += "I AM NOT A PRIMARY-"
            # Check the redshifts are nearly the same
            if np.all(has_redshift) and np.all(source_ok):
                log_message += "I AM OK-"
                if np.abs(galaxies[argsrt[-1]].redshift - galaxies[argsrt[-2]].redshift) > 0.003:
                    log_message += "I AM NOT CONSINSTENT-"
                    frb.status = TransientStatus.objects.get(name='AmbiguousHost') 
                    frb.save()
                    return log_message
                else:
                    log_message += "I AM CONSINSTENT-"
                    frb.status = TransientStatus.objects.get(name='Redshift')
                    frb.save()
                    return log_message
        '''

    # #########################################################
    # Too Faint?
    # #########################################################
    if frb.host is not None:
        mrs = frb_tags.values_from_tags(frb, 'max_mr')

        # Find mr_max (if it exists)
        if len(mrs) > 0:
            mr_max = np.max(mrs)
            # Use PATH host magnitudes
            if frb.mag_top_two_PATH > mr_max:
                frb.status = TransientStatus.objects.get(name='TooFaint')
                frb.save()
                return
    
    # #########################################################
    # Good Spectrum
    # #########################################################
    if FRBFollowUpObservation.objects.filter(
            transient=frb,
            success=True,
            mode__in=['longslit','mask']).exists():
        frb.status = TransientStatus.objects.get(name='GoodSpectrum')
        frb.save()
        return

    # #########################################################
    # Pending Spectrum
    # #########################################################

    if FRBFollowUpRequest.objects.filter(
            transient=frb,
            mode__in=['longslit','mask']).exists():
        frb.status = TransientStatus.objects.get(name='SpectrumPending')
        frb.save()
        return


    # #########################################################
    # Need Spectrum
    # #########################################################

    if frb.host is not None and (
        not FRBFollowUpRequest.objects.filter(
            transient=frb,
            mode__in=['longslit','mask']).exists()) and (
        not FRBFollowUpObservation.objects.filter(
            transient=frb,
            success=True,
            mode__in=['longslit','mask']).exists()): 

        # Require top 2 P(O|x) > min(P_Ox_min)
        print(f"Need spec :POx_mins = {POx_mins}, {frb.sum_top_two_PATH}")
        POx_mins = frb_tags.values_from_tags(frb, 'min_POx')
        if (len(POx_mins) == 0) or (
            frb.sum_top_two_PATH > np.min(POx_mins)):
            frb.status = TransientStatus.objects.get(name='NeedSpectrum') 
            frb.save()
            return

    # #########################################################
    # Unassigned
    # #########################################################

    # If you get to here, you are unassigned

    frb.status = TransientStatus.objects.get(name='Unassigned')
    frb.save()
    return