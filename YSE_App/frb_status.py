""" Code related to the Status of FRBs """

import numpy as np


from YSE_App import frb_tags
from YSE_App.chime import tags as chime_tags

from IPython import embed

# Create the FRB ones
all_status = [\
    'Unassigned', # Does not meet the criteria for FFFF FollowUp
    'RunPublicPATH', # Needs to be run through PATH with public data
        # P_Ux is None
        # At least one frb_tag is in the list of run_public_path entries below
    'NeedImage', # Needs deeper imaging
        # P_Ux > maximum of all P_Ux_max for the frb_tags
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
    'TooDusty', # Sightline exceeds E(B-V) threshold
    'TooFaint', # Host is too faint for spectroscopy
        # r-magnitude (or equivalent; we use the PATH band) of the top host candidate
        #   is fainter than the maximum(mr_max) for the sample/surveys
    'AmbiguousHost',  # Host is considered too ambiguous for further follow-up
        # At least one of the frb_tags has a min_POx value
        #  and the sum of the top two P(O|x) is less than the minimum of those
    'UnseenHost',  # Even with deep imaging, no compelling host was found
        # P(U|x) is set
        # Deep imaging must exist.  The list of telescope+intrument is below
        # P(U|x) > P_Ux_max
    'Redshift', # Redshift measured
        # P(O|x) of top 2 > P_Ox_min
        # If Primary does not exceed min_POx, then the top two redshifts must be nearly the same
        # Else, take primary
    'Ambiguousz',  
        # The top candidate has P(O|x) < P_Ox_min
        # but the top two candidates have Sum(P(O|x)) > P_Ox_min
        #  But their redshifts are not consistent at 0.003
]

# List of telescope+instruments that are considered Deep
deep_telinstr = []

# Good redshift sources
good_z_sources = ['FFFF', 'Keck', 'Lick', 'Gemini']


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

    # Run in reverse order of completion

    # Are top 2 P(O|x) > min(P_Ox_min)
    POx_satisfied_two = False  # Sum of top 2 exceed min_POx
    POx_satisfied_primary = False # Primary exceeds min_POx
    PATH_run = False
    if frb.host is not None:
        POx_mins = frb_tags.values_from_tags(frb, 'min_POx')
        if len(POx_mins) > 0: 
            PATH_run = True
            # Sum of two?
            if frb.sum_top_two_PATH > np.min(POx_mins):
                POx_satisfied_two = True
            # Primary?
            POx_values, galaxies, _ = frb.get_Path_values()
            primary_POx = np.max(POx_values)
            if primary_POx > np.min(POx_mins):
                POx_satisfied_primary = True

    # #########################################################
    # #########################################################
    # Too Dusty??
    # #########################################################
    if frb.mw_ebv is not None:
        ebv_maxs = frb_tags.values_from_tags(frb, 'max_EBV')
        if len(ebv_maxs) > 0:
            if frb.mw_ebv > np.min(ebv_maxs):
                frb.status = TransientStatus.objects.get(name='TooDusty')
                frb.save()
                return

    # #########################################################
    # Ambiguous host
    # #########################################################

    if frb.host is not None:
        # 
        if PATH_run and not POx_satisfied_two:
            frb.status = TransientStatus.objects.get(name='AmbiguousHost') 
            frb.save()
            return


    # #########################################################
    # Unseen host
    # #########################################################

    # In PATH table?
    path_qs = Path.objects.filter(transient=frb)
    if len(path_qs) > 0 and frb.P_Ux is not None:
        # Check on source of photometry (instrument) for all the galaxies
        all_telinstr = []
        for path in path_qs:
            # Grab a dict of the photometry
            phot_dict = path.galaxy.phot_dict
            all_telinstr += list(phot_dict.keys())

        # Too shallow?
        too_shallow = True
        for uni_telins in np.unique(all_telinstr):
            if uni_telins in deep_telinstr:
                too_shallow = False

        # P(U|x) too large?
        if not too_shallow:
            PUx_maxs = frb_tags.values_from_tags(frb, 'max_P_Ux')
            if len(PUx_maxs) > 0:
                # Use the max
                PUx_max = np.max(PUx_maxs)
                if frb.P_Ux > PUx_max:
                    frb.status = TransientStatus.objects.get(
                        name='UnseenHost')
                    frb.save()
                    return

    # #########################################################
    # Redshift?
    # #########################################################
    if frb.host is not None and frb.host.redshift is not None and POx_satisfied_two:
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
                    frb.status = TransientStatus.objects.get(name='Ambiguousz') 
                    frb.save()
                    return log_message
                else:
                    log_message += "I AM CONSINSTENT-"
                    frb.status = TransientStatus.objects.get(name='Redshift')
                    frb.save()
                    return log_message

        '''
        # Primary satisifies P(O|x) > min(P_Ox_min)
        if POx_satisfied_primary: 
            gal = galaxies[argsrt[-1]]
            if gal.redshift is None:
                has_redshift = False
            for gd_source in ['FFFF', 'Keck', 'Lick', 'Gemini']:
                if gd_source in frb.host.redshift_source:
                    source_ok = True
        else: # Primary does not satisify P(O|x) > min(P_Ox_min), but top 2 do
            for idx in argsrt[-2:]:
                gal = galaxies[idx]
                # TODO -- consider checking the redshift_quality
                if gal.redshift is None:
                    has_redshift = False
            # Check the redshifts are nearly the same
            if has_redshift:
                if np.abs(galaxies[argsrt[-1]].redshift - galaxies[argsrt[-2]].redshift) > 0.003:
                    #has_redshift = False
                    frb.status = TransientStatus.objects.get(name='Ambiguousz') 
                    frb.save()
                    return

        # Require redshift come from our measurement or was vetted
        source_ok = False
        for gd_source in ['FFFF', 'Keck', 'Lick', 'Gemini']:
            if gd_source in frb.host.redshift_source:
                source_ok = True

        # Do it?
        if has_redshift and source_ok:
            frb.status = TransientStatus.objects.get(name='Redshift')
            frb.save()
            return
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
        POx_mins = frb_tags.values_from_tags(frb, 'min_POx')
        print(f"Need spec :POx_mins = {POx_mins}, {frb.sum_top_two_PATH}")
        if (len(POx_mins) == 0) or (
            frb.sum_top_two_PATH > np.min(POx_mins)):
            frb.status = TransientStatus.objects.get(name='NeedSpectrum') 
            frb.save()
            return

    # #########################################################
    # Need Image
    # #########################################################

    if frb.P_Ux is not None and frb.host is not None and (
        not FRBFollowUpRequest.objects.filter(
            transient=frb,
            mode='image').exists()) and (
        not FRBFollowUpObservation.objects.filter(
            transient=frb,
            success=True,
            mode='image').exists()):

        # Require top P_Ux > min(P_Ux_max)
        PUx_maxs = frb_tags.values_from_tags(frb, 'max_P_Ux')
        if (len(PUx_maxs) == 0) or (
            frb.P_Ux > np.min(PUx_maxs)):
            frb.status = TransientStatus.objects.get(name='NeedImage')
            frb.save()
            return

    # #########################################################
    # Run Public PATH
    # #########################################################

    if frb.P_Ux is None:
        path_flags = frb_tags.values_from_tags(frb, 'run_public_path')
        if np.any(path_flags):
            frb.status = TransientStatus.objects.get(name='RunPublicPATH')
            frb.save()
            return


    # #########################################################
    # Unassigned
    # #########################################################

    # If you get to here, you are unassigned

    frb.status = TransientStatus.objects.get(name='Unassigned')
    frb.save()
    return