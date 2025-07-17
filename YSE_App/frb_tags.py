""" Code related to FRB tags """

import numpy as np

from YSE_App import frb_status
from YSE_App import frb_utils

def add_frb_tags(transient, user):
    from YSE_App.models import FRBTag

    # Add new ones
    for tag_name in transient['tags'].split(','):
        tag = frb_utils.add_or_grab_obj(
            FRBTag, dict(name=tag_name), {}, user)
        transient.frb_tags.add(tag)

    # Set status
    #frb_status.set_status(transient)

    # Save me!
    transient.save()

def values_from_tags(frb, key:str, debug:bool=False):
    """ Grab a list of values for a given key from the tags
      of a given FRB

    Args:
        frb (FRBTransient): FRBTransient instance
        key (str): key to grab
        debug (bool, optional): Debug flag. Defaults to False.

    Returns:
        list: list of values for the key;  can be empty
    """
    # Hiding here to avoid circular import (I hope)
    from YSE_App.models import FRBSampleCriteria

    # Prep
    tag_names = [frb_tag.name for frb_tag in frb.frb_tags.all()]
    if debug:
        print(f"tag_names = {tag_names} for {frb.name} and key {key}")

    # Get all samples with the frb survey
    samples = FRBSampleCriteria.objects.filter(
        frb_survey=frb.frb_survey)

    # Loop through 
    values = []
    for sample in samples:
        if sample.name in tag_names and getattr(sample,key) is not None:
            values.append(getattr(sample,key))

    return values

def chk_all_criteria(frb):
    from YSE_App.models import FRBSampleCriteria

    criteria = {}
    criteria['sample'] = []
    criteria['bright_star'] = [] # True if the FRB has a bright star and it is to be enforced
    criteria['EBV'] = []  # True if the FRB has a low enough MW E(B-V)
    criteria['POx'] = []  # True if the FRB has a high enough P(O|x)
    criteria['run_public_PATH'] = []  # True if we intend to run public PATH
    criteria['ran_deep_PATH'] = []  # True if we ran PATH on deeper imaging
    criteria['N_POx'] = []
    criteria['PUx'] = [] # True if P(U|x) > max_PUx

    # Grab the Sample object
    for frb_tag in frb.frb_tags.all():
        # Grab the criteria
        sample = frb_utils.add_or_grab_obj(
            FRBSampleCriteria, dict(name=frb_tag.name), {})
        criteria['sample'].append(sample.name)

        # Bright star?
        if sample.apply_bright_star and frb.bright_star is not None and frb.bright_star:
            criteria['bright_star'].append(True)
        else:
            criteria['bright_star'].append(False)

        # Dust
        if frb.mw_ebv < sample.max_EBV:
            criteria['EBV'].append(True)
        else:
            criteria['EBV'].append(False)

        # Run Public PATH?
        criteria['run_public_PATH'].append(sample.run_public_path)


        # PATH items
        if frb.host is not None:
            POx_values, galaxies, _ = frb.get_Path_values()
            argsrt = np.argsort(POx_values)
            pri_gal = galaxies[argsrt[-1]]  # Primary galaxy

            # P(O|x)
            if sample.use_top_two:
                criteria['N_POx'].append(2)
                if frb.sum_top_two_PATH > sample.min_POx:
                    criteria['POx'].append(True)
                else:
                    criteria['POx'].append(False)
            else:
                criteria['N_POx'].append(1)
                primary_POx = np.max(POx_values)
                if primary_POx > sample.min_POx:
                    criteria['POx'].append(True)
                else:  
                    criteria['POx'].append(False)
            # P(U|x)
            if sample.max_PUx is not None and frb.P_Ux > sample.max_PUx:
                criteria['PUx'].append(True)
            else:
                criteria['PUx'].append(False)

            # Ran deep PATH?
            rfilter = pri_gal.FilterMagString()
            if 'Blanco' in rfilter or 'DECam' in rfilter or 'Pan-STARRS' in rfilter: # Public
                criteria['ran_deep_PATH'].append(False)
            else:
                criteria['ran_deep_PATH'].append(True)
        else:
            criteria['POx'].append(False)
            criteria['N_POx'].append(0)
            criteria['PUx'].append(False)
            criteria['ran_deep_PATH'].append(False)


    # Convert to arrays
    for key in criteria:
        criteria[key] = np.array(criteria[key])

    # Return
    return criteria