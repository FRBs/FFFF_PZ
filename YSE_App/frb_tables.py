import numpy as np
import pandas 

from YSE_App.models import FRBTransient

from IPython import embed


def summary_table():
    """
    Generate a summary table of FRB transients.

    Returns:
        pandas.DataFrame: A DataFrame containing the summary information of FRB transients.
    """
    # Get it started
    all_frbs = FRBTransient.objects.all()
    all_tns = [frb.name for frb in all_frbs]
    frbs = pandas.DataFrame()
    frbs['TNS'] = all_tns

    # Add basic columns
    cols = ['ra', 'dec', 'a_err', 'b_err', 'theta', 'DM', 'DM_ISM', 'event_id', 'repeater', 'mw_ebv']
    for col in cols:
        frbs[col] = [getattr(frb, col) for frb in all_frbs]

    # Foreign keys
    fkeys = ['frb_survey', 'status']
    for key in fkeys:
        frbs[key] = [str(getattr(frb, key)) for frb in all_frbs]

    # Host and other Strings
    for col, key in zip(['Tags', 'Resources', 'Host'],
                        ['FRBTagsString', 
                         'FRBFollowUpResourcesString',
                         'HostString', 
                         ]):
        frbs[col] = [getattr(frb, key)() for frb in all_frbs]

    # Host 
    mags = [frb.host.path_mag if frb.host else np.nan for frb in all_frbs]
    frbs['Host_mag'] = mags
    POx = [frb.host.P_Ox if frb.host else np.nan for frb in all_frbs]
    frbs['POx'] = POx

    # Redshifts
    z = [frb.host.redshift if frb.host else np.nan for frb in all_frbs]
    frbs['z'] = z

    z_qual = [frb.host.redshift_quality if frb.host else -1 for frb in all_frbs]
    z_qual = [-1 if item is None else item for item in z_qual]
    frbs['z_qual'] = z_qual

    z_src = [frb.host.redshift_source if frb.host else '' for frb in all_frbs]
    z_src = ['' if item is None else item for item in z_src]
    frbs['z_src'] = np.array(z_src)


    # Top two candidates
    cand_poxs = [get_top_two_pox_gal_attr(frb,attr="P_Ox") for frb in all_frbs]
    frbs['cand_POx'] = cand_poxs

    cand_gal_names = [get_top_two_pox_gal_attr(frb,attr="name") for frb in all_frbs]
    frbs['cand_gal_names'] = cand_gal_names

    cand_gal_redshifts = [get_top_two_pox_gal_attr(frb,attr="redshift") for frb in all_frbs]
    frbs['cand_gal_redshifts'] = cand_gal_redshifts

    # Return
    return frbs


def get_gal_attr_from_qs(qs,attr="name"):
    """
    Given a QuerySet of galaxies, return a list of their attributes.

    Parameters:
    qs (QuerySet): A QuerySet of galaxy objects.

    Returns:
        list: A list of galaxy attribute.
    """
    if attr == "name":
        default_val = ""

    elif attr in ["redshift","P_Ox","path_mag"]:
        default_val = np.nan

    else:
        default_val = None


    return [getattr(gal,attr,default_val) for gal in qs]


def get_top_two_pox_gal_attr(frb_obj,attr="name"):
    """
    Given an FRBTransient object, return the names of the top two galaxies based on P_Ox.

    Parameters:
    frb_obj (FRBTransient): An FRBTransient object.

    Returns:
        list: A list of the names of the top two galaxies based on P_Ox.
    """
    if not frb_obj.host:
        return []
    
    path_values = frb_obj.get_Path_values()[0]
    # Get the indices of the top two P_Ox values
    top_two_indices = np.argsort(path_values)[-2:][::-1]


    if attr == "P_Ox":
        gal_attr = [path_values[i] for i in top_two_indices]

    else:        
        # Get the corresponding galaxy attributes
        gal_qs = [frb_obj.get_Path_values()[1][i] for i in top_two_indices]

        gal_attr = get_gal_attr_from_qs(gal_qs,attr=attr)

    return gal_attr
