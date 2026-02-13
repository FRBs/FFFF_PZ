""" Methods to init the DB """

import logging
import pandas

from django.db import IntegrityError
from django.db.models import ForeignKey

from YSE_App.models import TransientStatus
from YSE_App.models import ObservationGroup
from YSE_App.models import FRBSurvey
from YSE_App.frb_utils import add_or_grab_obj

from YSE_App import frb_status
from YSE_App import frb_utils
from YSE_App.models import FRBTransient
from YSE_App.models.frbtransient_models import FRBIngestionError
from YSE_App import frb_tags
from YSE_App.models import FRBSampleCriteria

from YSE_App.serializers import FRBSampleCriteriaSerializer

logger = logging.getLogger(__name__)


def init_obsgroups(user):
    """ Initialize the transient status table

    Args:
        user (): user
    """

    # Add into DB
    for obsgroup in ['DECaLS']:
        _ = add_or_grab_obj(ObservationGroup,
                        dict(name=obsgroup), {}, user)

def init_statuses(user):
    """ Initialize the transient status table

    Args:
        user (): user
    """

    # Delete all existing
    TransientStatus.objects.all().delete()

    # Add into DB
    for status in frb_status.all_status:
        _ = add_or_grab_obj(TransientStatus,
                        dict(name=status), {}, user)

def init_surveys(user):
    """ Initialize the FRBSurvey table 

    Args:
        user (): user
    """

    # Delete all existing
    FRBSurvey.objects.all().delete()

    # Add into DB
    for survey in ['CHIME/FRB']:
        _ = add_or_grab_obj(FRBSurvey,
                        dict(name=survey), {}, user)


def add_df_to_db(df_frbs:pandas.DataFrame, user,
                 delete_existing:bool=False):
    """ Add a pandas DataFrame of FRBs to the database

    Args:
        df_frbs (pandas.DataFrame): pandas DataFrame of FRBs
        user (_type_): autheticated user
        delete_existing (bool, optional): If True, delete any
            existing FRBs with the same TNS first. Defaults to False.

    Returns:
        tuple: (status_code, message) where status_code is 200 for success
            (even with partial failures) or 500 for complete failure.
            Message includes details of any failed FRBs.
    """

    # TNS names
    tns_names = [t.name for t in FRBTransient.objects.all()]

    # Track successes and failures
    dbtransients = []
    failed_frbs = []  # List of (name, error_message) tuples
    skipped_frbs = []  # List of names that already existed

    # Loop on me
    for ss in range(len(df_frbs)):

        transient = df_frbs.iloc[ss]
        frb_name = transient.get('name', f'Row_{ss}')

        try:
            # Nearly all of the following was taken from add_transient() in data_utils.py
            transientkeys = transient.keys()

            if transient['name'] in tns_names:
                t = FRBTransient.objects.get(name=transient['name'])
                if delete_existing:
                    t.delete()
                else:
                    logger.info(f"FRBTransient {transient['name']} already exists in the database. Skipping.")
                    skipped_frbs.append(transient['name'])
                    dbtransients.append(t)
                    continue

            transientdict = {'created_by_id':user.id,
                             'modified_by_id':user.id}

            for transientkey in transientkeys:
                if transientkey == 'transientphotometry' or \
                    transientkey == 'transientspectra' or \
                    transientkey == 'host' or \
                    transientkey == 'tags' or \
                    transientkey == 'gw' or \
                    transientkey == 'non_detect_instrument' or \
                    transientkey == 'internal_names': continue
                if not isinstance(FRBTransient._meta.get_field(transientkey), ForeignKey):
                    if transient[transientkey] is not None: transientdict[transientkey] = transient[transientkey]
                else:
                    fkmodel = FRBTransient._meta.get_field(transientkey).remote_field.model
                    fk = fkmodel.objects.filter(name=transient[transientkey])
                    try:
                        transientdict[transientkey] = fk[0]
                    except IndexError:
                        logger.warning(f'{frb_name}: Bad FK key "{transientkey}" - no matching object found')

            # Build it
            dbtransient = FRBTransient(**transientdict)

            # Set null status
            dbtransient.status = TransientStatus.objects.get(name='Unassigned')

            # Save me!
            dbtransient.save()

            # Tags
            if hasattr(transient, 'tags'):
                frb_tags.add_frb_tags(dbtransient, transient['tags'], user)

            # Add to list
            dbtransients.append(dbtransient)
            logger.info(f"Successfully ingested {frb_name}")

        except IntegrityError as e:
            # Handle duplicate event_id or other integrity errors
            error_msg = str(e)
            if 'Duplicate entry' in error_msg:
                # Try to extract the duplicate value (e.g., event_id)
                # Error format: "Duplicate entry 'VALUE' for key 'KEY'"
                import re
                match = re.search(r"Duplicate entry '([^']+)' for key '([^']+)'", error_msg)
                if match:
                    dup_value, dup_key = match.groups()
                    # Try to find the existing FRB with this duplicate value
                    if 'event_id' in dup_key:
                        try:
                            existing = FRBTransient.objects.get(event_id=int(dup_value))
                            failed_frbs.append((frb_name, f"Duplicate event_id ({dup_value}) - already used by {existing.name}"))
                        except:
                            failed_frbs.append((frb_name, f"Duplicate event_id ({dup_value})"))
                    else:
                        failed_frbs.append((frb_name, f"Duplicate {dup_key}: {dup_value}"))
                else:
                    failed_frbs.append((frb_name, f"Duplicate entry error"))
            else:
                failed_frbs.append((frb_name, f"Database integrity error: {error_msg}"))
            logger.error(f"IntegrityError for {frb_name}: {error_msg}")
            continue

        except FRBIngestionError as e:
            # Handle IRSA/DM_ISM failures - FRB was already deleted by the signal handler
            error_msg = str(e)
            failed_frbs.append((frb_name, error_msg))
            logger.error(f"FRBIngestionError for {frb_name}: {error_msg}")
            continue

        except Exception as e:
            # Catch any other unexpected exceptions
            error_msg = str(e)
            failed_frbs.append((frb_name, error_msg))
            logger.error(f"Unexpected error ingesting {frb_name}: {error_msg}")
            continue

    # Build structured response
    total = len(df_frbs)
    succeeded_count = len(dbtransients) - len(skipped_frbs)
    skipped_count = len(skipped_frbs)
    failed_count = len(failed_frbs)

    # Get list of successfully ingested FRB names
    succeeded_frbs = [t.name for t in dbtransients if t.name not in skipped_frbs]

    # Build response dictionary
    response_data = {
        "total": total,
        "succeeded": succeeded_count,
        "skipped": skipped_count,
        "failed": failed_count,
        "succeeded_frbs": succeeded_frbs,
        "skipped_frbs": skipped_frbs,
        "failed_frbs": [{"name": name, "error": err} for name, err in failed_frbs],
    }

    # Build human-readable message for backwards compatibility
    msg_parts = [f"Ingestion complete: {succeeded_count}/{total} succeeded"]
    if skipped_count > 0:
        msg_parts.append(f"{skipped_count} skipped (already existed)")
    if failed_count > 0:
        msg_parts.append(f"{failed_count} failed")
    response_data["message"] = ". ".join(msg_parts)

    # Return 200 even with partial failures so the client gets the detailed message
    # Return 500 only if ALL FRBs failed
    if succeeded_count == 0 and failed_count > 0:
        return 500, response_data

    return 200, response_data

