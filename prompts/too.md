# ToO code

## FFFF-PZ code

### Claude, please create a new special tag labeled "CHIME-ToO" which can be manually given to an FRB, which will override statuses such as "BrightStar", "TooDusty", "AmbiguousHost" so that it will enter usual follow-up status "NeedSpectrum" or "NeedImage". Perhaps this can be done by checking for this tag before any criteria are checked and proceeding to "NeedSpectrum" or "NeedImage" in frb_status.py. This tag will also have high priority, so it is associated with a very large weight when assigning probability for targeting.

## chime-ffff-pz code

### Generate a new script in chime-ffff-pz/scripts/update_too.py to add CHIME-ToO tags for a table of FRBs into FFFF-PZ.  It should read in a table of FRBs and add the CHIME-ToO tag to each FRB.  The script should be called from the command line witht an optional argument tbl_name.  If tbl_name is not provided, then the script will read from the file FRB_ToOs.csv in the chime_ffff_pz/data/ directory.  The script should use the update_tags method in data_utils.py to add the CHIME-ToO tag to each FRB.  The arguments to the script are:

- --tbl_name: The name of the table to read in.  Default is FRB_ToOs.csv.
- --keep_original: If True, the original tags will be kept.  Default is True.

#### Let me approve the edits as you go.