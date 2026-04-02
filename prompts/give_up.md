# Code related to the GiveUp status


## Please add a new model in frbtransient_models.py which will hold a Table of FRBs that have been given up on.

### The table should have the following columns:

- name
- reason
- date

### The model should be named FRBGiveUp and should inherit from BaseModel.

## Add code in frb_status.py to set the status to GiveUp if the FRB has been given up on.  Put it right before the Bright star check.

## Add a new method to data_utils.py to sync an input table of FRBs with the FRBGiveUp table.  Model it after the modify_frbs method.

## Add a script to chime-ffff-pz to sync the FRBGiveUp table with the input table.  Put it in chime_ffff_pz/scripts/sync_giveups.py and model it after the modify_frbs script.  You will need to add a new file to chime-ffff-pz/bin/ to call it.

## Have tbl_name be an optional argument to the script.  If it is not provided, then the script will read from a file named give_ups.csv in the chime_ffff_pz/data/ directory.

## Please expose the GiveUps table as a view in the API.  It should be accessible at the URL /api/giveups. 