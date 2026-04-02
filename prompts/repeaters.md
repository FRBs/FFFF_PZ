# Prompts for Repeaters

# FFFF-PZ code

## Add to FRBTransient model in frbtransient_models.py an optional column called "repeater_ids" which is a list of integers of the event IDs of the repeaters that the FRB is associated with.  Use a JSONField to store the list.  The default value should be an empty list.

# chime-ffff-pz code

## Generate a new script in chime-ffff-pz/scripts/update_repeaters.py to update information about repeaters in the FFFF-PZ database.  It should read in a table of repeaters with default name FRB_Repeaters.csv in the chime_ffff_pz/data/FRBs directory.  

### The table will have the following required columns:

- name (TNS name of FRB to )
- event_id (entry can be blank)
- ra
- dec
- a_err
- b_err
- theta
- DM
- frb_survey (always CHIME/FRB)
- repeater (True/False) – almost always True
- Rname 
- Comment (entry can be blank)

### It should also grab the full set of FRBs currently in FFFF-PZ.  I will then edit that file.

#### This failed because there was already an update_repeaters.py script.


