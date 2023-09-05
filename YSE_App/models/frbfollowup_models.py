""" Models for FRB Followup Resources """
from autoslug import AutoSlugField
from auditlog.registry import auditlog

from django.db import models
from django.dispatch import receiver

import numpy as np
import pandas

from YSE_App.models.base import BaseModel
from YSE_App.models.frbtransient_models import *
from YSE_App.models.instrument_models import *
from YSE_App.models.principal_investigator_models import *
from YSE_App import frb_targeting


class FRBFollowUpResource(BaseModel):
    """ FRBFollowUpResource model

    Defines a follow-up resource for observing FRBs
    """

    ### Entity relationships ###
    # Required
    name = models.CharField(max_length=64, unique=True)
    instrument = models.ForeignKey(Instrument, on_delete=models.CASCADE)
    valid_start = models.DateTimeField() # Start of observing run (UT)
    valid_stop = models.DateTimeField() # End of observing run (UT)

    # Number or targets for each observing mode [imaging, masks, longslit]
    num_targ_img = models.IntegerField()
    num_targ_mask = models.IntegerField()
    num_targ_longslit = models.IntegerField()

    # Maximum AM for the observing
    max_AM = models.FloatField()

    # Surveys to which this resource is applicable
    #  e.g. CHIME/FRB, CRAFT, etc.
    #  all is an allowed value
    frb_surveys = models.CharField(max_length=64)


    # Optional

    # Classical, ToO, Queue 
    obs_type = models.CharField(max_length=64, null=True, blank=True)
    # PI
    PI = models.ForeignKey(PrincipalInvestigator, on_delete=models.CASCADE, null=True, blank=True)
    # Program name
    Program = models.CharField(max_length=64, null=True, blank=True)
    # Semester, period, etc.,  e.g. 2024A
    Period = models.CharField(max_length=64, null=True, blank=True)

    # Optional items for selecting targets
    min_POx = models.FloatField(null=True, blank=True)
    frb_tags = models.CharField(max_length=64, null=True, blank=True)
    frb_statuses = models.CharField(max_length=64, null=True, blank=True)

    slug = AutoSlugField(null=True, default=None, unique=True, populate_from='name')

    def __str__(self):
        return f'Name: {self.name} Instrument: {self.instrument.name}, start={self.valid_start}, stop={self.valid_stop}'

    def natural_key(self):
        return f'Instrument: {self.instrument.name}, start={self.valid_start}, stop={self.valid_stop}'

    def NameString(self):
        return self.name

    def InstrString(self):
        return self.instrument.name

    def StartString(self):
        return self.valid_start.strftime('%Y-%b-%d')

    def StopString(self):
        return self.valid_stop.strftime('%Y-%b-%d')


    def get_frbs_by_mode(self):
        """ Generate a dict of valid FRBs for targeting by observing mode

        calls frb_targeting.grab_targets_by_mode()

        Returns:
            dict: dict of FRBTransients by observing mode
        """
        # Cut down transients by selection criteria
        #  magnitude
        #  FRB Survey
        #  P(O|x)
        #  E(B-V)
        #  Bright star?
        gd_frbs = frb_targeting.targetfrbs_for_fu(self)

        # Caculate minimum airmasses during the Resource period
        min_AM = frb_targeting.calc_airmasses(self, gd_frbs)

        # Parse on AM
        keep_frbs = []
        for ss in range(len(gd_frbs)):
            if min_AM[ss] <= self.max_AM:
                keep_frbs.append(gd_frbs[ss].id)
        gd_frbs = gd_frbs.filter(id__in=keep_frbs)

        # Now the various modes!
        frbs_by_mode = frb_targeting.grab_targets_by_mode(self, gd_frbs)

        return frbs_by_mode

    def generate_target_table(self):
        """ Generate a table of FRBTransients that are valid for this resource, 
        i.e. meet all the criteria including AM

        Returns:
            pandas.DataFrame: table of FRBTransients for observing
        """
        # All FRBs by mode satisfying criteria
        frbs_by_mode = self.get_frbs_by_mode()

        # Cut down by number requested and priority
        final_targs_by_mode = frb_targeting.select_with_priority(self, frbs_by_mode)

        # Table me
        tbls = []
        for key in final_targs_by_mode.keys():
            tbls.append(frb_targeting.target_table_from_frbs(final_targs_by_mode[key], key))

        # Combine
        target_table = pandas.concat(tbls, ignore_index=True)

        # Add in the resource
        target_table['Resource'] = self.name

        return target_table