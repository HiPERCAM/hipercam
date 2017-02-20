"""
constants to do with timing and various changes to ULTRCAM
"""

import datetime

# Some fixed dates needed by utimer. Put them here so that they are only
# computed once.
DSEC            = 86400
MJD0            = datetime.date(1858,11,17).toordinal()
UNIX            = datetime.date(1970,1,1).toordinal()  - MJD0
DEFDAT          = datetime.date(2000,1,1).toordinal()  - MJD0
FIRST           = datetime.date(2002,5,16).toordinal() - MJD0
MAY2002         = datetime.date(2002,5,12).toordinal() - MJD0
SEP2002         = datetime.date(2002,9,8).toordinal()  - MJD0
TSTAMP_CHANGE1  = datetime.date(2003,8,1).toordinal()  - MJD0
TSTAMP_CHANGE2  = datetime.date(2005,1,1).toordinal()  - MJD0
TSTAMP_CHANGE3  = datetime.date(2010,3,1).toordinal()  - MJD0
USPEC_CHANGE    = datetime.date(2011,9,21).toordinal() - MJD0

# ULTRACAM timing parameters from Vik
INVERSION_DELAY = 110.          # microseconds
HCLOCK          = 0.48          # microseconds
CDS_TIME_FDD    = 2.2           # microseconds
CDS_TIME_FBB    = 4.4           # microseconds
CDS_TIME_CDD    = 10.           # microseconds
SWITCH_TIME     = 1.2           # microseconds

# ULTRASPEC timing parameters
USPEC_FT_ROW    = 14.4e-6       # seconds
USPEC_FT_OFF    = 49.e-6        # seconds
USPEC_CLR_TIME  = 0.0309516     # seconds

# Bit masks needed for Meinberg GPS data.  See description in read_header.cc
# in the ULTRACAM pipeline for more
PCPS_FREER            = 0x01   # DCF77 clock running on xtal, GPS receiver has not verified its position
PCPS_DL_ENB           = 0x02   # daylight saving enabled
PCPS_SYNCD            = 0x04   # clock has sync'ed at least once after pwr up
PCPS_DL_ANN           = 0x08   # a change in daylight saving is announced
PCPS_UTC              = 0x10   # a special UTC firmware is installed
PCPS_LS_ANN           = 0x20   # leap second announced, (requires firmware rev. REV_PCPS_LS_ANN_...)
PCPS_IFTM             = 0x40   # the current time was set via PC, (requires firmware rev. REV_PCPS_IFTM_...)
PCPS_INVT             = 0x80   # invalid time because battery was disconn'd
PCPS_LS_ENB           = 0x0100 # current second is leap second
PCPS_ANT_FAIL         = 0x0200 # antenna failure
PCPS_UCAP_OVERRUN     = 0x2000 # events interval too short
PCPS_UCAP_BUFFER_FULL = 0x4000 # events read too slow
PCPS_IO_BLOCKED       = 0x8000 # access to microprocessor blocked
