#
# File of constants to do with timing and past runs.
#

import datetime

# Some fixed dates needed by utimer. Put them here so
# that they are only computed once.
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

# Bias level change dates and associated before and after levels
# 3 changes imply 4 bias levels.
BIAS_CHANGES    = (52450,52850,53900)

# for each readout mode there are len(BIAS_CHANGES) + 1 sets of values.
# Each set gives typical bias levels for left and right of each CCD. There
# are not values in all cases.
BIAS_LEVELS = {'cdd' : (
        ((1975,2050),(1825,1848),(1830,2020)),
        ((2125,2180),(1985,2040),(2290,2190)),
        ((2125,2180),(1985,2040),(2390,2190)),
        ((2245,2300),(2150,2210),(2470,2360)),
        ((1635,1695),(1532,1562),(1835,1685))
        ),
               'fbb' : (
        None,
        None,
        ((2660,2720),(2350,2330),(2900,2840)),
        ((1620,1690),(1280,1252),(1827,1685))
        ),
               'fdd' : (
        ((2125,2180),(1985,2040),(2390,2190)),
        None,
        None,
        None
        )}

# ULTRACAM Timing parameters from Vik
INVERSION_DELAY = 110.          # microseconds
HCLOCK          = 0.48          # microseconds
CDS_TIME_FDD    = 2.2           # microseconds
CDS_TIME_FBB    = 4.4           # microseconds
CDS_TIME_CDD    = 10.           # microseconds
SWITCH_TIME     = 1.2           # microseconds

USPEC_FT_ROW    = 14.4e-6       # seconds
USPEC_FT_OFF    = 49.e-6        # seconds
USPEC_CLR_TIME  = 0.0309516     # seconds

# Run dates (for checking for spurious times).
# Format: overall run id, start, stop.
RUN_DATES = (('2002-05','2002-05-16','2002-05-20'),
             ('2002-09','2002-09-09','2002-09-14'),
             ('2002-09','2002-09-16','2002-09-17'),
             ('2002-09','2002-09-19','2002-09-21'),
             ('2003-05','2003-05-19','2003-05-25'),
             ('2003-05','2003-06-06','2003-06-08'),
             ('2003-11','2003-10-29','2003-11-06'),
             ('2003-11','2003-11-10','2003-11-13'),
             ('2004-05','2004-04-29','2004-04-29'),
             ('2004-05','2004-05-02','2004-05-05'),
             ('2004-05','2004-05-17','2004-05-19'),
             ('2004-08','2004-08-20','2004-08-30'),
             ('2005-05','2005-05-04','2005-05-21'),
             ('2005-08','2005-08-09','2005-08-15'),
             ('2005-08','2005-08-25','2005-09-01'),
             ('2005-11','2005-11-23','2005-11-28'),
             ('2006-03','2006-03-01','2006-03-12'),
             ('2007-06','2007-06-08','2007-06-23'),
             ('2007-10','2007-10-16','2007-10-28'),
             ('2007-11','2007-11-19','2007-11-24'),
             ('2008-08','2008-08-04','2008-08-11'),
             ('2009-01','2008-12-31','2009-01-06'),
             ('2010-01','2010-01-05','2010-01-07'),
             ('2010-04','2010-04-20','2010-06-15'),
             ('2010-12','2010-11-09','2010-12-17'),
             ('2010-12','2011-01-06','2011-01-11'),
             ('2010-12','2011-01-14','2011-01-14'),
             ('2010-12','2011-01-16','2011-01-18'),
             ('2011-05','2011-04-21','2011-04-27'),
             ('2011-05','2011-05-18','2011-06-01'),
             ('2011-08','2011-08-15','2011-08-21'),
             ('2011-08','2011-08-25','2011-08-27'),
             ('2011-10','2011-10-29','2011-11-03'),
             ('2012-01','2012-01-09','2012-01-22'),
             ('2012-01','2012-01-28','2012-02-05'),
             ('2012-04','2012-04-24','2012-04-29'),
             ('2012-09','2012-08-31','2012-09-13'),
             ('2012-10','2012-10-08','2012-10-16'),
             ('2013-04','2013-04-19','2013-04-24'),
             ('2013-07','2013-07-10','2013-07-16'),
             ('2013-07','2013-07-20','2013-07-21'),
             ('2013-07','2013-07-25','2013-08-05'),
             ('2013-11','2013-11-05','2013-11-13'),
             )

RUN_TELS = {'2002-05' : 'WHT',
            '2002-09' : 'WHT',
            '2003-05' : 'WHT',
            '2003-11' : 'WHT',
            '2004-05' : 'WHT',
            '2004-08' : 'WHT',
            '2005-05' : 'VLT',
            '2005-08' : 'WHT',
            '2005-11' : 'VLT',
            '2006-03' : 'WHT',
            '2007-06' : 'VLT',
            '2007-10' : 'WHT',
            '2007-11' : 'WHT',
            '2008-08' : 'WHT',
            '2009-01' : 'WHT',
            '2010-01' : 'WHT',
            '2010-04' : 'NTT',
            '2010-12' : 'NTT',
            '2011-05' : 'NTT',
            '2011-08' : 'WHT',
            '2011-10' : 'WHT',
            '2012-01' : 'WHT',
            '2012-04' : 'WHT',
            '2012-09' : 'WHT',
            '2012-10' : 'WHT',
            '2013-04' : 'WHT',
            '2013-07' : 'WHT',
            '2013-11' : 'TNT',
            }



# Integer type numbers for ucm files. Commented out ones
# are yet to be implemented
ITYPE_DOUBLE    =  0
#ITYPE_CHAR      =  1
ITYPE_INT       =  2
ITYPE_UINT      =  3
#ITYPE_LINT      =  4
#ITYPE_ULINT     =  5
ITYPE_FLOAT     =  6
ITYPE_STRING    =  7
ITYPE_BOOL      =  8
ITYPE_DIR       =  9
#ITYPE_DATE      = 10
ITYPE_TIME      = 11
#ITYPE_POSITION  = 12
ITYPE_DVECTOR   = 13
ITYPE_UCHAR     = 14
#ITYPE_TELESCOPE = 15
ITYPE_USINT     = 16
ITYPE_IVECTOR   = 17
ITYPE_FVECTOR   = 18

TNAME = {ITYPE_DOUBLE : 'double', ITYPE_INT : 'int', ITYPE_UINT : 'uint',
         ITYPE_FLOAT : 'float', ITYPE_STRING : 'string', ITYPE_BOOL : 'bool',
         ITYPE_DIR : 'directory', ITYPE_TIME : 'time', ITYPE_DVECTOR : 'dvector',
         ITYPE_UCHAR : 'uchar', ITYPE_USINT : 'usint', ITYPE_IVECTOR : 'ivector',
         ITYPE_FVECTOR : 'fvector'}

# ucm magic number
MAGIC           = 47561009

# Bit masks needed for Meinberg GPS data.
# See description in read_header.cc in pipeline for more
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
