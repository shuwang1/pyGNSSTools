import calendar, time,  os, base64, socket, ssl, sys, datetime, GnssUtil, json
import urllib
import ftplib
from math import sin, cos, atan2, sqrt, acos, floor
import numpy as np
#import urllib2 

MAX_GPS_SATS        = 33
ICD200PI            = 3.1415926535898
PI_DEG              = ICD200PI/180
DAY_SECONDS         = 86400
WEEK_SECONDS        = 604800
HALF_WEEK_SECONDS   = 302400

C0          = 299792458  #GPS official speed of light.
L1_CODE_LEN = 1023
L1_FREQ     = 1575.42e6
L1_CHIP_HZ  = 1.023e6 #Hz C/A Chipping Frequency.
L1_CHIP_M   = C0 / L1_CHIP_HZ
L1_CODE_M   = L1_CODE_LEN * L1_CHIP_M
L1_BIT_M    = L1_CODE_M * 10

'''
#####################################################################
# As of June 30 2015, and until the leap second of December 31 2016
#        TAI is ahead of UTC   by 36 seconds.
#        TAI is ahead of GPS   by 19 seconds.
#        GPS is ahead of UTC   by 17 seconds.
# After December 2016,
#        TAI is ahead of UTC   by 37 seconds.
#        TAI is ahead of GPS   by 19 seconds.
#        GPS is ahead of UTC   by 18 seconds.
######################################################################
'''
LEAPSECONDS_2017 = 18

'''
#######################################################################################################
Predict satellite state ( position & velocity )
The parameters needed for computing the broadcast orbits are:

toe, TimeOfApplicability -- Time of the applicability
SqrtA -- Square root of the semi-major axis of the satellite orbit
M0, meanAnomaly -- mean anomaly at a reference time
OMEGA0, asceLon -- Longitude of ascending node of orbit plane at weekly epoch
dotOMEGA, dotAsceLon -- Rate of right ascension
i0, OrbitalInclination -- Inclination Angle at reference time
omega, ArgumentOfPerigee -- Argument of Perigee

dotInc,  -- Rate of inclination angle
dMotion, deltaMotion -- Mean Motion Difference
Cis -- Amplitude of the Sine Harmonic Correction term to the Angle of Inclination
Cic -- Amplitude of the Cosine Harmonic Correction term to the Angle of Inclination
Cus -- Amplitude of the Sine Harmonic Correction term to the argument of latitude
Cuc -- Amplitude of the Cosine Harmonic Correction term to the argument of latitude
Crs --  Amplitude of the Sine Harmonic Correction Term to the orbit radius
Crc --  Amplitude of the Cosine Harmonic Correction Term to the orbit radius

Reference:
Benjamin W. Remondi, Computing satellite velocity using the broadcast ephemeris, GPS Solution (2004) 8:181-183
Xiaofan Li, A Study on the GPS Satellites Orbit, ASEN 5050 Semester Project, http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
#################################################################################################################
'''

###################################
# Define some constant parameters
###################################
WGS84_a_m       = 6378137  # Semi-major axis radius of the Earth r = 6378.137 km
WGS84_a2_m2     = 6378137 * 6378137
WGS84_b_m       = 6356752.3142  # Semi-minor axis radius of the Earth r = 6356.7523142 km
WGS84_e         = 0.0818191908426  # Eccentricity 0.0818191908426
WGS84_e2        = 0.0818191908426 * 0.0818191908426
WGS84_u_km3s2   = 398600.4418  # Gravitational Parameter ? = 398600.4418 km3/s2
WGS84_u_m3s2    = 3.986004418e+14  # Gravitational Parameter ? = 398600.4418 km3/s2
WGS84_f         = 1/298.257223563  # Flattening f = 1/298.257223563
WGS84_w_rads    = 7.292158553e-5  # Rotation rate of the Earth ? =7.292158553e-5 rad/s

F = -4.442807633e-10   # What is this paramter?  I need find it out.

####################################################################
#/*                     CRC24 LOOKUP TABLE                        */
####################################################################
CRC24 =  (
 0x000000, 0x864CFB, 0x8AD50D, 0x0C99F6, 0x93E6E1, 0x15AA1A, 0x1933EC, 0x9F7F17,
 0xA18139, 0x27CDC2, 0x2B5434, 0xAD18CF, 0x3267D8, 0xB42B23, 0xB8B2D5, 0x3EFE2E,
 0xC54E89, 0x430272, 0x4F9B84, 0xC9D77F, 0x56A868, 0xD0E493, 0xDC7D65, 0x5A319E,
 0x64CFB0, 0xE2834B, 0xEE1ABD, 0x685646, 0xF72951, 0x7165AA, 0x7DFC5C, 0xFBB0A7,
 0x0CD1E9, 0x8A9D12, 0x8604E4, 0x00481F, 0x9F3708, 0x197BF3, 0x15E205, 0x93AEFE,
 0xAD50D0, 0x2B1C2B, 0x2785DD, 0xA1C926, 0x3EB631, 0xB8FACA, 0xB4633C, 0x322FC7,
 0xC99F60, 0x4FD39B, 0x434A6D, 0xC50696, 0x5A7981, 0xDC357A, 0xD0AC8C, 0x56E077,
 0x681E59, 0xEE52A2, 0xE2CB54, 0x6487AF, 0xFBF8B8, 0x7DB443, 0x712DB5, 0xF7614E,
 0x19A3D2, 0x9FEF29, 0x9376DF, 0x153A24, 0x8A4533, 0x0C09C8, 0x00903E, 0x86DCC5,
 0xB822EB, 0x3E6E10, 0x32F7E6, 0xB4BB1D, 0x2BC40A, 0xAD88F1, 0xA11107, 0x275DFC,
 0xDCED5B, 0x5AA1A0, 0x563856, 0xD074AD, 0x4F0BBA, 0xC94741, 0xC5DEB7, 0x43924C,
 0x7D6C62, 0xFB2099, 0xF7B96F, 0x71F594, 0xEE8A83, 0x68C678, 0x645F8E, 0xE21375,
 0x15723B, 0x933EC0, 0x9FA736, 0x19EBCD, 0x8694DA, 0x00D821, 0x0C41D7, 0x8A0D2C,
 0xB4F302, 0x32BFF9, 0x3E260F, 0xB86AF4, 0x2715E3, 0xA15918, 0xADC0EE, 0x2B8C15,
 0xD03CB2, 0x567049, 0x5AE9BF, 0xDCA544, 0x43DA53, 0xC596A8, 0xC90F5E, 0x4F43A5,
 0x71BD8B, 0xF7F170, 0xFB6886, 0x7D247D, 0xE25B6A, 0x641791, 0x688E67, 0xEEC29C,
 0x3347A4, 0xB50B5F, 0xB992A9, 0x3FDE52, 0xA0A145, 0x26EDBE, 0x2A7448, 0xAC38B3,
 0x92C69D, 0x148A66, 0x181390, 0x9E5F6B, 0x01207C, 0x876C87, 0x8BF571, 0x0DB98A,
 0xF6092D, 0x7045D6, 0x7CDC20, 0xFA90DB, 0x65EFCC, 0xE3A337, 0xEF3AC1, 0x69763A,
 0x578814, 0xD1C4EF, 0xDD5D19, 0x5B11E2, 0xC46EF5, 0x42220E, 0x4EBBF8, 0xC8F703,
 0x3F964D, 0xB9DAB6, 0xB54340, 0x330FBB, 0xAC70AC, 0x2A3C57, 0x26A5A1, 0xA0E95A,
 0x9E1774, 0x185B8F, 0x14C279, 0x928E82, 0x0DF195, 0x8BBD6E, 0x872498, 0x016863,
 0xFAD8C4, 0x7C943F, 0x700DC9, 0xF64132, 0x693E25, 0xEF72DE, 0xE3EB28, 0x65A7D3,
 0x5B59FD, 0xDD1506, 0xD18CF0, 0x57C00B, 0xC8BF1C, 0x4EF3E7, 0x426A11, 0xC426EA,
 0x2AE476, 0xACA88D, 0xA0317B, 0x267D80, 0xB90297, 0x3F4E6C, 0x33D79A, 0xB59B61,
 0x8B654F, 0x0D29B4, 0x01B042, 0x87FCB9, 0x1883AE, 0x9ECF55, 0x9256A3, 0x141A58,
 0xEFAAFF, 0x69E604, 0x657FF2, 0xE33309, 0x7C4C1E, 0xFA00E5, 0xF69913, 0x70D5E8,
 0x4E2BC6, 0xC8673D, 0xC4FECB, 0x42B230, 0xDDCD27, 0x5B81DC, 0x57182A, 0xD154D1,
 0x26359F, 0xA07964, 0xACE092, 0x2AAC69, 0xB5D37E, 0x339F85, 0x3F0673, 0xB94A88,
 0x87B4A6, 0x01F85D, 0x0D61AB, 0x8B2D50, 0x145247, 0x921EBC, 0x9E874A, 0x18CBB1,
 0xE37B16, 0x6537ED, 0x69AE1B, 0xEFE2E0, 0x709DF7, 0xF6D10C, 0xFA48FA, 0x7C0401,
 0x42FA2F, 0xC4B6D4, 0xC82F22, 0x4E63D9, 0xD11CCE, 0x575035, 0x5BC9C3, 0xDD8538 )

PWR2_5        = pow( 2, -5 )
PWR2_19       = pow( 2, -19 )
PWR2_29       = pow( 2, -29 )
PWR2_31       = pow( 2, -31 )
PWR2_33       = pow( 2, -33 )
PWR2_43       = pow( 2, -43 )
PWR2_55       = pow( 2, -55 )
PWR2_31_PI    = PWR2_31 * ICD200PI
PWR2_43_PI    = PWR2_43 * ICD200PI

class Ephemeris( object ) :

   
    def __init__( self ):
        self.DEBUG = 0

        self.alpha0, self.alpha1, self.alpha2, self.alpha3 = 0, 0, 0, 0
        self.beta0, self.beta1, self.beta2, self.beta3 = 0, 0, 0, 0
        self.A0, self.A1, self.tot, self.wnt, self.dtls = 0, 0, 0, 0, 0        
        #############################################################################################
        # Initialize Ephemeris Variables
        #
        # SVNs are "space vehicle numbers", the serial numbers assigned to each GPS satellite, and 
        # PRNs are the "pseudo-random noise" sequences, or Gold codes, that each satellite transmits 
        # to differentiate itself from other satellites in the active constellation.
        #############################################################################################
        self.svid       = [ int(0) for sat in range( MAX_GPS_SATS ) ]

        self.YY         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.MM         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.DD         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.hh         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.mm         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.sec        = [ float(0) for sat in range( MAX_GPS_SATS ) ]
        
        self.week       = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.acc        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.l2code     = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.idot       = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.iode       = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.toc        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.af2        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_55 
        self.af1        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43 
        self.af0        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31
        self.iodc       = [ int(0) for sat in range( MAX_GPS_SATS ) ]     
        self.Crs        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_5
        self.deltan     = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.M0         = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Cuc        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29
        self.e          = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_33
        self.Cus        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29  
        self.sqrtA      = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_19
        self.toe        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.Cic        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29
        self.Omega0     = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Cis        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.i0         = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Crc        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_5
        self.omega      = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Omegadot   = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.tgd        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31
        self.health     = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.l2p        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.fit        = [ int(0) for sat in range( MAX_GPS_SATS ) ]

        self.localDir = os.getcwd()
        if not self.localDir.endswith( '/' ) :
            self.localDir += '/' 

        directory = os.path.dirname( self.localDir )
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        self.ephemerisFilename  = None
        self.NtripSessionTime   = 300
        return
    
    def OLD__init__( self, ephemerisFilename = None, mountpoint = "abc", host = 'www.igs-ip.net',  NtripSSL = False ) :
        self.DEBUG = 1

        self.alpha0, self.alpha1, self.alpha2, self.alpha3 = 0, 0, 0, 0
        self.beta0, self.beta1, self.beta2, self.beta3 = 0, 0, 0, 0
        self.A0, self.A1, self.tot, self.wnt, self.dtls = 0, 0, 0, 0, 0        
        #############################################################################################
        # Initialize Ephemeris Variables
        #
        # SVNs are "space vehicle numbers", the serial numbers assigned to each GPS satellite, and 
        # PRNs are the "pseudo-random noise" sequences, or Gold codes, that each satellite transmits 
        # to differentiate itself from other satellites in the active constellation.
        #############################################################################################
        self.svid       = [ int(0) for sat in range( MAX_GPS_SATS ) ]

        self.YY         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.MM         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.DD         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.hh         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.mm         = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.sec        = [ float(0) for sat in range( MAX_GPS_SATS ) ]
        
        self.week       = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.acc        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.l2code     = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.idot       = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.iode       = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.toc        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.af2        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_55 
        self.af1        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43 
        self.af0        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31
        self.iodc       = [ int(0) for sat in range( MAX_GPS_SATS ) ]     
        self.Crs        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_5
        self.deltan     = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.M0         = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Cuc        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29
        self.e          = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_33
        self.Cus        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29  
        self.sqrtA      = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_19
        self.toe        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.Cic        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_29
        self.Omega0     = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Cis        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.i0         = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Crc        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_5
        self.omega      = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31_PI
        self.Omegadot   = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_43_PI
        self.tgd        = [ float(0) for sat in range( MAX_GPS_SATS ) ]     # PWR2_31
        self.health     = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.l2p        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        self.fit        = [ int(0) for sat in range( MAX_GPS_SATS ) ]
        
        self.NtripSessionTime = 300
        self.NtripHost, self.NtripMountpoint, self.NtripUsername, self.NtripPassword = host, mountpoint , 'swang', 'SigmaDesigns' #username & password for RTCM correction service        
        self.NtripSSL = NtripSSL
        if self.NtripSSL:
            self.NtripPort = 443 
        else:
            self.NtripPort = 2101
        
        
        if ephemerisFilename is None :
            self.ephemerisFilename = "hour" + time.strftime( "%j%H.%yO.eph", time.gmtime() )
        elif os.path.isfile( ephemerisFilename ) :
            self.ephemerisFilename = ephemerisFilename
            self.NtripSessionTime /= 2
            ephemerisFile = open( self.ephemerisFilename, "r" )
            jj = 0
            for line in ephemerisFile :
                param = line.split()
                if len( param ) < 2 :
                    continue
                if 'svid:' == param[0] :
                    jj = int( param[1] ) - 1
                    self.svid[ jj ] = jj + 1
                elif 'week:' == param[0] :
                    self.week[jj] = int( param[1] )
                elif 'acc:' == param[0] :
                    self.acc[jj] = int( param[1] )
                elif 'l2code:' == param[0] :
                    self.l2code[jj] = int( param[1] )
                elif 'idot:' == param[0] :
                    self.idot[jj] = float( param[1] )
                elif 'iode:' == param[0] :
                    self.iode[jj] = int( param[1] )
                elif 'toc:' == param[0] :
                    self.toc[jj] = int( param[1] )
                elif 'af2:' == param[0] :
                    self.af2[jj] = float( param[1] )
                elif 'af1:' == param[0] :
                    self.af1[jj] = float( param[1] )
                elif 'af0:' == param[0] :
                    self.af0[jj] = float( param[1] )
                elif 'iodc:' == param[0] :
                    self.iodc[jj] = int( param[1] )
                elif 'Crs:' == param[0] :
                    self.Crs[jj] = float( param[1] )
                elif 'deltan:' == param[0] :
                    self.deltan[jj] = float( param[1] )
                elif 'M0:' == param[0] :
                    self.M0[jj] = float( param[1] )
                elif 'Cuc:' == param[0] :
                    self.Cuc[jj] = float( param[1] )
                elif 'e:' == param[0] :
                    #########################################################################################################
                    # The eccentricity of the Earth's orbit is currently about 0.0167; the Earth's orbit is nearly circular. 
                    # Venus and Neptune have even lower eccentricity. Over hundreds of thousands of years, the eccentricity 
                    # of the Earth's orbit varies from nearly 0.0034 to almost 0.058 as a result of gravitational attractions 
                    # among the planets
                    ###########################################################################################################
                    self.e[jj] = float( param[1] )
                elif 'Cus:' == param[0] :
                    self.Cus[jj] = float( param[1] )
                elif 'sqrtA:' == param[0] :
                    self.sqrtA[jj] = float( param[1] )
                elif 'toe:' == param[0] :
                    self.toe[jj] = int( param[1] )
                elif 'Cic:' == param[0] :
                    self.Cic[jj] = float( param[1] )
                elif 'Omega0:' == param[0] :
                    self.Omega0[jj] = float( param[1] )
                elif 'Cis:' == param[0] :
                    self.Cis[jj] = float( param[1] )
                elif 'i0:' == param[0] :
                    self.i0[jj] = float( param[1] )
                elif 'Crc:' == param[0] :
                    self.Crc[jj] = float( param[1] )
                elif 'omega:' == param[0] :
                    self.omega[jj] = float( param[1] )
                elif 'Omegadot:' == param[0] :
                    self.Omegadot[jj] = float( param[1] )
                elif 'tgd:' == param[0] :
                    self.tgd[jj] = float( param[1] )
                elif 'health:' == param[0] :
                    self.health[jj] = int( param[1] )
                elif 'l2p:' == param[0] :
                    self.l2p[jj] = int( param[1] )
                elif 'fit:' == param[0] :
                    self.fit[jj] = int( param[1] )
            ephemerisFile.close()
        else:
            raise Exception( "Ephemeris.__init__():: The input file %s doesn't exist." % ephemerisFilename )
        
        self.localDir = os.getcwd() 
        if not self.localDir.endswith( '/' ) :
            self.localDir += '/' 

        directory = os.path.dirname( self.localDir )
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        return

    def Save_Ephemeris( self, ephemerisFilename = None ):
        localDir = self.localDir
        if ephemerisFilename is None :
            ephemerisFilename = "hour" + time.strftime( "%j%H.%yO.eph", time.gmtime() )
       
        self.ephemerisFilename = ephemerisFilename + '.json'
        ephemeris_ = []
        for svid in self.svid :
            if svid == 0 :
                continue
            jj = svid - 1
            ephemeris_.append({
                'svid':     self.svid[jj],
                'week':     self.week[jj],
                'acc':      self.acc[jj],
                'l2code':   self.l2code[jj],
                'idot':     self.idot[jj],
                'iode':     self.iode[jj],
                'toc':      self.toc[jj],
                'af2':      self.af2[jj],
                'af1':      self.af1[jj],
                'af0':      self.af0[jj],
                'iodc':     self.iodc[jj],
                'Crs':      self.Crs[jj],
                'deltan':   self.deltan[jj],
                'M0':       self.M0[jj],
                'Cuc':      self.Cuc[jj],
                'e':        self.e[jj],
                'Cus':      self.Cus[jj],
                'sqrtA':    self.sqrtA[jj],
                'toe':      self.toe[jj],
                'Cic':      self.Cic[jj],
                'Omega0':   self.Omega0[jj],
                'Cis':      self.Cis[jj],
                'i0':       self.i0[jj],
                'Crc':      self.Crc[jj],
                'omega':    self.omega[jj],
                'Omegadot': self.Omegadot[jj],
                'tgd':      self.tgd[jj],
                'health':   self.health[jj],
                'l2p':      self.l2p[jj],
                'fit':      self.fit[jj]
            })
        with open( localDir + self.ephemerisFilename, 'w') as f:
            json.dump( ephemeris_, f )

        return self.ephemerisFilename

    def Get_Ephemeris_From_NASA_Rinex( self, host = "NASA", backDays = 0, server = "cddis.gsfc.nasa.gov", username = "anonymous", password = "guest" ) :
        
        localDir = self.localDir + 'stations/' + host + "/"
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        dt = datetime.datetime.utcnow() - datetime.timedelta( days = backDays ) # some days back from now
        dailyFtpDir = "/gnss/data/daily/%04d/brdc/" % dt.year 
        dailyFilename = "brdc%03d0.%02dn.Z" % ( dt.timetuple().tm_yday, (dt.year-2000) )
        hourlyFtpDir = "/gnss/data/hourly/%04d/%03d/" % ( dt.year, dt.timetuple().tm_yday ) 
        hourlyFilename = "hour%03d0.%02dn" % ( dt.timetuple().tm_yday, (dt.year-2000)  )
        
        zFilename = localDir + hourlyFilename + '.Z'
        if ( not os.path.exists( zFilename ) ) or ( os.stat( zFilename ).st_size == 0 ) :
            if os.path.exists( zFilename ) :
                os.remove( zFilename )
            ftp = ftplib.FTP( server )  # Establish the connection
            ftp.login( username, password )
            ftp.cwd( hourlyFtpDir ) # Change to the proper directory        
            fObj = open( zFilename, 'wb')
            ftp.retrbinary('RETR ' + hourlyFilename + '.Z', fObj.write)
            fObj.close()
            ftp.close()
            sys.stdout.flush()

        rinexFilename = localDir + hourlyFilename + '.rinex'
        if ( not os.path.exists( rinexFilename ) ) or ( os.stat( rinexFilename ).st_size == 0 ) :
            if os.path.exists( rinexFilename ) :
                os.remove( rinexFilename )
            rinexFObj = open( rinexFilename, 'wb')
            rinexFObj.write( GnssUtil.unlzw( open( zFilename, 'rb' ).read() ) )
            rinexFObj.close()

        with open( rinexFilename, 'rt' ) as rinexFObj :
            while True :
                line = rinexFObj.readline()
                if line[60:73] == "END OF HEADER" :
                    break
                if line[60:69] == "ION ALPHA" :
                    self.alpha0 = float( line[2:14].replace('D', 'E') )
                    self.alpha1 = float( line[14:26].replace('D', 'E') )
                    self.alpha2 = float( line[26:38].replace('D', 'E') )
                    self.alpha3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "ION BETA" :
                    self.beta0 = float( line[2:14].replace('D', 'E') )
                    self.beta1 = float( line[14:26].replace('D', 'E') )
                    self.beta2 = float( line[26:38].replace('D', 'E') )
                    self.beta3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "DELTA-UTC" :
                    self.A0 = float( line[3:22].replace('D', 'E') )
                    self.A1 = float( line[22:41].replace('D', 'E') )
                    self.tot = int( line[41:50] )
                    self.wnt = int( line[50:59] )
                elif line[60:72] == "LEAP SECONDS" :
                    self.dtls = int( line[0:6].replace('D', 'E') )

            line = rinexFObj.readline( )
            while line:
                ii = int( line[0:2] ) - 1
                self.svid[ii] = ii + 1
                self.YY[ii] = int( line[3:5] ) + 2000
                self.MM[ii] = int( line[6:8] )
                self.DD[ii] = int( line[9:11] )
                self.hh[ii] = int( line[12:14] )
                self.mm[ii] = int( line[15:17] )
                self.sec[ii] = float( line[18:22] )
                self.af0[ii] = float( line[22:41].replace('D', 'E')  )
                self.af1[ii] = float( line[41:60].replace('D', 'E')  )
                self.af2[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.iode[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Crs[ii] = float( line[22:41].replace('D', 'E')  )
                self.deltan[ii] = float( line[41:60].replace('D', 'E')  )
                self.M0[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.Cuc[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.e[ii] = float( line[22:41].replace('D', 'E')  )
                self.Cus[ii] = float( line[41:60].replace('D', 'E')  )
                self.sqrtA[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                t_toe = int( float( line[3:22].replace('D', 'E')  ) )
                if t_toe < self.toe[ii] :
                    print ( "Ephemeris.Get_Ephemeris_From_NASA_Rinex():: The new Toe %d for Sat %02d is older than the current Toe %d." % (t_toe, self.svid[ii], self.toe[ii] ) )
                self.toe[ii] = t_toe
                self.Cic[ii] = float( line[22:41].replace('D', 'E')  )
                self.Omega0[ii] = float( line[41:60].replace('D', 'E')  )
                self.Cis[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.i0[ii] = float( line[3:22].replace('D', 'E')  )
                self.Crc[ii] = float( line[22:41].replace('D', 'E')  )
                self.omega[ii] = float( line[41:60].replace('D', 'E')  )
                self.Omegadot[ii] = float( line[60:79].replace('D', 'E')  )

                line = rinexFObj.readline( )
                self.idot[ii] = float( line[3:22].replace('D', 'E')  )
                self.l2code[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.week[ii] = int( float( line[41:60].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                self.health[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.health[ii] = ( self.health[ii] % 32 + 32 ) if ( self.health[ii] < 32 ) else self.health[ii] 
                self.tgd[ii] = float( line[41:60].replace('D', 'E')  ) 
                self.iodc[ii] = int( float( line[60:79].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                line = rinexFObj.readline( )

        rinexFObj.close()

        return
        
    def Get_Ephemeris_From_Ntrip( self, version = 'R', hostIP = '132.239.154.80', hostPort = 2101, mountpoint = 'DSME', username = "CRTNERLANGNE", password = 'ERLANGNESURV', SSL = False ) :
        #
        # https://igs.bkg.bund.de/root_ftp/NTRIP/streams/streamlist_world-wide.htm
        #
        # *
        # * This class implements the reading of some Message Types
        # * defined in the RTCM 3.2 Standard, plus some utilities to handle messages.
        # *
        # * Generation of the following Message Types:
        # *   1001, 1002, 1003, 1004, 1005, 1006, 1008, 1019, 1020, 1029, 1045
        # *
        # * Decoding of the following Message Types:
        # *   1019, 1045
        # *
        # * RTCM 3 message format (size in bits):
        # *   +----------+--------+-----------+--------------------+----------+
        # *   | preamble | 000000 |  length   |    data message    |  parity  |
        # *   +----------+--------+-----------+--------------------+----------+
        # *   |<-- 8 --->|<- 6 -->|<-- 10 --->|<--- length x 8 --->|<-- 24 -->|
        # *   +----------+--------+-----------+--------------------+----------+        
        # *
        # RTCM 1019 GPS Ephemerides
        # Sets of these messages (one per SV) are used to send the broadcast orbits for GPS in a Kepler format.
        
        self.NtripHost, self.NtripPort, self.NtripMountpoint, self.NtripUsername, self.NtripPassword, self.NtripSSL = hostIP, hostPort, mountpoint , username, password, SSL #username & password for RTCM correction service        

        localDir = self.localDir + "/stations/" + mountpoint + '/'
        directory = os.path.dirname( localDir )
        if not os.path.exists( directory ) :
            os.makedirs( directory )
    
            '''
            ################################################
            # Ntrip client for receiving RTCM corrections.
            ################################################
            '''
            ##
            # Generate an encoding of the username:password for the service.
            # The string must be first encoded in ascii to be correctly parsed by the
            # base64.b64encode function.
            #

        #'Authorization: Basic {}\r\n'.format( base64.b64encode( ( '%s:%s'%(self.NtripUsername,self.NtripPassword) ).encode("utf-8") ) ) +\

        NtripMountString = 'GET /%s HTTP/1.1\r\n' % self.NtripMountpoint +\
        'Host: %s\r\n' % self.NtripHost +\
        'Ntrip-Version: Ntrip/%s\r\n' % str( 'R' ) +\
        'User-Agent: NTRIP BNC/2.12.3(MAC)\r\n' +\
        'Authorization: Basic {}\r\n'.format( base64.b64encode( ( '%s:%s'%(self.NtripUsername,self.NtripPassword) ).encode('ascii') ).decode('ascii') ) +\
        'Connection: close\r\n' +\
        'Accept-Encoding: gzip\r\n' +\
        'Accept-Language: en-US,*\r\n\r\n'

        ntripConnectRetry = 0
        maxConnectRetry = 10

        ntripRecvBufLen =  1024 * 64
        ntripRecvBuf = bytearray( ntripRecvBufLen )
        ntripRecvBufView = memoryview( ntripRecvBuf )

        RTCM3_PREAMBLE = 211
        RTCM_MSG1019_LENGTH = 61

        now, Tstop = 0, time.time() + self.NtripSessionTime / 4
        while now < Tstop :
            foundNtripHeader = False
            while ( ntripConnectRetry < maxConnectRetry ) and ( not foundNtripHeader ):
                try:
                    sent, totalSent, msgLen, casterResponse = 0, 0, len( NtripMountString ), ''
                    ntripSocket = socket.socket( socket.AF_INET, socket.SOCK_STREAM )
                    if self.NtripSSL:
                        ntripSocket = ssl.wrap_socket( ntripSocket )
                    retErrNo = ntripSocket.connect_ex( ( self.NtripHost, int(self.NtripPort) ) )
                    if retErrNo != 0:
                        if self.DEBUG:
                            print( 'Ephemeris.Get_Ephemeris_From_Ntrip():: ntripSocket.connect_ex() is connecting %d of %d: retErrNo = %d' % ( ntripConnectRetry, maxConnectRetry, retErrNo ) )
                        continue
                    else:
                        ntripSocket.settimeout( 60 )
                        while totalSent < msgLen:
                            sent = ntripSocket.send( NtripMountString[totalSent:].encode("utf-8") )
                            totalSent = totalSent + sent
                            if sent == 0:
                                raise RuntimeError("socket connection broken")
                        casterResponse = ntripSocket.recv( 4096 )

                        if casterResponse.startswith("STREAMTABLE".encode("utf-8")):
                            raise Exception("Ephemeris.Get_Ephemeris_From_Ntrip():: Invalid or No Mountpoint")
                        elif casterResponse.startswith("HTTP/1.1 401 Unauthorized".encode("utf-8")):
                            sys.stderr.write('Username and/or password are not correct.\n')
                            sys.exit(2)
                        elif casterResponse.startswith("HTTP/1.1 404 Not Found".encode("utf-8")):
                            print ( "Ephemeris.Get_Ephemeris_From_Ntrip():: Mount Point %s does not exist\n" % self.NtripMountpoint )
                            ntripSocket.close()
                            return 0
                        elif casterResponse.startswith("HTTP/1.1 503 Service Unavailable".encode("utf-8")):
                            raise Exception("Ephemeris.Get_Ephemeris_From_Ntrip():: Invalid Server Response")
                        elif casterResponse.startswith("HTTP/1.1 200 OK".encode("utf-8")) or casterResponse.startswith("ICY 200 OK".encode("utf-8")):
                            foundNtripHeader = True
                            if self.DEBUG :
                                print ( 'foundNtripHeader = True: total %d sent and %d bytes received: %s\n' % ( totalSent, len(casterResponse), casterResponse ) )

                except socket.timeout:
                    if self.DEBUG :
                        print ( 'ntripSocket timeout %d of %d : retErrNo = %d; %d bytes sent and %d bytes received\n' % ( ntripConnectRetry, maxConnectRetry, retErrNo, totalSent, len(casterResponse) ) )
                    ntripConnectRetry = ntripConnectRetry + 1
                    ntripSocket.close()
                    ntripSocket = None
                    continue
                except socket.error as msg:
                    if self.DEBUG :
                        print ( msg + ': ntripSocket.send() sent %d bytes, %d in total of %d; retry %d\n' % ( sent, totalSent, msgLen, ntripConnectRetry ) )
                    ntripConnectRetry = ntripConnectRetry + 1
                    ntripSocket.close()
                    ntripSocket = None
                    continue
            
            remainNtripRecvBufView = ntripRecvBufView[0:]
            remainNtripRecvBufLen = ntripRecvBufLen
            while remainNtripRecvBufLen :
                numBytesRecv = ntripSocket.recv_into( remainNtripRecvBufView, remainNtripRecvBufLen )
                remainNtripRecvBufView = remainNtripRecvBufView[ numBytesRecv: ]
                remainNtripRecvBufLen -= numBytesRecv
                # if self.DEBUG:
                #     print( "total %d, received %d, remain %d, %X %X %X %X " % ( ntripRecvBufLen, numBytesRecv, remainNtripRecvBufLen, ntripRecvBuf[(ntripRecvBufLen-remainNtripRecvBufLen)], ntripRecvBuf[(ntripRecvBufLen-remainNtripRecvBufLen)+1], ntripRecvBuf[(ntripRecvBufLen-remainNtripRecvBufLen)+2], ntripRecvBuf[(ntripRecvBufLen-remainNtripRecvBufLen)+3] ) ) 
                # ret = ntripSocket.recv( remainNtripRecvBufLen )
                # numBytesRecv = len( ret )
                # if self.DEBUG:
                #     print( "total %d, received %d, remain %d" % ( ntripRecvBufLen, numBytesRecv, remainNtripRecvBufLen ) ) 
                #     print( ret ) 

            ii = 0
            while ii < ntripRecvBufLen - 70 :
                if int( ntripRecvBuf[ii] ) == RTCM3_PREAMBLE:
                    msgLen = ( int( ntripRecvBuf[ii+1] ) & 0x03 ) << 8
                    msgLen = msgLen | ( int( ntripRecvBuf[ii+2] ) & 0xFF )
                    frmEnd = ii + 3 + msgLen + 3
                    #if msgLen != RTCM_MSG1019_LENGTH or frmEnd >= ntripRecvBufLen :
                    if frmEnd >= ntripRecvBufLen :
                        ii += 1
                        continue
                    crc, msgCRC = int( 0 ), int( ntripRecvBuf[ii+3+msgLen] )
                    msgCRC = ( ( ( msgCRC << 8 ) | ntripRecvBuf[ii+3+msgLen+1] ) << 8 ) | ntripRecvBuf[ii+3+msgLen+2]
                    for msgByte in ntripRecvBuf[ ii : (ii+msgLen+3) ] :
                        crc = CRC24[ ( ( crc >> 16 )^msgByte ) & 0xFF ] ^ ( crc << 8 )
                    crc &= 0xFFFFFF

                    if crc != msgCRC:
                        msgType = int( ntripRecvBuf[ii+3] ) << 4
                        msgType = msgType | ( ( ntripRecvBuf[ii+4] >> 4 ) & 0x0F ) 
                        if self.DEBUG :
                            print ('# ii=%d Msg Type = %d, Message Length: %d, computed CRC: %d != CRC: %d ' % ( ii, msgType, msgLen, crc, msgCRC ) )
                        ii += 1
                        continue

                    msgType = int( ntripRecvBuf[ii+3] ) << 4
                    msgType = msgType | ( ( ntripRecvBuf[ii+4] >> 4 ) & 0x0F )
                    if msgType != 1019 :
                        ii += 1
                        if self.DEBUG :
                            print( 'ii=%d, Msg Length: %d, computed CRC: %d == CRC: %d, Msg Type: %d' % ( ii, msgLen, crc, msgCRC, msgType ) )
                        continue
                    else:
                        print( "ii = %04d , Msg Type = %d" % ( ii, msgType ) ) 

                    msg = ntripRecvBuf[(ii+3):(ii+msgLen+3)]

                    svid = ( ( msg[1] & 0x0F ) << 2 ) | ( msg[2] >> 6 ) #msg.read(6).uint
                    if svid > 32 or svid < 1  :
                        raise Exception( "Ephemeris.Get_Ephemeris_From_Ntrip():: The SVID number %02d is incorrect." % svid )
                    jj = svid - 1
                    self.svid[jj]       = svid
                    self.week[jj]       = ( ( ( msg[2] & 0x3F ) << 4 ) | ( msg[3] >> 4 ) )
                    self.acc[jj]        = msg[3] >> 4 #msg.read(4).uint
                    self.l2code[jj]     = msg[4] >> 6 #msg.read(2).uint
                    self.idot[jj]       = ( ( ( msg[4] & 0x3F ) << 2 ) | msg[5] ) * PWR2_43_PI #idot = msg.read(14).int
                    self.iode[jj]       = msg[6] #msg.read(8).uint
                    self.toc[jj]        = ( ( msg[7] << 8 ) | msg[8] ) << 4  #toc = msg.read(16).uint
                    self.af2[jj]        = msg[9] * PWR2_55 #af2 = msg.read(8).int
                    self.af1[jj]        = ( ( msg[10] << 8 ) | msg[11] ) * PWR2_43 #af1 = msg.read(16).int
                    self.af0[jj]        = ( ( ( ( msg[12] << 8 ) | msg[13] ) << 6 ) | ( msg[14] >> 2 ) ) * PWR2_31 #af0 = msg.read(22).int
                    self.iodc[jj]       = ( ( msg[14] & 0x03 )  << 8 ) | msg[15] #msg.read(10).uint
                    self.Crs[jj]        = ( ( msg[16] << 8 ) | msg[17] ) * PWR2_5 #crs = msg.read(16).int
                    self.deltan[jj]     = ( ( msg[18] << 8 ) | msg[19] ) * PWR2_43_PI #deltan  = ( msg[18] << 8 ) | msg[19] = msg.read(16).int
                    self.M0[jj]         = ( ( ( ( ( ( msg[20] << 8 ) | msg[21] ) << 8 ) | msg[22] ) << 8 ) | msg[23] ) * PWR2_31_PI #msg.read(32).int
                    self.Cuc[jj]        = ( ( msg[24] << 8 ) | msg[25] ) * PWR2_29 #msg.read(16).int
                    #########################################################################################################
                    # The eccentricity of the Earth's orbit is currently about 0.0167; the Earth's orbit is nearly circular. 
                    # Venus and Neptune have even lower eccentricity. Over hundreds of thousands of years, the eccentricity 
                    # of the Earth's orbit varies from nearly 0.0034 to almost 0.058 as a result of gravitational attractions 
                    # among the planets
                    ###########################################################################################################
                    self.e[jj]          = ( ( ( ( ( ( msg[26] << 8 ) | msg[27] ) << 8 ) | msg[28] ) << 8 ) | msg[29] ) * PWR2_33 #msg.read(32).uint
                    self.Cus[jj]        = ( ( msg[30] << 8 ) | msg[31] ) * PWR2_29 #msg.read(16).int
                    self.sqrtA[jj]      = ( ( ( ( ( ( msg[32] << 8 ) | msg[33] ) << 8 ) | msg[34] ) << 8 ) | msg[35] ) * PWR2_19 #msg.read(32).uint
                    t_toe               = ( ( msg[36] << 8 ) | msg[37] ) << 4 #msg.read(16).uint
                    if t_toe < self.toe[jj] :
                        print ( "Ephemeris.Get_Ephemeris_From_Ntrip():: The new Toe %d for Sat %02d is older than the current Toe %d." % (t_toe, self.svid[jj], self.toe[jj] ) )
                    if self.DEBUG:
                        print ( "self.svid[jj] = %d, self.week[jj] = %d, t_toe = %d "%(self.svid[jj], self.week[jj] , t_toe) )
                    self.toe[jj]        = t_toe                    
                    self.Cic[jj]        = ( ( msg[38] << 8 ) | msg[39] ) * PWR2_29 #msg.read(16).int
                    self.Omega0[jj]     = ( ( ( ( ( ( msg[40] << 8 ) | msg[41] ) << 8 ) | msg[42] ) << 8 ) | msg[43] ) * PWR2_31_PI #msg.read(32).int
                    self.Cis[jj]        = ( ( msg[44] << 8 ) | msg[45] ) * PWR2_29 #msg.read(16).int
                    self.i0[jj]         = ( ( ( ( ( ( msg[46] << 8 ) | msg[47] ) << 8 ) | msg[48] ) << 8 ) | msg[49] ) * PWR2_31_PI #msg.read(32).int
                    self.Crc[jj]        = ( ( msg[50] << 8 ) | msg[51] ) * PWR2_5 #msg.read(16).int
                    self.omega[jj]      = ( ( ( ( ( ( msg[52] << 8 ) | msg[53] ) << 8 ) | msg[54] ) << 8 ) | msg[55] ) * PWR2_31_PI #msg.read(32).int
                    self.Omegadot[jj]   = ( ( ( ( msg[56] << 8 ) | msg[57] ) << 8 ) | msg[58] ) * PWR2_43_PI #msg.read(24).int
                    self.tgd[jj]        = msg[59] * PWR2_31 #msg.read(8).int
                    self.health[jj]     = msg[60] >> 2 #msg.read(6).uint
                    self.l2p[jj]        = ( msg[60] & 0x02 ) >> 1 #msg.read(1).uint
                    self.fit[jj]        = msg[60] & 0x01 #msg.read(1).uint
                    ii = frmEnd
                else :
                    ii += 1
    
            if self.DEBUG :
                now = time.time()
                print ( 'Ephemeris.Get_Ephemeris_From_Ntrip():: %6.3f seconds are left.' % ( Tstop - now ) )

            ntripSocket.close()

        if ephemerisFilename is None :
            ephemerisFilename = self.NtripMountpoint + time.strftime( "%j%H.%yO.eph", time.gmtime() )

        self.ephemerisFilename = ephemerisFilename + '.txt'
        ephemerisFObj = open( localDir + self.ephemerisFilename, "w" )
        ii = 0
        for svid in self.svid :
            if svid == 0 :
                continue
            jj = svid - 1
            ephemerisFObj.write( 'svid: %d\n' % self.svid[jj] )
            ephemerisFObj.write( 'week: %d\n' % self.week[jj] )
            ephemerisFObj.write( 'acc: %d\n' % self.acc[jj] )
            ephemerisFObj.write( 'l2code: %d\n' % self.l2code[jj] )
            ephemerisFObj.write( 'idot: %E\n' % self.idot[jj] )
            ephemerisFObj.write( 'iode: %d\n' % self.iode[jj] )
            ephemerisFObj.write( 'toc: %d\n' % self.toc[jj] )
            ephemerisFObj.write( 'af2: %E\n' % self.af2[jj] )
            ephemerisFObj.write( 'af1: %E\n' % self.af1[jj] )
            ephemerisFObj.write( 'af0: %E\n' % self.af0[jj] )
            ephemerisFObj.write( 'iodc: %d\n' % self.iodc[jj] )
            ephemerisFObj.write( 'Crs: %E\n' % self.Crs[jj] )
            ephemerisFObj.write( 'deltan: %E\n' % self.deltan[jj] )
            ephemerisFObj.write( 'M0: %E\n' % self.M0[jj] )
            ephemerisFObj.write( 'Cuc: %E\n' % self.Cuc[jj] )
            ephemerisFObj.write( 'e: %E\n' % self.e[jj] )
            ephemerisFObj.write( 'Cus: %E\n' % self.Cus[jj] )
            ephemerisFObj.write( 'sqrtA: %E\n' % self.sqrtA[jj] )
            ephemerisFObj.write( 'toe: %d\n' % self.toe[jj] )
            ephemerisFObj.write( 'Cic: %E\n' % self.Cic[jj] )
            ephemerisFObj.write( 'Omega0: %E\n' % self.Omega0[jj] )
            ephemerisFObj.write( 'Cis: %E\n' % self.Cis[jj] )
            ephemerisFObj.write( 'i0: %E\n' % self.i0[jj] )
            ephemerisFObj.write( 'Crc: %E\n' % self.Crc[jj] )
            ephemerisFObj.write( 'omega: %E\n' % self.omega[jj] )
            ephemerisFObj.write( 'Omegadot: %E\n' % self.Omegadot[jj] )
            ephemerisFObj.write( 'tgd: %E\n' % self.tgd[jj] )
            ephemerisFObj.write( 'health: %d\n' % self.health[jj] )
            ephemerisFObj.write( 'l2p: %d\n' % self.l2p[jj] )
            ephemerisFObj.write( 'fit: %d\n\r\n\r\n' % self.fit[jj] )
            ii += 1
        ephemerisFObj.close()

        self.ephemerisFilename = ephemerisFilename + '.json'
        ephemeris_ = []
        for svid in self.svid :
            if svid == 0 :
                continue
            jj = svid - 1
            ephemeris_.append({
                'svid':     self.svid[jj],
                'week':     self.week[jj],
                'acc':      self.acc[jj],
                'l2code':   self.l2code[jj],
                'idot':     self.idot[jj],
                'iode':     self.iode[jj],
                'toc':      self.toc[jj],
                'af2':      self.af2[jj],
                'af1':      self.af1[jj],
                'af0':      self.af0[jj],
                'iodc':     self.iodc[jj],
                'Crs':      self.Crs[jj],
                'deltan':   self.deltan[jj],
                'M0':       self.M0[jj],
                'Cuc':      self.Cuc[jj],
                'e':        self.e[jj],
                'Cus':      self.Cus[jj],
                'sqrtA':    self.sqrtA[jj],
                'toe':      self.toe[jj],
                'Cic':      self.Cic[jj],
                'Omega0':   self.Omega0[jj],
                'Cis':      self.Cis[jj],
                'i0':       self.i0[jj],
                'Crc':      self.Crc[jj],
                'omega':    self.omega[jj],
                'Omegadot': self.Omegadot[jj],
                'tgd':      self.tgd[jj],
                'health':   self.health[jj],
                'l2p':      self.l2p[jj],
                'fit':      self.fit[jj]
            })
        with open( localDir + self.ephemerisFilename, 'w') as f:
            json.dump( ephemeris_, f )
        
        return ii, ephemerisFilename

    def Restore_Ephemeris( self, subdirectory = "IGS" , ephemerisFilename = None ) :
        localDir = self.localDir + "stations/" + subdirectory + '/'
        if ephemerisFilename :
            if os.path.isfile( localDir + ephemerisFilename ) :
                self.ephemerisFilename = ephemerisFilename
            else:
                raise Exception( "Ephemeris.Get_Ephemeris_From_Backup():: The input file %s doesn't exist." % ephemerisFilename )
        else:
            self.ephemerisFilename = 'hour' + time.strftime( "%j%H.%yO.eph", time.gmtime() )

        if os.path.isfile( localDir + self.ephemerisFilename ) :
            ephemerisFObj = open( localDir + self.ephemerisFilename, "r" )
            restored = list() 
            for line in ephemerisFObj :
                param = line.split()
                if len( param ) < 2 :
                    continue
                if 'svid:' == param[0] :
                    jj = int( param[1] ) - 1
                    self.svid[ jj ] = jj + 1
                elif 'week:' == param[0] :
                    self.week[jj] = int( param[1] )
                elif 'acc:' == param[0] :
                    self.acc[jj] = int( param[1] )
                elif 'l2code:' == param[0] :
                    self.l2code[jj] = int( param[1] )
                elif 'idot:' == param[0] :
                    self.idot[jj] = float( param[1] )
                elif 'iode:' == param[0] :
                    self.iode[jj] = int( param[1] )
                elif 'toc:' == param[0] :
                    self.toc[jj] = int( param[1] )
                elif 'af2:' == param[0] :
                    self.af2[jj] = float( param[1] )
                elif 'af1:' == param[0] :
                    self.af1[jj] = float( param[1] )
                elif 'af0:' == param[0] :
                    self.af0[jj] = float( param[1] )
                elif 'iodc:' == param[0] :
                    self.iodc[jj] = int( param[1] )
                elif 'Crs:' == param[0] :
                    self.Crs[jj] = float( param[1] )
                elif 'deltan:' == param[0] :
                    self.deltan[jj] = float( param[1] )
                elif 'M0:' == param[0] :
                    self.M0[jj] = float( param[1] )
                elif 'Cuc:' == param[0] :
                    self.Cuc[jj] = float( param[1] )
                elif 'e:' == param[0] :
                    #########################################################################################################
                    # The eccentricity of the Earth's orbit is currently about 0.0167; the Earth's orbit is nearly circular. 
                    # Venus and Neptune have even lower eccentricity. Over hundreds of thousands of years, the eccentricity 
                    # of the Earth's orbit varies from nearly 0.0034 to almost 0.058 as a result of gravitational attractions 
                    # among the planets
                    ###########################################################################################################
                    self.e[jj] = float( param[1] )
                elif 'Cus:' == param[0] :
                    self.Cus[jj] = float( param[1] )
                elif 'sqrtA:' == param[0] :
                    self.sqrtA[jj] = float( param[1] )
                elif 'toe:' == param[0] :
                    t_toe = int( param[1] )
                    if t_toe < self.toe[jj] :
                        print ( "Ephemeris.Get_Ephemeris_From_Ntrip():: The new Toe %d for Sat %02d is older than the current Toe %d." % (t_toe, self.svid[jj], self.toe[jj] ) )
                    self.toe[jj] = t_toe
                elif 'Cic:' == param[0] :
                    self.Cic[jj] = float( param[1] )
                elif 'Omega0:' == param[0] :
                    self.Omega0[jj] = float( param[1] )
                elif 'Cis:' == param[0] :
                    self.Cis[jj] = float( param[1] )
                elif 'i0:' == param[0] :
                    self.i0[jj] = float( param[1] )
                elif 'Crc:' == param[0] :
                    self.Crc[jj] = float( param[1] )
                elif 'omega:' == param[0] :
                    self.omega[jj] = float( param[1] )
                elif 'Omegadot:' == param[0] :
                    self.Omegadot[jj] = float( param[1] )
                elif 'tgd:' == param[0] :
                    self.tgd[jj] = float( param[1] )
                elif 'health:' == param[0] :
                    self.health[jj] = int( param[1] )
                elif 'l2p:' == param[0] :
                    self.l2p[jj] = int( param[1] )
                elif 'fit:' == param[0] :
                    self.fit[jj] = int( param[1] )
            ephemerisFObj.close()
            self.NtripSessionTime /= 2
        else:
            self.ephemerisFilename = None
        
        return self.ephemerisFilename

    def Load_EphemerisJsonFile( self, ephemerisFilename ):
        with open( ephemerisFilename, 'r' ) as f:
            ephData = json.load( f )
            for eph in ephData:
                ii = eph['svid'] - 1
                self.svid[ii]       = ii + 1
                self.week[ii]       = eph['week']
                self.idot[ii]       = eph['idot']
                self.af2[ii]        = eph['af2']
                self.af1[ii]        = eph['af1']
                self.Crs[ii]        = eph['Crs']
                self.deltan[ii]     = eph['deltan']
                self.M0[ii]         = eph['M0']
                self.Cuc[ii]        = eph['Cuc']
                self.e[ii]          = eph['e']
                self.Cus[ii]        = eph['Cus']
                self.sqrtA[ii]      = eph['sqrtA']
                self.toe[ii]        = eph['toe']
                self.Cic[ii]        = eph['Cic']
                self.Omega0[ii]     = eph['Omega0']
                self.Cis[ii]        = eph['Cis']
                self.i0[ii]         = eph['i0']
                self.Crc[ii]        = eph['Crc']
                self.omega[ii]      = eph['omega']
                self.Omegadot[ii]   = eph['Omegadot']
                self.tgd[ii]        = eph['tgd']

                if 'YY' in eph:
                    self.acc[ii]    = eph['YY']
                if 'MM' in eph:
                    self.acc[ii]    = eph['MM']
                if 'DD' in eph:
                    self.acc[ii]    = eph['DD']
                if 'hh' in eph:
                    self.acc[ii]    = eph['hh']
                if 'mm' in eph:
                    self.acc[ii]    = eph['mm']
                if 'sec' in eph:
                    self.acc[ii]    = eph['sec']

                if 'acc' in eph:
                    self.acc[ii]    = eph['acc']
                if 'l2code' in eph:
                    self.l2code[ii] = eph['l2code']
                if 'iode' in eph:
                    self.iode[ii]   = eph['iode']
                if 'iodc' in eph:
                    self.iodc[ii]   = eph['iodc']
                if 'toc' in eph:
                    self.toc[ii]    = eph['toc']
                if 'health' in eph:
                    self.health[ii] = eph['health']
                if 'l2p' in eph:
                    self.l2p[ii]    = eph['l2p']
                if 'fit' in eph:
                    self.fit[ii]    = eph['fit']

        return

    def Read_Rinex(self, rinexFilename ) :        
        localDir = self.localDir
        if not os.path.isfile( localDir + rinexFilename ) :
            raise Exception( "Ephemeris.Read_Rinex():: Ths input file %s doesn't exist." % (localDir + rinexFilename) )
        
        with open( localDir + rinexFilename, 'rt' ) as rinexFile :
            while True :
                line = rinexFile.readline()
                if line[60:73] == "END OF HEADER" :
                    break
                if line[60:69] == "ION ALPHA" :
                    self.alpha0 = float( line[2:14].replace('D', 'E') )
                    self.alpha1 = float( line[14:26].replace('D', 'E') )
                    self.alpha2 = float( line[26:38].replace('D', 'E') )
                    self.alpha3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "ION BETA" :
                    self.beta0 = float( line[2:14].replace('D', 'E') )
                    self.beta1 = float( line[14:26].replace('D', 'E') )
                    self.beta2 = float( line[26:38].replace('D', 'E') )
                    self.beta3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "DELTA-UTC" :
                    self.A0 = float( line[3:22].replace('D', 'E') )
                    self.A1 = float( line[22:41].replace('D', 'E') )
                    self.tot = int( line[41:50] )
                    self.wnt = int( line[50:59] )
                elif line[60:72] == "LEAP SECONDS" :
                    self.dtls = int( line[0:6].replace('D', 'E') )

            line = rinexFile.readline( )
            while line:
                ii = int( line[0:2] ) - 1
                self.svid[ii] = ii + 1
                self.YY[ii] = int( line[3:5] ) + 2000
                self.MM[ii] = int( line[6:8] )
                self.DD[ii] = int( line[9:11] )
                self.hh[ii] = int( line[12:14] )
                self.mm[ii] = int( line[15:17] )
                self.sec[ii] = float( line[18:22] )
                self.af0[ii] = float( line[22:41].replace('D', 'E')  )
                self.af1[ii] = float( line[41:60].replace('D', 'E')  )
                self.af2[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFile.readline( )
                self.iode[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Crs[ii] = float( line[22:41].replace('D', 'E')  )
                self.deltan[ii] = float( line[41:60].replace('D', 'E')  )
                self.M0[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFile.readline( )
                self.Cuc[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.e[ii] = float( line[22:41].replace('D', 'E')  )
                self.Cus[ii] = float( line[41:60].replace('D', 'E')  )
                self.sqrtA[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFile.readline( )
                self.toe[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Cic[ii] = float( line[22:41].replace('D', 'E')  )
                self.Omega0[ii] = float( line[41:60].replace('D', 'E')  )
                self.Cis[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFile.readline( )
                self.i0[ii] = float( line[3:22].replace('D', 'E')  )
                self.Crc[ii] = float( line[22:41].replace('D', 'E')  )
                self.omega[ii] = float( line[41:60].replace('D', 'E')  )
                self.Omegadot[ii] = float( line[60:79].replace('D', 'E')  )

                line = rinexFile.readline( )
                self.idot[ii] = float( line[3:22].replace('D', 'E')  )
                self.l2code[ii] = int( float( line[22:41].replace('D', 'E')  ) )  #### Is self.l2code  an integer?
                self.week[ii] = int( float( line[41:60].replace('D', 'E')  ) )  #### Is self.week  an integer?

                line = rinexFile.readline( )
                self.health[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.health[ii] = ( self.health[ii] % 32 + 32 ) if ( self.health[ii] < 32 ) else self.health[ii] 
                self.tgd[ii] = float( line[41:60].replace('D', 'E')  ) 
                self.iodc[ii] = int( float( line[60:79].replace('D', 'E')  ) )

                line = rinexFile.readline( )
                line = rinexFile.readline( )

        rinexFile.close()   
        return        
    
    
    def Set_WeekNumber( self, weekNum, svid = None ) :
        if svid is None :
            svid = list( self.svid )

        for ii in svid:
            self.week[ii-1] = weekNum
    
    def Get_Almanac( self, svid = None ):
        if svid is None :
            svid = list()
            for ii in range( MAX_GPS_SATS ):
                if self.svid[ii] == 0 :
                    continue
                svid.append( self.svid[ii] )
        NUM_SAT = len( svid )

        alm  = [ list([ float(0.0)  for col in range(13) ]) for row in range( NUM_SAT ) ]
        for ii in range( NUM_SAT ) :
            sat = svid[ii] - 1
            alm[ii][0] = self.svid[sat]
            alm[ii][1] = self.health[sat]
            alm[ii][2] = self.e[sat]
            alm[ii][3] = self.toe[sat]
            alm[ii][4] = self.i0[sat]
            alm[ii][5] = self.Omegadot[sat]
            alm[ii][6] = self.sqrtA[sat]
            alm[ii][7] = self.Omega0[sat]
            alm[ii][8] = self.omega[sat]
            alm[ii][9] = self.M0[sat]
            alm[ii][10] = self.af0[sat]
            alm[ii][11] = self.af1[sat]
            alm[ii][12] = self.week[sat]

        return alm
    
    def Get_SatState( self, epoch = None, svid = None  ) :
        #
        # http://www.navipedia.net/index.php/GPS_and_Galileo_Satellite_Coordinates_Computation
        #
        if ( epoch is None ) or ( epoch == 0 ) :
            now = time.gmtime()
            epoch_now = calendar.timegm( now ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00', '%Y-%m-%d %H:%M:%S' ) ) + LEAPSECONDS_2017
            epoch = [ epoch_now // WEEK_SECONDS, epoch_now % WEEK_SECONDS ]

        if svid is None :
            svid = list()
            for ii in range( MAX_GPS_SATS ):
                if self.svid[ii] == 0 :
                    continue
                svid.append( self.svid[ii] )
        
        NUM_SAT = len( svid )

        lari        = [ float(0.0) for col in range(6) ]
        xyOrb       = [ float(0.0) for col in range(4) ]
        satState    = [ list([ float(0.0) for col in range( 10 ) ]) for row in range( NUM_SAT ) ]

        for ii in range( NUM_SAT ) :
            if svid[ii] == 0 :  # The Sat or its ephemeris is unavailable.
                continue
            sat = svid[ii] - 1
            if self.svid[sat] == 0 :
                continue
            satState[ii][0] = svid[ii]
            tk =  epoch[1] - self.toe[sat] # + ( epoch[0]-self.week[sat] ) * WEEK_SECONDS
            if tk > 302400 :
                tk -= WEEK_SECONDS
            if tk < -302400 :
                tk += WEEK_SECONDS
            
            #############################################################################
            # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
            # The time of transmissionis is modified from TOW from navigation data and
            # the correction of is performed as
            #############################################################################
            A = self.sqrtA[sat]*self.sqrtA[sat]
            meanMotion = sqrt(WGS84_u_m3s2) / (self.sqrtA[sat]*A) + self.deltan[sat]# This value should be caliberated with adding Mean Motion Difference.
            meanAnomaly  = self.M0[sat] + meanMotion * tk

                #######################################################################################################################
                # Kepler equation for calculating Eccentric Anomaly
                #
                # In orbital mechanics, Kepler's equation relates various geometric properties of the orbit
                # of a body subject to a central force.
                #
                # ecceAnomaly = Kepler( meanAnomaly, Eccentricity[sat] )
                # := ecceAnomaly - (ecceAnomaly-Eccentricity[sat]*sin(ecceAnomaly)-meanAnomaly)/(1-Eccentricity[sat]*cos(ecceAnomaly))
                #
                # Reference: https://en.wikipedia.org/wiki/Kepler's_equation
                ########################################################################################################################
            E0 = meanAnomaly
            ecceAnomaly = meanAnomaly + self.e[sat]*sin(meanAnomaly) / ( 1 - self.e[sat]*cos(meanAnomaly) )
            while( abs(ecceAnomaly-E0) > 1.0e-14 ):
                E0 = ecceAnomaly
                ecceAnomaly = ecceAnomaly - ( ecceAnomaly - self.e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - self.e[sat]*cos(ecceAnomaly) )
            sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
            e_cosEcceAnom_1 = 1 -  self.e[sat]*cosEcceAnom

            lari[0] = atan2( sqrt(1-self.e[sat]*self.e[sat]) * sinEcceAnom,  cosEcceAnom-self.e[sat] ) #True anomaly. Later, Argument of latitude
            lari[3] = sinEcceAnom * meanMotion * (1+self.e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

            lari[0] += self.omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            lari[0:3] = [ lari[0] + self.Cus[sat]*sin2Lat + self.Cuc[sat]*cos2Lat ,
                A*e_cosEcceAnom_1 + self.Crs[sat]*sin2Lat + self.Crc[sat]*cos2Lat ,
                     self.i0[sat] + self.Cis[sat]*sin2Lat + self.Cic[sat]*cos2Lat + self.idot[sat]*tk ]

            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            ddLat, ddRad, ddInc = 2 * ( self.Cus[sat] * cos2Lat - self.Cuc[sat] * sin2Lat ) * lari[3], 2 * ( self.Crs[sat] * cos2Lat - self.Crc[sat] * sin2Lat ) * lari[3], 2 * ( self.Cis[sat] * cos2Lat - self.Cic[sat] * sin2Lat ) * lari[3]

            lari[3:6] = [ lari[3] + ddLat, A*self.e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, self.idot[sat] + ddInc ]

            cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
            satState[ii][8:10] = [ self.af0[sat] + F*self.e[sat]*self.sqrtA[sat]*sinEcceAnom - self.tgd[sat], self.af1[sat] + self.af2[sat]*tk ] 
            satState[ii][7] = satState[ii][8] + satState[ii][9] * tk 

            dotAsceLon = self.Omegadot[sat] - WGS84_w_rads
            asceLon = self.Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * self.toe[sat] # satState[ii][7] ## satState[ii][7] is the clock bias
            cosLon, sinLon = cos( asceLon ), sin( asceLon )

            xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
            xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

            satState[ii][1:4] = [ +cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1] ,
                                  +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                                  +sinInc*xyOrb[1] ]
            satState[ii][4:7] = [ +cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satState[ii][2] ,
                                  +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satState[ii][1] ,
                                  +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

        return satState

    def Get_PositionAssist( self, rcvXyz, epoch = None , propagationDelayCompensation = False  ) :
        #
        # The assistance information timed when the receiver received the GPS signals.
        #
        #############################################################################################
        # Convert the Lat-Lon-Height of the GPS signal receiver to X-Y-Z
        #
        # %Reference Position, the truth
        # rcvLalh = [ +(32+56.4806/60) * (ICD200PI/180), -( 117+14.4715/60 ) * (ICD200PI/180), 50 ] 
        ############################################################################################# 

        if epoch is None :
            now = time.gmtime()
            epoch_now = calendar.timegm( now ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00', '%Y-%m-%d %H:%M:%S' ) ) + LEAPSECONDS_2017
            epoch = [ epoch_now // WEEK_SECONDS, epoch_now % WEEK_SECONDS ]

        svid = list()
        for ii in range( MAX_GPS_SATS ):
            if self.svid[ii] == 0 :
                continue
            svid.append( self.svid[ii] )

        NUM_SAT = len( svid )
        visibleSat  = list()
        xyzdxyz     = list()
        lari        = [ float(0.0) for col in range(6) ]
        xyOrb       = [ float(0.0) for col in range(4) ]
        satState    = [ list([ float(0.0) for col in range(10) ]) for row in range( MAX_GPS_SATS ) ]
        dxyz        = [ list([ float(0.0) for col in range(3) ]) for row in range( MAX_GPS_SATS ) ]

        for ii in range( NUM_SAT ) :
            if svid[ii] == 0 :
                continue
            sat = svid[ii] - 1
            satState[sat][0] = svid[ii]
            
            tk = epoch[1] - self.toe[sat] # + ( epoch[0]-self.week[sat] ) * WEEK_SECONDS
            if tk > 302400 :
                tk -= WEEK_SECONDS
            if tk < -302400 :
                tk += WEEK_SECONDS
            #############################################################################
            # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
            # The time of transmissionis is modified from TOW from navigation data and
            # the correction of is performed as
            #############################################################################
            A = self.sqrtA[sat]*self.sqrtA[sat]
            meanMotion = sqrt(WGS84_u_m3s2) / (self.sqrtA[sat]*A) + self.deltan[sat]# This value should be caliberated with adding Mean Motion Difference.
            meanAnomaly  = self.M0[sat] + meanMotion * tk

                #######################################################################################################################
                # Kepler equation for calculating Eccentric Anomaly
                #
                # In orbital mechanics, Kepler's equation relates various geometric properties of the orbit
                # of a body subject to a central force.
                #
                # ecceAnomaly = Kepler( meanAnomaly, Eccentricity[sat] )
                # := ecceAnomaly - (ecceAnomaly-Eccentricity[sat]*sin(ecceAnomaly)-meanAnomaly)/(1-Eccentricity[sat]*cos(ecceAnomaly))
                #
                # Reference: https://en.wikipedia.org/wiki/Kepler's_equation
                ########################################################################################################################
            E0 = meanAnomaly
            ecceAnomaly = meanAnomaly + self.e[sat]*sin(meanAnomaly) / ( 1 - self.e[sat]*cos(meanAnomaly) )
            while( abs(ecceAnomaly-E0) > 1.0e-14 ):
                E0 = ecceAnomaly
                ecceAnomaly = ecceAnomaly - ( ecceAnomaly - self.e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - self.e[sat]*cos(ecceAnomaly) )
            sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
            e_cosEcceAnom_1 = 1 -  self.e[sat]*cosEcceAnom

            lari[0] = atan2( sqrt(1-self.e[sat]*self.e[sat]) * sinEcceAnom,  cosEcceAnom-self.e[sat] ) #True anomaly. Later, Argument of latitude
            lari[3] = sinEcceAnom * meanMotion * (1+self.e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

            lari[0] += self.omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            lari[0:3] = [ lari[0] + self.Cus[sat]*sin2Lat + self.Cuc[sat]*cos2Lat ,
                A*e_cosEcceAnom_1 + self.Crs[sat]*sin2Lat + self.Crc[sat]*cos2Lat ,
                     self.i0[sat] + self.Cis[sat]*sin2Lat + self.Cic[sat]*cos2Lat + self.idot[sat]*tk ]

            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            ddLat, ddRad, ddInc = 2 * ( self.Cus[sat] * cos2Lat - self.Cuc[sat] * sin2Lat ) * lari[3], 2 * ( self.Crs[sat] * cos2Lat - self.Crc[sat] * sin2Lat ) * lari[3], 2 * ( self.Cis[sat] * cos2Lat - self.Cic[sat] * sin2Lat ) * lari[3]

            lari[3:6] = [ lari[3] + ddLat, A*self.e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, self.idot[sat] + ddInc ]

            cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
            
            dotAsceLon = self.Omegadot[sat] - WGS84_w_rads
            asceLon = self.Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * self.toe[sat] # satState[sat][7] ## satState[sat][7] is the clock bias 
            cosLon, sinLon = cos( asceLon ), sin( asceLon )

            xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
            xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

            satState[sat][1:4] = [  +cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1] ,
                                    +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                                    +sinInc*xyOrb[1] ]
            satState[sat][4:7] = [  +cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satState[sat][2] ,
                                    +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satState[sat][1] ,
                                    +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

            satState[sat][8:10] = [ self.af0[sat] + F*self.e[sat]*self.sqrtA[sat]*sinEcceAnom - self.tgd[sat], self.af1[sat] + self.af2[sat]*tk ] 
            satState[sat][7] = satState[sat][8] + satState[sat][9] * tk

            dxyz[sat] = [ satState[sat][1] - rcvXyz[0], satState[sat][2] - rcvXyz[1], satState[sat][3] - rcvXyz[2]  ]

            if propagationDelayCompensation : # propagationDelayCompensation = True
                for ii in range( 1 ) :
                    tau = sqrt( dxyz[sat][0]*dxyz[sat][0] + dxyz[sat][1]*dxyz[sat][1] + dxyz[sat][2]*dxyz[sat][2] ) / C0  #The light propogation time between the receiver and each satellite.                                
                    dSatState = [ satState[sat][4]*tau, satState[sat][5]*tau, satState[sat][6]*tau ] #The satellite state Difference during the light propogaion 

                    satState[sat][1:4] = [ satState[sat][1] - dSatState[0], satState[sat][2] - dSatState[1], satState[sat][3] - dSatState[2] ] #extrapolate the satellite position backwards to the transmission time.
                    ## The change in satellite velocity may be neglected.
                    satState[sat][7] -= satState[sat][9] * tau #Caliberate the GPS clock bias with subtrating the extra clock drifting.
                
                    dxyz[sat] = [ dxyz[sat][0] - dSatState[0], dxyz[sat][1] - dSatState[1], dxyz[sat][2] - dSatState[2] ] #recalculate the spatial difference between the receiver and each satellite.
                
            xyzdxyz1 =  dxyz[sat][0]*rcvXyz[0] + dxyz[sat][1]*rcvXyz[1] + dxyz[sat][2]*rcvXyz[2]
            if xyzdxyz1 > 0 :
                visibleSat.append( svid[ii] )
                xyzdxyz.append( xyzdxyz1 ) 
        
        numVisibleSat   = len( visibleSat )
        satXyzvt        = [ list([ float(0.0) for col in range(8) ]) for row in range( numVisibleSat ) ]
        assist          = [ list([ float(0.0) for col in range(4) ]) for row in range( numVisibleSat ) ]
        for ii in range( numVisibleSat ) :
            jj = visibleSat[ii] - 1
            satXyzvt[ii]    = satState[jj][0:8]
            pseudorange     = sqrt( dxyz[jj][0]*dxyz[jj][0] + dxyz[jj][1]*dxyz[jj][1] + dxyz[jj][2]*dxyz[jj][2] ) ## - satState[jj][7] * C0 
            pseudorangedot  = ( satState[jj][4]*dxyz[jj][0] + satState[jj][5]*dxyz[jj][1] + satState[jj][6]*dxyz[jj][2]  ) / pseudorange
            
            #PRN, L1 Doppler (Hz), C/A codephase(chips), pseudorange?
            assist[ii] = [ jj+1, ( pseudorangedot/C0 - satState[jj][9] ) * L1_FREQ, ( pseudorange % L1_CODE_M ) / L1_CHIP_M, pseudorange  ]
        
        return assist, satXyzvt
    
    #
    # 03/27/2019  Shu Wang, shuwang1@outlook.com
    #   Review Get_Receiver_XyzPosition() and make sure the calibration works and improved.
    #
    def Get_Receiver_XyzPosition( self, epoch, visibleSats, codeShifts, rcvXyz, propagationDelayCompensation = False ) :
        # WEEK_SECONDS      = 604800
        # HALF_WEEK_SECONDS = 302400

        lari            = [ float(0.0) for col in range(6) ]
        xyOrb           = [ float(0.0) for col in range(4) ]

        numVisibleSat   = len( visibleSats )
        satXyzvt        = [ list([ float(0.0) for col in range(8) ]) for row in range( numVisibleSat ) ]
        rcvXyzt         = [ rcvXyz[0], rcvXyz[1], rcvXyz[2], 0 ]

        for ii in range( numVisibleSat ) :
            sat = visibleSats[ii] - 1
            if self.svid[sat] == 0 :
                if self.DEBUG:
                    print("No ephemeris data for SV %d" % (sat+1))
                continue
            satXyzvt[ii][0] = sat
            tk = epoch[1] - self.toe[sat] #+ ( epoch[0] - self.week[sat] ) * WEEK_SECONDS
            if tk > 302400 :
                tk -= WEEK_SECONDS
            if tk < -302400 :
                tk += WEEK_SECONDS
            #############################################################################
            # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
            # The time of transmissionis is modified from TOW from navigation data and
            # the correction of is performed as
            #############################################################################
            A = self.sqrtA[sat]*self.sqrtA[sat]
            meanMotion = sqrt(WGS84_u_m3s2) / (self.sqrtA[sat]*A) + self.deltan[sat] # This value should be caliberated with adding Mean Motion Difference.
            meanAnomaly  = self.M0[sat] + meanMotion * tk

                #######################################################################################################################
                # Kepler equation for calculating Eccentric Anomaly
                #
                # In orbital mechanics, Kepler's equation relates various geometric properties of the orbit
                # of a body subject to a central force.
                #
                # ecceAnomaly = Kepler( meanAnomaly, Eccentricity[sat] )
                # := ecceAnomaly - (ecceAnomaly-Eccentricity[sat]*sin(ecceAnomaly)-meanAnomaly)/(1-Eccentricity[sat]*cos(ecceAnomaly))
                #
                # Reference: https://en.wikipedia.org/wiki/Kepler's_equation
                ########################################################################################################################
            E0 = meanAnomaly
            ecceAnomaly = meanAnomaly + self.e[sat]*sin(meanAnomaly) / ( 1 - self.e[sat]*cos(meanAnomaly) )
            watchDog = 0
            while( ( abs(ecceAnomaly-E0) > 1.0e-14 ) and ( watchDog < 100 ) ) :
                E0 = ecceAnomaly
                ecceAnomaly = ecceAnomaly - ( ecceAnomaly - self.e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - self.e[sat]*cos(ecceAnomaly) )
                watchDog += 1
            sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
            e_cosEcceAnom_1 = 1 -  self.e[sat]*cosEcceAnom

            lari[0] = atan2( sqrt(1-self.e[sat]*self.e[sat]) * sinEcceAnom,  cosEcceAnom-self.e[sat] ) #True anomaly. Later, Argument of latitude
            lari[3] = sinEcceAnom * meanMotion * (1+self.e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

            lari[0] += self.omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            lari[0:3] = [ lari[0] + self.Cus[sat]*sin2Lat + self.Cuc[sat]*cos2Lat ,
                A*e_cosEcceAnom_1 + self.Crs[sat]*sin2Lat + self.Crc[sat]*cos2Lat ,
                     self.i0[sat] + self.Cis[sat]*sin2Lat + self.Cic[sat]*cos2Lat + self.idot[sat]*tk ]

            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            ddLat, ddRad, ddInc = 2 * ( self.Cus[sat] * cos2Lat - self.Cuc[sat] * sin2Lat ) * lari[3], 2 * ( self.Crs[sat] * cos2Lat - self.Crc[sat] * sin2Lat ) * lari[3], 2 * ( self.Cis[sat] * cos2Lat - self.Cic[sat] * sin2Lat ) * lari[3]

            lari[3:6] = [ lari[3] + ddLat, A*self.e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, self.idot[sat] + ddInc ]

            cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
            satXyzvt[ii][7] =  self.af0[sat] + ( self.af1[sat] + self.af2[sat]*tk ) * tk + F*self.e[sat]*self.sqrtA[sat]*sinEcceAnom - self.tgd[sat]

            dotAsceLon = self.Omegadot[sat] - WGS84_w_rads
            asceLon = self.Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * self.toe[sat]
            cosLon, sinLon = cos( asceLon ), sin( asceLon )

            xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
            xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

            satXyzvt[ii][1:4] = [+cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1], 
                                +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                                +sinInc*xyOrb[1] ]

            satXyzvt[ii][4:7] = [+cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satXyzvt[ii][2] ,
                                +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satXyzvt[ii][1] ,
                                +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

        ########################################################################################################
        # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
        ########################################################################################################
        dXyz    = [ float(0.0), float(0.0), float(0.0) ]
        b       = [ float(0.0) for col in range( numVisibleSat ) ]
        A       = [ [float(0.0), float(0.0), float(0.0), float(0.0) ] for row in range( numVisibleSat ) ]
        AA_     = [ [float(0.0), float(0.0), float(0.0), float(0.0) ] for row in range( 4 ) ]

        if propagationDelayCompensation :
            
            for ii in range( 100 ) :
                Ab = [ float(0.0), float(0.0), float(0.0), float(0.0) ]
                for ll in range( numVisibleSat ) :
                    dXyz    = [ rcvXyzt[0] - satXyzvt[ll][1], rcvXyzt[1] - satXyzvt[ll][2], rcvXyzt[2] - satXyzvt[ll][3] ]
                    tau     = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] ) / C0
                    dXyz    = [ dXyz[0]+satXyzvt[ll][4]*tau, dXyz[1]+satXyzvt[ll][5]*tau, dXyz[2]+satXyzvt[ll][6]*tau ]
                    predict = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] )
                    observ  = codeShifts[ll] * L1_CHIP_M
                    b[ll], observ = observ - predict % L1_CODE_M, observ + floor( predict / L1_CODE_M ) * L1_CODE_M
                    A[ll]   = [ dXyz[0] / observ, dXyz[1] / observ, dXyz[2] / observ, 1  ]
                    Ab = [ Ab[0] + A[ll][0]*b[ll], Ab[1] + A[ll][1]*b[ll], Ab[2] + A[ll][2]*b[ll], Ab[3] + A[ll][3]*b[ll] ]
                
                AA = [ [float(0.0), float(0.0), float(0.0), float(0.0) ] for row in range( 4 ) ]
                for jj in range( 4 ) :
                    for kk in range( 4 ) :
                        for ll in range( numVisibleSat ) :
                            AA[jj][kk] += ( A[ll][jj] * A[ll][kk] )
                ##########################################################################
                # http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
                # http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
                ###########################################################################
                AA_[0][0] = AA[1][1]*AA[2][2]*AA[3][3]+AA[1][2]*AA[2][3]*AA[3][1]+AA[1][3]*AA[2][1]*AA[3][2]-AA[1][1]*AA[2][3]*AA[3][2]-AA[1][2]*AA[2][1]*AA[3][3]-AA[1][3]*AA[2][2]*AA[3][1]
                AA_[0][1] = AA[0][1]*AA[2][3]*AA[3][2]+AA[0][2]*AA[2][1]*AA[3][3]+AA[0][3]*AA[2][2]*AA[3][1]-AA[0][1]*AA[2][2]*AA[3][3]-AA[0][2]*AA[2][3]*AA[3][1]-AA[0][3]*AA[2][1]*AA[3][2]
                AA_[0][2] = AA[0][1]*AA[1][2]*AA[3][3]+AA[0][2]*AA[1][3]*AA[3][1]+AA[0][3]*AA[1][1]*AA[3][2]-AA[0][1]*AA[1][3]*AA[3][2]-AA[0][2]*AA[1][1]*AA[3][3]-AA[0][3]*AA[1][2]*AA[3][1]
                AA_[0][3] = AA[0][1]*AA[1][3]*AA[2][2]+AA[0][2]*AA[1][1]*AA[2][3]+AA[0][3]*AA[1][2]*AA[2][1]-AA[0][1]*AA[1][2]*AA[2][3]-AA[0][2]*AA[1][3]*AA[2][1]-AA[0][3]*AA[1][1]*AA[2][2]
                AA_[1][0] = AA[1][0]*AA[2][3]*AA[3][2]+AA[1][2]*AA[2][0]*AA[3][3]+AA[1][3]*AA[2][2]*AA[3][0]-AA[1][0]*AA[2][2]*AA[3][3]-AA[1][2]*AA[2][3]*AA[3][0]-AA[1][3]*AA[2][0]*AA[3][2]
                AA_[1][1] = AA[0][0]*AA[2][2]*AA[3][3]+AA[0][2]*AA[2][3]*AA[3][0]+AA[0][3]*AA[2][0]*AA[3][2]-AA[0][0]*AA[2][3]*AA[3][2]-AA[0][2]*AA[2][0]*AA[3][3]-AA[0][3]*AA[2][2]*AA[3][0]
                AA_[1][2] = AA[0][0]*AA[1][3]*AA[3][2]+AA[0][2]*AA[1][0]*AA[3][3]+AA[0][3]*AA[1][2]*AA[3][0]-AA[0][0]*AA[1][2]*AA[3][3]-AA[0][2]*AA[1][3]*AA[3][0]-AA[0][3]*AA[1][0]*AA[3][2]
                AA_[1][3] = AA[0][0]*AA[1][2]*AA[2][3]+AA[0][2]*AA[1][3]*AA[2][0]+AA[0][3]*AA[1][0]*AA[2][2]-AA[0][0]*AA[1][3]*AA[2][2]-AA[0][2]*AA[1][0]*AA[2][3]-AA[0][3]*AA[1][2]*AA[2][0]
                AA_[2][0] = AA[1][0]*AA[2][1]*AA[3][3]+AA[1][1]*AA[2][3]*AA[3][0]+AA[1][3]*AA[2][0]*AA[3][1]-AA[1][0]*AA[2][3]*AA[3][1]-AA[1][1]*AA[2][0]*AA[3][3]-AA[1][3]*AA[2][1]*AA[3][0]
                AA_[2][1] = AA[0][0]*AA[2][3]*AA[3][1]+AA[0][1]*AA[2][0]*AA[3][3]+AA[0][3]*AA[2][1]*AA[3][0]-AA[0][0]*AA[2][1]*AA[3][3]-AA[0][1]*AA[2][3]*AA[3][0]-AA[0][3]*AA[2][0]*AA[3][1]
                AA_[2][2] = AA[0][0]*AA[1][1]*AA[3][3]+AA[0][1]*AA[1][3]*AA[3][0]+AA[0][3]*AA[1][0]*AA[3][1]-AA[0][0]*AA[1][3]*AA[3][1]-AA[0][1]*AA[1][0]*AA[3][3]-AA[0][3]*AA[1][1]*AA[3][0]
                AA_[2][3] = AA[0][0]*AA[1][3]*AA[2][1]+AA[0][1]*AA[1][0]*AA[2][3]+AA[0][3]*AA[1][1]*AA[2][0]-AA[0][0]*AA[1][1]*AA[2][3]-AA[0][1]*AA[1][3]*AA[2][0]-AA[0][3]*AA[1][0]*AA[2][1]
                AA_[3][0] = AA[1][0]*AA[2][2]*AA[3][1]+AA[1][1]*AA[2][0]*AA[3][2]+AA[1][2]*AA[2][1]*AA[3][0]-AA[1][0]*AA[2][1]*AA[3][2]-AA[1][1]*AA[2][2]*AA[3][0]-AA[1][2]*AA[2][0]*AA[3][1]
                AA_[3][1] = AA[0][0]*AA[2][1]*AA[3][2]+AA[0][1]*AA[2][2]*AA[3][0]+AA[0][2]*AA[2][0]*AA[3][1]-AA[0][0]*AA[2][2]*AA[3][1]-AA[0][1]*AA[2][0]*AA[3][2]-AA[0][2]*AA[2][1]*AA[3][0]
                AA_[3][2] = AA[0][0]*AA[1][2]*AA[3][1]+AA[0][1]*AA[1][0]*AA[3][2]+AA[0][2]*AA[1][1]*AA[3][0]-AA[0][0]*AA[1][1]*AA[3][2]-AA[0][1]*AA[1][2]*AA[3][0]-AA[0][2]*AA[1][0]*AA[3][1]
                AA_[3][3] = AA[0][0]*AA[1][1]*AA[2][2]+AA[0][1]*AA[1][2]*AA[2][0]+AA[0][2]*AA[1][0]*AA[2][1]-AA[0][0]*AA[1][2]*AA[2][1]-AA[0][1]*AA[1][0]*AA[2][2]-AA[0][2]*AA[1][1]*AA[2][0]

                detA = AA[0][0]*AA_[0][0]  + AA[0][1]*AA_[1][0] + AA[0][2]*AA_[2][0] + AA[0][3]*AA_[3][0]
                if detA :                
                    rcvXyzt = [ rcvXyzt[0] + (AA_[0][0]*Ab[0]+AA_[0][1]*Ab[1]+AA_[0][2]*Ab[2]+AA_[0][3]*Ab[3])/detA ,
                                rcvXyzt[1] + (AA_[1][0]*Ab[0]+AA_[1][1]*Ab[1]+AA_[1][2]*Ab[2]+AA_[1][3]*Ab[3])/detA , 
                                rcvXyzt[2] + (AA_[2][0]*Ab[0]+AA_[2][1]*Ab[1]+AA_[2][2]*Ab[2]+AA_[2][3]*Ab[3])/detA ,
                                rcvXyzt[3] + (AA_[3][0]*Ab[0]+AA_[3][1]*Ab[1]+AA_[3][2]*Ab[2]+AA_[3][3]*Ab[3])/detA ]
                else :
                    print ( "Ephemeris.Get_XyzPosition():: detA = %f " % detA )

        else :  # propagationDelayCompensation == False
            ##
            # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
            #
            for ii in range( 100 ) :
                Ab = [ float(0.0), float(0.0), float(0.0), float(0.0) ]
                for ll in range( numVisibleSat ) :
                    dXyz    = [ rcvXyzt[0] - satXyzvt[ll][1], rcvXyzt[1] - satXyzvt[ll][2], rcvXyzt[2] - satXyzvt[ll][3] ]
                    predict = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] )
                    observ  = codeShifts[ll] * L1_CHIP_M
                    b[ll], observ = observ - predict % L1_CODE_M, observ + floor( predict / L1_CODE_M ) * L1_CODE_M
                    A[ll]   = [ dXyz[0] / observ, dXyz[1] / observ, dXyz[2] / observ, 1  ]
                    Ab = [ Ab[0] + A[ll][0]*b[ll], Ab[1] + A[ll][1]*b[ll], Ab[2] + A[ll][2]*b[ll], Ab[3] + A[ll][3]*b[ll] ]
                
                AA = [ [float(0.0), float(0.0), float(0.0), float(0.0) ] for row in range( 4 ) ]
                for jj in range( 4 ) :
                    for kk in range( 4 ) :
                        for ll in range( numVisibleSat ) :
                            AA[jj][kk] += ( A[ll][jj] * A[ll][kk] )
                ##########################################################################
                # http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
                # http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
                ###########################################################################
                AA_[0][0] = AA[1][1]*AA[2][2]*AA[3][3]+AA[1][2]*AA[2][3]*AA[3][1]+AA[1][3]*AA[2][1]*AA[3][2]-AA[1][1]*AA[2][3]*AA[3][2]-AA[1][2]*AA[2][1]*AA[3][3]-AA[1][3]*AA[2][2]*AA[3][1]
                AA_[0][1] = AA[0][1]*AA[2][3]*AA[3][2]+AA[0][2]*AA[2][1]*AA[3][3]+AA[0][3]*AA[2][2]*AA[3][1]-AA[0][1]*AA[2][2]*AA[3][3]-AA[0][2]*AA[2][3]*AA[3][1]-AA[0][3]*AA[2][1]*AA[3][2]
                AA_[0][2] = AA[0][1]*AA[1][2]*AA[3][3]+AA[0][2]*AA[1][3]*AA[3][1]+AA[0][3]*AA[1][1]*AA[3][2]-AA[0][1]*AA[1][3]*AA[3][2]-AA[0][2]*AA[1][1]*AA[3][3]-AA[0][3]*AA[1][2]*AA[3][1]
                AA_[0][3] = AA[0][1]*AA[1][3]*AA[2][2]+AA[0][2]*AA[1][1]*AA[2][3]+AA[0][3]*AA[1][2]*AA[2][1]-AA[0][1]*AA[1][2]*AA[2][3]-AA[0][2]*AA[1][3]*AA[2][1]-AA[0][3]*AA[1][1]*AA[2][2]
                AA_[1][0] = AA[1][0]*AA[2][3]*AA[3][2]+AA[1][2]*AA[2][0]*AA[3][3]+AA[1][3]*AA[2][2]*AA[3][0]-AA[1][0]*AA[2][2]*AA[3][3]-AA[1][2]*AA[2][3]*AA[3][0]-AA[1][3]*AA[2][0]*AA[3][2]
                AA_[1][1] = AA[0][0]*AA[2][2]*AA[3][3]+AA[0][2]*AA[2][3]*AA[3][0]+AA[0][3]*AA[2][0]*AA[3][2]-AA[0][0]*AA[2][3]*AA[3][2]-AA[0][2]*AA[2][0]*AA[3][3]-AA[0][3]*AA[2][2]*AA[3][0]
                AA_[1][2] = AA[0][0]*AA[1][3]*AA[3][2]+AA[0][2]*AA[1][0]*AA[3][3]+AA[0][3]*AA[1][2]*AA[3][0]-AA[0][0]*AA[1][2]*AA[3][3]-AA[0][2]*AA[1][3]*AA[3][0]-AA[0][3]*AA[1][0]*AA[3][2]
                AA_[1][3] = AA[0][0]*AA[1][2]*AA[2][3]+AA[0][2]*AA[1][3]*AA[2][0]+AA[0][3]*AA[1][0]*AA[2][2]-AA[0][0]*AA[1][3]*AA[2][2]-AA[0][2]*AA[1][0]*AA[2][3]-AA[0][3]*AA[1][2]*AA[2][0]
                AA_[2][0] = AA[1][0]*AA[2][1]*AA[3][3]+AA[1][1]*AA[2][3]*AA[3][0]+AA[1][3]*AA[2][0]*AA[3][1]-AA[1][0]*AA[2][3]*AA[3][1]-AA[1][1]*AA[2][0]*AA[3][3]-AA[1][3]*AA[2][1]*AA[3][0]
                AA_[2][1] = AA[0][0]*AA[2][3]*AA[3][1]+AA[0][1]*AA[2][0]*AA[3][3]+AA[0][3]*AA[2][1]*AA[3][0]-AA[0][0]*AA[2][1]*AA[3][3]-AA[0][1]*AA[2][3]*AA[3][0]-AA[0][3]*AA[2][0]*AA[3][1]
                AA_[2][2] = AA[0][0]*AA[1][1]*AA[3][3]+AA[0][1]*AA[1][3]*AA[3][0]+AA[0][3]*AA[1][0]*AA[3][1]-AA[0][0]*AA[1][3]*AA[3][1]-AA[0][1]*AA[1][0]*AA[3][3]-AA[0][3]*AA[1][1]*AA[3][0]
                AA_[2][3] = AA[0][0]*AA[1][3]*AA[2][1]+AA[0][1]*AA[1][0]*AA[2][3]+AA[0][3]*AA[1][1]*AA[2][0]-AA[0][0]*AA[1][1]*AA[2][3]-AA[0][1]*AA[1][3]*AA[2][0]-AA[0][3]*AA[1][0]*AA[2][1]
                AA_[3][0] = AA[1][0]*AA[2][2]*AA[3][1]+AA[1][1]*AA[2][0]*AA[3][2]+AA[1][2]*AA[2][1]*AA[3][0]-AA[1][0]*AA[2][1]*AA[3][2]-AA[1][1]*AA[2][2]*AA[3][0]-AA[1][2]*AA[2][0]*AA[3][1]
                AA_[3][1] = AA[0][0]*AA[2][1]*AA[3][2]+AA[0][1]*AA[2][2]*AA[3][0]+AA[0][2]*AA[2][0]*AA[3][1]-AA[0][0]*AA[2][2]*AA[3][1]-AA[0][1]*AA[2][0]*AA[3][2]-AA[0][2]*AA[2][1]*AA[3][0]
                AA_[3][2] = AA[0][0]*AA[1][2]*AA[3][1]+AA[0][1]*AA[1][0]*AA[3][2]+AA[0][2]*AA[1][1]*AA[3][0]-AA[0][0]*AA[1][1]*AA[3][2]-AA[0][1]*AA[1][2]*AA[3][0]-AA[0][2]*AA[1][0]*AA[3][1]
                AA_[3][3] = AA[0][0]*AA[1][1]*AA[2][2]+AA[0][1]*AA[1][2]*AA[2][0]+AA[0][2]*AA[1][0]*AA[2][1]-AA[0][0]*AA[1][2]*AA[2][1]-AA[0][1]*AA[1][0]*AA[2][2]-AA[0][2]*AA[1][1]*AA[2][0]

                detA = AA[0][0]*AA_[0][0] + AA[0][1]*AA_[1][0] + AA[0][2]*AA_[2][0] + AA[0][3]*AA_[3][0]
                if detA :                
                    rcvXyzt = [ rcvXyzt[0] + (AA_[0][0]*Ab[0]+AA_[0][1]*Ab[1]+AA_[0][2]*Ab[2]+AA_[0][3]*Ab[3])/detA ,
                                rcvXyzt[1] + (AA_[1][0]*Ab[0]+AA_[1][1]*Ab[1]+AA_[1][2]*Ab[2]+AA_[1][3]*Ab[3])/detA , 
                                rcvXyzt[2] + (AA_[2][0]*Ab[0]+AA_[2][1]*Ab[1]+AA_[2][2]*Ab[2]+AA_[2][3]*Ab[3])/detA ,
                                rcvXyzt[3] + (AA_[3][0]*Ab[0]+AA_[3][1]*Ab[1]+AA_[3][2]*Ab[2]+AA_[3][3]*Ab[3])/detA ]
                else :
                    print ( "Ephemeris.Get_XyzPosition():: detA = %f " % detA )

        return rcvXyzt

    def Get_Receiver_XyzPosition_Numpy( self, epoch, visibleSats, codeShifts, rcvXyz = [0, 0, 0], propagationDelayCompensation = False ) :
        lari            = [ float(0.0) for col in range(6) ]
        xyOrb           = [ float(0.0) for col in range(4) ]

        numVisibleSat   = len( visibleSats )
        satXyzvt        = [ list([ float(0.0) for col in range(8) ]) for row in range( numVisibleSat ) ]
        rcvXyzt         = [ rcvXyz[0], rcvXyz[1], rcvXyz[2], 0 ]

        for ii in range( numVisibleSat ) :
            sat = visibleSats[ii] - 1
            if self.svid[sat] == 0 :
                continue
            satXyzvt[ii][0] = sat
            tk = epoch[1] - self.toe[sat] #+ ( epoch[0] - self.week[sat] ) * WEEK_SECONDS
            if tk > 302400 :
                tk -= WEEK_SECONDS
            if tk < -302400 :
                tk += WEEK_SECONDS
            #############################################################################
            # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
            # The time of transmissionis is modified from TOW from navigation data and
            # the correction of is performed as
            #############################################################################
            A = self.sqrtA[sat]*self.sqrtA[sat]
            meanMotion = sqrt(WGS84_u_m3s2) / (self.sqrtA[sat]*A) + self.deltan[sat]# This value should be caliberated with adding Mean Motion Difference.
            meanAnomaly  = self.M0[sat] + meanMotion * tk

                #######################################################################################################################
                # Kepler equation for calculating Eccentric Anomaly
                #
                # In orbital mechanics, Kepler's equation relates various geometric properties of the orbit
                # of a body subject to a central force.
                #
                # ecceAnomaly = Kepler( meanAnomaly, Eccentricity[sat] )
                # := ecceAnomaly - (ecceAnomaly-Eccentricity[sat]*sin(ecceAnomaly)-meanAnomaly)/(1-Eccentricity[sat]*cos(ecceAnomaly))
                #
                # Reference: https://en.wikipedia.org/wiki/Kepler's_equation
                ########################################################################################################################
            E0 = meanAnomaly
            ecceAnomaly = meanAnomaly + self.e[sat]*sin(meanAnomaly) / ( 1 - self.e[sat]*cos(meanAnomaly) )
            while( abs(ecceAnomaly-E0) > 1.0e-14 ) :
                E0 = ecceAnomaly
                ecceAnomaly = ecceAnomaly - ( ecceAnomaly - self.e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - self.e[sat]*cos(ecceAnomaly) )
            sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
            e_cosEcceAnom_1 = 1 -  self.e[sat]*cosEcceAnom

            lari[0] = atan2( sqrt(1-self.e[sat]*self.e[sat]) * sinEcceAnom,  cosEcceAnom-self.e[sat] ) #True anomaly. Later, Argument of latitude
            lari[3] = sinEcceAnom * meanMotion * (1+self.e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

            lari[0] += self.omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            lari[0:3] = [ lari[0] + self.Cus[sat]*sin2Lat + self.Cuc[sat]*cos2Lat ,
                A*e_cosEcceAnom_1 + self.Crs[sat]*sin2Lat + self.Crc[sat]*cos2Lat ,
                     self.i0[sat] + self.Cis[sat]*sin2Lat + self.Cic[sat]*cos2Lat + self.idot[sat]*tk ]

            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            ddLat, ddRad, ddInc = 2 * ( self.Cus[sat] * cos2Lat - self.Cuc[sat] * sin2Lat ) * lari[3], 2 * ( self.Crs[sat] * cos2Lat - self.Crc[sat] * sin2Lat ) * lari[3], 2 * ( self.Cis[sat] * cos2Lat - self.Cic[sat] * sin2Lat ) * lari[3]

            lari[3:6] = [ lari[3] + ddLat, A*self.e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, self.idot[sat] + ddInc ]

            cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
            satXyzvt[ii][7] =  self.af0[sat] + ( self.af1[sat] + self.af2[sat]*tk ) * tk + F*self.e[sat]*self.sqrtA[sat]*sinEcceAnom - self.tgd[sat] #Calculate clock bias error of satellites

            dotAsceLon = self.Omegadot[sat] - WGS84_w_rads
            asceLon = self.Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * self.toe[sat]
            cosLon, sinLon = cos( asceLon ), sin( asceLon )

            xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
            xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

            satXyzvt[ii][1:4] = [+cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1], 
                                +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                                +sinInc*xyOrb[1] ]

            satXyzvt[ii][4:7] = [+cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satXyzvt[ii][2] ,
                                +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satXyzvt[ii][1] ,
                                +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

        ########################################################################################################
        # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
        ########################################################################################################
        dXyz    = [ float(0.0) for col in range( 3 ) ]
        b       = [ float(0.0) for col in range( numVisibleSat ) ]
        A       = [ list([ float(0.0) for col in range( 4 ) ]) for row in range( numVisibleSat ) ]
        AA_     = [ list([ float(0.0) for col in range( 4 ) ]) for row in range( 4 ) ]

        if propagationDelayCompensation :
            
            for ii in range( 1000 ) :
                Ab = [ float(0.0) for row in range( 4 ) ]
                for ll in range( numVisibleSat ) :
                    dXyz    = [ rcvXyzt[0] - satXyzvt[ll][1], rcvXyzt[1] - satXyzvt[ll][2], rcvXyzt[2] - satXyzvt[ll][3] ]
                    tau = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] ) / C0
                    dXyz    = [ dXyz[0]+satXyzvt[ll][4]*tau, dXyz[1]+satXyzvt[ll][5]*tau, dXyz[2]+satXyzvt[ll][6]*tau ]
                    predict = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] )
                    observ  = codeShifts[ll] * L1_CHIP_M
                    b[ll], observ = observ - predict % L1_CODE_M, observ + floor( predict / L1_CODE_M ) * L1_CODE_M
                    A[ll]   = [ dXyz[0] / observ, dXyz[1] / observ, dXyz[2] / observ, 1  ]
                
                dRcvXyzt = np.linalg.solve( np.array( A ), np.array( b ) )[0]
                rcvXyzt  = [rcvXyzt[0] + dRcvXyzt[0], rcvXyzt[1] + dRcvXyzt[1], rcvXyzt[2] + dRcvXyzt[2], rcvXyzt[3] + dRcvXyzt[3] ]
                

        else :    
            ## This is 
            # 
            # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
            #
            for ii in range( 1000 ) :
                Ab = [ float(0.0) for row in range( 4 ) ]
                for ll in range( numVisibleSat ) :
                    dXyz    = [ rcvXyzt[0] - satXyzvt[ll][0], rcvXyzt[1] - satXyzvt[ll][1], rcvXyzt[2] - satXyzvt[ll][2] ]
                    predict = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] )
                    observ  = codeShifts[ll] * L1_CHIP_M
                    b[ll], observ = observ - predict % L1_CODE_M, observ + floor( predict / L1_CODE_M ) * L1_CODE_M
                    A[ll]   = [ dXyz[0] / observ, dXyz[1] / observ, dXyz[2] / observ, 1  ]

                dRcvXyzt = np.linalg.solve( np.array( A ), np.array( b ) )
                rcvXyzt  = [rcvXyzt[0] + dRcvXyzt[0], rcvXyzt[1] + dRcvXyzt[1], rcvXyzt[2] + dRcvXyzt[2], rcvXyzt[3] + dRcvXyzt[3] ]

        return rcvXyzt


    def Generate_Subframes( self, sat, sfn = None ):
        PWR2_31_    = 2**31
        PWR2_43_    = 2**43
        PWR2_55_    = 2**55
        PWR2_5_     = 2**5
        PWR2_29_    = 2**29
        PWR2_19_    = 2**19
        PWR2_31_PI_ = PWR2_31_ / ICD200PI
        PWR2_43_PI_ = PWR2_43_ / ICD200PI
        ## 
        # https://github.com/osqzss/gps-sdr-sim/blob/master/gpssim.c
        #
        WN, URA, DATA_ID, SF4_PAGE25_SVID, SF5_PAGE25_SVID = 0, 0, 1, 0x3F, 0x1B
        if sfn == None :
            sfLen = 200
            sf = bytearray( sfLen )
            sfView = memoryview( sf )

            # Subframe 1
            sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
            sf[4], sf[5], sf[6], sf[7] = 0x00, 0x01, 0x00, 0x00
            value = ( ( WN&0x3FF ) << 20 ) | ( ( self.l2code[sat]&0x3 ) << 18 ) | ( ( URA&0xF ) << 14 ) | ( ( self.health[sat]&0x3 ) << 8 ) | ( ( self.iodc[sat]&0x3 ) << 6 )
            sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            sf[12], sf[13], sf[14], sf[15] = 0x00, 0x00, 0x00, 0x00
            sf[16], sf[17], sf[18], sf[19] = 0x00, 0x00, 0x00, 0x00
            sf[20], sf[21], sf[22], sf[23] = 0x00, 0x00, 0x00, 0x00
            value = ( int(self.tgd[sat]*PWR2_31_) & 0xFF ) << 6
            sf[24], sf[25], sf[26], sf[27] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.iodc[sat] & 0xFF ) << 22 ) | ( ( self.toc[sat] & 0xFFFF ) << 6 )
            sf[28], sf[29], sf[30], sf[31] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( int(self.af2[sat]*PWR2_55_) & 0xFF ) << 22 ) | ( ( int(self.af1[sat]*PWR2_43_) & 0xFFFF ) << 6 )
            sf[32], sf[33], sf[34], sf[35] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( int(self.af0[sat]*PWR2_31_) & 0x3FFFFF ) << 8
            sf[36], sf[37], sf[38], sf[39] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF

            # Subframe 2
            sf[40], sf[41], sf[42], sf[43] = 0x22, 0xC0, 0x00, 0x00
            sf[44], sf[45], sf[46], sf[47] = 0x00, 0x02, 0x00, 0x00
            value = ( ( self.iode[sat] & 0xFF ) << 22 ) | ( ( int(self.Crs[sat]*PWR2_5_) & 0xFFFF ) << 6 ) 
            sf[48], sf[49], sf[50], sf[51] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            M0 = int( self.M0[sat] * PWR2_31_PI_ )
            value = ( ( int(self.deltan[sat]*PWR2_43_PI_) & 0xFFFF ) << 14 ) | ( ( ( M0 >> 24 ) & 0xFF ) << 6 ) 
            sf[52], sf[53], sf[54], sf[55] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( M0&0xFFFFFF ) << 6         
            sf[56], sf[57], sf[58], sf[59] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            e = int( self.e[sat] * ( 2**33 ) )
            value = ( ( int(self.Cuc[sat])&0xFFFF ) << 14 ) | ( ( ( e >> 24 )&0xFF ) << 6 ) 
            sf[60], sf[61], sf[62], sf[63] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( e & 0xFFFFFF ) << 6
            sf[64], sf[65], sf[66], sf[67] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            sqrtA = int( self.sqrtA[sat] * PWR2_19_ )
            value = ( ( int( self.Cus[sat] * PWR2_29_ ) & 0xFFFF ) << 14 ) | ( ( ( sqrtA >> 24 ) & 0xFFFF ) << 6 )
            sf[68], sf[69], sf[70], sf[71] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( sqrtA & 0xFFFFFF ) << 6
            sf[72], sf[73], sf[74], sf[75] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.toe[sat] & 0xFFFF ) << 14
            sf[76], sf[77], sf[78], sf[79] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF

            # Subframe 3
            sf[80], sf[81], sf[82], sf[83] = 0x22, 0xC0, 0x00, 0x00
            sf[84], sf[85], sf[86], sf[87] = 0x00, 0x03, 0x00, 0x00
            Omega0 = int( self.Omega0[sat] * PWR2_31_PI_ )
            value = ( ( int( self.Cic[sat] * PWR2_29_ ) & 0xFFFF ) << 14 ) | ( ( ( Omega0 >> 24 )&0xFF ) << 6 ) 
            sf[88], sf[89], sf[90], sf[91] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( Omega0 & 0xFFFFFF ) << 6         
            sf[92], sf[93], sf[94], sf[95] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            i0 = int( self.i0[sat] * PWR2_31_PI_ )
            value = ( ( int( self.Cis[sat] * PWR2_31_PI_ ) & 0xFFFF ) << 14 ) | ( ( ( i0 >> 24 ) & 0xFF ) << 6 ) 
            sf[96], sf[97], sf[98], sf[99] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( i0 & 0xFFFFFF ) <<  6  
            sf[100], sf[101], sf[102], sf[103] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            omega = self.omega[sat] * PWR2_31_PI_
            value = ( ( int( self.Crc[sat] * PWR2_5_ ) & 0xFFFF ) << 14 ) | ( ( ( int( omega ) >> 24 ) & 0xFF ) << 6 )
            sf[104], sf[105], sf[106], sf[107] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( int( omega ) & 0xFFFFFF ) << 6
            sf[108], sf[109], sf[110], sf[111] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( int( self.Omegadot[sat] * PWR2_43_PI_ ) & 0xFFFFFF ) << 6
            sf[112], sf[113], sf[114], sf[115] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.iode[sat] & 0xFF ) << 22 ) | ( ( int( self.idot[sat] * PWR2_43_PI_ ) & 0x3FFF ) << 8 )
            sf[116], sf[117], sf[118], sf[119] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF

            # Subframe 4
            if True :
                sf[120], sf[121], sf[122], sf[123] = 0x22, 0xC0, 0x00, 0x00
                sf[124], sf[125], sf[126], sf[127] = 0x00, 0x04, 0x00, 0x00
                value = ( DATA_ID << 28 ) | ( SF4_PAGE18_SVID << 22 ) | ( ( SF4_PAGE18_SVID & 0xFF ) << 14 ) | ( ( SF4_PAGE18_SVID & 0xFF ) << 6 )
                sf[128], sf[129], sf[130], sf[131] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( Omega0 & 0xFFFFFF ) << 6
                sf[132], sf[133], sf[134], sf[135] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.Cis[sat] & 0xFFFF ) << 14 ) | ( ( ( i0 >> 24 ) & 0xFF ) << 6 )
                sf[136], sf[137], sf[138], sf[139] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( i0 & 0xFFFFFF ) <<  6
                sf[140], sf[141], sf[142], sf[143] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.Crc[sat] & 0xFFFF ) << 14 ) | ( ( ( omega >> 24 ) & 0xFF ) << 6 )
                sf[144], sf[145], sf[146], sf[147] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( omega & 0xFFFFFF ) << 6
                sf[148], sf[149], sf[150], sf[151] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( self.Omegadot[sat] & 0xFFFFFF ) << 6
                sf[152], sf[153], sf[154], sf[155] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.iode[sat] & 0xFF ) << 22 ) | ( ( int(self.idot[sat]) & 0x3FFF ) << 8 )
                sf[156], sf[157], sf[158], sf[159] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            else :
                sf[120], sf[121], sf[122], sf[123] = 0x22, 0xC0, 0x00, 0x00
                sf[124], sf[125], sf[126], sf[127] = 0x00, 0x04, 0x00, 0x00
                value = ( DATA_ID << 28 ) | ( SF4_PAGE25_SVID << 22 )
                sf[128], sf[129], sf[130], sf[131] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                sf[132], sf[133], sf[134], sf[135] = 0x00, 0x04, 0x00, 0x00
                sf[136], sf[137], sf[138], sf[139] = 0x00, 0x04, 0x00, 0x00
                sf[140], sf[141], sf[142], sf[143] = 0x00, 0x04, 0x00, 0x00
                sf[144], sf[145], sf[146], sf[147] = 0x00, 0x04, 0x00, 0x00
                sf[148], sf[149], sf[150], sf[151] = 0x00, 0x04, 0x00, 0x00
                sf[152], sf[153], sf[154], sf[155] = 0x00, 0x04, 0x00, 0x00
                sf[156], sf[157], sf[158], sf[159] = 0x00, 0x04, 0x00, 0x00
            
            # Subframe 5
            sf[160], sf[161], sf[162], sf[163] = 0x22, 0xC0, 0x00, 0x00
            sf[164], sf[165], sf[166], sf[167] = 0x00, 0x05, 0x00, 0x00
            value = ( DATA_ID << 28 ) | ( SF5_PAGE25_SVID << 22 ) | ( ( ( self.toe[sat] / 4096 ) &0xFF ) << 14 ) | ( ( ( self.week[sat] % 256 ) &0xFF ) << 6 )
            sf[168], sf[169], sf[170], sf[171] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            sf[172], sf[173], sf[174], sf[175] = 0x00, 0x04, 0x00, 0x00
            sf[176], sf[177], sf[178], sf[179] = 0x00, 0x04, 0x00, 0x00
            sf[180], sf[181], sf[182], sf[183] = 0x00, 0x04, 0x00, 0x00
            sf[184], sf[185], sf[186], sf[187] = 0x00, 0x04, 0x00, 0x00
            sf[188], sf[189], sf[190], sf[191] = 0x00, 0x04, 0x00, 0x00
            sf[192], sf[193], sf[194], sf[195] = 0x00, 0x04, 0x00, 0x00
            sf[196], sf[197], sf[198], sf[199] = 0x00, 0x04, 0x00, 0x00
        elif sfn == 1:  # Subframe 1
            sfLen = 40
            sf = bytearray( sfLen )
            sfView = memoryview( sf )
        
            sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
            sf[4], sf[5], sf[6], sf[7] = 0x00, 0x01, 0x00, 0x00
            value = ( ( WN&0x3FF ) << 20 ) | ( ( self.l2code[sat]&0x3 ) << 18 ) | ( ( URA&0xF ) << 14 ) | ( ( self.health[sat]&0x3 ) << 8 ) | ( ( self.iodc[sat]&0x3 ) << 6 )
            sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            sf[12], sf[13], sf[14], sf[15] = 0x00, 0x00, 0x00, 0x00
            sf[16], sf[17], sf[18], sf[19] = 0x00, 0x00, 0x00, 0x00
            sf[20], sf[21], sf[22], sf[23] = 0x00, 0x00, 0x00, 0x00
            value = ( int(self.tgd[sat]*PWR2_31_)& 0xFF ) << 6
            sf[24], sf[25], sf[26], sf[27] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.iodc[sat] & 0xFF ) << 22 ) | ( ( self.toc[sat] & 0xFFFF ) << 6 )
            sf[28], sf[29], sf[30], sf[31] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.af2[sat] & 0xFF ) << 22 ) | ( ( self.af1[sat] & 0xFFFF ) << 6 )
            sf[32], sf[33], sf[34], sf[35] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.af0[sat] & 0x3FFFFF ) << 8
            sf[36], sf[37], sf[38], sf[39] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            
        elif sfn == 2: # Subframe 2
            sfLen = 40
            sf = bytearray( sfLen )
            sfView = memoryview( sf )

            sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
            sf[4], sf[5], sf[6], sf[7] = 0x00, 0x02, 0x00, 0x00
            value = ( ( self.iode[sat]&0xFF ) << 22 ) | ( ( self.Crs[sat]&0xFFFF ) << 6 ) 
            sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.deltan[sat]&0xFFFF ) << 14 ) | ( ( ( self.M0[sat] >> 24 )&0xFF ) << 6 ) 
            sf[12], sf[13], sf[14], sf[15] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.M0[sat]&0xFFFFFF ) << 6         
            sf[16], sf[17], sf[18], sf[19] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.Cuc[sat]&0xFFFF ) << 14 ) | ( ( ( self.e[sat] >> 24 )&0xFF ) << 6 ) 
            sf[20], sf[21], sf[22], sf[23] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.e[sat] & 0xFFFFFF ) << 6
            sf[24], sf[25], sf[26], sf[27] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( ( self.Cus[sat] & 0xFFFF ) << 14 ) | ( ( ( self.sqrtA[sat] >> 24 ) & 0xFFFF ) << 6 )
            sf[28], sf[29], sf[30], sf[31] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.sqrtA[sat] & 0xFFFFFF ) << 6
            sf[32], sf[33], sf[34], sf[35] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            value = ( self.toe[sat] & 0xFFFF ) << 14
            sf[36], sf[37], sf[38], sf[39] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                
        elif sfn == 3: # Subframe 3
            sfLen = 40
            sf = bytearray( sfLen )
            sfView = memoryview( sf )

            if True :
                sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
                sf[4], sf[5], sf[6], sf[7] = 0x00, 0x04, 0x00, 0x00
                value = ( DATA_ID << 28 ) | ( SF4_PAGE18_SVID << 22 ) | ( ( SF4_PAGE18_SVID & 0xFF ) << 14 ) | ( ( SF4_PAGE18_SVID & 0xFF ) << 6 )
                sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( self.Omega0[sat]&0xFFFFFF ) << 6
                sf[12], sf[13], sf[14], sf[15] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.Cis[sat]&0xFFFF ) << 14 ) | ( ( ( self.i0[sat] >> 24 )&0xFF ) << 6 )
                sf[16], sf[17], sf[18], sf[19] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( self.i0[sat]&0xFFFFFF ) <<  6
                sf[20], sf[21], sf[22], sf[23] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.Crc[sat] & 0xFFFF ) << 14 ) | ( ( ( self.omega[sat] >> 24 ) & 0xFF ) << 6 )
                sf[24], sf[25], sf[26], sf[27] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( self.omega[sat] & 0xFFFFFF ) << 6
                sf[28], sf[29], sf[30], sf[31] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( self.Omegadot[sat] & 0xFFFFFF ) << 6
                sf[32], sf[33], sf[34], sf[35] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                value = ( ( self.iode[sat] & 0xFF ) << 22 ) | ( ( int(self.idot[sat]) & 0x3FFF ) << 8 )
                sf[36], sf[37], sf[38], sf[39] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            else :
                sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
                sf[4], sf[5], sf[6], sf[7] = 0x00, 0x04, 0x00, 0x00
                value = ( DATA_ID << 28 ) | ( SF4_PAGE25_SVID << 22 )
                sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
                sf[12], sf[13], sf[14], sf[15] = 0x00, 0x04, 0x00, 0x00
                sf[16], sf[17], sf[18], sf[19] = 0x00, 0x04, 0x00, 0x00
                sf[20], sf[21], sf[22], sf[23] = 0x00, 0x04, 0x00, 0x00
                sf[24], sf[25], sf[26], sf[27] = 0x00, 0x04, 0x00, 0x00
                sf[28], sf[29], sf[30], sf[31] = 0x00, 0x04, 0x00, 0x00
                sf[32], sf[33], sf[34], sf[35] = 0x00, 0x04, 0x00, 0x00
                sf[36], sf[37], sf[38], sf[39] = 0x00, 0x04, 0x00, 0x00

        elif sfn == 5: # Subframe 5
            sfLen = 40
            sf = bytearray( sfLen )
            sfView = memoryview( sf )

            sf[0], sf[1], sf[2], sf[3] = 0x22, 0xC0, 0x00, 0x00
            sf[4], sf[5], sf[6], sf[7] = 0x00, 0x05, 0x00, 0x00
            value = ( DATA_ID << 28 ) | ( SF5_PAGE25_SVID << 22 ) | ( ( ( self.toe[sat] / 4096 ) &0xFF ) << 14 ) | ( ( ( self.week[sat] % 256 ) &0xFF ) << 6 )
            sf[8], sf[9], sf[10], sf[11] = ( value & 0xFF000000 ) >> 24, ( value & 0xFF0000 ) >> 16, ( value & 0xFF00 ) >> 8, value & 0xFF
            sf[12], sf[13], sf[14], sf[15] = 0x00, 0x04, 0x00, 0x00
            sf[16], sf[17], sf[18], sf[19] = 0x00, 0x04, 0x00, 0x00
            sf[20], sf[21], sf[22], sf[23] = 0x00, 0x04, 0x00, 0x00
            sf[24], sf[25], sf[26], sf[27] = 0x00, 0x04, 0x00, 0x00
            sf[28], sf[29], sf[30], sf[31] = 0x00, 0x04, 0x00, 0x00
            sf[32], sf[33], sf[34], sf[35] = 0x00, 0x04, 0x00, 0x00
            sf[36], sf[37], sf[38], sf[39] = 0x00, 0x04, 0x00, 0x00
        
        
        return sf

    def unlzw( self, zData ):
    # This function was adapted for Python from Mark Adler's C implementation
    # https://github.com/umeat/unlzw
    # ftp.stu.edu.tw/pub/FreeBSD/branches/4.0-stable/src/gnu/usr.bin/gzip/unlzw.c
    # https://opensource.apple.com/source/gnuzip/gnuzip-22.2/gzip/unlzw.c.auto.html
    # http://aluigi.altervista.org/mytoolz/unlzw.c
    #
    # Decompress compressed data generated by the Unix compress utility (LZW
    # compression, files with .Z suffix). Input can be given as any type which
    # can be 'converted' to a bytearray (e.g. string, or bytearray). Returns 
    # decompressed data as string, or raises error.
    
    # Convert input data stream to byte array, and get length of that array
        try:
            ba_in = bytearray( zData )
        except ValueError:
            raise TypeError("Unable to convert inputted data to bytearray")
            
        inlen = len(ba_in)
        prefix = [None] * 65536         # index to LZW prefix string
        suffix = [None] * 65536         # one-character LZW suffix
        
        # Process header
        if inlen < 3:
            raise ValueError("Invalid Input: Length of input too short for processing")
        
        if (ba_in[0] != 0x1f) or (ba_in[1] != 0x9d):
            raise ValueError("Invalid Header Flags Byte: Incorrect magic bytes")
        
        flags = ba_in[2]
        if flags & 0x60:
            raise ValueError("Invalid Header Flags Byte: Flag byte contains invalid data")
            
        max_ = flags & 0x1f
        if (max_ < 9) or (max_ > 16):
            raise ValueError("Invalid Header Flags Byte: Max code size bits out of range")
            
        if (max_ == 9): max_ = 10       # 9 doesn't really mean 9 
        flags &= 0x80                   # true if block compressed
        
        # Clear table, start at nine bits per symbol
        bits = 9
        mask = 0x1ff
        end = 256 if flags else 255
        
        # Ensure stream is initially valid
        if inlen == 3: return 0         # zero-length input is permitted
        if inlen == 4:                  # a partial code is not okay
            raise ValueError("Invalid Data: Stream ended in the middle of a code")
        
        # Set up: get the first 9-bit code, which is the first decompressed byte,
        # but don't create a table entry until the next code
        buf = ba_in[3]
        buf += ba_in[4] << 8
        final = prev = buf & mask       # code
        buf >>= bits
        left = 16 - bits
        if prev > 255: 
            raise ValueError("Invalid Data: First code must be a literal")
        
        # We have output - allocate and set up an output buffer with first byte
        put = [final]
        
        # Decode codes
        mark = 3                        # start of compressed data
        nxt = 5                         # consumed five bytes so far
        while nxt < inlen:
            # If the table will be full after this, increment the code size
            if (end >= mask) and (bits < max_):
                # Flush unused input bits and bytes to next 8*bits bit boundary
                # (this is a vestigial aspect of the compressed data format
                # derived from an implementation that made use of a special VAX
                # machine instruction!)
                rem = (nxt - mark) % bits
                
                if (rem):
                    rem = bits - rem
                    if rem >= inlen - nxt: 
                        break
                    nxt += rem
                
                buf = 0
                left = 0
                
                # mark this new location for computing the next flush
                mark = nxt
                
                # increment the number of bits per symbol
                bits += 1
                mask <<= 1
                mask += 1
            
            # Get a code of bits bits
            buf += ba_in[nxt] << left
            nxt += 1
            left += 8
            if left < bits: 
                if nxt == inlen:
                    raise ValueError("Invalid Data: Stream ended in the middle of a code")
                buf += ba_in[nxt] << left
                nxt += 1
                left += 8
            code = buf & mask
            buf >>= bits
            left -= bits
            
            # process clear code (256)
            if (code == 256) and flags:
                # Flush unused input bits and bytes to next 8*bits bit boundary
                rem = (nxt - mark) % bits
                if rem:
                    rem = bits - rem
                    if rem > inlen - nxt:
                        break
                    nxt += rem
                buf = 0
                left = 0
                
                # Mark this location for computing the next flush
                mark = nxt
                
                # Go back to nine bits per symbol
                bits = 9                    # initialize bits and mask
                mask = 0x1ff
                end = 255                   # empty table
                continue                    # get next code
            
            # Process LZW code
            temp = code                     # save the current code
            stack = []                      # buffer for reversed match - empty stack
            
            # Special code to reuse last match
            if code > end:
                # Be picky on the allowed code here, and make sure that the
                # code we drop through (prev) will be a valid index so that
                # random input does not cause an exception
                if (code != end + 1) or (prev > end):
                    raise ValueError("Invalid Data: Invalid code detected")
                stack.append(final)
                code = prev

            # Walk through linked list to generate output in reverse order
            while code >= 256:
                stack.append(suffix[code])
                code = prefix[code]

            stack.append(code)
            final = code
            
            # Link new table entry
            if end < mask:
                end += 1
                prefix[end] = prev
                suffix[end] = final
            
            # Set previous code for next iteration
            prev = temp
            
            # Write stack to output in forward order
            put += stack[::-1]

        # Return the decompressed data as string
        return bytes( bytearray( put ) )
        
    def Get_Latlons_Distance( self, latlon1, latlon2 ) :
        dLat, dLon = latlon2[0]-latlon1[0], latlon2[1]-latlon1[1]
        a = sin( dLat/2 )**2 + cos( latlon1[0] ) * cos( latlon2[0] ) * sin( dLon/2 )**2
        c = 2 * atan2( sqrt(a), sqrt(1 - a) )
        
        return 6373.0 * c
    
    """
    # #########################################################################################
    # Daily GPS Broadcast Ephemeris Files
    # The daily GPS broadcast ephemeris file is a merge of the individual site navigation files
    # into one, non-redundant file that can be utilized by users instead of the many individual 
    # navigation files.
    #
    # The daily file created at BKG each day contains unique navigation messages from sites in 
    # Europe.
    #
    # The starting directory for the daily files is
    #
    # ftp://cddis.nasa.gov/gnss/data/daily/.
    #
    # Append the following directory and file names to the starting directory:
    #
    # YYYY/DDD/YYn/brdcDDD0.YYn.Z   (merged GPS broadcast ephemeris file)
    #
    # OR
    #
    # YYYY/brdc/brdcDDD0.YYn.Z   (merged GPS broadcast ephemeris file)
    #
    # YYYY/DDD/YYn/ifagDDD0.YYn.Z   (daily file historically created at BKG)
    #==========================================================================================
    # Hourly GPS Broadcast Ephemeris Files
    #
    # The combined broadcast ephemeris file is generated on an hourly basis from all hourly 
    # navigation files archived at the CDDIS. The hourly navigation file contains all broadcast 
    # messages with the TOE of the day that are available when the file is created at the top 
    # of the hour. The file is updated each hour with new navigation messages.
    #
    # At the end of the UTC day, when the final version of the file is generated, the file 
    # is copied to the daily directories and becomes the "daily" broadcast ephemeris file
    #
    # The starting directory for the hourly files is
    #
    # ftp://cddis.nasa.gov/gnss/data/hourly/
    #
    # Append the following directory and file names to the starting directory:
    # YYYY/DDD/hourDDDm.YYn.Z 
    """
    def Get_GPS_Ephemeris_From_NASA( self, forceUpdate = False, isHourly = True , host = "CDDIS" ):

        # try:
        #     urllib.request.urlopen('http://cddis.nasa.gov', timeout=1)
        # except urllib.request.URLError as err: 
        #     print( err )
        
        localDir = self.localDir + 'stations/' +  host + "/"
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        # Gets the day of year
        DoY = time.localtime().tm_yday
        fourdyear   = time.localtime().tm_year
        if isHourly:
            rinexFilename = "hour%03d0.%dn.Z" % ( DoY, fourdyear % 100 )
        else:
            rinexFilename = "brdc%03d0.%dn.Z" % ( DoY, fourdyear % 100 )
        localRinexFilename = localDir + rinexFilename
        unzippedRinexFilename = rinexFilename + '.rinex'
        localUnzippedRinexFilename = localDir + unzippedRinexFilename

        if os.path.exists( localRinexFilename )  and ( os.stat( localRinexFilename ).st_size <= 1024 ) :
            os.remove( localRinexFilename )
        elif os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size > 1024 ) and ( False == forceUpdate ):
            if os.path.exists( localUnzippedRinexFilename ) and ( os.stat( localUnzippedRinexFilename ).st_size < 1024 ) :
                os.remove( localUnzippedRinexFilename )
            elif not os.path.exists( localUnzippedRinexFilename ):
                rinexFObj = open( localUnzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
                if self.DEBUG:
                    print( "\tA local copy of %s was found and we unzipped it." % rinexFilename )
        elif ( not os.path.exists( localRinexFilename ) ) or ( True == forceUpdate ) :

            ### Generate URL from ftp://cddis.gsfc.nasa.gov/gnss/data/
            if isHourly:
                remoteDir = "ftp://cddis.gsfc.nasa.gov/gnss/data/hourly/%4d/%03d/" % ( fourdyear, DoY )
            else:
                remoteDir = "ftp://cddis.nasa.gov/gnss/data/daily/%4d/%03d/%2dn/" % ( fourdyear, DoY, fourdyear % 100 )

            parsedRemoteRinexFilename = urllib.parse.urlparse( remoteDir + rinexFilename )
            #print ("URL Generated " + pathfileurl + ephemerisfilename + '\n')

            ###Download File from URL, ftp://cddis.gsfc.nasa.gov/gnss/data/
            try:
                #print ( "Connecting to FTP " + parsedurl.netloc )
                ftp = ftplib.FTP( parsedRemoteRinexFilename.netloc )
                #print ( "FTP Login..." )
                ftp.login()
                #print ( "FTP Changing directory..." )
                ftp.cwd( os.path.dirname( parsedRemoteRinexFilename.path )[1:] )
                #ftp.dir()
                #print ("Downloading hourly ephemeris data to " + hourlyEphemerisfilename)
                ftp.retrbinary('RETR %s' % rinexFilename, open( localRinexFilename, 'wb' ).write)
                #print ("Disconnecting from " + parsedurl.netloc)
                ftp.quit()
            except ftplib.all_errors as err:
                print (err)

            if( ( os.path.exists( localRinexFilename ) ) or ( os.stat( localRinexFilename ).st_size > 0 ) ) :
                rinexFObj = open( localUnzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()

        with open( localUnzippedRinexFilename, 'rt' ) as rinexFObj :
            while True :
                line = rinexFObj.readline()
                if line[60:73] == "END OF HEADER" :
                    break
                if line[60:69] == "ION ALPHA" :
                    self.alpha0 = float( line[2:14].replace('D', 'E') )
                    self.alpha1 = float( line[14:26].replace('D', 'E') )
                    self.alpha2 = float( line[26:38].replace('D', 'E') )
                    self.alpha3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "ION BETA" :
                    self.beta0 = float( line[2:14].replace('D', 'E') )
                    self.beta1 = float( line[14:26].replace('D', 'E') )
                    self.beta2 = float( line[26:38].replace('D', 'E') )
                    self.beta3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "DELTA-UTC" :
                    self.A0 = float( line[3:22].replace('D', 'E') )
                    self.A1 = float( line[22:41].replace('D', 'E') )
                    self.tot = int( line[41:50] )
                    self.wnt = int( line[50:59] )
                elif line[60:72] == "LEAP SECONDS" :
                    self.dtls = int( line[0:6].replace('D', 'E') )

            line = rinexFObj.readline( )
            while line:
                ii = int( line[0:2] ) - 1
                self.svid[ii] = ii + 1
                self.YY[ii] = int( line[3:5] ) + 2000
                self.MM[ii] = int( line[6:8] )
                self.DD[ii] = int( line[9:11] )
                self.hh[ii] = int( line[12:14] )
                self.mm[ii] = int( line[15:17] )
                self.sec[ii] = float( line[18:22] )
                self.af0[ii] = float( line[22:41].replace('D', 'E')  )
                self.af1[ii] = float( line[41:60].replace('D', 'E')  )
                self.af2[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.iode[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Crs[ii] = float( line[22:41].replace('D', 'E')  )
                self.deltan[ii] = float( line[41:60].replace('D', 'E')  )
                self.M0[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.Cuc[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.e[ii] = float( line[22:41].replace('D', 'E')  )
                self.Cus[ii] = float( line[41:60].replace('D', 'E')  )
                self.sqrtA[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.toe[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Cic[ii] = float( line[22:41].replace('D', 'E')  )
                self.Omega0[ii] = float( line[41:60].replace('D', 'E')  )
                self.Cis[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.i0[ii] = float( line[3:22].replace('D', 'E')  )
                self.Crc[ii] = float( line[22:41].replace('D', 'E')  )
                self.omega[ii] = float( line[41:60].replace('D', 'E')  )
                self.Omegadot[ii] = float( line[60:79].replace('D', 'E')  )

                line = rinexFObj.readline( )
                self.idot[ii] = float( line[3:22].replace('D', 'E')  )
                self.l2code[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                t_week = int( float( line[41:60].replace('D', 'E')  ) )
                if ( t_week != self.week[ii] ) and ( self.week[ii] != 0 ) :
                    print( "t_week = %d and self.week = %d " % ( t_week, self.week[ii] ) )
                self.week[ii] = t_week

                line = rinexFObj.readline( )
                self.health[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.health[ii] = ( self.health[ii] % 32 + 32 ) if ( self.health[ii] < 32 ) else self.health[ii] 
                self.tgd[ii] = float( line[41:60].replace('D', 'E')  ) 
                self.iodc[ii] = int( float( line[60:79].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                line = rinexFObj.readline( )

        rinexFObj.close()
        
        return ii, unzippedRinexFilename

    #############################################################################################
    # broadcast nav files for sites are available from the UNAVCO ftp site in this directory:
    # ftp://data-out.unavco.org/pub/rinex/nav/
    # 
    # Files are organized by year, then day of year, then station. For example:
    # ftp://data-out.unavco.org/pub/rinex/nav/2019/070/azu10700.19n.Z
    ###############################################################################################
    def Get_GPS_Ephemeris_From_UNAVCO( self, forceUpdate = False, host = "azu1" ):
        localDir = self.localDir + 'stations/' + host + "/"
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        #Get the Ephemeris filename.
        DoY = time.localtime().tm_yday
        fourdyear   = time.localtime().tm_year
        ephemerisFilename = "%s%03d0.%dn.Z" % ( host, DoY, fourdyear % 100 )
        localRinexFilename = localDir + ephemerisFilename
        rinexFilename = localRinexFilename + '.rinex'

        if ( not os.path.exists( rinexFilename ) ) or ( os.stat( rinexFilename ).st_size == 0 ) or ( forceUpdate == True ):
            if os.path.exists( rinexFilename ) :
                os.remove( rinexFilename )

            remoteDir = "ftp://anonymous:anonymous@data-out.unavco.org/pub/rinex/nav/%4d/%03d/" % ( fourdyear, DoY )
 
            try:
                # Download Ephemeris file from URL, ftp://data-out.unavco.org/pub/rinex/nav/
                # for example, 
                # urllib.request.urlretrieve("ftp://data-out.unavco.org/pub/rinex/nav/2019/070/azu10700.19n.Z", localDir + "azu10700.19n.Z")
                urllib.request.urlretrieve( remoteDir + ephemerisFilename, localRinexFilename )
            except urllib.error.URLError as err:
                print (err)

            rinexFObj = open( rinexFilename, 'wb')
            rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
            rinexFObj.close()

        with open( rinexFilename, 'rt' ) as rinexFObj :
            while True :
                line = rinexFObj.readline()
                if line[60:73] == "END OF HEADER" :
                    break
                if line[60:69] == "ION ALPHA" :
                    self.alpha0 = float( line[2:14].replace('D', 'E') )
                    self.alpha1 = float( line[14:26].replace('D', 'E') )
                    self.alpha2 = float( line[26:38].replace('D', 'E') )
                    self.alpha3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "ION BETA" :
                    self.beta0 = float( line[2:14].replace('D', 'E') )
                    self.beta1 = float( line[14:26].replace('D', 'E') )
                    self.beta2 = float( line[26:38].replace('D', 'E') )
                    self.beta3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "DELTA-UTC" :
                    self.A0 = float( line[3:22].replace('D', 'E') )
                    self.A1 = float( line[22:41].replace('D', 'E') )
                    self.tot = int( line[41:50] )
                    self.wnt = int( line[50:59] )
                elif line[60:72] == "LEAP SECONDS" :
                    self.dtls = int( line[0:6].replace('D', 'E') )

            line = rinexFObj.readline( )
            while line:
                ii = int( line[0:2] ) - 1
                self.svid[ii] = ii + 1
                self.YY[ii] = int( line[3:5] ) + 2000
                self.MM[ii] = int( line[6:8] )
                self.DD[ii] = int( line[9:11] )
                self.hh[ii] = int( line[12:14] )
                self.mm[ii] = int( line[15:17] )
                self.sec[ii] = float( line[18:22] )
                self.af0[ii] = float( line[22:41].replace('D', 'E')  )
                self.af1[ii] = float( line[41:60].replace('D', 'E')  )
                self.af2[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.iode[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Crs[ii] = float( line[22:41].replace('D', 'E')  )
                self.deltan[ii] = float( line[41:60].replace('D', 'E')  )
                self.M0[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.Cuc[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.e[ii] = float( line[22:41].replace('D', 'E')  )
                self.Cus[ii] = float( line[41:60].replace('D', 'E')  )
                self.sqrtA[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                t_toe = int( float( line[3:22].replace('D', 'E')  ) )
                if t_toe < self.toe[ii] :
                    print ( "Ephemeris.Get_GPS_Ephemeris_From_UNAVCO():: The new Toe %d for Sat %02d is older than the current Toe %d." % (t_toe, self.svid[ii], self.toe[ii] ) )
                self.toe[ii] = t_toe
                self.Cic[ii] = float( line[22:41].replace('D', 'E')  )
                self.Omega0[ii] = float( line[41:60].replace('D', 'E')  )
                self.Cis[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.i0[ii] = float( line[3:22].replace('D', 'E')  )
                self.Crc[ii] = float( line[22:41].replace('D', 'E')  )
                self.omega[ii] = float( line[41:60].replace('D', 'E')  )
                self.Omegadot[ii] = float( line[60:79].replace('D', 'E')  )

                line = rinexFObj.readline( )
                self.idot[ii] = float( line[3:22].replace('D', 'E')  )
                self.l2code[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.week[ii] = int( float( line[41:60].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                self.health[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.health[ii] = ( self.health[ii] % 32 + 32 ) if ( self.health[ii] < 32 ) else self.health[ii] 
                self.tgd[ii] = float( line[41:60].replace('D', 'E')  ) 
                self.iodc[ii] = int( float( line[60:79].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                line = rinexFObj.readline( )

        rinexFObj.close()
        
        return ii, rinexFilename

    # ################################################################################
    # Data Download Information, http://sopac-csrc.ucsd.edu/index.php/data-download/
    #
    # To access to the SOPAC/CSRC GARNER archive, click on the links below. Depending 
    # on your browser, a dialog box may appear asking for a username and password. 
    # Enter the username "anonymous" and your email address as the password.
    #
    # Access the archive through ftp, ftp://garner.ucsd.edu/pub/
    #
    # Access the archive through http, http://garner.ucsd.edu/pub/
    #
    # You can also access the archive through the GPS Seamless Archive (GSAC).
    # http://geogsac.ucsd.edu:8080/gsacws/gsacapi/site/form
    #
    # You can include anonymous authentication in URL access to GARNER:
    # Format http://username:password@garner.ucsd.edu/
    # To login as an anonymous user, please use your email address as the password.
    # To enter your email address as the password, replace the '@' symbol with '%40'.
    # Example http://anonymous:user%40host.edu@garner.ucsd.edu/pub$ wget
    #
    # http://anonymous:jason%40ucsd.edu@garner.ucsd.edu/pub/rinex/2002/001/apex0010.02d.Z
    #######################################################################################
    def Get_GPS_Ephemeris_From_UCSD( self, host = "Garner", forceUpdate = False ):
        localDir = self.localDir + 'stations/' + host + "/"
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        #Gets the day of year
        DoY = time.localtime().tm_yday
        fourdyear   = time.localtime().tm_year
        rinexFilename = "p472%03d0.%dn.Z" % ( DoY, fourdyear % 100 )
        localRinexFilename = localDir + rinexFilename
        unzippedRinexFilename = rinexFilename + '.rinex'
        localUnzippedRinexFilename = localDir + unzippedRinexFilename

        if os.path.exists( localRinexFilename )  and ( os.stat( localRinexFilename ).st_size <= 1024 ) :
            os.remove( localRinexFilename )
        elif os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size > 1024 ) and ( False == forceUpdate ):
            if os.path.exists( localUnzippedRinexFilename ) and ( os.stat( localUnzippedRinexFilename ).st_size < 1024 ) :
                os.remove( localUnzippedRinexFilename )
            elif not os.path.exists( localUnzippedRinexFilename ):
                rinexFObj = open( localUnzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
                if self.DEBUG:
                    print( "\tA local copy of %s was found and we unzipped it." % rinexFilename )
        elif ( not os.path.exists( localRinexFilename ) ) or ( True == forceUpdate ) :

            ### Generate URL from ftp://garner.ucsd.edu/pub/
            remoteDir = "ftp://garner.ucsd.edu/pub/nav/%4d/%03d/" % ( fourdyear, DoY )
            parsedRemoteRinexFilename = urllib.parse.urlparse( remoteDir + rinexFilename )
            
            ###Download File from URL, ftp://garner.ucsd.edu/pub/
            try:
                ftp = ftplib.FTP( parsedRemoteRinexFilename.netloc )
                ftp.login( "anonymous", "shuwang@sandiego.edu" )
                ftp.cwd( os.path.dirname( parsedRemoteRinexFilename.path )[1:] )
                #ftp.dir()
                ftp.retrbinary('RETR %s' % rinexFilename, open( localRinexFilename, 'wb' ).write)
                ftp.quit()
            except ftplib.all_errors as err:
                print (err)

            if ( os.path.exists( localRinexFilename ) ) and ( os.stat( localRinexFilename ).st_size > 0 ):
                rinexFObj = open( unzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
            else:
                unzippedRinexFilename = "auto%03d0.%dn" % ( DoY, fourdyear % 100 )
                localUnzippedRinexFilename = localDir + unzippedRinexFilename

                if os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size <= 1024 ) :
                    os.remove( localRinexFilename )

                if ( not os.path.exists( localRinexFilename ) ) or ( True == forceUpdate ):
                    remoteDir = "ftp://anonymous:shuwang@sandiego.edu@garner.ucsd.edu/pub/rinex/%4d/%03d/" % ( fourdyear, DoY )
                    try:
                        urllib.request.urlretrieve( remoteDir + unzippedRinexFilename, localUnzippedRinexFilename )                    
                    except urllib.error.URLError as err:
                        print (err)

        with open( localUnzippedRinexFilename, 'rt' ) as rinexFObj :
            while True :
                line = rinexFObj.readline()
                if line[60:73] == "END OF HEADER" :
                    break
                if line[60:69] == "ION ALPHA" :
                    self.alpha0 = float( line[2:14].replace('D', 'E') )
                    self.alpha1 = float( line[14:26].replace('D', 'E') )
                    self.alpha2 = float( line[26:38].replace('D', 'E') )
                    self.alpha3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "ION BETA" :
                    self.beta0 = float( line[2:14].replace('D', 'E') )
                    self.beta1 = float( line[14:26].replace('D', 'E') )
                    self.beta2 = float( line[26:38].replace('D', 'E') )
                    self.beta3 = float( line[38:50].replace('D', 'E') )
                elif line[60:69] == "DELTA-UTC" :
                    self.A0 = float( line[3:22].replace('D', 'E') )
                    self.A1 = float( line[22:41].replace('D', 'E') )
                    self.tot = int( line[41:50] )
                    self.wnt = int( line[50:59] )
                elif line[60:72] == "LEAP SECONDS" :
                    self.dtls = int( line[0:6].replace('D', 'E') )

            line = rinexFObj.readline( )
            while line:
                ii = int( line[0:2] ) - 1
                self.svid[ii] = ii + 1
                self.YY[ii] = int( line[3:5] ) + 2000
                self.MM[ii] = int( line[6:8] )
                self.DD[ii] = int( line[9:11] )
                self.hh[ii] = int( line[12:14] )
                self.mm[ii] = int( line[15:17] )
                self.sec[ii] = float( line[18:22] )
                self.af0[ii] = float( line[22:41].replace('D', 'E')  )
                self.af1[ii] = float( line[41:60].replace('D', 'E')  )
                self.af2[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.iode[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.Crs[ii] = float( line[22:41].replace('D', 'E')  )
                self.deltan[ii] = float( line[41:60].replace('D', 'E')  )
                self.M0[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.Cuc[ii] = int( float( line[3:22].replace('D', 'E')  ) )
                self.e[ii] = float( line[22:41].replace('D', 'E')  )
                self.Cus[ii] = float( line[41:60].replace('D', 'E')  )
                self.sqrtA[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                t_toe = int( float( line[3:22].replace('D', 'E')  ) )
                if t_toe < self.toe[ii] :
                    print ( "Ephemeris.Get_GPS_Ephemeris_From_UCSD():: The new Toe %d for Sat %02d is older than the current Toe %d." % (t_toe, self.svid[ii], self.toe[ii] ) )
                self.toe[ii] = t_toe
                self.Cic[ii] = float( line[22:41].replace('D', 'E')  )
                self.Omega0[ii] = float( line[41:60].replace('D', 'E')  )
                self.Cis[ii] = float( line[60:79].replace('D', 'E')  )
                
                line = rinexFObj.readline( )
                self.i0[ii] = float( line[3:22].replace('D', 'E')  )
                self.Crc[ii] = float( line[22:41].replace('D', 'E')  )
                self.omega[ii] = float( line[41:60].replace('D', 'E')  )
                self.Omegadot[ii] = float( line[60:79].replace('D', 'E')  )

                line = rinexFObj.readline( )
                self.idot[ii] = float( line[3:22].replace('D', 'E')  )
                self.l2code[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.week[ii] = int( float( line[41:60].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                self.health[ii] = int( float( line[22:41].replace('D', 'E')  ) )
                self.health[ii] = ( self.health[ii] % 32 + 32 ) if ( self.health[ii] < 32 ) else self.health[ii] 
                self.tgd[ii] = float( line[41:60].replace('D', 'E')  ) 
                self.iodc[ii] = int( float( line[60:79].replace('D', 'E')  ) )

                line = rinexFObj.readline( )
                line = rinexFObj.readline( )

        rinexFObj.close()
        
        return ii, unzippedRinexFilename


    def Search_GPS_Ephemeris( self, epoch, ephemerisFilename = None, svids = None, forceUpdate = False ):
        
        gpsDateTime = GnssUtil.Epoch2Gps( epoch )
        if None == ephemerisFilename:   
            ephemerisFilename = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( 
                gpsDateTime.tm_year, 
                gpsDateTime.tm_yday, 
                gpsDateTime.tm_hour, 
                gpsDateTime.tm_min 
                )

        if None == svids:
            svids = range(1, MAX_GPS_SATS , 1 )
    
        DoY         = gpsDateTime.tm_yday
        fourdyear   = gpsDateTime.tm_year
        rinexFilename           = "hour%03d0.%dn.Z" % ( DoY, fourdyear % 100 )
        remoteDir               = "ftp://cddis.gsfc.nasa.gov/gnss/data/hourly/%4d/%03d/" % ( fourdyear, DoY )
        unzippedRinexFilename   = self.Download_RinexFile( rinexFilename, remoteDir, type = 1, forceUpdate = True )
        if unzippedRinexFilename :
            self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
        else:
            if self.DEBUG:
                print( "Ephemeris.Get_Ephemeris(): ", remoteDir )

        rinexFilename   = "brdc%03d0.%dn.Z" % ( DoY, fourdyear % 100 )
        remoteDir       = "ftp://cddis.nasa.gov/gnss/data/daily/%4d/%03d/%2dn/" % ( fourdyear, DoY, fourdyear % 100 )
        unzippedRinexFilename = self.Download_RinexFile( rinexFilename, remoteDir, type = 1 )
        if unzippedRinexFilename :
            self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
        else:
            if self.DEBUG:
                print( "Ephemeris.Get_Ephemeris(): ", remoteDir + rinexFilename )

        host = "azu1"
        rinexFilename   = "%s%03d0.%dn.Z" % ( host, DoY, fourdyear % 100 )
        remoteDir       = "ftp://anonymous:anonymous@data-out.unavco.org/pub/rinex/nav/%4d/%03d/" % ( fourdyear, DoY )
        unzippedRinexFilename = self.Download_RinexFile( rinexFilename, remoteDir, type = 2 )
        if unzippedRinexFilename :
            self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
        else:
            if self.DEBUG:
                print( "Ephemeris.Get_Ephemeris(): ", remoteDir + rinexFilename )
            remoteDir       = "ftp://data-out.unavco.org/pub/rinex/nav/%4d/%03d/" % ( fourdyear, DoY )
            unzippedRinexFilename = self.Download_RinexFile( rinexFilename, remoteDir, username = "anonymous", password = "shuwang@sandiego.edu", type = 1 )
            if unzippedRinexFilename :
                self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
            else:
                if self.DEBUG:
                    print( "Ephemeris.Get_Ephemeris(): ", remoteDir + rinexFilename )

        host = "p472"
        rinexFilename   = "%s%03d0.%dn.Z" % ( host, DoY, fourdyear % 100 )
        remoteDir       = "ftp://garner.ucsd.edu/pub/nav/%4d/%03d/" % ( fourdyear, DoY )
        unzippedRinexFilename = self.Download_RinexFile( rinexFilename, remoteDir, username = "anonymous", password = "shuwang@sandiego.edu", type = 1 )
        if unzippedRinexFilename :
            self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
        else:
            if self.DEBUG:
                print( "Ephemeris.Get_Ephemeris(): ", remoteDir + rinexFilename )

        host = "auto"
        rinexFilename   = "%s%03d0.%dn" % ( host, DoY, fourdyear % 100 )
        remoteDir       = "ftp://anonymous:shuwang@sandiego.edu@garner.ucsd.edu/pub/rinex/%4d/%03d/" % ( fourdyear, DoY )
        unzippedRinexFilename = self.Download_RinexFile( rinexFilename, remoteDir, type = 2, forceUpdate = True)
        if unzippedRinexFilename :
            self.Update_GPS_Ephemeris( unzippedRinexFilename, ephemerisFilename, epoch )
        else:
            if self.DEBUG:
                print( "Ephemeris.Get_Ephemeris(): ", remoteDir + rinexFilename )

        return ephemerisFilename

    def Download_RinexFile( self, rinexFilename, remoteDir, username = None, password = None, type = 1, forceUpdate = False ):
        localDir = os.getcwd()
        if not localDir.endswith('/') :
            localDir  += '/'
        
        localRinexFilename = localDir + rinexFilename
        if os.path.exists( localRinexFilename )  and ( os.stat( localRinexFilename ).st_size <= 1024 ) :
            os.remove( localRinexFilename )
        elif os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size > 1024 ) and ( False == forceUpdate ):
            if localRinexFilename.endswith('.Z') or localRinexFilename.endswith('.z'):
                unzippedRinexFilename = rinexFilename + '.rinex'
                localUnzippedRinexFilename = localDir + unzippedRinexFilename
                if os.path.exists( localUnzippedRinexFilename )  and ( os.stat( localUnzippedRinexFilename ).st_size > 1024 ) :
                    return unzippedRinexFilename
                rinexFObj = open( localDir + unzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
                return unzippedRinexFilename
            else:
                return rinexFilename

        if type == 1:
            parsedRemoteRinexFilename = urllib.parse.urlparse( remoteDir + rinexFilename )
            try:
                ftp = ftplib.FTP( parsedRemoteRinexFilename.netloc )
                if username:
                    ftp.login( username, password )
                else:
                    ftp.login()
                ftp.cwd( os.path.dirname( parsedRemoteRinexFilename.path )[1:] )
                # if self.DEBUG:
                #     ftp.dir()
                ftp.retrbinary('RETR %s' % rinexFilename, open( localRinexFilename, 'wb' ).write)
                ftp.quit()
            except ftplib.all_errors as err:
                if os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size == 0 ) :
                    if self.DEBUG:
                        print("\tdownloaded a broken copy of %s" % rinexFilename )
                    os.remove( localRinexFilename )
                print (err)
                return None

        elif type == 2:
            try:
                urllib.request.urlretrieve( remoteDir + rinexFilename, localRinexFilename )
            except urllib.error.URLError as err:
                if os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size == 0 ) :
                    if self.DEBUG:
                        print("\tdownloaded a broken copy of %s" % rinexFilename )
                    os.remove( localRinexFilename )
                print (err)
                return None

        if localRinexFilename.endswith('.Z') or localRinexFilename.endswith('.z'):
            unzippedRinexFilename = rinexFilename + '.rinex'
            if( ( os.path.exists( localRinexFilename ) ) and ( os.stat( localRinexFilename ).st_size > 0 ) ) :
                rinexFObj = open( localDir + unzippedRinexFilename, 'wb')
                rinexFObj.write( GnssUtil.unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
        else:
            unzippedRinexFilename = rinexFilename

        return unzippedRinexFilename

    def Update_GPS_Ephemeris( self, unzippedRinexFilename, ephemerisFilename, epochs = None, svids = None ):
        if None == unzippedRinexFilename:
            print( "Update_GPS_Ephemeris(): unzippedRinexFilename cannot be None" )
            return None

        localDir = os.getcwd()
        if not localDir.endswith('/') :
            localDir  += '/'

        localRinexFile  = localDir + unzippedRinexFilename
        if ( not os.path.exists( localRinexFile ) ) or ( os.stat( localRinexFile ).st_size == 0 ) :
            print( "Update_GPS_Ephemeris(): %s cannot be empty or broken." % unzippedRinexFilename )
            return None

        if None == svids:
            svids = range( 1, MAX_GPS_SATS + 1, 1 )

        if None == epochs:
            if self.DEBUG:
                print( "Update_GPS_Ephemeris(): epochs = None")
            epochs = [ [0.0, 0.0] for ii in svids ]
            nowEpoch = GnssUtil.Get_NowEpoch()
            for ii in range( len(svids) ):
                if svids[ ii ] > 0:
                    epochs[ ii ] = nowEpoch
        elif len( epochs ) <= 2:
            epoch = epochs
            epochs = [ [0.0, 0.0] for ii in svids ]
            for ii in range( len(svids) ):
                if svids[ ii ] > 0:
                    epochs[ ii ] = epoch

        ll          = [ 0 for ii in svids ]
        depoch      = [ 14400 for ii in svids ]
        updatedEph  = []
        alpha0, alpha1, alpha2, alpha3 = 0, 0, 0, 0
        beta0, beta1, beta2, beta3 = 0, 0, 0, 0
        A0, A1, tot, wnt = 0, 0, 0, 0
        dtls = LEAPSECONDS_2017
        with open( localRinexFile, 'rb' ) as f:
            lines = f.readlines()
            for end_of_header_line_number in range( len( lines ) ) :
                if b'END OF HEADER' in lines[ end_of_header_line_number ]:
                    break
                elif b"ION ALPHA" in lines[ end_of_header_line_number ] :
                    alpha0 = float( lines[ end_of_header_line_number ][2:14].replace(b'D', b'E') )
                    alpha1 = float( lines[ end_of_header_line_number ][14:26].replace(b'D', b'E') )
                    alpha2 = float( lines[ end_of_header_line_number ][26:38].replace(b'D', b'E') )
                    alpha3 = float( lines[ end_of_header_line_number ][38:50].replace(b'D', b'E') )
                elif b"ION BETA" in lines[ end_of_header_line_number ]:
                    beta0 = float( lines[ end_of_header_line_number ][2:14].replace(b'D', b'E') )
                    beta1 = float( lines[ end_of_header_line_number ][14:26].replace(b'D', b'E') )
                    beta2 = float( lines[ end_of_header_line_number ][26:38].replace(b'D', b'E') )
                    beta3 = float( lines[ end_of_header_line_number ][38:50].replace(b'D', b'E') )
                elif b"DELTA-UTC" in lines[ end_of_header_line_number ] :
                    A0 = float( lines[ end_of_header_line_number ][3:22].replace(b'D', b'E') )
                    A1 = float( lines[ end_of_header_line_number ][22:41].replace(b'D', b'E') )
                    tot = int( lines[ end_of_header_line_number ][41:50] )
                    wnt = int( lines[ end_of_header_line_number ][50:59] )
                elif b"LEAP SECONDS" in lines[ end_of_header_line_number ] :
                    dtls = int( lines[ end_of_header_line_number ][0:6].replace(b'D', b'E') )

            for jj in range( end_of_header_line_number + 1, len( lines ) , 8 ):
                ii = int( lines[jj][0:2] ) - 1
                if 1 > svids[ ii ]:
                    continue 
                toe = GnssUtil.Gps2Epoch( "%4d-%02d-%02d %02d:%02d:%02d" % ( 
                    int( lines[jj][3:5] ) + 2000, 
                    int( lines[jj][6:8] ), 
                    int( lines[jj][9:11] ), 
                    int( lines[jj][12:14] ), 
                    int( lines[jj][15:17] ), 
                    float( lines[jj][18:22] ) ) )[1]
                depoch1 = abs( toe - epochs[ ii ][1] )
                if depoch1 < depoch[ ii ]:
                    depoch[ii]  = depoch1
                    ll[ii] = jj
            #print( ll, unzippedRinexFilename )

            localEphemerisFile = localDir + ephemerisFilename
            if os.path.exists( localEphemerisFile ) :
                with open( ephemerisFilename, 'rb' ) as f:
                    old_ephemeris = json.load( f )
                    mm = 0
                    for kk in range( len( old_ephemeris ) ):
                        eph     = old_ephemeris[kk]
                        ii      = eph['svid'] - 1
                        toe     = eph['toe']
                        depoch1 = abs( toe - epochs[ ii ][1] )
                        jj = ll[ ii ]
                        if ( depoch1 > depoch[ ii ] ) and ( jj > 0 ) :
                            updatedEph.append({
                                'alpha0': alpha0,
                                'alpha1': alpha1,
                                'alpha2': alpha2,
                                'alpha3': alpha3,

                                'beta0': beta0,
                                'beta1': beta1,
                                'beta2': beta2,
                                'beta3': beta3,

                                'A0': A0,
                                'A1': A1,
                                'tot': tot,
                                'wnt': wnt,

                                'dtls': dtls,

                                'svid': ii+1,
                                'YY': int( lines[jj][3:5] ) + 2000,
                                'MM': int( lines[jj][6:8] ),
                                'DD': int( lines[jj][9:11] ),
                                'hh': int( lines[jj][12:14] ),
                                'mm': int( lines[jj][15:17] ),
                                'sec': float( lines[jj][18:22] ),
                                'af0': float( lines[jj][22:41].replace(b'D', b'E')  ),
                                'af1': float( lines[jj][41:60].replace(b'D', b'E')  ),
                                'af2': float( lines[jj][60:79].replace(b'D', b'E')  ),

                                'iode': int( float( lines[jj+1][3:22].replace(b'D', b'E')  ) ),
                                'Crs': float( lines[jj+1][22:41].replace(b'D', b'E')  ),
                                'deltan': float( lines[jj+1][41:60].replace(b'D', b'E')  ),
                                'M0': float( lines[jj+1][60:79].replace(b'D', b'E')  ),

                                'Cuc': int( float( lines[jj+2][3:22].replace(b'D', b'E')  ) ),
                                'e': float( lines[jj+2][22:41].replace(b'D', b'E')  ),
                                'Cus': float( lines[jj+2][41:60].replace(b'D', b'E')  ),
                                'sqrtA': float( lines[jj+2][60:79].replace(b'D', b'E')  ),
                        
                                'toe': int( float( lines[jj+3][3:22].replace(b'D', b'E')  ) ),
                                'Cic': float( lines[jj+3][22:41].replace(b'D', b'E')  ),
                                'Omega0': float( lines[jj+3][41:60].replace(b'D', b'E')  ),
                                'Cis': float( lines[jj+3][60:79].replace(b'D', b'E')  ),

                                'i0': float( lines[jj+4][3:22].replace(b'D', b'E')  ),
                                'Crc': float( lines[jj+4][22:41].replace(b'D', b'E')  ),
                                'omega': float( lines[jj+4][41:60].replace(b'D', b'E')  ),
                                'Omegadot': float( lines[jj+4][60:79].replace(b'D', b'E')  ),

                                'idot': float( lines[jj+5][3:22].replace(b'D', b'E')  ),
                                'l2code': int( float( lines[jj+5][22:41].replace(b'D', b'E')  ) ),
                                'week': int( float( lines[jj+5][41:60].replace(b'D', b'E')  ) ),

                                'health': int( float( lines[jj+6][22:41].replace(b'D', b'E')  ) ),
                                'tgd': float( lines[jj+6][41:60].replace(b'D', b'E')  ),
                                'iodc': int( float( lines[jj+6][60:79].replace(b'D', b'E')  ) )
                            })
                            mm += 1
                            ll[ii] = -mm
                            
                        elif ( depoch1 <= depoch[ ii ] ) and ( jj > 0 ) :
                            updatedEph.append( eph )
                            mm += 1
                            ll[ii] = -mm
                            depoch[ ii ] = depoch1
                        elif ( depoch1 < depoch[ ii ] ) and ( jj < 0 ) :
                            updatedEph[ -jj ] = eph
                        elif jj == 0:
                            updatedEph.append( eph )
                            mm += 1
                            ll[ii] = -mm
                            depoch[ ii ] = depoch1

            for ii in range( len(ll) ):
                if ll[ii] > 0:
                    jj = ll[ ii ]
                    updatedEph.append({
                        'alpha0': alpha0,
                        'alpha1': alpha1,
                        'alpha2': alpha2,
                        'alpha3': alpha3,

                        'beta0': beta0,
                        'beta1': beta1,
                        'beta2': beta2,
                        'beta3': beta3,

                        'A0': A0,
                        'A1': A1,
                        'tot': tot,
                        'wnt': wnt,

                        'dtls': dtls,

                        'svid': ii+1,
                        'YY': int( lines[jj][3:5] ) + 2000,
                        'MM': int( lines[jj][6:8] ),
                        'DD': int( lines[jj][9:11] ),
                        'hh': int( lines[jj][12:14] ),
                        'mm': int( lines[jj][15:17] ),
                        'sec': float( lines[jj][18:22] ),
                        'af0': float( lines[jj][22:41].replace(b'D', b'E')  ),
                        'af1': float( lines[jj][41:60].replace(b'D', b'E')  ),
                        'af2': float( lines[jj][60:79].replace(b'D', b'E')  ),

                        'iode': int( float( lines[jj+1][3:22].replace(b'D', b'E')  ) ),
                        'Crs': float( lines[jj+1][22:41].replace(b'D', b'E')  ),
                        'deltan': float( lines[jj+1][41:60].replace(b'D', b'E')  ),
                        'M0': float( lines[jj+1][60:79].replace(b'D', b'E')  ),

                        'Cuc': int( float( lines[jj+2][3:22].replace(b'D', b'E')  ) ),
                        'e': float( lines[jj+2][22:41].replace(b'D', b'E')  ),
                        'Cus': float( lines[jj+2][41:60].replace(b'D', b'E')  ),
                        'sqrtA': float( lines[jj+2][60:79].replace(b'D', b'E')  ),
                
                        'toe': int( float( lines[jj+3][3:22].replace(b'D', b'E')  ) ),
                        'Cic': float( lines[jj+3][22:41].replace(b'D', b'E')  ),
                        'Omega0': float( lines[jj+3][41:60].replace(b'D', b'E')  ),
                        'Cis': float( lines[jj+3][60:79].replace(b'D', b'E')  ),

                        'i0': float( lines[jj+4][3:22].replace(b'D', b'E')  ),
                        'Crc': float( lines[jj+4][22:41].replace(b'D', b'E')  ),
                        'omega': float( lines[jj+4][41:60].replace(b'D', b'E')  ),
                        'Omegadot': float( lines[jj+4][60:79].replace(b'D', b'E')  ),

                        'idot': float( lines[jj+5][3:22].replace(b'D', b'E')  ),
                        'l2code': int( float( lines[jj+5][22:41].replace(b'D', b'E')  ) ),
                        'week': int( float( lines[jj+5][41:60].replace(b'D', b'E')  ) ),

                        'health': int( float( lines[jj+6][22:41].replace(b'D', b'E')  ) ),
                        'tgd': float( lines[jj+6][41:60].replace(b'D', b'E')  ),
                        'iodc': int( float( lines[jj+6][60:79].replace(b'D', b'E')  ) )
                    })
        #print( ll, unzippedRinexFilename )
        with open( ephemerisFilename, 'wt+' ) as f:
            json.dump( updatedEph, f )
        
        return ephemerisFilename

    def Identify_GPS_Ephemeris( self, unzippedRinexFilename, ephemerisFilename = None, svids = None, gpsDateTime = None ):
        
        if None == unzippedRinexFilename:
            print( "Identify_GPS_Ephemeris(): unzippedRinexFilename cannot be None" )
            return None

        localDir = self.localDir
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        localRinexFile = localDir + unzippedRinexFilename
        if ( not os.path.exists( localRinexFile ) ) or ( os.stat( localRinexFile ).st_size == 0 ) :
            print( "Identify_GPS_Ephemeris(): %s cannot be empty or broken." % unzippedRinexFilename )
            return None

        if None == ephemerisFilename:
            ephemerisFilename = unzippedRinexFilename + '.json'

        if None == svids:
            svids = range(1, MAX_GPS_SATS , 1 )
        
        epoch = [ [0.0, 0.0] for ii in range( len(svids) ) ]
        if None == gpsDateTime:
            if self.DEBUG:
                print( "Identify_GPS_Ephemeris(): gpsDateTime = None")
            now = GnssUtil.Get_NowEpoch()
            for ii in range( len(svids) ):
                if svids[ ii ] > 0:
                    epoch[ ii ] = now
        elif isinstance( gpsDateTime, list ):
            for ii in range( len(svids) ):
                if svids[ ii ] > 0:
                    epoch[ ii ] = GnssUtil.Gps2Epoch( gpsDateTime[ ii ]  ) 
        else:
            for ii in range( len(svids) ):
                if svids[ ii ] > 0:
                    epoch[ ii ] = GnssUtil.Gps2Epoch( gpsDateTime ) 

        ll      = [ 0.0 for ii in range( len(svids) ) ]
        depoch  = [ 43200 for ii in range( len(svids) ) ]

        with open( localRinexFile, 'rb' ) as f:
            lines = f.readlines()

        for kk in range( len( lines ) ) :
            if b'END OF HEADER' in lines[kk]:
                break

        for jj in range( kk+1, len( lines ), 8 ):
            ii      = int( lines[jj][0:2] ) - 1
            if 0 == svids[ ii ]:
                continue 
            YY      = int( lines[jj][3:5] ) + 2000
            MM      = int( lines[jj][6:8] )
            DD      = int( lines[jj][9:11] )
            hh      = int( lines[jj][12:14] )
            mm      = int( lines[jj][15:17] )
            sec     = float( lines[jj][18:22] )
            epoch1  = GnssUtil.Gps2Epoch( "%4d-%02d-%02d %02d:%02d:%02d" % ( YY, MM, DD, hh, mm, sec ) )
            depoch1 = abs( epoch1[1] - epoch[ ii ][1] )

            if depoch1 < depoch[ ii ]:
                depoch[ ii ] = depoch1
                ll[ii] = jj
            #     if self.DEBUG:
            #         print( "Identify_GPS_Ephemeris(): %s, svids = %02d, old depoch = %6d, new depoch = %6d, jj = %d " % ( unzippedRinexFilename, svids[ii], depoch[ ii ], depoch1, jj) )
            #     continue

            # if depoch1 >= depoch[ ii ]:
            #     if self.DEBUG:
            #         print( "Identify_GPS_Ephemeris(): %s, svids = %02d, old depoch = %6d, new depoch = %6d" % ( unzippedRinexFilename, svids[ii], depoch[ ii ], depoch1) )
            #     continue
                
        eph = []
        if kk == 7:
            for ii in range( len( ll ) ) :
                if ll[ ii ] :
                    jj = ll[ ii ]
                    eph.append({
                        'alpha0': float( lines[3][2:14].replace(b'D', b'E') ),
                        'alpha1': float( lines[3][14:26].replace(b'D', b'E') ),
                        'alpha2': float( lines[3][26:38].replace(b'D', b'E') ),
                        'alpha3': float( lines[3][38:50].replace(b'D', b'E') ),

                        'beta0': float( lines[4][2:14].replace(b'D', b'E') ),
                        'beta1': float( lines[4][14:26].replace(b'D', b'E') ),
                        'beta2': float( lines[4][26:38].replace(b'D', b'E') ),
                        'beta3': float( lines[4][38:50].replace(b'D', b'E') ),

                        'A0': float( lines[5][3:22].replace(b'D', b'E') ),
                        'A1': float( lines[5][22:41].replace(b'D', b'E') ),
                        'tot': int( lines[5][41:50] ),
                        'wnt': int( lines[5][50:59] ),

                        'LEAP SECONDS': int( lines[6][0:6].replace(b'D', b'E') ),

                        'svid': ii+1,
                        'YY': int( lines[jj][3:5] ) + 2000,
                        'MM': int( lines[jj][6:8] ),
                        'DD': int( lines[jj][9:11] ),
                        'hh': int( lines[jj][12:14] ),
                        'mm': int( lines[jj][15:17] ),
                        'sec': float( lines[jj][18:22] ),
                        'af0': float( lines[jj][22:41].replace(b'D', b'E')  ),
                        'af1': float( lines[jj][41:60].replace(b'D', b'E')  ),
                        'af2': float( lines[jj][60:79].replace(b'D', b'E')  ),

                        'iode': int( float( lines[jj+1][3:22].replace(b'D', b'E')  ) ),
                        'Crs': float( lines[jj+1][22:41].replace(b'D', b'E')  ),
                        'deltan': float( lines[jj+1][41:60].replace(b'D', b'E')  ),
                        'M0': float( lines[jj+1][60:79].replace(b'D', b'E')  ),

                        'Cuc': int( float( lines[jj+2][3:22].replace(b'D', b'E')  ) ),
                        'e': float( lines[jj+2][22:41].replace(b'D', b'E')  ),
                        'Cus': float( lines[jj+2][41:60].replace(b'D', b'E')  ),
                        'sqrtA': float( lines[jj+2][60:79].replace(b'D', b'E')  ),
                
                        'toe': int( float( lines[jj+3][3:22].replace(b'D', b'E')  ) ),
                        'Cic': float( lines[jj+3][22:41].replace(b'D', b'E')  ),
                        'Omega0': float( lines[jj+3][41:60].replace(b'D', b'E')  ),
                        'Cis': float( lines[jj+3][60:79].replace(b'D', b'E')  ),

                        'i0': float( lines[jj+4][3:22].replace(b'D', b'E')  ),
                        'Crc': float( lines[jj+4][22:41].replace(b'D', b'E')  ),
                        'omega': float( lines[jj+4][41:60].replace(b'D', b'E')  ),
                        'Omegadot': float( lines[jj+4][60:79].replace(b'D', b'E')  ),

                        'idot': float( lines[jj+5][3:22].replace(b'D', b'E')  ),
                        'l2code': int( float( lines[jj+5][22:41].replace(b'D', b'E')  ) ),
                        'week': int( float( lines[jj+5][41:60].replace(b'D', b'E')  ) ),

                        'health': int( float( lines[jj+6][22:41].replace(b'D', b'E')  ) ),
                        'tgd': float( lines[jj+6][41:60].replace(b'D', b'E')  ),
                        'iodc': int( float( lines[jj+6][60:79].replace(b'D', b'E')  ) )
                    })
        else:
            for ii in range( len( ll ) ) :
                if ll[ ii ] :
                    jj = ll[ ii ]
                    eph.append({
                        'svid': ii+1,
                        'YY': int( lines[jj][3:5] ) + 2000,
                        'MM': int( lines[jj][6:8] ),
                        'DD': int( lines[jj][9:11] ),
                        'hh': int( lines[jj][12:14] ),
                        'mm': int( lines[jj][15:17] ),
                        'sec': float( lines[jj][18:22] ),
                        'af0': float( lines[jj][22:41].replace(b'D', b'E')  ),
                        'af1': float( lines[jj][41:60].replace(b'D', b'E')  ),
                        'af2': float( lines[jj][60:79].replace(b'D', b'E')  ),

                        'iode': int( float( lines[jj+1][3:22].replace(b'D', b'E')  ) ),
                        'Crs': float( lines[jj+1][22:41].replace(b'D', b'E')  ),
                        'deltan': float( lines[jj+1][41:60].replace(b'D', b'E')  ),
                        'M0': float( lines[jj+1][60:79].replace(b'D', b'E')  ),

                        'Cuc': int( float( lines[jj+2][3:22].replace(b'D', b'E')  ) ),
                        'e': float( lines[jj+2][22:41].replace(b'D', b'E')  ),
                        'Cus': float( lines[jj+2][41:60].replace(b'D', b'E')  ),
                        'sqrtA': float( lines[jj+2][60:79].replace(b'D', b'E')  ),
                
                        'toe': int( float( lines[jj+3][3:22].replace(b'D', b'E')  ) ),
                        'Cic': float( lines[jj+3][22:41].replace(b'D', b'E')  ),
                        'Omega0': float( lines[jj+3][41:60].replace(b'D', b'E')  ),
                        'Cis': float( lines[jj+3][60:79].replace(b'D', b'E')  ),

                        'i0': float( lines[jj+4][3:22].replace(b'D', b'E')  ),
                        'Crc': float( lines[jj+4][22:41].replace(b'D', b'E')  ),
                        'omega': float( lines[jj+4][41:60].replace(b'D', b'E')  ),
                        'Omegadot': float( lines[jj+4][60:79].replace(b'D', b'E')  ),

                        'idot': float( lines[jj+5][3:22].replace(b'D', b'E')  ),
                        'l2code': int( float( lines[jj+5][22:41].replace(b'D', b'E')  ) ),
                        'week': int( float( lines[jj+5][41:60].replace(b'D', b'E')  ) ),

                        'health': int( float( lines[jj+6][22:41].replace(b'D', b'E')  ) ),
                        'tgd': float( lines[jj+6][41:60].replace(b'D', b'E')  ),
                        'iodc': int( float( lines[jj+6][60:79].replace(b'D', b'E')  ) )
                    })

        with open( ephemerisFilename, 'wt' ) as f:
            json.dump( eph, f )
        
        return ephemerisFilename

    def Old_Update_GPS_Ephemeris( self, searchedRinexFilename, ephemerisFilename = None, gpsDateTime = None, svid = None  ):
        
        if None == searchedRinexFilename:
            print( "The input parameter unzippedRinexFilename cannot be None" )
            return None

        localDir = self.localDir
        directory = os.path.dirname( localDir )        
        if not os.path.exists( directory ) :
            os.makedirs( directory )

        localRinexFile = localDir + searchedRinexFilename
        if ( not os.path.exists( localRinexFile ) ) or ( os.stat( localRinexFile ).st_size == 0 ) :
            print( "%s cannot be empty or broken." % localRinexFile  )
            return None
            
        if ( None == ephemerisFilename ) and ( None == gpsDateTime ) :
            now = GnssUtil.Epoch2Gps( GnssUtil.Get_NowEpoch() )
            ephemerisFilename = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( now.tm_year, now.tm_yday, now.tm_hour, now.tm_min )
        elif ( None == ephemerisFilename ) and ( None != gpsDateTime ):
            gpsDateTime_ = time.strptime( gpsDateTime, '%Y-%m-%d %H:%M:%S' )
            ephemerisFilename = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( 
                gpsDateTime_.tm_year, 
                gpsDateTime_.tm_yday, 
                gpsDateTime_.tm_hour, 
                gpsDateTime_.tm_min 
                )

        if None == svid:
            svid = range(1, MAX_GPS_SATS , 1 )
        
        epoch = [ [0.0, 0.0] for ii in range( len(svid) ) ]
        if None == gpsDateTime:
            now = GnssUtil.Get_NowEpoch()
            for ii in range( len(svid) ):
                if svid[ ii ] > 0:
                    epoch[ ii ] = now
        else:
            if isinstance( gpsDateTime, str ) :
                for ii in range( len(svid) ):
                    if svid[ ii ] > 0:
                        epoch[ ii ] = GnssUtil.Gps2Epoch( gpsDateTime  )                 
            else:
                for ii in range( len(svid) ):
                    if svid[ ii ] > 0:
                        epoch[ ii ] = GnssUtil.Gps2Epoch( gpsDateTime[ ii ]  ) 
        
        with open( localRinexFile, 'r' ) as f:
            rin = json.load( f )
        
        if os.path.exists( ephemerisFilename ) and ( os.stat( ephemerisFilename ).st_size > 1024 ):

            with open( ephemerisFilename, 'r' ) as f:
                eph = json.load( f )

            depoch1 = [ 14400 for ii in range( len(svid) ) ]
            depoch2 = [ 14400 for ii in range( len(svid) ) ]
            if isinstance( rin, list):
                for ii in range( len( svid ) ):
                    if 0 == svid[ ii ]:
                        continue
                    for jj in range( len( rin ) ):
                        if rin[jj]['svid'] == svid[ii] :
                            depoch1[svid[ii]-1] = abs( epoch[ svid[ii]-1 ][1] - int( rin[jj]['toe'] ) )
                            jj_ = jj
                    for kk in range( len(eph) ):
                        if eph[kk]['svid'] == svid[ii] :
                            depoch2[svid[ii]-1] = abs( epoch[ svid[ii]-1 ][1] - eph[kk]['toe'] )
                            kk_ = kk

                    if depoch1[svid[ii]-1] < depoch2[svid[ii]-1] :
                        eph[kk_] = rin[jj_]
                
                with open( ephemerisFilename, 'w' ) as f:
                    json.dump( eph, f )

            else:
                for ii in range( len( svid ) ):
                    if 0 == svid[ ii ]:
                        continue
                    for jj in range( len(rin['data']) ):
                        if rin['data'][jj]['svid'] == svid[ii] :
                            depoch1[svid[ii]-1] = abs( epoch[ svid[ii]-1 ][1] - rin['data'][jj]['toe'] )
                            jj_ = jj
                    for kk in range( len(eph) ):
                        if eph[kk]['svid'] == svid[ii] :
                            depoch2[svid[ii]-1] = abs( epoch[ svid[ii]-1 ][1] - eph[kk]['toe'] )
                            kk_ = kk

                    if depoch1[svid[ii]-1] < depoch2[svid[ii]-1] :
                        eph[kk_] = rin['data'][jj_]
                
                with open( ephemerisFilename, 'w' ) as f:
                    json.dump( eph, f )
        else:
            eph = []
            if isinstance( rin, list ):
                with open( ephemerisFilename, 'w+' ) as f:
                    json.dump( rin, f )
            else:
                with open( ephemerisFilename, 'w+' ) as f:
                    json.dump( rin['data'], f )

        return ephemerisFilename

    def UT( self ) :
        print ( "Ephemeris.UT():: start unit testing .... " )
        for ii in self.svid :
            if self.svid[ii] == 0 :
                continue
            if self.week[ii] < 0 :
                print ( '\tEphemeris.UT():: svid: %d, week: %d' % ( self.svid[ii], self.week[ii] ) )            
            if self.e[ii] < 0.0034 or 0.058 < self.e[ii] :
                print ( '\tEphemeris.UT():: ID: %d, Eccentricity: %f' % ( self.svid[ii], self.e[ii] ) )
            if self.i0[ii] < (50*PI_DEG) or (60*PI_DEG) < self.i0[ii] :
                print ( '\tEphemeris.UT():: ID: %d, Orbital Inclination(rad): %f' % ( self.svid[ii], self.i0[ii] ) )
        return

    def UT_Get_XyzPosition( self, lalh_rad, epoch ) :
        print ( "Ephemeris.UT_Get_XyzPosition():: start unit testing .... " )
        rcvXyz = GnssUtil.Lalh2Xyz( lalh_rad )
        self.Get_Ephemeris_From_NASA_Rinex()
        self.Get_Ephemeris_From_Ntrip( mountpoint = "JOG20" ) ### 'OSLS0', 'AZU10', 'GMSD0', 'KERG1', 'GRAS0'
        
        assist, satXyzvt = self.Get_PositionAssist( lalh_rad, epoch, propagationDelayCompensation = False )
        visibles = len( assist )
        visibleSats, pseudorange = [ 0.0 for ii in range( visibles ) ], [ 0.0 for ii in range( visibles ) ]
        for ii in range( visibles ) :
            visibleSats[ii] = assist[ii][0]
            pseudorange[ii] = assist[ii][3]
        print ( "\tGet_PositionAssist():: The %d visibles are " % visibles, visibleSats )
        
        eRcvXyzt = self.Get_Receiver_XyzPosition( epoch, visibleSats, pseudorange, rcvXyz, propagationDelayCompensation = False )
        print ( "\tGet_XyzPosition():: The original XYZ position is", rcvXyz )
        print ( "\tGet_XyzPosition():: The XYZ position estimate is", eRcvXyzt )

        return eRcvXyzt