'''
# def Get_NowEpoch( leapSeconds = 18 )

# def Utc2Epoch( utcDateTime = '2019-03-24 17:00:00 UTC', leapSeconds = 18 ) 

# def Download_Broadcast_RinexFiles( epoch, ephemeris_Server_Configuration = "Ephemeris_Servers_Configuration.json", forceUpdate = False, __DEBUG = False )

# def Get_GPS_Receiver_XyzPosition( ephemerisFilename, epoch, visibleSats, codeShifts, rcvXyz, __DEBUG = True ) 

# def Xyz2AzimElev( satXyz, rcvXyz )

# def Ionosphere_Delay( azel, lalh, tow, alpha, beta )

# def Ionosphere_DelayEx( satXyz, rcvXyz, tow, alpha, beta )

# def UT_Ionosphere_DelayEx( ephemerisFilename, epoch, visibleSats, codeShifts, rcvXyz, __DEBUG = True )

# def Get_GPS_Positioning_Assistance( ephemerisFilename, rcvXyz, epoch = None, __DEBUG = False )

'''
import calendar, time, usb, serial, threading, os, base64, socket, ssl, sys, Ephemeris, json
from datetime import datetime
from math import sin, cos, atan2, sqrt, acos, asin
from sys import platform
from pynmea import nmea
from matplotlib import pyplot
#import urllib2

MAX_SATS        = 33
MAX_GPS_SATS    = 32
ICD200PI        = 3.1415926535898
PI_DEG          = ICD200PI/180
DAY_SECONDS     = 86400
WEEK_SECONDS    = 604800

C0          = 299792458  #GPS official speed of light.
L1_CODE_LEN = 1023
L1_FREQ     = 1575.42e6
L1_CHIP_HZ  = 1.023e6 #Hz C/A Chipping Frequency.
L1_CHIP_M   = 293.052256 # C0 / L1_CHIP_HZ
L1_CODE_M   = 299792.458 # L1_CODE_LEN * L1_CHIP_M
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
LEAPSECONDS_1980_2017 = 18

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
WGS84_e         = 0.081819190842622  # Eccentricity 0.0818191908426
WGS84_e2        = 0.081819190842622 * 0.081819190842622
WGS84_u_km3s2   = 398600.4418  # Gravitational Parameter ? = 398600.4418 km3/s2
WGS84_u_m3s2    = 3.986004418e+14  # Gravitational Parameter ? = 398600.4418 km3/s2
WGS84_f         = 1/298.257223563  # Flattening f = 1/298.257223563
WGS84_w_rads    = 7.292158553e-5  # Rotation rate of the Earth ? =7.292158553e-5 rad/s

F = -4.442807633e-10   # What is this paramter?  I need find it out.

PWR2_5        = pow( 2, -5 )
PWR2_19       = pow( 2, -19 )
PWR2_29       = pow( 2, -29 )
PWR2_31       = pow( 2, -31 )
PWR2_33       = pow( 2, -33 )
PWR2_43       = pow( 2, -43 )
PWR2_55       = pow( 2, -55 )
PWR2_31_PI    = PWR2_31 * ICD200PI
PWR2_43_PI    = PWR2_43 * ICD200PI

'''
# GPS Epoch is a continuous time system for all satellites and observation systems. 
# It is represented in seconds since Jan 5, 1980 (presumably when GNSS went online). 
# It is measured as the average of the atomic clocks for the observations and satellites
'''
def Get_NowEpoch( leapSeconds = 18 ) :
    #import calendar, time
    #LEAPSECONDS_1980_2017 = 18
    #WEEK_SECONDS = 604800

    epoch_second = calendar.timegm( time.gmtime() ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00 UTC', '%Y-%m-%d %H:%M:%S UTC' ) ) + leapSeconds
    epoch = [ epoch_second // 604800, epoch_second % 604800 ]
    return epoch

def Utc2Epoch( utcDateTime = '2019-03-24 17:00:00 UTC', leapSeconds = 18 ) :
    #import calendar, time
    ##LEAPSECONDS_1980_2017 = 18
    #WEEK_SECONDS = 604800

    epoch_second = calendar.timegm( time.strptime( utcDateTime, '%Y-%m-%d %H:%M:%S UTC' ) ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00 UTC', '%Y-%m-%d %H:%M:%S UTC' ) ) + leapSeconds
    epoch = [ epoch_second // WEEK_SECONDS, epoch_second % WEEK_SECONDS ]
    return epoch

def Gps2Epoch( gpsDateTime = '2019-03-24 17:00:00' ) :
    epoch_second = calendar.timegm( time.strptime( gpsDateTime, '%Y-%m-%d %H:%M:%S' ) ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00', '%Y-%m-%d %H:%M:%S' ) ) 
    epoch = [ epoch_second // WEEK_SECONDS, epoch_second % WEEK_SECONDS ]
    return epoch

def Epoch2Utc( epoch, leapSeconds = 18 ):
    utcDateTime     = time.gmtime( epoch[0] * WEEK_SECONDS + epoch[1] - leapSeconds + calendar.timegm( time.strptime( '1980-01-06 00:00:00 UTC', '%Y-%m-%d %H:%M:%S UTC' ) ) )
    return utcDateTime

def Epoch2Gps( epoch ):
    gpsDateTime     = time.gmtime( epoch[0] * WEEK_SECONDS + epoch[1] + calendar.timegm( time.strptime( '1980-01-06 00:00:00 UTC', '%Y-%m-%d %H:%M:%S UTC' ) ) )
    return gpsDateTime

def Kepler( M, e ) :
    EPS, E0 = 0.000001, M
    ecceAnomaly = M - ( M - e*sin(M) - M ) / ( 1 - e*cos(M) )
    while( abs(ecceAnomaly-E0) > EPS ):
        E0 = ecceAnomaly
        ecceAnomaly = ecceAnomaly - ( ecceAnomaly - e*sin(ecceAnomaly) - M ) / ( 1 - e*cos(ecceAnomaly) )
    return ecceAnomaly

def Dms2Deg( dms ) :
    return ( dms // 100 ) + ( dms % 100 ) / 60

def Dms2Rad( dms ) :
    return ( ( dms // 100 ) + ( dms % 100 ) / 60 ) * PI_DEG

def Lalh2Xyz( lalh, a = 6378137, ecc2 = 0.0818191908426*0.0818191908426 ) :
##############################################################################
# Convert the Lat-Lon-Height of the GPS signal receiver to ECEF WGS-84 X-Y-Z
# 
# http://www.oc.nps.edu/oc2902w/coord/
###############################################################################
    sinLat, cosLat, sinLon, cosLon = sin( lalh[0] ), cos( lalh[0] ), sin( lalh[1] ), cos( lalh[1] )
    N = a / sqrt( 1.0 - ecc2 * sinLat * sinLat  ) 
    N_h = N + lalh[2]
    return [ cosLat * cosLon * N_h, cosLat * sinLon * N_h, ( (1.0-ecc2)*N + lalh[2] ) * sinLat ]

def Xyz2Lal( xyz, ecc2 = 0.0818191908426*0.0818191908426 ) :
###################################################################################
# Convert ECEF WGS-84 X-Y-Z to the Lat-Lon-Height of the GPS signal receiver
# 
# http://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html
# http://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf
###################################################################################
    p2 = xyz[0]*xyz[0] + xyz[1]*xyz[1]
    p1 = sqrt( p2 )
    r2 = p2 + xyz[2]*xyz[2]
    r1 = sqrt( r2 )
    b2 = r2 * ( 1 - ecc2 )
    b1 = sqrt( b2 )
    theta = atan2( r1 * xyz[2], b1 * p1  )

    return [ atan2( xyz[2] + (r2-b2)/b2 * b1 * sin(theta)**3, p1 - ecc2 * r1 * cos(theta)**2 ), atan2( xyz[1], xyz[0] ) ]

def Xyz2Lalh( xyz, a = 6378137, ecc2 = 0.0818191908426*0.0818191908426 ) :
###################################################################################
# Convert ECEF WGS-84 X-Y-Z to the Lat-Lon-Height of the GPS signal receiver
# 
# http://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html
# http://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf
# http://clynchg3c.com/Technote/geodesy/coordcvt.pdf
# https://github.com/aewallin/ppp-tools/blob/master/ppp_common.py
###################################################################################
    lalh = [0.0, 0.0, 0.0]
    p2 = xyz[0]*xyz[0] + xyz[1]*xyz[1]
    p1 = sqrt( p2 )
    r2 = p2 + xyz[2]*xyz[2]
    r1 = sqrt( r2 )
    b2 = r2 * ( 1 - ecc2 )
    b1 = sqrt( b2 )

    #phi = atan2( p1, xyz[2] )
    phi = atan2( xyz[2], p1 * ( 1 - ecc2 ) )
    sinPhi = sin(phi)
    for ii in range( 1000 ) :
        RN  = a / sqrt( 1.0 - ecc2 * sinPhi * sinPhi )
        h   = p1/cos(phi) - RN
        phi = atan2( xyz[2], p1 * ( 1 - ecc2 * RN / ( RN + h ) ) )
        sinPhi = sin(phi)

    return  [ phi, atan2( xyz[1], xyz[0] ), h ]
    #return  [ phi, atan2( xyz[1], xyz[0] + p1 ) * 2, h ]

def Lallh2EnuRotationMatrix( lalh ):
    #############################################################################
    # Calculates the ENU rotation matrix for a given Lat-Lon Height in radians.
    # phi = Lallh2EnuRotMatrix( lalh ) = 3x3 rotation matrix.
    #
    # Use the following equations to achieve to rotation
    # - Vector Rotations (E) enu system (X) xyz sytem (Phi) rotation matrix
    # E = Phi*X
    # X = phi^t * E 
    # - Matrix Rotations (E) enu System (X) xyz system (Phi) rotation Matrix
    # E = Phi * X *  Phi^t
    # X = Phi^t * X * Phi
    ##############################################################################
    sinLat, cosLat, sinLon, cosLon = sin( lalh[0] ), cos( lalh[0] ), sin( lalh[1] ), cos( lalh[1] )  
    return [ [-sinLat*cosLon,  -sinLat*sinLon,   cosLat], [-sinLon,          cosLon,          0], [cosLat*cosLon,   cosLat*sinLon,   sinLat] ]

class BkgdColors :
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def pInverse3( A ) :
    '''##########################################################################
    # The inverse of a 3x3 matrix: https://www.dr-lex.be/random/matrix-inv.html
    # #############################################################################
    | a11 a12 a13 |-1
    | a21 a22 a23 |    =  1/DET * A
    | a31 a32 a33 |
    
    with A  =

    |  a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13 |
    |-(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13)|
    |  a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12 |

    and DET  =  a11(a33a22-a32a23) - a21(a33a12-a32a13) + a31(a23a12-a22a13)
    ############################################################################# '''

    rows = len( A )
    AA  = [ list([ 0.0 for col in range( 3 ) ]) for row in range( 3 ) ]
    A_  = [ list([ 0.0 for col in range( rows ) ]) for row in range( 3 ) ]
    for jj in range( 3 ) :
        for kk in range( 3 ) :
            for ll in range( rows ) :
                AA[jj][kk] += ( A[ll][jj] * A[ll][kk] ) 

    AA_ = [ [ AA[2][2]*AA[1][1]-AA[2][1]*AA[1][2], -AA[2][2]*AA[0][1]+AA[2][1]*AA[0][2], AA[1][2]*AA[0][1]-AA[1][1]*AA[0][2] ] ,
    [ -AA[2][2]*AA[1][0]+AA[2][0]*AA[1][2], AA[2][2]*AA[0][0]-AA[2][0]*AA[0][2], -AA[1][2]*AA[0][0]+AA[1][0]*AA[0][2] ] ,
    [ AA[2][1]*AA[1][0]-AA[2][0]*AA[1][1], -AA[2][1]*AA[0][0]+AA[2][0]*AA[0][1], AA[1][1]*AA[0][0]-AA[1][0]*AA[0][1] ] ]

    det = AA[0][0]*AA_[0][0] + AA[1][0]*AA_[0][1] + AA[2][0]*AA_[0][2]

    for kk in range( 3 ) :
        for ll in range( rows ) :
            for jj in range( 3 ) :
                A_[kk][ll] += ( AA_[kk][jj]*A[ll][jj] )
            A_[kk][ll] /= det

    return A_

def pInverse4( A ) :
    rows = len( A )
    AA  = [ list([ 0.0 for col in range( 4 ) ]) for row in range( 4 ) ]
    AA_ = [ list([ 0.0 for col in range( 4 ) ]) for row in range( 4 ) ]
    A_  = [ list([ 0.0 for col in range( rows ) ]) for row in range( 4 ) ]
    for jj in range( 4 ) :
        for kk in range( 4 ) :
            for ll in range( rows ) :
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
                
    if not detA :
        detA = 1
    for kk in range( 4 ) :
        for ll in range( rows ) :
            for jj in range( 4 ) :
                A_[kk][ll] += ( AA_[kk][jj]*A[ll][jj] )
            A_[kk][ll] /= detA

    return A_

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

def Generate_CRC24( msg ) :
    for msgByte in msg :
        crc = CRC24[ ( ( crc >> 16 ) ^ msgByte ) & 0xFF ] ^ ( crc << 8 )
    
    return crc & 0xFFFFFF

def Plot_Positions( lalPositions, htmlfile = 'SanDiegoMap.html' ) :
    zoom = 4
    coloricon = os.path.join( os.path.dirname(__file__), 'markers/%s.png' )

    f = open( htmlfile, 'w' )
    f.write( '<html>\n' )
    f.write( '<head>\n' )
    f.write( '<meta name="viewport" content="initial-scale=1.0, user-scalable=no" />\n' )
    f.write( '<meta http-equiv="content-type" content="text/html; charset=UTF-8"/>\n' )
    f.write( '<title> Google Maps - pygmaps </title>\n' )
    f.write( '<script type="text/javascript" src="https://maps.googleapis.com/maps/api/js?libraries=visualization&sensor=true_or_false"></script>\n' )
    f.write( '<script type="text/javascript">\n' )
    f.write( '\tfunction initialize() {\n' )
    f.write( '\t\tvar centerlatlng = new google.maps.LatLng(%f, %f);\n' % ( float( lalPositions[0][0] ), float( lalPositions[0][1] ) ) )
    f.write( '\t\tvar myOptions = {\n')
    f.write( '\t\t\tzoom: %d,\n' % zoom )
    f.write( '\t\t\tcenter: centerlatlng,\n')
    f.write( '\t\t\tmapTypeId: google.maps.MapTypeId.ROADMAP\n')
    f.write( '\t\t};\n')
    f.write( '\t\tvar map = new google.maps.Map(document.getElementById("map_canvas"), myOptions);\n')
    f.write( '\n' )

    for ii in range( 1, len(lalPositions) ) :
        f.write( '\t\tvar latlng = new google.maps.LatLng(%f, %f);\n' % ( lalPositions[ii][0], lalPositions[ii][1] ) )
        f.write( '\t\tvar img = new google.maps.MarkerImage(\'%s\');\n' % ( coloricon % 'FF0000' ) )
        f.write( '\t\tvar marker = new google.maps.Marker({\n' )
        f.write( '\t\ttitle: "%s",\n' % 'No Implementation' )
        f.write( '\t\ticon: img,\n')
        f.write( '\t\tposition: latlng\n' )
        f.write( '\t\t});\n' )
        f.write( '\t\tmarker.setMap(map);\n' )
        f.write( '\n' )        

    f.write( '\t}\n' )
    f.write( '</script>\n' )
    f.write( '</head>\n' )
    f.write( '<body style="margin:0px; padding:0px;" onload="initialize()">\n' )
    f.write( '\t<div id="map_canvas" style="width: 100%; height: 100%;"></div>\n' )
    f.write( '</body>\n' )
    f.write( '</html>\n' )
    f.close()
    return

SIN_LUT_SIZE = 512
SIN_LUT_65535 = [0, 201, 402, 603, 804, 1005, 1206, 1407, 1608, 1809, 2010, 2211, 2412, 2613, 2814, 3015,
                3216, 3416, 3617, 3818, 4019, 4219, 4420, 4621, 4821, 5022, 5222, 5422, 5623, 5823, 6023, 6223,
                6424, 6624, 6824, 7024, 7223, 7423, 7623, 7823, 8022, 8222, 8421, 8620, 8820, 9019, 9218, 9417,
                9616, 9815, 10014, 10212, 10411, 10609, 10808, 11006, 11204, 11402, 11600, 11798, 11996, 12193, 12391, 12588,
                12785, 12982, 13179, 13376, 13573, 13770, 13966, 14163, 14359, 14555, 14751, 14947, 15142, 15338, 15533, 15729,
                15924, 16119, 16313, 16508, 16703, 16897, 17091, 17285, 17479, 17673, 17866, 18060, 18253, 18446, 18639, 18831,
                19024, 19216, 19408, 19600, 19792, 19984, 20175, 20366, 20557, 20748, 20939, 21129, 21319, 21509, 21699, 21889,
                22078, 22267, 22456, 22645, 22834, 23022, 23210, 23398, 23586, 23773, 23960, 24147, 24334, 24521, 24707, 24893,
                25079, 25265, 25450, 25635, 25820, 26005, 26189, 26374, 26557, 26741, 26925, 27108, 27291, 27473, 27656, 27838,
                28020, 28201, 28383, 28564, 28745, 28925, 29106, 29286, 29465, 29645, 29824, 30003, 30181, 30360, 30538, 30716,
                30893, 31070, 31247, 31424, 31600, 31776, 31952, 32127, 32302, 32477, 32651, 32826, 32999, 33173, 33346, 33519,
                33692, 33864, 34036, 34208, 34379, 34550, 34721, 34891, 35061, 35231, 35400, 35569, 35738, 35906, 36074, 36242,
                36409, 36576, 36743, 36909, 37075, 37241, 37406, 37571, 37736, 37900, 38064, 38227, 38390, 38553, 38715, 38877, 
                39039, 39200, 39361, 39522, 39682, 39842, 40001, 40161, 40319, 40478, 40635, 40793, 40950, 41107, 41263, 41419,
                41575, 41730, 41885, 42039, 42194, 42347, 42500, 42653, 42806, 42958, 43109, 43261, 43411, 43562, 43712, 43861,
                44011, 44159, 44308, 44456, 44603, 44750, 44897, 45043, 45189, 45334, 45479, 45624, 45768, 45912, 46055, 46198,
                46340, 46482, 46624, 46765, 46905, 47046, 47185, 47325, 47464, 47602, 47740, 47877, 48014, 48151, 48287, 48423,
                48558, 48693, 48827, 48961, 49095, 49228, 49360, 49492, 49624, 49755, 49885, 50016, 50145, 50274, 50403, 50531,
                50659, 50787, 50913, 51040, 51166, 51291, 51416, 51540, 51664, 51788, 51911, 52033, 52155, 52277, 52398, 52518,
                52638, 52758, 52877, 52995, 53113, 53231, 53348, 53464, 53580, 53696, 53811, 53925, 54039, 54153, 54266, 54378,
                54490, 54602, 54713, 54823, 54933, 55042, 55151, 55260, 55367, 55475, 55582, 55688, 55794, 55899, 56003, 56108,
                56211, 56314, 56417, 56519, 56620, 56721, 56822, 56922, 57021, 57120, 57218, 57316, 57413, 57510, 57606, 57702,
                57797, 57891, 57985, 58079, 58171, 58264, 58356, 58447, 58537, 58628, 58717, 58806, 58895, 58983, 59070, 59157,
                59243, 59329, 59414, 59498, 59582, 59666, 59749, 59831, 59913, 59994, 60075, 60155, 60234, 60313, 60391, 60469, 
                60546, 60623, 60699, 60775, 60850, 60924, 60998, 61071, 61144, 61216, 61287, 61358, 61429, 61498, 61567, 61636,
                61704, 61772, 61838, 61905, 61970, 62035, 62100, 62164, 62227, 62290, 62352, 62414, 62475, 62535, 62595, 62654,
                62713, 62771, 62829, 62886, 62942, 62998, 63053, 63107, 63161, 63214, 63267, 63319, 63371, 63422, 63472, 63522,
                63571, 63620, 63668, 63715, 63762, 63808, 63853, 63898, 63943, 63986, 64030, 64072, 64114, 64155, 64196, 64236,
                64276, 64315, 64353, 64391, 64428, 64464, 64500, 64535, 64570, 64604, 64638, 64671, 64703, 64734, 64765, 64796,
                64826, 64855, 64883, 64911, 64939, 64966, 64992, 65017, 65042, 65066, 65090, 65113, 65136, 65158, 65179, 65199,
                65219, 65239, 65258, 65276, 65293, 65310, 65327, 65342, 65357, 65372, 65386, 65399, 65412, 65424, 65435, 65446,
                65456, 65466, 65475, 65483, 65491, 65498, 65504, 65510, 65515, 65520, 65524, 65527, 65530, 65532, 65534, 65535]

def Get_SineCosine( theta ) :
    ## A Look-Up Table based simple sine & cosine function.
    lutSize2 = SIN_LUT_SIZE << 1 
    dPhi = ICD200PI / ( 2 * SIN_LUT_SIZE )

    ii = int( round( theta / dPhi ) ) % ( lutSize2 << 1 )
    sinTheta = ( -1.0 if ii > lutSize2  else +1.0 ) * SIN_LUT_65535[ ( SIN_LUT_SIZE - ii%SIN_LUT_SIZE ) if ( ii % lutSize2 ) > SIN_LUT_SIZE else ( ii % SIN_LUT_SIZE ) ]/65535.0

    ii = ( SIN_LUT_SIZE - ii ) % ( lutSize2 << 1 )
    cosTheta= ( -1.0 if ii > lutSize2  else +1.0 ) * SIN_LUT_65535[ ( SIN_LUT_SIZE - ii%SIN_LUT_SIZE ) if ( ii % lutSize2 ) > SIN_LUT_SIZE else ( ii % SIN_LUT_SIZE ) ]/65535.0

    return [ sinTheta, cosTheta ]

PRN_LFSR_TAPS = [ [ 1, 5 ], [ 2, 6 ], [ 3, 7 ], [ 4, 8 ], [ 0, 8 ], [ 1, 9 ], 
                [ 0, 7 ], [ 1, 8 ], [ 2, 9 ], [ 1, 2 ], [ 2, 3 ], [ 4, 5 ],
                [ 5, 6 ], [ 6, 7 ], [ 7, 8 ], [ 8, 9 ], [ 0, 3 ], [ 1, 4 ],
                [ 2, 5 ], [ 3, 6 ], [ 4, 7 ], [ 5, 8 ], [ 0, 2 ], [ 3, 5 ],
                [ 4, 6 ], [ 5, 7 ], [ 6, 8 ], [ 7, 9 ], [ 0, 5 ], [ 1, 6 ],
                [ 2, 7 ], [ 3, 8 ], [ 4, 9 ], [ 3, 9 ], [ 0, 6 ], [ 1, 7 ], [ 3, 9 ] ]

def Get_CACode( prn ) :
    # G1 LFSR: x^10+x^3+1
    s   = [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 ]
    n   = len( s )
    g1  = [ 1 for ii in range( n ) ]	#initialization vector for G1
    
    # G2j LFSR: x^10+x^9+x^8+x^6+x^3+x^2+1
    t   = [ 0, 1, 1, 0, 0, 1, 0, 1, 1, 1 ]
    q   = [ 1 for ii in range( n ) ]	#initialization vector for G2
    g2  = [ 1 for ii in range( L1_CODE_LEN ) ]
    
    # generate C/A Code sequences:
    tap = PRN_LFSR_TAPS[ prn ]
    code = [ 0 for ii in range( L1_CODE_LEN ) ]
    for ii in range( L1_CODE_LEN ) :
        g2[ii] = ( q[tap[0]] + q[tap[1]] ) % 2
        code[ii] = +1 if ( ( g1[n-1] + g2[ii] ) % 2 ) > 0 else -1
        t_g, t_q = 0, 0
        for jj in range( n ) :
            t_g += g1[jj] * s[jj]
            t_q += q[jj] * t[jj] 
        for jj in range( n-1, 0, -1 ) :
            g1[jj] = g1[jj-1]
            q[jj] = q[jj-1]
        g1[0] = t_g % 2
        q[0] = t_q % 2

    return code

def unlzw( zData ):
    # This function was adapted for Python from Mark Adler's C implementation
    # https://github.com/umeat/unlzw
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

def Get_Latlons_Distance( latlon1, latlon2 ) :
    dLat, dLon = latlon2[0]-latlon1[0], latlon2[1]-latlon1[1]
    a = sin( dLat/ 2 )**2  + cos( latlon1[0] ) * cos( latlon2[0] ) * ( sin( dLon / 2 )**2 )
    c = 2.0 * atan2( sqrt(a), sqrt(1 - a) )
    
    return 6371009 * c

def Download_Broadcast_RinexFiles( epoch, ephemeris_Server_Configuration = "Ephemeris_Servers_Configuration.json", forceUpdate = False, __DEBUG = False ):
    import urllib
    import ftplib
    WEEK_SECONDS        = 604800

    localDir = os.getcwd()
    if not localDir.endswith('/') :
        localDir  += '/'

    gpsDateTime     = time.gmtime( epoch[0] * WEEK_SECONDS + epoch[1] + calendar.timegm( time.strptime( '1980-01-06 00:00:00 UTC', '%Y-%m-%d %H:%M:%S UTC' ) ) )
    DoY             = gpsDateTime.tm_yday
    fourdyear       = gpsDateTime.tm_year
    
    with open( ephemeris_Server_Configuration, 'r' ) as f:
        server_Configuration = json.load( f )

    unzippedRinexFilenameList = []
    for config in server_Configuration:
        rinexFilename   = config["Rinex_Filename_Template"] % ( DoY, fourdyear % 100 )
        localRinexFilename = localDir + rinexFilename
        if os.path.exists( localRinexFilename )  and ( ( os.stat( localRinexFilename ).st_size <= 1024 )  or config["Forced_Update"] ) :
            os.remove( localRinexFilename )
        elif os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size > 1024 ) and ( False == forceUpdate ):
            if localRinexFilename.endswith('.Z') or localRinexFilename.endswith('.z'):
                unzippedRinexFilename = rinexFilename + '.rinex'
                localUnzippedRinexFilename = localDir + unzippedRinexFilename
                if os.path.exists( localUnzippedRinexFilename )  and ( os.stat( localUnzippedRinexFilename ).st_size > 1024 ) :
                    unzippedRinexFilenameList.append( unzippedRinexFilename )
                rinexFObj = open( localDir + unzippedRinexFilename, 'wb')
                rinexFObj.write( unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()
                unzippedRinexFilenameList.append( unzippedRinexFilename )
            else:
                unzippedRinexFilenameList.append( rinexFilename )
            continue

        if config["Download_Method"] == "FTP":
            remoteDir = config["Remote_Directory_Template"] 
            if 2 == remoteDir.count("/%"):
                remoteDir = remoteDir % ( fourdyear, DoY )
            elif 3 == remoteDir.count("/%"):
                remoteDir = remoteDir % ( fourdyear, DoY, fourdyear % 100 )
            parsedRemoteRinexFilename = urllib.parse.urlparse( remoteDir + rinexFilename )
            try:
                ftp = ftplib.FTP( parsedRemoteRinexFilename.netloc )
                username = config["Username"]
                if username:
                    ftp.login( username, config["Password"] )
                else:
                    ftp.login()
                ftp.cwd( os.path.dirname( parsedRemoteRinexFilename.path )[1:] )
                if __DEBUG:
                    ftp.dir()
                ftp.retrbinary('RETR %s' % rinexFilename, open( localRinexFilename, 'wb' ).write)
                ftp.quit()
            except ftplib.all_errors as err:
                if os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size == 0 ) :
                    if __DEBUG:
                        print("\tdownloaded a broken copy of %s" % rinexFilename )
                        print (err)
                    os.remove( localRinexFilename )
                continue

        elif config["Download_Method"] == "HTTP":
            try:
                remoteDir = config["Remote_Directory_Template"] 
                if 2 == remoteDir.count("/%"):
                    remoteDir = remoteDir % ( fourdyear, DoY )
                elif 3 == remoteDir.count("/%"):
                    remoteDir = remoteDir % ( fourdyear, DoY, fourdyear % 100 )

                urllib.request.urlretrieve( remoteDir + rinexFilename, localRinexFilename )
            except urllib.error.URLError as err:
                if os.path.exists( localRinexFilename ) and ( os.stat( localRinexFilename ).st_size == 0 ) :
                    if __DEBUG:
                        print("\tdownloaded a broken copy of %s" % rinexFilename )
                        print (err)
                    os.remove( localRinexFilename )
                continue

        if localRinexFilename.endswith('.Z') or localRinexFilename.endswith('.z'):
            unzippedRinexFilename = rinexFilename + '.rinex'
            if( ( os.path.exists( localRinexFilename ) ) and ( os.stat( localRinexFilename ).st_size > 0 ) ) :
                rinexFObj = open( localDir + unzippedRinexFilename, 'wb')
                rinexFObj.write( unlzw( open( localRinexFilename, 'rb' ).read() ) )
                rinexFObj.close()

            unzippedRinexFilenameList.append( unzippedRinexFilename )
        else:
            unzippedRinexFilenameList.append( rinexFilename )

    return unzippedRinexFilenameList


def Merge_GPS_Ephemeris( unzippedRinexFilenameList, ephemerisFilename = None, epochs = None, svids = None, __DEBUG = True ):
    LEAPSECONDS_2017 = 18

    localDir = os.getcwd()
    if not localDir.endswith('/') :
        localDir  += '/'

    gpsDateTime = time.gmtime()
    if None == ephemerisFilename:   
        ephemerisFilename = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( 
            gpsDateTime.tm_year, 
            gpsDateTime.tm_yday, 
            gpsDateTime.tm_hour if gpsDateTime.tm_min < 30 else ( gpsDateTime.tm_hour + 1) , 
            30 if gpsDateTime.tm_min < 30 else 0
            )
   

    if None == svids:
        svids = range( 1, MAX_GPS_SATS + 1, 1 )

    if None == epochs:
        epochs = [ [0.0, 0.0] for ii in svids ]
        nowEpoch = Get_NowEpoch()
        for ii in range( len(svids) ):
            if svids[ ii ] > 0:
                epochs[ ii ] = nowEpoch
    elif len( epochs ) <= 2:
        epoch = epochs
        epochs = [ [0.0, 0.0] for ii in svids ]
        for ii in range( len(svids) ):
            if svids[ ii ] > 0:
                epochs[ ii ] = epoch

    depoch      = [ 14400 for ii in svids ]
    updatedEph  = []
    alpha0, alpha1, alpha2, alpha3 = 0, 0, 0, 0
    beta0, beta1, beta2, beta3 = 0, 0, 0, 0
    A0, A1, tot, wnt = 0, 0, 0, 0
    dtls = LEAPSECONDS_2017
    for unzippedRinexFilename in unzippedRinexFilenameList:
        localRinexFile  = localDir + unzippedRinexFilename
        if ( not os.path.exists( localRinexFile ) ) or ( os.stat( localRinexFile ).st_size == 0 ) :
            continue
        ll = [ 0 for ii in svids ]
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
            
            for jj in range( end_of_header_line_number + 1, len( lines ) , 1 ):
                try: 
                    ii = int( lines[jj][0:2] ) - 1
                    if 1 > svids[ ii ]:
                        continue 
                except ValueError:
                    continue
                ######################################################################
                toe = Gps2Epoch( "%4d-%02d-%02d %02d:%02d:%02d" % ( 
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

#
# 03/27/2019  Shu Wang, shuwang1@outlook.com
#   Review Get_Receiver_XyzPosition() and make sure the calibration works and improved.
#
def Get_GPS_Receiver_XyzPosition( ephemerisFilename, epoch, visibleSats, codeShifts, rcvXyz, __DEBUG = True ) :
    from math import floor, sqrt
    MAX_GPS_SATS        = 32
    WEEK_SECONDS        = 604800
    HALF_WEEK_SECONDS   = 302400
    #WGS84_u_m3s2        = 3.986004418e+14  # Gravitational Parameter ? = 398600.4418 km3/s2
    SQRT_WGS84_u_m3s2   = 1.99649803857e+7  # Gravitational Parameter ? = 398600.4418 km3/s2
    WGS84_w_rads        = 7.292158553e-5  # Rotation rate of the Earth ? =7.292158553e-5 rad/s  
    L1_CHIP_M           = 293.052256109
    L1_CODE_M           = 299792.458
    HALF_L1_CODE_M      = 149896.229

    alpha   = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    beta    = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    svid    = [ 0 for ii in range(MAX_GPS_SATS) ]
    week    = [ 0 for ii in range(MAX_GPS_SATS) ]
    toe     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    sqrtA   = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    deltan  = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    M0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    e       = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    omega   = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cus     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cuc     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Crs     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Crc     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cis     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cic     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    i0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    idot    = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af0     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af1     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af2     = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    Omega0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Omegadot    = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    tgd     = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    DOP = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
    with open( ephemerisFilename, 'r' ) as f:
        ephData = json.load( f )
        for eph in ephData:

            ii = eph['svid'] - 1
            svid[ii]    = ii + 1

            if 'alpha0' in eph:
                alpha[ii][0]    = eph['alpha0']
                alpha[ii][1]    = eph['alpha1']
                alpha[ii][2]    = eph['alpha2']
                alpha[ii][3]    = eph['alpha3']

            if 'beta0' in eph:
                beta[ii][0]     = eph['beta0']
                beta[ii][1]     = eph['beta1']
                beta[ii][2]     = eph['beta2']
                beta[ii][3]     = eph['beta3']

            week[ii]    = eph['week']
            # if 'acc' in eph:
            #     self.acc[ii]          = eph['acc']
            # if 'l2code' in eph:
            #     self.l2code[ii]       = eph['l2code']
            idot[ii]       = eph['idot']
            # if 'iode' in eph:
            #      self.iode[ii]       = eph['iode']
            af0[ii]        = eph['af0']
            af1[ii]        = eph['af1']
            af2[ii]        = eph['af2']
            # if 'iodc' in eph:
            #      self.iodc[ii]       = eph['iodc']

            Crs[ii]        = eph['Crs']
            Crc[ii]        = eph['Crc']
            deltan[ii]     = eph['deltan']
            M0[ii]         = eph['M0']
            Cus[ii]        = eph['Cus']
            Cuc[ii]        = eph['Cuc']
            e[ii]          = eph['e']
            sqrtA[ii]      = eph['sqrtA']
            toe[ii]        = eph['toe']
            Cis[ii]        = eph['Cis']
            Cic[ii]        = eph['Cic']
            i0[ii]         = eph['i0']
            omega[ii]      = eph['omega']
            Omegadot[ii]   = eph['Omegadot']
            Omega0[ii]     = eph['Omega0']
            tgd[ii]        = eph['tgd']
            # if 'health' in eph:
            #     self.health[ii] = eph['health']
            # if 'toc' in eph:
            #     self.toc[ii]    = eph['toc']
            # if 'l2p' in eph:
            #     self.l2p[ii]    = eph['l2p']
            # if 'fit' in eph:
            #     self.fit[ii]    = eph['fit']

    lari            = [ 0.0 for col in range(6) ]
    xyOrb           = [ 0.0 for col in range(4) ]

    numVisibleSat   = len( visibleSats )
    satXyzvt        = [ list([ 0.0 for col in range(8) ]) for row in range( numVisibleSat ) ]
    tau             = [ 0.0 for row in range( numVisibleSat ) ]  
    rcvXyzt         = [ rcvXyz[0], rcvXyz[1], rcvXyz[2], 0 ]
    b0              = 0.0   #initial receiver xlock biase
    dRcvXyzt        = [ 1.0, 1.0, 1.0, 1.0 ] 

    watchDog1 = 0
    while ( watchDog1 < 10000 ) and ( (dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2]) > 0.01 ) and ( abs( dRcvXyzt[3] ) > 0.001 ) :
        for ii in range( numVisibleSat ) :
            sat = visibleSats[ii] - 1
            if svid[sat] == 0 :                
                del visibleSats[ii]
                del satXyzvt[ii][:]
                del tau[ii]
                numberVisibleSat -= 1
                if __DEBUG:
                    print("No ephemeris data for SV %d" % (sat+1))
                continue
            satXyzvt[ii][0] = sat
            tk = epoch[1] - tau[ii] - toe[sat] + ( epoch[0] - week[sat] ) * WEEK_SECONDS
            if tk > HALF_WEEK_SECONDS :
                tk -= WEEK_SECONDS
            if tk < -HALF_WEEK_SECONDS :
                tk += WEEK_SECONDS
            #############################################################################
            # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
            # The time of transmissionis is modified from TOW from navigation data and
            # the correction of is performed as
            #############################################################################
            A = sqrtA[sat] * sqrtA[sat]
            meanMotion = SQRT_WGS84_u_m3s2 / ( sqrtA[sat]*A ) + deltan[sat] # This value should be caliberated with adding Mean Motion Difference.
            meanAnomaly  = M0[sat] + meanMotion * tk

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
            ecceAnomaly = meanAnomaly + e[sat]*sin(meanAnomaly) / ( 1 - e[sat]*cos(meanAnomaly) )
            watchDog2 = 0
            while( ( abs(ecceAnomaly-E0) > 1.0e-14 ) and ( watchDog2 < 1000 ) ) :
                E0 = ecceAnomaly
                ecceAnomaly = ecceAnomaly - ( ecceAnomaly - e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - e[sat]*cos(ecceAnomaly) )
                watchDog2 += 1
            sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
            e_cosEcceAnom_1 = 1 -  e[sat]*cosEcceAnom

            lari[0] = atan2( sqrt(1-e[sat]*e[sat]) * sinEcceAnom,  cosEcceAnom-e[sat] ) #True anomaly. Later, Argument of latitude
            lari[3] = sinEcceAnom * meanMotion * (1+e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

            lari[0] += omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            lari[0:3] = [ lari[0] + Cus[sat]*sin2Lat + Cuc[sat]*cos2Lat ,
                A*e_cosEcceAnom_1 + Crs[sat]*sin2Lat + Crc[sat]*cos2Lat ,
                        i0[sat] + Cis[sat]*sin2Lat + Cic[sat]*cos2Lat + idot[sat]*tk ]

            sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
            ddLat, ddRad, ddInc = 2 * ( Cus[sat] * cos2Lat - Cuc[sat] * sin2Lat ) * lari[3], 2 * ( Crs[sat] * cos2Lat - Crc[sat] * sin2Lat ) * lari[3], 2 * ( Cis[sat] * cos2Lat - Cic[sat] * sin2Lat ) * lari[3]

            lari[3:6] = [ lari[3] + ddLat, A*e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, idot[sat] + ddInc ]

            cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
            satXyzvt[ii][7] =  af0[sat] + ( af1[sat] + af2[sat]*tk ) * tk + F*e[sat]*sqrtA[sat]*sinEcceAnom - tgd[sat]  #Satellite clock bias errors

            dotAsceLon = Omegadot[sat] - WGS84_w_rads
            asceLon = Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * toe[sat]
            cosLon, sinLon = cos( asceLon ), sin( asceLon )

            xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
            xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

            satXyzvt[ii][1:4] = [+cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1], 
                                +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                                +sinInc*xyOrb[1] ]

            satXyzvt[ii][4:7] = [+cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satXyzvt[ii][2] ,
                                +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satXyzvt[ii][1] ,
                                +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

            theta = WGS84_w_rads * tau[ii]
            cosTheta, sinTheta  = cos( theta ), sin( theta )
            satXyzvt[ii][1:4]   = [ +cosTheta * satXyzvt[ii][1] + sinTheta * satXyzvt[ii][2],
                                    -sinTheta * satXyzvt[ii][1] + cosTheta * satXyzvt[ii][2],
                                    satXyzvt[ii][3] ]

        ########################################################################################################
        # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
        ########################################################################################################
        dXyz    = [ 0.0, 0.0, 0.0 ]
        b       = [ 0.0 for col in range( numVisibleSat ) ]
        A       = [ [0.0, 0.0, 0.0, 0.0 ] for row in range( numVisibleSat ) ]
        AA_     = [ [0.0, 0.0, 0.0, 0.0 ] for row in range( 4 ) ]

        ##
        # REFERENCE: http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/index_files/comp_pos.m
        #
        watchDog3 = 0
        while ( watchDog3 < 10 ) and ( (dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2]) > 0.01 ) and ( abs( dRcvXyzt[3] ) > 0.001 ) :
            Ab = [ 0.0, 0.0, 0.0, 0.0 ]
            for ll in range( numVisibleSat ) :
                dXyz    = [ satXyzvt[ll][1]-rcvXyzt[0], satXyzvt[ll][2]-rcvXyzt[1], satXyzvt[ll][3]-rcvXyzt[2] ]
                predict = sqrt( dXyz[0]*dXyz[0] + dXyz[1]*dXyz[1] + dXyz[2]*dXyz[2] )
                observ  = codeShifts[ll] * L1_CHIP_M
                A[ll]   = [ -dXyz[0] / observ, -dXyz[1] / observ, -dXyz[2] / observ, 1 ]

                sat = satXyzvt[ll][0] - 1
                #satAzel, rcvLalh = Xyz2AzimElev( satXyzvt[ll][1:4], rcvXyzt[0:3] )
                satAzel, ionoDelay = Ionosphere_DelayEx( satXyzvt[ll][1:4], rcvXyzt[0:3], epoch[1], alpha[sat], beta[sat] )
                E2      = satAzel[1] * satAzel[1]
                calib   = 0
                calib   = rcvXyzt[3]
                calib   -= C0*satXyzvt[ll][7]
                calib   += ( 2.312 / sin( sqrt( E2 + 0.001904 ) ) + 0.084 / sin( sqrt( E2 + 0.0006854 ) ) ) 
                #calib   += C0*ionoDelay

                b[ll]   = ( observ - predict - calib ) % ( L1_CODE_M )
                b[ll]   = b[ll] if b[ll] < HALF_L1_CODE_M else b[ll] - L1_CODE_M
                Ab      = [ Ab[0] + A[ll][0]*b[ll], Ab[1] + A[ll][1]*b[ll], Ab[2] + A[ll][2]*b[ll], Ab[3] + A[ll][3]*b[ll] ]
                tau[ll] +=  - b[ll] / C0

            #print(b)
            AA = [ [0.0, 0.0, 0.0, 0.0 ] for row in range( 4 ) ]
            for jj in range( 4 ) :
                for kk in range( 4 ) :
                    for ll in range( numVisibleSat ) :
                        AA[jj][kk] += ( A[ll][jj] * A[ll][kk] )
            ##########################################################################
            # http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
            # http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
            ###########################################################################
            AA_[0][0] = AA[1][1]*AA[2][2]*AA[3][3]+AA[1][2]*AA[2][3]*AA[3][1]+AA[1][3]*AA[2][1]*AA[3][2] -AA[1][1]*AA[2][3]*AA[3][2]-AA[1][2]*AA[2][1]*AA[3][3]-AA[1][3]*AA[2][2]*AA[3][1]
            AA_[0][1] = AA[0][1]*AA[2][3]*AA[3][2]+AA[0][2]*AA[2][1]*AA[3][3]+AA[0][3]*AA[2][2]*AA[3][1] -AA[0][1]*AA[2][2]*AA[3][3]-AA[0][2]*AA[2][3]*AA[3][1]-AA[0][3]*AA[2][1]*AA[3][2]
            AA_[0][2] = AA[0][1]*AA[1][2]*AA[3][3]+AA[0][2]*AA[1][3]*AA[3][1]+AA[0][3]*AA[1][1]*AA[3][2] -AA[0][1]*AA[1][3]*AA[3][2]-AA[0][2]*AA[1][1]*AA[3][3]-AA[0][3]*AA[1][2]*AA[3][1]
            AA_[0][3] = AA[0][1]*AA[1][3]*AA[2][2]+AA[0][2]*AA[1][1]*AA[2][3]+AA[0][3]*AA[1][2]*AA[2][1] -AA[0][1]*AA[1][2]*AA[2][3]-AA[0][2]*AA[1][3]*AA[2][1]-AA[0][3]*AA[1][1]*AA[2][2]
            AA_[1][0] = AA[1][0]*AA[2][3]*AA[3][2]+AA[1][2]*AA[2][0]*AA[3][3]+AA[1][3]*AA[2][2]*AA[3][0] -AA[1][0]*AA[2][2]*AA[3][3]-AA[1][2]*AA[2][3]*AA[3][0]-AA[1][3]*AA[2][0]*AA[3][2]
            AA_[1][1] = AA[0][0]*AA[2][2]*AA[3][3]+AA[0][2]*AA[2][3]*AA[3][0]+AA[0][3]*AA[2][0]*AA[3][2] -AA[0][0]*AA[2][3]*AA[3][2]-AA[0][2]*AA[2][0]*AA[3][3]-AA[0][3]*AA[2][2]*AA[3][0]
            AA_[1][2] = AA[0][0]*AA[1][3]*AA[3][2]+AA[0][2]*AA[1][0]*AA[3][3]+AA[0][3]*AA[1][2]*AA[3][0] -AA[0][0]*AA[1][2]*AA[3][3]-AA[0][2]*AA[1][3]*AA[3][0]-AA[0][3]*AA[1][0]*AA[3][2]
            AA_[1][3] = AA[0][0]*AA[1][2]*AA[2][3]+AA[0][2]*AA[1][3]*AA[2][0]+AA[0][3]*AA[1][0]*AA[2][2] -AA[0][0]*AA[1][3]*AA[2][2]-AA[0][2]*AA[1][0]*AA[2][3]-AA[0][3]*AA[1][2]*AA[2][0]
            AA_[2][0] = AA[1][0]*AA[2][1]*AA[3][3]+AA[1][1]*AA[2][3]*AA[3][0]+AA[1][3]*AA[2][0]*AA[3][1] -AA[1][0]*AA[2][3]*AA[3][1]-AA[1][1]*AA[2][0]*AA[3][3]-AA[1][3]*AA[2][1]*AA[3][0]
            AA_[2][1] = AA[0][0]*AA[2][3]*AA[3][1]+AA[0][1]*AA[2][0]*AA[3][3]+AA[0][3]*AA[2][1]*AA[3][0] -AA[0][0]*AA[2][1]*AA[3][3]-AA[0][1]*AA[2][3]*AA[3][0]-AA[0][3]*AA[2][0]*AA[3][1]
            AA_[2][2] = AA[0][0]*AA[1][1]*AA[3][3]+AA[0][1]*AA[1][3]*AA[3][0]+AA[0][3]*AA[1][0]*AA[3][1] -AA[0][0]*AA[1][3]*AA[3][1]-AA[0][1]*AA[1][0]*AA[3][3]-AA[0][3]*AA[1][1]*AA[3][0]
            AA_[2][3] = AA[0][0]*AA[1][3]*AA[2][1]+AA[0][1]*AA[1][0]*AA[2][3]+AA[0][3]*AA[1][1]*AA[2][0] -AA[0][0]*AA[1][1]*AA[2][3]-AA[0][1]*AA[1][3]*AA[2][0]-AA[0][3]*AA[1][0]*AA[2][1]
            AA_[3][0] = AA[1][0]*AA[2][2]*AA[3][1]+AA[1][1]*AA[2][0]*AA[3][2]+AA[1][2]*AA[2][1]*AA[3][0] -AA[1][0]*AA[2][1]*AA[3][2]-AA[1][1]*AA[2][2]*AA[3][0]-AA[1][2]*AA[2][0]*AA[3][1]
            AA_[3][1] = AA[0][0]*AA[2][1]*AA[3][2]+AA[0][1]*AA[2][2]*AA[3][0]+AA[0][2]*AA[2][0]*AA[3][1] -AA[0][0]*AA[2][2]*AA[3][1]-AA[0][1]*AA[2][0]*AA[3][2]-AA[0][2]*AA[2][1]*AA[3][0]
            AA_[3][2] = AA[0][0]*AA[1][2]*AA[3][1]+AA[0][1]*AA[1][0]*AA[3][2]+AA[0][2]*AA[1][1]*AA[3][0] -AA[0][0]*AA[1][1]*AA[3][2]-AA[0][1]*AA[1][2]*AA[3][0]-AA[0][2]*AA[1][0]*AA[3][1]
            AA_[3][3] = AA[0][0]*AA[1][1]*AA[2][2]+AA[0][1]*AA[1][2]*AA[2][0]+AA[0][2]*AA[1][0]*AA[2][1] -AA[0][0]*AA[1][2]*AA[2][1]-AA[0][1]*AA[1][0]*AA[2][2]-AA[0][2]*AA[1][1]*AA[2][0]

            detA = AA[0][0]*AA_[0][0] + AA[0][1]*AA_[1][0] + AA[0][2]*AA_[2][0] + AA[0][3]*AA_[3][0]
            if detA :
                dRcvXyzt= [ (AA_[0][0]*Ab[0]+AA_[0][1]*Ab[1]+AA_[0][2]*Ab[2]+AA_[0][3]*Ab[3])/detA ,
                            (AA_[1][0]*Ab[0]+AA_[1][1]*Ab[1]+AA_[1][2]*Ab[2]+AA_[1][3]*Ab[3])/detA , 
                            (AA_[2][0]*Ab[0]+AA_[2][1]*Ab[1]+AA_[2][2]*Ab[2]+AA_[2][3]*Ab[3])/detA ,
                            (AA_[3][0]*Ab[0]+AA_[3][1]*Ab[1]+AA_[3][2]*Ab[2]+AA_[3][3]*Ab[3])/detA ]               
                rcvXyzt = [ rcvXyzt[0] + dRcvXyzt[0] ,
                            rcvXyzt[1] + dRcvXyzt[1] , 
                            rcvXyzt[2] + dRcvXyzt[2] ,
                            rcvXyzt[3] + dRcvXyzt[3] ]

                # print ( detA )
                # #[ HDOP, VHDP, TDOP, PDOP, GDOP ]
                # DOP[0] = sqrt( AA_[0][0] + AA_[1][1] )
                # DOP[1] = sqrt( AA_[2][2] )
                # DOP[2] = sqrt( AA_[3][3] )
                # DOP[3] = sqrt( AA_[0][0] + AA_[1][1] + AA_[2][2] )
                # DOP[4] = sqrt( AA_[0][0] + AA_[1][1] + AA_[2][2] + AA_[3][3] )

                # # DOP = [ sqrt( AA_[0][0] + AA_[1][1] ), sqrt( AA_[2][2] ), sqrt( AA_[3][3] ), 
                # #         sqrt( AA_[0][0] + AA_[1][1] + AA_[2][2] ), sqrt( AA_[0][0] + AA_[1][1] + AA_[2][2] + AA_[3][3] ) ]
            elif ( 0 == detA ) and __DEBUG :
                print ( "Get_XyzPosition():: numVisibleSat = %d, detA = %f, [ HDOP, VHDP, TDOP, PDOP, GDOP ] = " % ( numVisibleSat, detA ) )
            watchDog3 += 1
        watchDog1 += 1
    return rcvXyzt, dRcvXyzt

def Xyz2AzimElev( satXyz, rcvXyz ):
    a, ecc2, lalh, satAzel = 6378137, 0.006694379990138, [0.0, 0.0, 0.0], [0.0, 0.0]

    p2 = rcvXyz[0]*rcvXyz[0] + rcvXyz[1]*rcvXyz[1]
    p1 = sqrt( p2 )
    r2 = p2 + rcvXyz[2]*rcvXyz[2]
    r1 = sqrt( r2 )
    b2 = r2 * ( 1 - ecc2 )
    b1 = sqrt( b2 )
    
    #phi = atan2( p1, xyz[2] )
    phi = atan2( rcvXyz[2], p1 * ( 1 - ecc2 ) )
    sinPhi = sin(phi)
    for ii in range( 1000 ) :
        RN  = a / sqrt( 1.0 - ecc2 * sinPhi * sinPhi )
        h   = p1/cos(phi) - RN
        phi = atan2( rcvXyz[2], p1 * ( 1 - ecc2 * RN / ( RN + h ) ) )
        sinPhi = sin(phi)
    lalh = [ phi, atan2( rcvXyz[1], rcvXyz[0] ), h ]

    #Find enu coordinates of position vector from user to satellite
    dXyz = [ satXyz[0] - rcvXyz[0], satXyz[1] - rcvXyz[1], satXyz[2] - rcvXyz[2] ]
    sinlat, coslat, sinlam, coslam = sin( lalh[0] ), cos( lalh[0] ), sin( lalh[1] ), cos( lalh[1] )
    enu = [ -sinlam*dXyz[0] + coslam*dXyz[1] ,
            -sinlat*coslam*dXyz[0] - sinlat*sinlam*dXyz[1] + coslat*dXyz[2] ,
            +coslat*coslam*dXyz[0] + coslat*sinlam*dXyz[1] + sinlat*dXyz[2]  ]
    satAzel = [ atan2( enu[0], enu[1] ), asin( enu[2]/sqrt( enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2] ) ) ]

    return satAzel, lalh

def Ionosphere_Delay( azel, lalh, tow, alpha, beta ):
    ele = azel[1] / ICD200PI

    #Calculate the earth-centred angle (elevation  in semicircles)
    psi = 0.0137 / ( ele + 0.11 ) -0.022
    #Compute Subionospheric latitude, the latitude of the Ionospheric Pierce Point
    phi_i = lalh[0] / ICD200PI + psi * cos( azel[0] )
    phi_i = 0.416 if phi_i > 0.416 else ( -0.416 if phi_i < -0.416  else phi_i )

    #Compute Subionospheric longitude, the longitude of the IPP
    lambda_i = lalh[1] / ICD200PI + psi * sin( azel[0] ) / cos( phi_i * ICD200PI )

    #Find the geomagnetic latitude of the IPP
    phi_m = phi_i + 0.064 * cos( (lambda_i -1.617) * ICD200PI )

    #Find the local time at the IPP
    t = ( 43200 * lambda_i + tow ) % 86400
    t = (t+86400 ) if t < 0 else ( (t-86400) if t >= 86400 else t )

    #Compute the amplitude of ionospheric delay
    Ai = ( ( alpha[3] * phi_m + alpha[2] ) * phi_m + alpha[1] ) * phi_m + alpha[0]
    
    if Ai < 0:
        Ai = 0.0
    
    #Compute the period of ionospheric delay
    Per = ( ( beta[3] * phi_m + beta[2] ) * phi_m + beta[1] ) * phi_m + beta[0]
    if Per < 72000:
        Per = 72000
    
    #Compute the phase of ionospheric delay
    Xi = 2 * ICD200PI * ( t - 50400 ) / Per
    Xi2 = Xi*Xi

    #Compute the slant factor (elevation  in semicircles)
    F = 1.0 + 16.0 * ( 0.53 - ele ) ** 3

    #Compute the ionospheric time delay
    ionoDelay = ( 5.0e-9 + Ai * ( 1.0 - Xi2 * ( 0.5 - Xi2/24.0 ) ) ) * F if abs( Xi ) < 1.57 else 5.0e-9 * F
    return ionoDelay


def Ionosphere_DelayEx( satXyz, rcvXyz, tow, alpha, beta ):
    a, ecc2, lalh, satAzel = 6378137, 0.0818191908426*0.0818191908426, [0.0, 0.0, 0.0], [0.0, 0.0]

    p2 = rcvXyz[0]*rcvXyz[0] + rcvXyz[1]*rcvXyz[1]
    p1 = sqrt( p2 )
    r2 = p2 + rcvXyz[2]*rcvXyz[2]
    r1 = sqrt( r2 )
    b2 = r2 * ( 1 - ecc2 )
    b1 = sqrt( b2 )
    
    #phi = atan2( p1, xyz[2] )
    phi = atan2( rcvXyz[2], p1 * ( 1 - ecc2 ) )
    sinPhi = sin(phi)
    for ii in range( 1000 ) :
        RN  = a / sqrt( 1.0 - ecc2 * sinPhi * sinPhi )
        h   = p1/cos(phi) - RN
        phi = atan2( rcvXyz[2], p1 * ( 1 - ecc2 * RN / ( RN + h ) ) )
        sinPhi = sin(phi)
    lalh = [ phi, atan2( rcvXyz[1], rcvXyz[0] ), h ]

    #Find enu coordinates of position vector from user to satellite
    dXyz = [ satXyz[0] - rcvXyz[0], satXyz[1] - rcvXyz[1], satXyz[2] - rcvXyz[2] ]
    sinlat, coslat, sinlam, coslam = sin( lalh[0] ), cos( lalh[0] ), sin( lalh[1] ), cos( lalh[1] )
    enu = [ -sinlam*dXyz[0] + coslam*dXyz[1] ,
            -sinlat*coslam*dXyz[0] - sinlat*sinlam*dXyz[1] + coslat*dXyz[2] ,
            +coslat*coslam*dXyz[0] + coslat*sinlam*dXyz[1] + sinlat*dXyz[2]  ]
    satAzel = [ atan2( enu[0], enu[1] ), asin( enu[2]/sqrt( enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2] ) ) ]


    # % function for computing an Ionospheric range correction for the   *
    # % GPS L1 frequency from the parameters broadcasted in the GPS      *
    # % Navigation Message.                                              *
    # % ==================================================================
    # % References:                                                      *
    # % Klobuchar, J.A., (1996) "Ionosphercic Effects on GPS", in        *
    # %   Parkinson, Spilker (ed), "Global Positioning System Theory and *
    # %   Applications, pp.513-514.                                      *
    # % ICD-GPS-200, Rev. C, (1997), pp. 125-128                         *
    # % NATO, (1991), "Technical Characteristics of the NAVSTAR GPS",    *
    # %   pp. A-6-31   -   A-6-33                                        *
    # % ==================================================================

    ele = satAzel[1] / ICD200PI

    #Calculate the earth-centred angle (elevation  in semicircles)
    psi = 0.0137 / ( ele + 0.11 ) -0.022
    #Compute Subionospheric latitude, the latitude of the Ionospheric Pierce Point
    phi_i = lalh[0] / ICD200PI + psi * cos( satAzel[0] )
    phi_i = 0.416 if phi_i > 0.416 else ( -0.416 if phi_i < -0.416  else phi_i )

    #Compute Subionospheric longitude, the longitude of the IPP
    lambda_i = lalh[1] / ICD200PI + psi * sin( satAzel[0] ) / cos( phi_i * ICD200PI )

    #Find the geomagnetic latitude of the IPP
    phi_m = phi_i + 0.064 * cos( (lambda_i -1.617) * ICD200PI )

    #Find the local time at the IPP
    t = ( 43200 * lambda_i + tow ) % 86400
    t = (t+86400 ) if t < 0 else ( (t-86400) if t >= 86400 else t )

    #Compute the amplitude of ionospheric delay
    Ai = ( ( alpha[3] * phi_m + alpha[2] ) * phi_m + alpha[1] ) * phi_m + alpha[0]
    
    if Ai < 0:
        Ai = 0.0
    
    #Compute the period of ionospheric delay
    Per = ( ( beta[3] * phi_m + beta[2] ) * phi_m + beta[1] ) * phi_m + beta[0]
    if Per < 72000:
        Per = 72000
    
    #Compute the phase of ionospheric delay
    Xi = 2 * ICD200PI * ( t - 50400 ) / Per
    Xi2 = Xi*Xi

    #Compute the slant factor (elevation  in semicircles)
    F = 1.0 + 16.0 * ( 0.53 - ele ) ** 3

    #Compute the ionospheric time delay
    ionoDelay = ( 5.0e-9 + Ai * ( 1.0 - Xi2 * ( 0.5 - Xi2/24.0 ) ) ) * F if abs( Xi ) < 1.57 else 5.0e-9 * F

    return satAzel, ionoDelay


def UT_Ionosphere_DelayEx( ephemerisFilename, epoch, visibleSats, codeShifts, rcvXyz, __DEBUG = True ):
    from math import floor, sqrt
    import gps3

    MAX_GPS_SATS        = 32
    WEEK_SECONDS        = 604800
    HALF_WEEK_SECONDS   = 302400
    #WGS84_u_m3s2        = 3.986004418e+14  # Gravitational Parameter ? = 398600.4418 km3/s2
    SQRT_WGS84_u_m3s2   = 1.99649803857e+7  # Gravitational Parameter ? = 398600.4418 km3/s2
    WGS84_w_rads        = 7.292158553e-5  # Rotation rate of the Earth ? =7.292158553e-5 rad/s  
    L1_CHIP_M           = 293.052256109
    L1_CODE_M           = 299792.458

    alpha   = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    beta    = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    svid    = [0 for ii in range(MAX_GPS_SATS) ]
    week    = [0 for ii in range(MAX_GPS_SATS) ]
    toe     = [0.0 for ii in range(MAX_GPS_SATS) ]
    sqrtA   = [0.0 for ii in range(MAX_GPS_SATS) ]
    deltan  = [0.0 for ii in range(MAX_GPS_SATS) ]
    M0      = [0.0 for ii in range(MAX_GPS_SATS) ]
    e       = [0.0 for ii in range(MAX_GPS_SATS) ]
    omega   = [0.0 for ii in range(MAX_GPS_SATS) ]
    Cus     = [0.0 for ii in range(MAX_GPS_SATS) ]
    Cuc     = [0.0 for ii in range(MAX_GPS_SATS) ]
    Crs     = [0.0 for ii in range(MAX_GPS_SATS) ]
    Crc     = [0.0 for ii in range(MAX_GPS_SATS) ]
    Cis     = [0.0 for ii in range(MAX_GPS_SATS) ]
    Cic     = [0.0 for ii in range(MAX_GPS_SATS) ]
    i0      = [0.0 for ii in range(MAX_GPS_SATS) ]
    idot    = [0.0 for ii in range(MAX_GPS_SATS) ]
    af0     = [0.0 for ii in range(MAX_GPS_SATS) ]
    af1     = [0.0 for ii in range(MAX_GPS_SATS) ]
    af2     = [0.0 for ii in range(MAX_GPS_SATS) ]

    Omega0      = [0.0 for ii in range(MAX_GPS_SATS) ]
    Omegadot    = [0.0 for ii in range(MAX_GPS_SATS) ]

    tgd     = [0.0 for ii in range(MAX_GPS_SATS) ]

    with open( ephemerisFilename, 'r' ) as f:
        ephData = json.load( f )
        for eph in ephData:

            ii = eph['svid'] - 1
            svid[ii]    = ii + 1

            if 'alpha0' in eph:
                alpha[ii][0]    = eph['alpha0']
                alpha[ii][1]    = eph['alpha1']
                alpha[ii][2]    = eph['alpha2']
                alpha[ii][3]    = eph['alpha3']

            if 'beta0' in eph:
                beta[ii][0]     = eph['beta0']
                beta[ii][1]     = eph['beta1']
                beta[ii][2]     = eph['beta2']
                beta[ii][3]     = eph['beta3']

            week[ii]    = eph['week']
            # if 'acc' in eph:
            #     self.acc[ii]          = eph['acc']
            # if 'l2code' in eph:
            #     self.l2code[ii]       = eph['l2code']
            idot[ii]       = eph['idot']
            # if 'iode' in eph:
            #      self.iode[ii]       = eph['iode']
            af0[ii]        = eph['af0']
            af1[ii]        = eph['af1']
            af2[ii]        = eph['af2']
            # if 'iodc' in eph:
            #      self.iodc[ii]       = eph['iodc']

            Crs[ii]        = eph['Crs']
            Crc[ii]        = eph['Crc']
            deltan[ii]     = eph['deltan']
            M0[ii]         = eph['M0']
            Cus[ii]        = eph['Cus']
            Cuc[ii]        = eph['Cuc']
            e[ii]          = eph['e']
            sqrtA[ii]      = eph['sqrtA']
            toe[ii]        = eph['toe']
            Cis[ii]        = eph['Cis']
            Cic[ii]        = eph['Cic']
            i0[ii]         = eph['i0']
            omega[ii]      = eph['omega']
            Omegadot[ii]   = eph['Omegadot']
            Omega0[ii]     = eph['Omega0']
            tgd[ii]        = eph['tgd']
            # if 'health' in eph:
            #     self.health[ii] = eph['health']
            # if 'toc' in eph:
            #     self.toc[ii]    = eph['toc']
            # if 'l2p' in eph:
            #     self.l2p[ii]    = eph['l2p']
            # if 'fit' in eph:
            #     self.fit[ii]    = eph['fit']


def Get_GPS_Positioning_Assistance( ephemerisFilename, rcvXyz, epoch = None, __DEBUG = False ) :
    from math import floor, sqrt
    MAX_GPS_SATS        = 32
    WEEK_SECONDS        = 604800
    HALF_WEEK_SECONDS   = 302400
    #WGS84_u_m3s2        = 3.986004418e+14  # Gravitational Parameter ? = 398600.4418 km3/s2
    SQRT_WGS84_u_m3s2   = 1.99649803857e+7  # Gravitational Parameter ? = 398600.4418 km3/s2
    WGS84_w_rads        = 7.292158553e-5  # Rotation rate of the Earth ? =7.292158553e-5 rad/s  
    L1_CHIP_M           = 293.052256109
    L1_CODE_M           = 299792.458

    alpha   = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    beta    = [ [0.0, 0.0, 0.0, 0.0] for ii in range(MAX_GPS_SATS) ]
    svid    = [ 0 for ii in range(MAX_GPS_SATS) ]
    week    = [ 0 for ii in range(MAX_GPS_SATS) ]
    toe     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    sqrtA   = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    deltan  = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    M0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    e       = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    omega   = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cus     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cuc     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Crs     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Crc     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cis     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Cic     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    i0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    idot    = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af0     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af1     = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    af2     = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    Omega0      = [ 0.0 for ii in range(MAX_GPS_SATS) ]
    Omegadot    = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    tgd     = [ 0.0 for ii in range(MAX_GPS_SATS) ]

    DOP = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
    with open( ephemerisFilename, 'r' ) as f:
        ephData = json.load( f )
        for eph in ephData:

            ii = eph['svid'] - 1
            svid[ii]    = ii + 1

            if 'alpha0' in eph:
                alpha[ii][0]    = eph['alpha0']
                alpha[ii][1]    = eph['alpha1']
                alpha[ii][2]    = eph['alpha2']
                alpha[ii][3]    = eph['alpha3']

            if 'beta0' in eph:
                beta[ii][0]     = eph['beta0']
                beta[ii][1]     = eph['beta1']
                beta[ii][2]     = eph['beta2']
                beta[ii][3]     = eph['beta3']

            week[ii]    = eph['week']
            # if 'acc' in eph:
            #     self.acc[ii]          = eph['acc']
            # if 'l2code' in eph:
            #     self.l2code[ii]       = eph['l2code']
            idot[ii]       = eph['idot']
            # if 'iode' in eph:
            #      self.iode[ii]       = eph['iode']
            af0[ii]        = eph['af0']
            af1[ii]        = eph['af1']
            af2[ii]        = eph['af2']
            # if 'iodc' in eph:
            #      self.iodc[ii]       = eph['iodc']

            Crs[ii]        = eph['Crs']
            Crc[ii]        = eph['Crc']
            deltan[ii]     = eph['deltan']
            M0[ii]         = eph['M0']
            Cus[ii]        = eph['Cus']
            Cuc[ii]        = eph['Cuc']
            e[ii]          = eph['e']
            sqrtA[ii]      = eph['sqrtA']
            toe[ii]        = eph['toe']
            Cis[ii]        = eph['Cis']
            Cic[ii]        = eph['Cic']
            i0[ii]         = eph['i0']
            omega[ii]      = eph['omega']
            Omegadot[ii]   = eph['Omegadot']
            Omega0[ii]     = eph['Omega0']
            tgd[ii]        = eph['tgd']
            # if 'health' in eph:
            #     self.health[ii] = eph['health']
            # if 'toc' in eph:
            #     self.toc[ii]    = eph['toc']
            # if 'l2p' in eph:
            #     self.l2p[ii]    = eph['l2p']
            # if 'fit' in eph:
            #     self.fit[ii]    = eph['fit']


    if epoch is None :
        now = time.gmtime()
        epoch_now = calendar.timegm( now ) - calendar.timegm( time.strptime( '1980-01-06 00:00:00', '%Y-%m-%d %H:%M:%S' ) ) + LEAPSECONDS_2017
        epoch = [ epoch_now // WEEK_SECONDS, epoch_now % WEEK_SECONDS ]


    visibleSat  = []
    xyzdxyz     = []
    lari        = [ 0.0 for col in range(6) ]
    xyOrb       = [ 0.0 for col in range(4) ]
    satXyzvt    = [ list([ float(0.0) for col in range(10) ]) for row in range( MAX_GPS_SATS ) ]
    dxyz        = [ list([ float(0.0) for col in range(3) ]) for row in range( MAX_GPS_SATS ) ]

    satXyzvt        = [ list([ 0.0 for col in range(8) ]) for row in range( MAX_GPS_SATS ) ]
    rcvXyzt         = [ rcvXyz[0], rcvXyz[1], rcvXyz[2], 0 ]
    b0              = 0.0   #initial receiver xlock biase
    dRcvXyzt        = [ 1.0, 1.0, 1.0, 1.0 ] 

    for ii in range( MAX_GPS_SATS ) :
        sat = svid[ii] - 1
        if svid[sat] == 0 :
            if __DEBUG:
                print("No ephemeris data for SV %d" % (sat+1))
            continue
        satXyzvt[ii][0] = svid[ii]
        tk = epoch[1] - toe[sat] #+ ( epoch[0] - week[sat] ) * WEEK_SECONDS
        if tk > HALF_WEEK_SECONDS :
            tk -= WEEK_SECONDS
        if tk < -HALF_WEEK_SECONDS :
            tk += WEEK_SECONDS
        #############################################################################
        # http://ccar.colorado.edu/asen5050/projects/projects_2008/xiaofanli/
        # The time of transmissionis is modified from TOW from navigation data and
        # the correction of is performed as
        #############################################################################
        A = sqrtA[sat] * sqrtA[sat]
        meanMotion = SQRT_WGS84_u_m3s2 / ( sqrtA[sat]*A ) + deltan[sat] # This value should be caliberated with adding Mean Motion Difference.
        meanAnomaly  = M0[sat] + meanMotion * tk

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
        ecceAnomaly = meanAnomaly + e[sat]*sin(meanAnomaly) / ( 1 - e[sat]*cos(meanAnomaly) )
        watchDog2 = 0
        while( ( abs(ecceAnomaly-E0) > 1.0e-14 ) and ( watchDog2 < 1000 ) ) :
            E0 = ecceAnomaly
            ecceAnomaly = ecceAnomaly - ( ecceAnomaly - e[sat]*sin(ecceAnomaly) - meanAnomaly ) / ( 1 - e[sat]*cos(ecceAnomaly) )
            watchDog2 += 1
        sinEcceAnom, cosEcceAnom = sin( ecceAnomaly ), cos( ecceAnomaly )
        e_cosEcceAnom_1 = 1 -  e[sat]*cosEcceAnom

        lari[0] = atan2( sqrt(1-e[sat]*e[sat]) * sinEcceAnom,  cosEcceAnom-e[sat] ) #True anomaly. Later, Argument of latitude
        lari[3] = sinEcceAnom * meanMotion * (1+e[sat]*cos(lari[0])) / ( sin(lari[0]) * e_cosEcceAnom_1 * e_cosEcceAnom_1 ) #Rate of true anomaly, later, latitude

        lari[0] += omega[sat]  #Latitude = True anomaly + Perigee : Previously True anomaly; Now, Argument of latitude
        sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
        lari[0:3] = [ lari[0] + Cus[sat]*sin2Lat + Cuc[sat]*cos2Lat ,
            A*e_cosEcceAnom_1 + Crs[sat]*sin2Lat + Crc[sat]*cos2Lat ,
                    i0[sat] + Cis[sat]*sin2Lat + Cic[sat]*cos2Lat + idot[sat]*tk ]

        sin2Lat, cos2Lat = sin( 2*lari[0] ), cos( 2*lari[0] )
        ddLat, ddRad, ddInc = 2 * ( Cus[sat] * cos2Lat - Cuc[sat] * sin2Lat ) * lari[3], 2 * ( Crs[sat] * cos2Lat - Crc[sat] * sin2Lat ) * lari[3], 2 * ( Cis[sat] * cos2Lat - Cic[sat] * sin2Lat ) * lari[3]

        lari[3:6] = [ lari[3] + ddLat, A*e[sat]*sinEcceAnom*meanMotion/e_cosEcceAnom_1 + ddRad, idot[sat] + ddInc ]

        cosLat, sinLat, cosInc, sinInc = cos( lari[0] ), sin( lari[0] ), cos( lari[2] ), sin( lari[2] )
        satXyzvt[ii][7] =  af0[sat] + ( af1[sat] + af2[sat]*tk ) * tk + F*e[sat]*sqrtA[sat]*sinEcceAnom - tgd[sat]  #Satellite clock bias errors

        dotAsceLon = Omegadot[sat] - WGS84_w_rads
        asceLon = Omega0[sat] + dotAsceLon*tk - WGS84_w_rads * toe[sat]
        cosLon, sinLon = cos( asceLon ), sin( asceLon )

        xyOrb[0:2] = [ cosLat*lari[1], sinLat*lari[1] ]
        xyOrb[2:4] = [ cosLat*lari[4] - xyOrb[1]*lari[3], sinLat*lari[4] + xyOrb[0]*lari[3] ]

        satXyzvt[sat][1:4] = [+cosLon*xyOrb[0] - sinLon*cosInc*xyOrb[1], 
                            +sinLon*xyOrb[0] + cosLon*cosInc*xyOrb[1] ,
                            +sinInc*xyOrb[1] ]

        satXyzvt[sat][4:7] = [+cosLon*xyOrb[2] - sinLon*cosInc*xyOrb[3] + sinLon*sinInc*lari[5]*xyOrb[1] - dotAsceLon*satXyzvt[ii][2] ,
                            +sinLon*xyOrb[2] + cosLon*cosInc*xyOrb[3] - cosLon*sinInc*lari[5]*xyOrb[1] + dotAsceLon*satXyzvt[ii][1] ,
                            +sinInc*xyOrb[3] + cosInc*lari[5]*xyOrb[1] ]

        satXyzvt[sat][8:10] = [ af0[sat] + F*e[sat]*sqrtA[sat]*sinEcceAnom - tgd[sat], af1[sat] + af2[sat]*tk ] 
        satXyzvt[sat][7] = satXyzvt[sat][8] + satXyzvt[sat][9] * tk

        dxyz[sat] = [ satXyzvt[sat][1] - rcvXyz[0], satXyzvt[sat][2] - rcvXyz[1], satXyzvt[sat][3] - rcvXyz[2]  ]

        xyzdxyz1 =  dxyz[sat][0]*rcvXyz[0] + dxyz[sat][1]*rcvXyz[1] + dxyz[sat][2]*rcvXyz[2]
        if xyzdxyz1 > 0 :
            visibleSat.append( svid[ii] )
            xyzdxyz.append( xyzdxyz1 ) 
    
    numVisibleSat   = len( visibleSat )
    if numVisibleSat > 0:
        satXyzvt_       = [ list([ float(0.0) for col in range(8) ]) for row in range( numVisibleSat ) ]
        assist          = [ list([ float(0.0) for col in range(4) ]) for row in range( numVisibleSat ) ]
        for ii in range( numVisibleSat ) :
            jj = visibleSat[ii] - 1
            satXyzvt_[ii]    = satXyzvt[jj][0:8]
            pseudorange     = sqrt( dxyz[jj][0]*dxyz[jj][0] + dxyz[jj][1]*dxyz[jj][1] + dxyz[jj][2]*dxyz[jj][2] ) ## - satXyzvt[jj][7] * C0 
            pseudorangedot  = ( satXyzvt[jj][4]*dxyz[jj][0] + satXyzvt[jj][5]*dxyz[jj][1] + satXyzvt[jj][6]*dxyz[jj][2]  ) / pseudorange
            
            #PRN, L1 Doppler (Hz), C/A codephase(chips), pseudorange?
            assist[ii] = [ jj+1, ( pseudorangedot/C0 - satXyzvt[jj][9] ) * L1_FREQ, ( pseudorange % L1_CODE_M ) / L1_CHIP_M, pseudorange  ]
    else:
        assist, satXyzvt_ = None, None
    
    return assist, satXyzvt_
