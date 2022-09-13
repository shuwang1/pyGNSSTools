import GnssUtil, Ephemeris, datetime, json, os, subprocess
from math import sqrt, floor

# print ( "\n[STEP 1] ", GnssUtil.Get_NowEpoch(), " Download Ephemeris Data From A GPS Worldwide Datacenter ... " )
# gpsEph_nasa = Ephemeris.Ephemeris( )
# numEph, ephemerisFilename = gpsEph_nasa.Get_GPS_Ephemeris_From_NASA( forceUpdate = True )
# gpsEph_nasa.Save_Ephemeris( ephemerisFilename )
# print ( "\tgpsEph_nasa.Get_GPS_Ephemeris_From_NASA():: The number of ephemeris is %d and they are saved as %s" % ( numEph, ephemerisFilename ) )

Test_Case = 9
if 9 == Test_Case :
    satelliteSignalEstimatesFilename = "rtl192034_20190922021_246104_L1C_GO.json"
    epoch = [ 2047, 246105 ]
    #epoch = [ 2047, 243777 ]

elif 8 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190860647_880MS_L1C_GO.json"
    #epoch = [ 2046, 283636 ]
    #epoch = [ 2046, 265636 ]
    #epoch = [ 2046, 265063 ]
    epoch = [ 2046, 265057 ]
elif 7 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190860414_880MS_L1C_GO.json"
    #epoch = [ 2046, 274486 ]
    #epoch = [ 2046, 256486 ]
    epoch = [ 2046, 257066 ]
elif 6 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190860147_880MS_L1C_GO.json"
    epoch = [2046, 265633]
elif 5 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190852321_880MS_L1C_GE.json"
    epoch = [2046, 256862]
elif 4 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190852045_880MS_L1C_GE.json"
    epoch = [2046, 247543]
elif 3 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190851820_880MS_L1C_GE.json"
    epoch = [2046, 238831]
elif 2 == Test_Case :
    satelliteSignalEstimatesFilename = "home00USA_20190851629_880MS_L1C_GE.json"
    epoch = [2046, 232158]
elif 1 == Test_Case :
    satelliteSignalEstimatesFilename = "Synapse_140423_1433_Search_Output.json"
    epoch = GnssUtil.Utc2Epoch( '2014-04-23 19:33:30 UTC' )

print("\n[STEP 0] ", GnssUtil.Get_NowEpoch(), " Start Initialization ... ")
MAX_NUM_EXPERIMENTS = 4
lats_deg, lons_deg, latlons_deg, latlons = [ 0.0 for ii in range( MAX_NUM_EXPERIMENTS ) ], [ 0.0 for ii in range( MAX_NUM_EXPERIMENTS ) ], [ [0.0, 0.0] for ii in range( MAX_NUM_EXPERIMENTS ) ], [ [0.0, 0.0] for ii in range( MAX_NUM_EXPERIMENTS ) ]
xyzt = [ [0.0, 0.0, 0.0, 0.0] for ii in range( MAX_NUM_EXPERIMENTS ) ]

#lalh_deg    = [ +(32+56.4806/60), -( 117+14.4715/60 ), 100.0 ] ## 12526 High Bluff Drive, San Diego CA 92130, United States
lalh_deg    = [ +32.945280, -117.218690, 100.0 ] ## 4548 Campobello St, San Diego CA 92130, United States
lalh        = [ lalh_deg[0] * GnssUtil.PI_DEG, lalh_deg[1] * GnssUtil.PI_DEG, 100.0 ]
rcvXyz      = GnssUtil.Lalh2Xyz( lalh )
rcvXyzt     = [ rcvXyz[0], rcvXyz[1], rcvXyz[2], 0 ]

mm = 0
lats_deg[mm], lons_deg[mm], latlons_deg[mm], latlons[mm] = lalh_deg[0], lalh_deg[1], [ lalh_deg[0], lalh_deg[1] ], [ lalh[0], lalh[1] ]
xyzt[mm] = rcvXyzt
print("\tThe original XYZ position is\n\txyzt[0] = rcvXyzt = ", xyzt[mm], "\n\tlalh_deg = ", lalh_deg)

print ( "\n[STEP 1] ", GnssUtil.Get_NowEpoch(), " Download Ephemeris Data From A GPS Worldwide Datacenter ... " )
gpsEph_nasa = Ephemeris.Ephemeris( )

print( "\n[STEP 3] ", GnssUtil.Get_NowEpoch() , " Download RF Signal Samples From A Compatible RF Front-End And Start Digital Signal Processing.")
localDir = os.getcwd() 
if not localDir.endswith( '/' ) :
    localDir += '/' 

workingSubdirectory = 'RTL_SDR_RF_Signal_Samples/'

print( "\n[STEP 5] ", GnssUtil.Get_NowEpoch(), " Parse and Preprocess The Satellite Signal Detection Results " )

print( "\tparse the GPS signal estimates from ", satelliteSignalEstimatesFilename )
with open( localDir + workingSubdirectory + satelliteSignalEstimatesFilename, 'r' ) as f:
    satelliteSignalEstimates = json.load( f )

obsNumVisibles  = 0
obsVisibleSats  = []
obsCodeShifts   = []
obsSNRs         = []       
for estimate in satelliteSignalEstimates:
    if ( estimate[ "isDetected" ] == True ) :
        obsVisibleSats.append( estimate["svid"] )
        obsCodeShifts.append( estimate["codeShift_chip"] )
        obsSNRs.append( estimate["snr_dB"] )
        obsNumVisibles += 1

status = os.chdir( localDir )
print( "\tstatus = os.chdir( localDir ) = ", status )

print( "\tThe number of visible SVs detected from the RF signal samples is: ", obsNumVisibles )
print( "\tsvids = ", obsVisibleSats )
print( "\tObserved codeShifts = ", obsCodeShifts )
print( "\tObserved SNRs = ", obsSNRs )   

if( obsNumVisibles < 4 ):
    print ( "\n[Step 7] ", GnssUtil.Get_NowEpoch(), " Need More RF Signals Of Visible Satellites And Start Searching Assistant RF Signals .... " )
    exit( 0 )

print( "\n[Step 8] ", GnssUtil.Get_NowEpoch(), "Start PostProcessing And Positioning From The RF Measurements On The %d Visible Satellites" % obsNumVisibles )
print ( "\tThe original position is ", rcvXyzt, ", the total number of visibles is ", obsNumVisibles, "\n\tthe epoch [Week Number, ToA] = ", epoch )

kk = 0
epoch_ = [ epoch[0], epoch[1] + kk ]
rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
print("\tThe broadcast ephemeris file is ", updatedEphemeris )
jj = 5
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris, GnssUtil.Epoch2Utc( epoch_ ) )
eRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
print ( "\t[Initial1] With ", epoch_ , " %d SVs, the estimate is " % jj, eRcvXyzt, " and the error is %.2f m. " % obsError )

from scipy.optimize import minimize
def XyzPositioningError( kk_ ):
    global updatedEphemeris, jj, epoch, obsVisibleSats, obsCodeShifts, rcvXyzt
    rcvXyzt_    = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris,  [ epoch[0], epoch[1]+kk_ ], obsVisibleSats[0:jj ], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
    obsErr      = sqrt( (rcvXyzt[0]-rcvXyzt_[0])**2 + (rcvXyzt[1]-rcvXyzt_[1])**2 + (rcvXyzt[2]-rcvXyzt_[2])**2 )
    return obsErr

ii = 0
res = minimize( XyzPositioningError, ii, method='nelder-mead', options={'xtol': 1e-10, 'disp': False})
print( "\t[Initial2] jj = %d, res.x = %f, res.fun = %f " % ( jj, res.x, res.fun ) )

obsError_, ii_, jj_ = 10000, 0, obsNumVisibles

# ii0 = 0
# for ii in range( -3600, 1201, 600 ):
#     kk = ii0 + ii
#     epoch_ = [ epoch[0], epoch[1]+kk ]
#     rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
#     updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
#     if not os.path.exists( updatedEphemeris ) :
#         updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
#         print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
#     gpsEph_nasa.Load_EphemerisJsonFile( updatedEphemeris )
#     for jj in range( 4 , obsNumVisibles+1, 1 ):
#         eRcvXyzt = gpsEph_nasa.Get_Receiver_XyzPosition( epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], propagationDelayCompensation = False )
#         obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
#         if( obsError < obsError_ ):
#             obsError_, ii_, jj_ = obsError, kk, jj
#         if( obsError < 5000 ):
#             print ( "\t\tWith ToA = %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m. " % obsError )
# print ( "\t[Rudimentary] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

# obsError_, ii_, jj_ = 10000, 0, obsNumVisibles
# ii0 = 0
# for ii in range( -3600, 1201, 600 ):
#     kk = ii0 + ii
#     epoch_ = [ epoch[0], epoch[1]+kk ]
#     rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
#     updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
#     if not os.path.exists( updatedEphemeris ) :
#         updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
#         print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
#     for jj in range( 4 , obsNumVisibles, 1 ) :
#         eRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
#         obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
#         if( obsError < obsError_ ):
#             obsError_, ii_, jj_ = obsError, kk, jj
#         if( obsError < 5000 ):
#             print ( "\t\tWith ToA = %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m. " % obsError )
# print ( "\t[Rudimentary] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

ii1 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
print( "\t\t", GnssUtil.Epoch2Utc( epoch_ )  )
for ii in range( - 100, 101, 2 ):
    kk = ii1 + ii
    epoch_ = [ epoch[0], epoch[1]+kk ]
    rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
    updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
    if not os.path.exists( updatedEphemeris ) :
        updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
        print("\tThe updated broadcast ephemeris file is ", updatedEphemeris, GnssUtil.Epoch2Utc( epoch_ ) )
    for jj in range( max( 4, jj_-3 ), min( jj_+3, obsNumVisibles ), 1 ):
        eRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 1000 ):
            print ( "\t\tWith ToA = %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m. " % obsError )
print ( "\t[Intermediate] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

ii2 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris, GnssUtil.Epoch2Utc( epoch_ ) )
for ii in range( - 2000, 2000, 1 ):
    kk = ii2 + ii*0.001
    epoch_ = [ epoch[0], epoch[1]+kk ]
    for jj in range( max( 4, jj_-2 ), min( jj_+2, obsNumVisibles ), 1 ):
        eRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_ , obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 100 ):
            print ( "\t\tWith %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m. " % obsError )
print ( "\t[Precious] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

ii3 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris, GnssUtil.Epoch2Utc( epoch_ ) )
for ii in range( - 1000, 1000, 1 ):
    kk = ii3 + ii*0.000001
    epoch_ = [ epoch[0], epoch[1]+kk ]
    for jj in range( max( 4, jj_-1 ), min( jj_+1, obsNumVisibles ), 1 ):
        eRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 10 ):
            print ( "\t\tWith %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m. " % obsError )
print ( "\t[Final1] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

exit(0)


mm = 1
xyzt[mm]    = eRcvXyzt
lalh_       = GnssUtil.Xyz2Lalh( xyzt[mm][0:3] )
lalh_deg_   = [ lalh_[0] / GnssUtil.PI_DEG, lalh_[1] / GnssUtil.PI_DEG, 100.0 ]
lats_deg[mm], lons_deg[mm], latlons_deg[mm], latlons[mm] = lalh_deg_[0], lalh_deg_[1], [ lalh_deg_[0], lalh_deg_[1] ], [ lalh_[0], lalh_[1] ]



obsError_, ii_, ll_= 5515.10, 0, 0
dEpoch = range( -18000, 18001, 180 )
for ll in dEpoch:
    rcvEpoch = GnssUtil.Epoch2Gps( [ epoch[0], epoch[1]+ll ] )
    updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
    for jj in range( max( 4, jj_-3 ), min( jj_+3, obsNumVisibles ), 1 ):
        res = minimize( XyzPositioningError, ii_, method='nelder-mead', options={'xtol': 1e-10, 'disp': False})
        #print( "ll = %d, jj = %d, res.x = %f, res.fun = %f " % ( ll, jj, res.x, res.fun ) )
        if( res.fun < obsError_ ):
            obsError_, ii_, jj_, ll_ = res.fun, res.x, jj, ll
            print ( "\tobsError_ = %f, ii_ = %f, jj_ = %d, ll_ = %d" % (obsError_, ii_, jj_, ll_) )
            print( res )
        if( res.fun < 10 ):
            print ( "\t\tWith %d s & %d visibles, the estimate is " % (epoch[1]+res.x, jj), eRcvXyzt, " and the error is %.2f m. " % res.fun )
print ( "\t[Final2] With ", [ epoch[0], epoch[1] + ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m." % obsError_ )

mm = 2
xyzt[mm]    = eRcvXyzt
lalh_       = GnssUtil.Xyz2Lalh( xyzt[mm][0:3] )
lalh_deg_   = [ lalh_[0] / GnssUtil.PI_DEG, lalh_[1] / GnssUtil.PI_DEG, 100.0 ]
lats_deg[mm], lons_deg[mm], latlons_deg[mm], latlons[mm] = lalh_deg_[0], lalh_deg_[1], [ lalh_deg_[0], lalh_deg_[1] ], [ lalh_[0], lalh_[1] ]

print ( latlons_deg )

import gmplot, webbrowser

localDir = os.getcwd()
if not localDir.endswith( '/' ) :
    localDir += '/'
htmlMapFile = localDir + 'SanDiegoMap.html'


gmap        = gmplot.GoogleMapPlotter( lats_deg[0], lons_deg[0], 14 )
#gmap.apikey = "AIzaSyDq7hcNIAssb7DboRsln6qhBsbcuYnmhy0"
gmap.apikey = "AIzaSyA8i39YHp23B_ayfm9Yg8bunfAqlJQRDRE"

#gmap.plot( lat, lat, 'cornflowerblue', edge_width=10)
gmap.scatter( lats_deg, lons_deg, 'r', 0.1, marker = True )
gmap.draw( htmlMapFile )

GnssUtil.Plot_Positions( latlons_deg, htmlMapFile  )
webbrowser.open( htmlMapFile, new = 2, autoraise = True )

UTC_OFFSET_TIMEDELTA = datetime.datetime.utcnow() - datetime.datetime.now()
local_datetime = datetime.datetime.strptime("2008-09-17 14:04:00", "%Y-%m-%d %H:%M:%S")
result_utc_datetime = local_datetime + UTC_OFFSET_TIMEDELTA
result_utc_datetime.strftime("%Y-%m-%d %H:%M:%S")