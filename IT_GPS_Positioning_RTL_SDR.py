import GnssUtil, Ephemeris, datetime, time, json, os, subprocess
from math import sqrt, floor

# print ( "\n[STEP 1] ", GnssUtil.Get_NowEpoch(), " Download Ephemeris Data From A GPS Worldwide Datacenter ... " )
# gpsEph_nasa = Ephemeris.Ephemeris( )
# numEph, ephemerisFilename = gpsEph_nasa.Get_GPS_Ephemeris_From_NASA( forceUpdate = True )
# gpsEph_nasa.Save_Ephemeris( ephemerisFilename )
# print ( "\tgpsEph_nasa.Get_GPS_Ephemeris_From_NASA():: The number of ephemeris is %d and they are saved as %s" % ( numEph, ephemerisFilename ) )


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

localDir = os.getcwd() 
if not localDir.endswith( '/' ) :
    localDir += '/' 

workingSubdirectory = 'RTL_SDR_RF_Signal_Samples/'
if not localDir.endswith( workingSubdirectory ) :
    localDir = localDir + workingSubdirectory

directory = os.path.dirname( localDir )
if not os.path.exists( directory ) :
    os.makedirs( directory )

status = os.chdir( localDir  )
print( "\tThe current working directory is %s; os.chdir() = " % localDir , status )

epoch = GnssUtil.Get_NowEpoch()
print ( "\n[STEP 1] ", epoch, " Download Ephemeris Data From GPS Worldwide Datacenters ... " )

Broadcast_Ephemeris_Servers = "Ephemeris_Servers_Configuration.json"
if not os.path.exists( Broadcast_Ephemeris_Servers ):
    print( "\t The configuration file %s for downloading broadcast Ephemeris data is unavailable." % Broadcast_Ephemeris_Servers )
unzippedRinexFilenameList = GnssUtil.Download_Broadcast_RinexFiles( epoch, ephemeris_Server_Configuration = Broadcast_Ephemeris_Servers, forceUpdate = False, __DEBUG = False )
print( "\tThe downloaded Rinex files are ", unzippedRinexFilenameList  )

rcvGpsTime = time.gmtime()
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvGpsTime.tm_year, rcvGpsTime.tm_yday, rcvGpsTime.tm_hour if rcvGpsTime.tm_min < 30 else ( rcvGpsTime.tm_hour + 1) , 30 if rcvGpsTime.tm_min < 30 else 0 )
updatedEphemeris = GnssUtil.Merge_GPS_Ephemeris( unzippedRinexFilenameList, updatedEphemeris, epochs = None, svids = None )
print ( "\tGnssUtil.Merge_GPS_Ephemeris():: The ephemeris data are saved in %s" % ( updatedEphemeris ) )

gpsEph_nasa = Ephemeris.Ephemeris( )
numEph, ephemerisFilename = gpsEph_nasa.Get_GPS_Ephemeris_From_NASA( forceUpdate = True )
gpsEph_nasa.Save_Ephemeris( ephemerisFilename )
print ( "\tgpsEph_nasa.Get_GPS_Ephemeris_From_NASA():: The number of ephemeris is %d and they are saved as %s" % ( numEph, ephemerisFilename ) )

print("\n[STEP 2.1] ", GnssUtil.Get_NowEpoch(), " Predict The Satellite Signal Parameters For All Healthy GPS Satellites And Identify The Visibles ... ")
epoch = GnssUtil.Get_NowEpoch()
ephAssist, satXyzvt = gpsEph_nasa.Get_PositionAssist( rcvXyz, epoch, propagationDelayCompensation = False )
ephVisibles     = len( ephAssist )
ephVisibleSats  = [ 0.0 for ii in range(ephVisibles) ]
preCodeShifts   = [ 0.0 for ii in range(ephVisibles) ]
obsPseudoranges = [ 0.0 for ii in range(ephVisibles) ]
dPseudoranges   = [ 0.0 for ii in range(ephVisibles) ]
est             = [ 0.0 for ii in range(ephVisibles) ]
for ii in range( ephVisibles ) :
    ephVisibleSats[ii]  = ephAssist[ii][0]
    preCodeShifts[ii]   = round( ephAssist[ii][2] * 32 ) / 32.0
    obsPseudoranges[ii] = round( ephAssist[ii][3] * 1000 ) / 1000
    dPseudoranges[ii]   = round( ( ephAssist[ii][3] - ( floor( ephAssist[ii][3] / GnssUtil.L1_CODE_M ) * GnssUtil.L1_CODE_M )  - preCodeShifts[ii] * GnssUtil.L1_CHIP_M ) * 100 ) / 100
    print( '\t', ephAssist[ii][0], '\t', ephAssist[ii][1], '\t', preCodeShifts[ii], '\t',  satXyzvt[ii][0] )
print( "\tgpsEph.Get_PositionAssist( lalh, epoch ):: From the location rcvXyz = ", rcvXyz, ", the %d visibles are " % ephVisibles )
print( "\tpreCodeShift = ", preCodeShifts )
print( "\tobsPseudorange = ", obsPseudoranges )
print( "\tdPseudorange = ", dPseudoranges )

print( "\n[STEP 2.2] ", GnssUtil.Get_NowEpoch(), " Generate Satellite Positioning Assistance File, which includes the satellite signal predicts for %d visible satellites." % ephVisibles )
satellitePositioningAssistance = []
for ii in range( ephVisibles ) :
    satellitePositioningAssistance.append({
        'svid': ephAssist[ii][0],
        'dopplerShift_Hz': round( ephAssist[ii][1] * 1000 ) / 1000 ,
        'codeShift_chip': round( ephAssist[ii][2] * 32 ) / 32 ,
        'pseudorange_m': round( ephAssist[ii][3] * 1000 ) / 1000 ,
        'weekNumber': epoch[0],
        'timeOfApplicability_s': epoch[1]
    }) 
    print( "\tPRN = %02d,  dopplerShift = %+8.2f,  codeShift = %4d,  pseudorange = %f" % ( ephAssist[ii][0], ephAssist[ii][1], ephAssist[ii][2], ephAssist[ii][3] ) )

gpsDateTime = GnssUtil.Epoch2Gps( epoch )
satellitePositioningAssistanceFilename = "assist_%4d%03d%02d%02d_L1C_GE.json" % ( 
                                        gpsDateTime.tm_year, 
                                        gpsDateTime.tm_yday, 
                                        gpsDateTime.tm_hour, 
                                        gpsDateTime.tm_min 
                                        )
with open( satellitePositioningAssistanceFilename, 'w') as f:  
    json.dump( satellitePositioningAssistance, f )

with open( satellitePositioningAssistanceFilename ) as f:
    satellitePositioningAssistance = json.load( f )

print( "\t save into a Json file: ", satellitePositioningAssistanceFilename )

print( "\n[STEP 3] ", GnssUtil.Get_NowEpoch() , " Download RF Signal Samples From A Compatible RF Front-End And Start Digital Signal Processing.")

status = subprocess.call( "cp  %s  ./" % ( satellitePositioningAssistanceFilename ) , shell=True )
print( "\tcp  %s  ./; status=" % ( satellitePositioningAssistanceFilename ), status )

status = subprocess.call( "cp  %s  ~/Documents/Cloud\ Based\ GPS/Code_Phase_Search/DerivedData/Code_Phase_Search/Build/Products/Debug/" % ( satellitePositioningAssistanceFilename ) , shell=True )
print( "\tcp  %s  ~/Documents/Cloud\ Based\ GPS/Code_Phase_Search/DerivedData/Code_Phase_Search/Build/Products/Debug/; status=" % ( satellitePositioningAssistanceFilename ), status )

status = os.getcwd()
print( "\tstatus = os.getcwd() = ", status )

status = subprocess.call( "cp  ~/Documents/Cloud\ Based\ GPS/Code_Phase_Search/DerivedData/Code_Phase_Search/Build/Products/Debug/Code_Phase_Search  ./", shell=True )
print( "\tcp  ~/Documents/Cloud\ Based\ GPS/Code_Phase_Search/DerivedData/Code_Phase_Search/Build/Products/Debug/Code_Phase_Search  ./; status=", status )


gpsDateTime = GnssUtil.Epoch2Gps( epoch )
satelliteSignalSamplesFilename      = "rtl192034_%4d%03d%02d%02d_%d_L1C_GO.piq" % ( 
                                        gpsDateTime.tm_year, 
                                        gpsDateTime.tm_yday, 
                                        gpsDateTime.tm_hour, 
                                        gpsDateTime.tm_min,
                                        epoch[1]
                                        )
satelliteSignalEstimatesFilename    = "rtl192034_%4d%03d%02d%02d_%d_L1C_GO.json" % ( 
                                        gpsDateTime.tm_year, 
                                        gpsDateTime.tm_yday, 
                                        gpsDateTime.tm_hour, 
                                        gpsDateTime.tm_min,
                                        epoch[1]
                                        )

CodePhaseSearchConfigFilename       = "RTL_SDR_GPS_Positioning_Engine_Configuration.json"
print( "\n\tthe configuration file, %s \n\tthe positioning assistance file, %s" % ( CodePhaseSearchConfigFilename, satellitePositioningAssistanceFilename ) )
with open( CodePhaseSearchConfigFilename, 'rt' ) as f:
    CodePhaseSearchConfig = json.load( f )
    CodePhaseSearchConfig['Default_Carrier_Frequency_Hz']                   = 1575420000
    CodePhaseSearchConfig['Default_Sample_Rate_Hz']                         = 2046000
    CodePhaseSearchConfig['Default_Oversampling_Ratio']                     = 1
    CodePhaseSearchConfig['Default_Upsampling_Ratio']                       = 2

    CodePhaseSearchConfig['Default_Initial_Offset_MS']                      = 0.5
    CodePhaseSearchConfig['Default_Initial_Offset_Sample']                  = 2046
    CodePhaseSearchConfig['Default_Sample_Duration_MS']                     = 900

    CodePhaseSearchConfig['Default_Coherent_Integration_Time_MS']           = 11
    CodePhaseSearchConfig['Default_Coherence_Interval_MS']                  = 11
    CodePhaseSearchConfig['Default_Noncoherence_Integration']               = 40
    CodePhaseSearchConfig['Default_Doppler_Shift_Search_Window_Low_Bin']    = -150
    CodePhaseSearchConfig['Default_Doppler_Shift_Search_Window_High_Bin']   = 150
    
    CodePhaseSearchConfig['Default_Sample_Input_Filename']                  = satelliteSignalSamplesFilename
    
    CodePhaseSearchConfig['NeedAssistance']                                 = True
    CodePhaseSearchConfig['Default_Assistance_Filename']                    = satellitePositioningAssistanceFilename

    CodePhaseSearchConfig['Default_Use_Remote_Host']                        = True
    CodePhaseSearchConfig['Default_Remote_Host_IP_Address']                 = "70.95.74.204"
    #CodePhaseSearchConfig['Default_Remote_Host_IP_Address']                 = "192.168.0.110"
    #CodePhaseSearchConfig['Default_Remote_Host_IP_Address']                 = "192.168.0.59"
    #CodePhaseSearchConfig['Default_Remote_Host_IP_Address']                 = "192.168.0.13"
    CodePhaseSearchConfig['Default_Remote_Host_IP_Port']                    = 1234
    
    CodePhaseSearchConfig['Default_Sample_Output_Filename']                 = satelliteSignalSamplesFilename
    CodePhaseSearchConfig['Default_Search_Output_Filename']                 = satelliteSignalEstimatesFilename

    
with open( CodePhaseSearchConfigFilename, 'wt' ) as f:
    json.dump( CodePhaseSearchConfig, f )

code_phase_search_parameter = "--tc %d --tci %d --quant 0 --ot %s" % ( CodePhaseSearchConfig['Default_Coherent_Integration_Time_MS'], CodePhaseSearchConfig['Default_Coherence_Interval_MS'], CodePhaseSearchConfig['Default_Search_Output_Filename'] )
try:
    print("\t", GnssUtil.Get_NowEpoch(), " Try to execute 'Code_Phase_Search %s --assist %s %s' " % ( CodePhaseSearchConfig['Default_Remote_Host_IP_Address'], CodePhaseSearchConfig['Default_Assistance_Filename'], code_phase_search_parameter ) )
    epoch   = GnssUtil.Get_NowEpoch()
    status  = subprocess.call( "./Code_Phase_Search", shell=True )
except subprocess.CalledProcessError as err:
    print (err)

print( "\n[STEP 5] ", GnssUtil.Get_NowEpoch(), " Parse and Preprocess The Satellite Signal Detection Results " )

print( "\tparse the GPS signal estimates from ", satelliteSignalEstimatesFilename )
with open( satelliteSignalEstimatesFilename, 'r' ) as f:
    satelliteSignalEstimates = json.load( f )
    for estimate in satelliteSignalEstimates:
        estimate["weekNumber"]  = epoch[0]
        estimate["ToA_s"]       = epoch[1]

with open( satelliteSignalEstimatesFilename, 'w' ) as f:
    json.dump( satelliteSignalEstimates, f )

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
jj = obsNumVisibles - 1 if obsNumVisibles > 7 else obsNumVisibles
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
eRcvXyzt, dRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
print ( "\t[Initial1] With ", epoch_ , " %d SVs, the estimate is " % jj, eRcvXyzt, " and the error is %.2f m, norm(dRcvXyzt) = %f " %( obsError, sqrt( dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2] ) ) ) 

from scipy.optimize import minimize
def XyzPositioningError( kk_ ):
    global updatedEphemeris, jj, epoch, obsVisibleSats, obsCodeShifts, rcvXyzt
    rcvXyzt_, dRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris,  [ epoch[0], epoch[1]+kk_ ], obsVisibleSats[0:jj ], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
    obsErr      = sqrt( (rcvXyzt[0]-rcvXyzt_[0])**2 + (rcvXyzt[1]-rcvXyzt_[1])**2 + (rcvXyzt[2]-rcvXyzt_[2])**2 )
    return obsErr

ii = 0
res = minimize( XyzPositioningError, ii, method='nelder-mead', options={'xtol': 1e-10, 'disp': False})
print( "\t[Initial2] jj = %d, res.x = %f, res.fun = %f " % ( jj, res.x, res.fun ) )

obsError_, ii_, jj_ = 10000, 0, obsNumVisibles

ii1 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
for ii in range( - 100, 101, 1 ):
    kk = ii1 + ii/10
    epoch_ = [ epoch[0], epoch[1]+kk ]
    rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
    updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
    if not os.path.exists( updatedEphemeris ) :
        updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
        print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
    for jj in range( max( 6, jj_-3 ), min( jj_+3, obsNumVisibles ), 1 ):
        eRcvXyzt, dRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 1000 ):
            print ( "\t\tWith ToA = %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m, norm(dRcvXyzt) = %f" % (obsError, sqrt( dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2] ) ) )
print ( "\t[Intermediate] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

ii2 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
for ii in range( - 100, 101, 1 ):
    kk = ii2 + ii*0.001
    epoch_ = [ epoch[0], epoch[1]+kk ]
    for jj in range( max( 6, jj_-2 ), min( jj_+2, obsNumVisibles ), 1 ):
        eRcvXyzt, dRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_ , obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 100 ):
            print ( "\t\tWith %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m, norm(dRcvXyzt) = %f" % ( obsError, sqrt( dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2] ) ) )
print ( "\t[Precious] With ", [ epoch[0], epoch[1]+ii_ ], " %d SVs, the estimate is " % jj_, eRcvXyzt, " and the error is %.2f m. " % obsError_ )

ii3 = ii_
epoch_ = [ epoch[0], epoch[1]+ii_ ]
rcvEpoch = GnssUtil.Epoch2Gps( epoch_ )
updatedEphemeris = "brdc_%4d%03d%02d%02d_L1C_GE.json" % ( rcvEpoch.tm_year, rcvEpoch.tm_yday, rcvEpoch.tm_hour, ( rcvEpoch.tm_min // 30 ) * 30 )
if not os.path.exists( updatedEphemeris ) :
    updatedEphemeris = gpsEph_nasa.Search_GPS_Ephemeris( epoch_, ephemerisFilename = updatedEphemeris )
    print("\tThe updated broadcast ephemeris file is ", updatedEphemeris )
for ii in range( - 100, 101, 1 ):
    kk = ii3 + ii*0.00001
    epoch_ = [ epoch[0], epoch[1]+kk ]
    for jj in range( max( 6, jj_-1 ), min( jj_+1, obsNumVisibles ), 1 ):
        eRcvXyzt, dRcvXyzt = GnssUtil.Get_GPS_Receiver_XyzPosition( updatedEphemeris, epoch_, obsVisibleSats[0:jj], obsCodeShifts[0:jj], rcvXyzt[0:3], __DEBUG = True )
        obsError = sqrt( (rcvXyzt[0]-eRcvXyzt[0])**2 + (rcvXyzt[1]-eRcvXyzt[1])**2 + (rcvXyzt[2]-eRcvXyzt[2])**2 )
        if( obsError < obsError_ ):
            obsError_, ii_, jj_ = obsError, kk, jj
        if( obsError < 10 ):
            print ( "\t\tWith %d s & %d visibles, the estimate is " % (epoch[1]+kk, jj), eRcvXyzt, " and the error is %.2f m, norm(dRcvXyzt) = %f" % ( obsError, sqrt( dRcvXyzt[0]*dRcvXyzt[0] + dRcvXyzt[1]*dRcvXyzt[1] + dRcvXyzt[2]*dRcvXyzt[2] ) ) )
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