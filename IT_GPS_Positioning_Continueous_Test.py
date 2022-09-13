import GnssUtil, Ephemeris, datetime, time, json, os, subprocess

ii = 0
while (ii < 1000):
    try:
        print("\n[ii=%04d]"%ii, datetime.datetime.now(), "Execute 'python3 IT_GPS_Positioning_RTL_SDR_11ms.py' "  )
        status  = subprocess.call( "python3 IT_GPS_Positioning_RTL_SDR_11ms.py", shell=True )
    except subprocess.CalledProcessError as err:
        print (err)

    try:
        print("\n[ii=%04d]"%ii, datetime.datetime.now(), "Execute 'python3 IT_GPS_Positioning_RTL_SDR_16ms.py' "  )
        status  = subprocess.call( "python3 IT_GPS_Positioning_RTL_SDR_16ms.py", shell=True )
    except subprocess.CalledProcessError as err:
        print (err)

    ii += 1