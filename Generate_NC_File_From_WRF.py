#! /usr/bin/env python

"""

--------
Details:
--------

    This code post processes WRF files as followings:

    1. Converts all staggered variables to normal variables on the theta grid. These variables include:

       U['t','alt','lat','stagerred-lon']          to    U['t','alt','lat','lon']
       V['t','alt','staggered-lat','lon']          to    V['t','alt','lat','lon']
       W['t','stagerred-alt','lat','lon']          to    W['t','alt','lat','lon']
       Longitude['t','lat','staggered-lon']        to    Longitude['lat','lon']
       Latitude['t','staggered-lat','lon']         to    Latitude['lat','lon']
       Altitude['t','staggered-alt','lat','lon']   to    Altitude['t','alt','lat','lon']

       Note: Longitude and Latitude does not change over time. So we only output the Latitude and Longitude at t=0
             and squeeze the array so that they are only a function of ['lat','lon']
             However, Altitudes change over time hence we output them as Altitude['t','alt,','lat','lon']

    2. The Times variable is Nx19 array where each element is a single character. This code parses them are converts 
       them to netCDF num times, with respect to 'gregorian' calendar starting from 'days since 1970-01-01 00:00:00 UTC'.

    3. The output file is saved according to CF convension v37 available at: http://cfconventions.org/
       The CF convension can be checked online at: http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl

-----------
Input File:
-----------

    The input file should be netCDF nc file from WRF model that has the following dimensions and variabels:

    Dimensions:
    
    t
    alt
    staggered-alt         The length is alt+1
    lat
    staggered-lat         The length is lat+1
    lon
    staggered-lon         The length is lon+1

    Variables:

    Times     [NumTimes,19]                        Each row are 19 characters
    XLONG_U   ['t','lat','staggered-lon']          Longitudes on staggered grid of U               (not a function of 't' indeed)
    XLAT_U    ['t','lat','staggered-lon']          Latitudes on staggered grid of U                (not a function of 't' indeed)
    XLONG_V   ['t','staggered-lat','lon']          Longitudes on staggered grid of V               (not a function of 't' indeed)
    XLAT_V    ['t','staggered-lat',lon]            Latitudes on staggered grid of V                (not a function of 't' indeed)
    PH        ['t','staggered-alt','lat','lon']    Potential energy
    PHB       ['t','staggered-alt','lat','lon']    Potential energy base                           (not a function of 't' indeed)
    U         ['t','alt','lat','staggered-lon']    U velocity where lon direction is staggered
    V         ['t','alt','ataggered-lat','lon']    V velocity where lat direction is staggered
    W         ['t','staggered-alt','lat','lon']    W velocity where alt direction is staggered

------------
Output File:
------------

    The output file is netCDF nc file that has the following dimensions and variables:

    Dimensions

    t
    alt
    lat
    lon

    Variables:

    Datetimes     [NumTimes]                      (Each row is a float number)
    Longitude     ['lat','lon']
    Latitude      ['lat','lon']
    Altitude      ['t','alt','lat','lon']         (function of 't')
    U             ['t','alt','lat','lon']
    V             ['t','alt','lat','lon']
    W             ['t','alt','lat','lon']

------
Usage:
------

    Preprocess -i VirginiaTech_Real_WRF_Data_Fall2016.nc -o VirginiaTech_Real_WRF_Data_Fall2016_preprocessed.nc


-------
Author:
-------

    Siavash Ameli
    University of California, Berkeley
    December 21, 2016

"""

# =======
# Imports
# =======

import netCDF4
import numpy
import datetime
import time
import sys
import getopt
import warnings
import os

# =============
# Compute Times
# =============

def ComputeTimes(InputFile,TimeUnit,TimeCalendar):
    """
    Converts date char format to datetime numeric format.
    This parses the times chars and converts them to date times.
    """

    Times = InputFile.variables['time']
    DaysList = []

    for i in range(Times.shape[0]):

        # Parse chars to integers
        #Year = int(Times[i,0] + Times[i,1] + Times[i,2] + Times[i,3])
        #Month = int(Times[i,5] + Times[i,6])
        ##Day = int(Times[i,8] + Times[i,9])
        #Hour = int(Times[i,11] + Times[i,12])
        #Minute = int(Times[i,14] + Times[i,15])
        #Second = int(Times[i,17] + Times[i,18])
        day = float(Times[i])/(60*60*24)
        # Create Day object
        #DaysList.append(datetime.datetime(Year,Month,Day,Hour,Minute,Second))
     #   DaysList.append(day)
    # Convert dates to numbers
    #DateTimes = netCDF4.date2num(DaysList,units=TimeUnit,calendar=TimeCalendar)

    return day

# ==================
# Compute Longitudes
# ==================

def ComputeLongitudes(InputFile):
    """
    Converts the staggered longitudes to longitudes on the middle grid.

    Note: The Longitudes in the WRF dataset are Longitude[t,y,x], however the Longitudes are
    not a function of t. Here we output Longitude[y,x] for when t=0.
    """

    StaggeredLongitudes = InputFile.variables['XLONG_U']
    TimeIndex = 0
    Longitudes = (StaggeredLongitudes[TimeIndex,:,0:-1] + StaggeredLongitudes[TimeIndex,:,1:]) / 2.0
    return Longitudes

# =================
# Compute Latitudes
# =================

def ComputeLatitudes(InputFile):
    """
    Converts the staggered latitudes to latitudes on the middle grid.

    Note: The Latitudes in the WRF dataset are Latitude[t,y,x], however the Latitudes are
    not a function of t. Here we output Latitutde[y,x] for when t=0.
    """

    StaggeredLatitudes = InputFile.variables['XLAT_V']
    Latitudes = (StaggeredLatitudes[0,0:-1,:] + StaggeredLatitudes[0,1:,:]) / 2.0
    return Latitudes

# =================
# Compute Altitudes
# =================

def ComputeAltitudes(InputFile):
    """
    Converts potential energy to height.

    Note: 
        - PHB[t,z,y,x] does not vary with time, but it is a function of [z,y,x]
        - PH[t,z,y,x] varies with time, hence is a function of time.

    As a result, the altitudes arrays change over time, hence Altitude[t,z,y,x]
    """

    # Stagerred potential energy (PH) and potential energy base (PHB)
    StagerredPH = InputFile.variables['PH']
    StagerredPHB = InputFile.variables['PHB']

    # Convert stagerred variables to normal variables
    PH = (StagerredPH[:,0:-1,:,:] + StagerredPH[:,1:,:,:]) / 2.0
    PHB = (StagerredPHB[:,0:-1,:,:] + StagerredPHB[:,1:,:,:]) / 2.0

    # Compute altitudes from potential energies
    g = 9.81
    Altitudes = (PH + PHB)/g
    
    # Check ascending order
    # for t in range(Altitudes.shape[0]):
    #     for j in range(Altitudes.shape[2]):
    #         for i in range(Altitudes.shape[3]):
    #             for k in range(Altitudes.shape[1]-1):
    #                 if Altitudes[t,k,j,i] > Altitudes[t,k+1,j,i]:
    #                     print("not ascending order: k:%d, %f,%f"%(k,Altitudes[t,k,j,i],Altitudes[t,k+1,j,i]))

    return Altitudes

# =========
# Compute U
# =========

def ComputeU(InputFile):
    """
    Averages U in the middle points of x axis.
    """

    StaggeredU = InputFile.variables['U']
    U = (StaggeredU[:,:,:,0:-1] + StaggeredU[:,:,:,1:]) / 2.0
    return U

# =========
# Compute V
# =========

def ComputeV(InputFile):
    """
    Averages V in the middle points of y axis.
    """

    StaggeredV = InputFile.variables['V']
    V = (StaggeredV[:,:,0:-1,:] + StaggeredV[:,:,1:,:]) / 2.0
    return V

# =========
# Compute W
# =========

def ComputeW(InputFile):
    """
    Averages W in the middle points of z axis.
    """

    StaggeredW = InputFile.variables['W']
    W = (StaggeredW[:,0:-1,:,:] + StaggeredW[:,1:,:,:]) / 2.0
    return W

# ===============
# Parse Arguments
# ===============

def ParseArguments(argv):
    """
    Parses the arguments and set the input/output filenames.
    """

    def PrintUsage(ExecName):
        print("Usage: " + ExecName + " -i <InputFile.nc> -o <OutputFile.nc>")

    # variables
    InputFilename = ''
    OutputFilename = ''

    # Read arguments
    try:
        opts,args = getopt.getopt(argv[1:],"hvi:o:",["input=","output="])
    except getopt.GetoptError:
        PrintUsage(argv[0])
        sys.exit(2)

    # Assign options
    for opt,arg in opts:
        
        if opt == '-h':
            PrintUsage(argv[0])
            sys.exit()
        elif opt == "-v":
            print("version 0.0.1")
            print("Author: Siavash Ameli")
            sys.exit()
        elif opt in ("-i","--input"):
            InputFilename = arg
        elif opt in ("-o","--output"):
            OutputFilename = arg

    # Check
    if (InputFilename == '') or (OutputFilename == ''):
            PrintUsage(argv[0])
            sys.exit(2)

    return InputFilename,OutputFilename

def fixfill(n2Darray, fillvalue):
    dim = n2Darray.shape
    print dim
    for i in range(dim[1]):
        for j in range(dim[2]):
            if n2Darray[0,i,j] > fillvalue:
                n2Darray[0,i,j] = fillvalue
                
    return n2Darray
# ====
# Main
# ====

#def main(argv):
def main(InputFilename,OutputFilename):
    
    """
    Main function.
    This function makes a wrf file to be cf complaint, and also converts the staggered grid/variables to
    the middle grid/variables.
    """
    TimeSize = 0
    'Initialize an INFINITE loop to read all of the files in the sequence'
    while 1: 
        filename = '%s%02d%s' % (InputFilename[:21], TimeSize,InputFilename[-9:])
        'Check that the file actually exists, if it doesnt break the loop'
        if not(os.path.isfile(filename)):
            break        
        TimeSize+=1
        
    InputFile = netCDF4.Dataset(InputFilename)   
    # Dimensions size
    LatDimSize = len(InputFile.dimensions['y'])
    LonDimSize = len(InputFile.dimensions['x'])
    #Coordinates
    Lat = InputFile.variables['latitude'][:]
    Latitude = Lat[:,1]
    del Lat
    Long = InputFile.variables['longitude'][:]
    Longitude = Long[1,:]
    InputFile.close()
    
    Eastvel = numpy.zeros([TimeSize,LatDimSize,LonDimSize])    
    Northvel = numpy.zeros([TimeSize,LatDimSize,LonDimSize])
    DateTime = numpy.zeros([TimeSize])
    TimeUnit = 'days since 1970-01-01 00:00:00 UTC'
    TimeCalendar ='gregorian'
    for tt in range(0,TimeSize):
        filename = '%s%02d%s' % (InputFilename[:21], tt,InputFilename[-9:])
        #print filename
        InputFile = netCDF4.Dataset(filename)
        DateTime[tt] = ComputeTimes(InputFile,TimeUnit,TimeCalendar)
        #Eastvel[tt,:,:] = fixfill(InputFile.variables['UGRD_10maboveground'][:],9.999e+20)
        #Northvel[tt,:,:] = fixfill(InputFile.variables['VGRD_10maboveground'][:],9.999e+20)
        Eastvel[tt,:,:] = InputFile.variables['UGRD_10maboveground'][:,:]
        Northvel[tt,:,:] = InputFile.variables['VGRD_10maboveground'][:,:]
        InputFile.close()
        
    # Filenames
    #InputFilename,OutputFilename = ParseArguments(argv)
    
    # Files
    #InputFile = netCDF4.Dataset(InputFilename)
    OutputFile = netCDF4.Dataset(OutputFilename,'w',format='NETCDF4_CLASSIC')
    # Dimensions size
#    LatDimSize = len(InputFile.dimensions['y'])
#    LonDimSize = len(InputFile.dimensions['x'])

    # Dimensions
    OutputFile.createDimension('time',None)
    OutputFile.createDimension('lat',LatDimSize)
    OutputFile.createDimension('lon',LonDimSize)

    # Datetime
    OutputDatetime = OutputFile.createVariable('time',numpy.dtype('float64').char,('time',))
    OutputDatetime[:] = DateTime
    #OutputDatetime[:] = InputFile.variables['time']86400
    OutputDatetime.units = TimeUnit
    OutputDatetime.calendar = TimeCalendar
    OutputDatetime.standard_name = 'time'
    OutputDatetime._CoordinateAxisType = "Time"
    OutputDatetime.axis = "T"

    # Longitude
    OutputLongitude = OutputFile.createVariable('lon',numpy.dtype('float64').char,('lon',),zlib=True)
    #OutputLongitude[:] = InputFile.variables['longitude'][:]
    OutputLongitude[:] = Longitude
    OutputLongitude.units = 'degree_east'
    OutputLongitude.standard_name = 'longitude'
    OutputLongitude.positive = 'east'
    OutputLongitude._CoordinateAxisType = "Lon"
    OutputLongitude.axis = "X"
    OutputLongitude.coordsys = "geographic"

    # Latitude
    OutputLatitude = OutputFile.createVariable('lat',numpy.dtype('float64').char,('lat',),zlib=True)
    #OutputLatitude[:] = InputFile.variables['latitude'][:]
    OutputLatitude[:] = Latitude
    OutputLatitude.units = 'degree_north'
    OutputLatitude.standard_name = 'latitude'
    OutputLatitude.positive = 'up'
    OutputLatitude._CoordinateAxisType = "Lat"
    OutputLatitude.axis = "Y"
    OutputLatitude.coordsys = "geographic"

    # Altitude
    """
    OutputAltitude = OutputFile.createVariable('alt',numpy.dtype('float64').char,('time','alt','lat','lon',))
    OutputAltitude[:] = ComputeAltitudes(InputFile)
    OutputAltitude.units = 'm'
    OutputAltitude.standard_name = 'altitude'
    OutputAltitude.positive = 'up'
    OutputAltitude._CoordinateAxisType = "Height"
    OutputAltitude.axis = "Z"
    OutputAltitude._CoordinateZisPositive = "up"
    OutputAltitude.coordsys = "geographic"
    """
    # Velocity U
    OutputU = OutputFile.createVariable('U',numpy.dtype('float64').char,('time','lat','lon',),zlib=True,fill_value=9.999e+20)
    OutputU[:] = Eastvel
    OutputU.units = 'm s-1'
    OutputU.standard_name = 'surface_eastward_sea_water_velocity'#'eastward_wind'
    OutputU.coordinates = 'Longitude Latitude'
    OutputU.positive = 'toward east'
    OutputU.coordsys = "geographic"

    # Velocity V
    OutputV = OutputFile.createVariable('V',numpy.dtype('float64').char,('time','lat','lon',),zlib=True,fill_value=9.999e+20)
    OutputV[:] = Northvel
    OutputV.units = 'm s-1'
    OutputV.standard_name = 'surface_northward_sea_water_velocity'#'northward_wind'
    OutputV.coordinates = 'Longitude Latitude'
    OutputV.positive = 'toward north'
    OutputV.coordsys = "geographic"
   
    # Velocity W
    """
    OutputW = OutputFile.createVariable('W',numpy.dtype('float64').char,('time','alt','lat','lon',),zlib=True)
    OutputW[:] = ComputeW(InputFile)
    OutputW.units = 'm s-1'
    OutputW.standard_name = 'upward_air_velocity'
    OutputW.coordinates = 'Longitude Latitude Altitude Datetime'
    OutputW.positive = 'upward'
    OutputW.coordsys = "geographic"
    """
    # Global Attributes
    OutputFile.Conventions = 'CF-1.6'
    OutputFile.COORD_SYSTEM = 'GEOGRAPHIC'
    OutputFile.title = 'cf-compliance wrf file'
    OutputFile.creator_name = "Peter Nolan"
    OutputFile.creator_email = "pnolan86@vt.edu"
    OutputFile.institution = 'Virginia Polytechnical Institute and State University'
    OutputFile.creation_date = time.strftime("%x")

#    OutputFile.contributor_name = 'Siavash Ameli'
#    OutputFile.contributor_email = 'sameli@berkeley.edu'
#    OutputFile.contributor_role = 'Post process data to fill missing points.'
#    OutputFile.institution = 'University of California, Berkeley'
#    OutputFile.date_modified = time.strftime("%x")
    OutputFile.title = 'NCEP Atmospheric Data.'
    OutputFile.source = 'NCEP model data'
    OutputFile.references = 'not available'
    OutputFile.history = 'Original data: The data was generated by WRF software. Modified data: processed by Siavash Ameli using python script.'
#    OutputFile.summary = """The WRF data is not CF complaince. This dataset is converted to cf-compliance dataset, the datetime is parsed and
#            converted to num times. The staggered grid is converted to theta grid. The vertical axis is computed and added to variabes. Since the
#            vertical axis is time variable. the velocity data are probled at fixed grid locations using interplation."""
    OutputFile.project = 'Advanced Lagrangian Predictions for Hazards Assessments (NSF-ALPHA)'
    OutputFile.acknowledgement = 'This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant No. 1520825.'
    OutputFile.geospatial_lat_min = "%f"%(numpy.min(OutputLatitude[:]))
    OutputFile.geospatial_lat_max = "%f"%(numpy.max(OutputLatitude[:]))
    OutputFile.geospatial_lat_units = 'degree_north'
    OutputFile.geospatial_lon_min = "%f"%(numpy.min(OutputLongitude[:]))
    OutputFile.geospatial_lon_max = "%f"%(numpy.max(OutputLongitude[:]))
    OutputFile.geospatial_lon_units = 'degree_east'
    #OutputFile.geospatial_vertical_min = "%f"%(numpy.min(OutputAltitude[:]))
    #OutputFile.geospatial_vertical_max = "%f"%(numpy.max(OutputAltitude[:]))
    #OutputFile.geospatial_vertial_units = 'm'
    #OutputFile.geospatial_vertial_positive = 'up'
    OutputFile.time_coverage_start = "%s"%(netCDF4.num2date(OutputDatetime[0],units=OutputDatetime.units,calendar=OutputDatetime.calendar))
    OutputFile.time_coverage_end = "%s"%(netCDF4.num2date(OutputDatetime[-1],units=OutputDatetime.units,calendar=OutputDatetime.calendar))
    OutputFile.cdm_data_type = 'grid'
    
    # Close streams
    #InputFile.close()
    OutputFile.close()

# ===========
# System Main
# ===========

#if __name__ == '__main__':
#    main(sys.argv)
main('hiresw.t00z.arw_5km.f00.conus.nc','Ready_4_Thredds.nc')
