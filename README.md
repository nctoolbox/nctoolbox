NCTOOLBOX
=========

Download the latest release:
[You can download the latest stable release at by clicking here.](https://github.com/nctoolbox/nctoolbox/releases)


##Brief summary:

nctoolbox is a [Matlab](http://www.mathworks.com/) toolbox that provides read-only access to [common data model](http://www.unidata.ucar.edu/software/netcdf-java/CDM/index.html) datasets. Under the hood, nctoolbox uses [NetCDF-Java](http://www.unidata.ucar.edu/software/netcdf-java/) as the data access layer. This allows nctoolbox to access [NetCDF](http://www.unidata.ucar.edu/software/netcdf/), [OPeNDAP](http://opendap.org/), [HDF5](http://www.hdfgroup.org/HDF5/), GRIB, GRIB2, HDF4 and many (15+) other file formats and services using the same API. It works with Matlab 2008a and later.

###Prequisites

Matlab R2008a+.  You can verify the version of Matlab by typing:
    
    
    version


Java 6+.  You can verify the version of Java used by Matlab by typing: 
    

    version('-java'). 


The version returned should start with 'Java 1.6'. If it starts with 'Java 1.5' you can try updating the Matlab JVM: http://www.mathworks.com/support/solutions/en/data/1-1812J/ or use the older nctoolbox-20091112 version of this toolbox.

###Setup

In Matlab, change to the nctoolbox directory. For example,
 

    cd ~/Documents/MATLAB/nctoolbox


Run the setup_nctoolbox.m function
 

    setup_nctoolbox


This sets up nctoolbox for the current Matlab session only. You will need to add the follwing lines to your startup.m file if you would like nctoolbox automatically when you start Matlab:


    % Edit '/Path/to/nctoolbox' to correct nctoolbox directory
    addpath('/Path/To/nctoolbox')
    setup_nctoolbox
      
###Demos

  * Demos that display basic functionality are in the 'demos' subdirectory.  These demos
     can fail if the machines serving the remote data access URLs are unavailable. A gallery
     of these demos is [visible on github](http://nctoolbox.github.io/nctoolbox/demos.html).
  * Contributed demos that display additional or specialized functionality are in 
     the demos/contrib directory.  Some of these depend on accessing remote sites for
     data that can be less reliable than the data URLs in the 'demos' directory.

###Documentation

  * We are migrating the documentation from the [GoogleCode nctoolbox documentation](http://code.google.com/p/nctoolbox/wiki/Documentation) to the [Github nctoolbox documentation](https://github.com/nctoolbox/nctoolbox/wiki/Documentation).
  * Descriptions of the nctoolbox classes are available as [NctoolboxClasses](https://github.com/nctoolbox/nctoolbox/wiki/NctoolboxClasses).
  * This repository includes a number of useful [utility functions](https://github.com/nctoolbox/nctoolbox/wiki/NctoolboxUtilityFunctions).
  * The [wiki](https://github.com/nctoolbox/nctoolbox/wiki) includes [documentation](https://github.com/nctoolbox/nctoolbox/wiki/Documentation) on [using ncdataset](https://github.com/nctoolbox/nctoolbox/wiki/UsingNcdataset) and [using ncgeodataset](https://github.com/nctoolbox/nctoolbox/wiki/UsingNcgeodataset).
  * Several of the [demos](https://github.com/nctoolbox/nctoolbox/tree/master/demos) are rendered in the [demos gallery](http://nctoolbox.github.io/nctoolbox/demos.html).


