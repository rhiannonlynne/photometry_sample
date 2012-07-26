# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
from copy import deepcopy
import numpy
import pylab
from Sed import Sed
from Bandpass import Bandpass

# These are some convenience functions to access Bandpass and Sed with some of the 'bulk'
# functionality that seemed like it would be useful in engineering tests.
# Please feel free to add to this or to send me suggestions. 

# Note that the way these functions are set up is not optimal for large numbers of
# SEDs, but for the few which are stored within the photometry_sample directories, this works well.

figformat = 'png'

EXPTIME = 15                      # Default exposure time. (option for method calls).                   
NEXP = 2                          # Default number of exposures. (option for methods).
EFFAREA = numpy.pi*(6.5*100/2.0)**2   # Default effective area of primary mirror. (option for methods).
GAIN = 2.3                        # Default gain. (option for method call).
RDNOISE = 5                       # Default value - readnoise electrons or adu per pixel (per exposure)
DARKCURRENT = 0.2                 # Default value - dark current electrons or adu per pixel per second
OTHERNOISE = 4.69                 # Default value - other noise electrons or adu per pixel per exposure
PLATESCALE = 0.2                  # Default value - "/pixel
SEEING = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}  # Default seeing values (in ")     


def readBandpasses(filterlist=('u', 'g', 'r', 'i', 'z', 'y'),
                   filterDir = None,
                   prefix = 'total_', suffix='.dat'):
    """Read a series of bandpass files, such as might be provided by
    the 'throughputs' package, creating a dictionary of Bandpass objects."""
    # If rootdir is "None", then check for some likely alternatives
    # (this just makes this function more convenient to use).
    # First check for environment variable LSST_THROUGHPUTS_DEFAULT
    #  (which is set by the throughputs package, or by 'phot_ref.csh'. 
    if filterDir == None:
        filterDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    # If this failed, then environment variable is not set.
    # Try defaulting to current directory (otherwise user should set this value).
    if filterDir == None:
        filterDir  = '.'
    # Now let's read the bandpass files, and store in a dictionary. 
    bandpassDict = {}
    for f in filterlist:
        # Instantiate the bandpass object in the dictionary.
        bandpassDict[f] = Bandpass()
        # Read the throughput file.
        filename = os.path.join(filterDir, prefix + f + suffix)
        bandpassDict[f].readThroughput(filename)
    # Return bandpass dictionary.
    return bandpassDict

def buildBandpasses(componentList_common, componentList_filter=None,
                    filterlist=('u', 'g', 'r', 'i', 'z', 'y'),
                    filterDir = None):
    """Build a set of bandpasses from individual files - combining
    each component by multiplying the transmission curves together.
    componentList_common contains elements which are common to all
    filters (such as the detector, lenses, mirrors) while
    componentList_filter contains elements which are specific to
    individual filters (and must specify the full name of the item,
    one for each filter in filterlist).
    If filterDir is specified, it is used as the root directory
    for all components. """
    if filterDir == None:
        filterDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    if filterDir == None:
        filterDir = '.'
    bandpassDict = {}
    # If componentList_filter has been set to None, then not using any filter-specific
    # components, but still want to generate bandpasses in all filters (using only common components).
    if componentList_filter == None:
        if componentList_common == None:
            raise Exception('Cannot set both componentList_filter and componentList_common to None.')
        inputfiles = componentList_common
        for f in filterlist:
            bandpassDict[f] = Bandpass()
            bandpassDict[f].readThroughputList(componentList = inputfiles,
                                               rootDir = filterDir)
    # Otherwise, have some filter specific values.
    # Loop over the filters and the matching items in the componentList_filter. 
    else:
        for f, c in zip(filterlist, componentList_filter):
            if componentList_common == None:
                inputfiles = [c,]
            else:
                inputfiles = componentList_common + [c,]
            bandpassDict[f] = Bandpass()
            bandpassDict[f].readThroughputList(componentList = inputfiles,
                                               rootDir = filterDir)
    return bandpassDict


def multiplyBandpassDict(bandpassDict1, bandpassDict2):
    """Multiply two bandpass dictionaries together, filter by filter, returning a new dictionary -
    that is multiply bandpassDict1[key] * bandpassDict2[key]. """
    # This could be useful if you have read in one set of bandpass information,
    # then done some changes to that data in python (for example, shifting the filter curve in wavelength)
    # and then multiplying them all back together again to get the total bandpass throughput curves.
    # First, check the dictionary keys match (if they have the same filters):
    if set(bandpassDict1.keys()) != set(bandpassDict2.keys()):
        raise Exception('The two dictionaries must have the same keys')
    bandpassDict_new = {}
    for f in bandpassDict1.keys():
        wavelen, sb = bandpassDict1[f].multiplyThroughputs(bandpassDict2[f].wavelen, bandpassDict2[f].sb)
        bandpassDict_new[f] = Bandpass(wavelen=wavelen, sb=sb)
    return bandpassDict_new


def plotBandpasses(bandpassDict, titletext=None, newfig=True, savefig=False, addlegend=True,
                   linestyle='-', linewidth=2):
    """Plot the bandpass throughput curves. """
    # Generate a new figure, if desired. 
    if newfig: 
        pylab.figure()
    # Plot the bandpass curves. 
    for f in bandpassDict.keys():
        pylab.plot(bandpassDict[f].wavelen, bandpassDict[f].sb, marker="", linestyle=linestyle, 
                   linewidth=linewidth, label=f)
    # Only draw the legend if desired (many bandpassDicts plotted together could make the legend unwieldy).
    if addlegend:
        pylab.legend(numpoints=1, fancybox=True, shadow=True)
    # Limit wavelengths to the LSST range. 
    pylab.xlim(300, 1150)
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Throughput')
    if newfig:
        # Only add the grid if it's a new figure (otherwise, it toggles on/off).
        pylab.grid()
    # Add a plot title.
    if titletext != None:
        pylab.title(titletext)
    # Save the figure, if desired.
    if savefig:
        if titletext!=None:
            pylab.savefig('%s.%s' %(titletext, figformat), format=figformat)
        else:
            pylab.savefig('throughputs.%s' %(titletext, figformat), format=figformat)
    return
        

def readPhotSeds(sedDir=None):
    """Read all the seds provided by this package, storing them in a
    dictionary and saving lists of the SEDs of each type (so that
    they can be separated later if desired). """
    # The environment variable PHOT_REF_SEDS_DIR is set if the file
    #  'phot_ref.csh' is sourced. 
    if sedDir == None:
        sedDir = os.getenv('PHOT_REF_SEDS_DIR')
    # Again, if this was not set .. try defaulting to current directory.
    # Otherwise, the user should specify the base directory for the SEDs.
    if sedDir == None:
        sedDir = '.'
    # Set up the lists for the photometry reference SEDs.
    #  (I'm doing this by hand to keep these lists well-known over time.)
    sedlists = {}
    sedlists['quasar'] = ['quasar.dat',]
    sedlists['stars'] = ['km10_7250.fits_g45', 'km10_6500.fits_g45',
                         'km10_6000.fits_g45', 'km10_5250.fits_g45',
                         'km10_4500.fits_g45', 'm3.0Full.dat']
    sedlists['sn'] = ['sn1a_15.0.dat', 'sn1a_20.0.dat', 'sn1a_10.0.dat']
    sedlists['galaxies'] = ['Sa_template_norm.sed.dat', 'Sdm_template_norm.sed0.dat',
                            'Ell2_template_norm.sed.dat']
    sedlists['photoZ_outliers'] = ['xspec_172.sed.dat', 'xspec_173.sed.dat',
                                   'xspec_175.sed.dat', 'xspec_176.sed.dat',
                                   'xspec_90.sed.dat', 'xspec_91.sed.dat']
    # Let's go read the files.
    sedDict = {}
    # Loop through quasar, stars, galaxies, sn in the sedlist dictionary.
    for objtype in sedlists.keys():
        # And for each type of object, loop through and read each SED.
        for s in sedlists[objtype]:
            sedDict[s] = Sed()
            sedDict[s].readSED_flambda(os.path.join(sedDir, objtype, s))
    # Return the sed dictionary and the dictionary containing the lists of each
    #  type of object. 
    return sedDict, sedlists

def readAnySeds(inputfileList, sedDir=None):
    """Read the seds in a list and store them in a format appropriate for use with the 
    rest of these routines. sedDir (if set) can be the root directory for the files in the list."""
    # Set up the sedlists dictionary to hold the sed names. We'll store them all keyed under 'any'. 
    sedlists = {}
    sedlists['any'] = deepcopy(inputfileList)
    # If the root directory is set, add it to the input file names. 
    if sedDir != None:
        ifiles = []    
        for i in inputfileList:
            ifiles.append(os.path.join(sedDir, i))
        inputfileList = ifiles
    # Read the files.
    sedDict = {} 
    for filename, s in zip(inputfileList, sedlists['any']):
        sedDict[s] = Sed()
        sedDict[s].readSED_flambda(filename)
    # Return the sed dictionary and the dictionary containing the lists of each type of object.
    return sedDict, sedlists


def makeRedshiftedSeds(sedDict, sedlists, redshifts):
    """Redshift the quasar, galaxies and SN by the amounts given in redshifts.
    If redshift is a numpy array, all objects except stars are redshifted to the same values.
    If redshift is a dictionary containing numpy arrays in the 'quasar', 'sn' or 'galaxies' keys, then
    those redshifts are applied to those objects only. """
    # Check what kind of redshift object we received.
    if isinstance(redshifts, dict):
        # Loop over all the object types in the redshift dictionary. 
        for objtype in redshifts.keys():
            # Then loop over all redshifts for this type. 
            newSedsForList = []
            for z in redshifts[objtype]:
                # And add new SEDs redshifted to this z for this type of object. 
                for s in sedlists[objtype]:
                    # Make the new name for the new SED at this redshift. 
                    newsedname = s + '_Z_%.3f' %(z)
                    # Add it to the SED dictionary. 
                    sedDict[newsedname] = redshiftSingleSED(sedDict[s], z)
                    newSedsForList.append(newsedname)
            # Done redshifting all objects of this type; updated sedlist.
            sedlists[objtype] += newSedsForList
    else:
        # Using the same redshift list for everything (except stars).
        for z in redshifts:
            for objtype in sedlists.keys():
                newSedsForList = []
                if objtype == 'stars':
                    # Skip stars.
                    continue
                for s in sedlists[objtype]:
                    newsedname = s + '_Z_%.3f' %(z)
                    sedDicts[newsedname] = redshiftSingleSED(sedDict[s], z)
                    newSedsForList.append(newsedname)
                sedlists[objtype] += newSedsForList
    return sedDict, sedlists


def redshiftSingleSED(sed_in, z):
    # Make a copy of the input SED, as otherwise we will overwrite original arrays.
    sed_out = deepcopy(sed_in)
    sed_out.redshiftSED(z, dimming=True)
    return sed_out


def matchSedsBp(sedDict, bpDict, refFilter=None):
    """Match the wavelength ranges for all the Seds and the bandpass dictionary.
    This will speed up later calculation of magnitudes (if you're doing a lot of them),
     but is not strictly necessary. """
    if refFilter == None:
        refFilter = bpDict.keys()[0]
    wavelen_match = bpDict[refFilter].wavelen
    # Check all filters in bpDict match in wavelength space (note bandpasses must be regular grid).
    for f in bpDict.keys():
        if numpy.any(bpDict[f].wavelen != wavelen_match):
            bpDict[f].resampleBandpass(wavelen_min=wavelen_match.min(), wavelen_max=wavelen_match.max(),
                                       wavelen_step = wavelen_match[1] - wavelen_match[0])            
    # Check all seds in sedDict match in wavelength space. 
    for s in sedDict.keys():    
        if sedDict[s].needResample(wavelen_match=wavelen_match):
            sedDict[s].resampleSED(wavelen_match=wavelen_match)
    return sedDict, bpDict


def calcNatMags(bandpassDict, sedDict, sedlists):
    """Calculate (zeropoint-calibrated) natural magnitudes for each SED,
    in all bandpasses. Changes in this magnitude includes only color/wavelength-dependent effects. """
    # Create a dictionary to hold all of the magnitude information, for all filters.
    mags = {}
    # Loop over SEDs:
    for o in sedlists.keys():
        for s in sedlists[o]:
            # Create a dictionary for each object to hold the multiple filter information.
            mags[s] = {}
            for f in bandpassDict.keys():
                # Calculate the magnitudes. 
                mags[s][f] = sedDict[s].calcMag(bandpassDict[f])
    return mags

def calcInstMags(bandpassDict, sedDict, sedlists, 
                 expTime=EXPTIME, effarea=EFFAREA, gain=GAIN):
    """Calculate instrumental magnitudes for each SED, in all bandpasses. 
    Changes in this magnitude includes gray-scale effects as well as color/wavelength dependent effects. """
    mags = {}
    for o in sedlists.keys():
        for s in sedlists[o]:
            mags[s] = {}
            for f in bandpassDict.keys():
                mags[s][f] = sedDict[s].calcADU(bandpassDict[f], expTime=expTime, gain=gain, effarea=effarea)
                mags[s][f] = -2.5*numpy.log10(mags[s][f])
    return mags


def calcDeltaMags(mags1, mags2, mmags=True, matchBlue=False):
    """Calculate the difference in magnitudes between two sets of magnitudes, mags calculated as above.
    If 'mmags' is True, then returns delta mags in mmags. 
    If 'matchBlue' is True, then scales change in magnitude between mags1 and mags2 so that the blue 
     star 'km10_7250.fits_g45' has zero magnitude change."""
    dmags_seds = list(set(mags1) & set(mags2))
    s = mags1.keys()[0]
    dmags_filters = list(set(mags1[s]) & set(mags2[s]))
    dmags= {}
    for s in dmags_seds:
        dmags[s] = {}
        for f in dmags_filters:
            dmags[s][f] = mags1[s][f] - mags2[s][f]
            # Convert to millimags if desired.
            if mmags:
                dmags[s][f] *= 1000.0
    # Apply scaling so that bluest star remains constant, if desired. 
    if matchBlue:
        # Calculate offset to apply.
        offset = {}
        for f in dmags_filters:
            offset[f] = dmags['km10_7250.fits_g45'][f]
        # Apply offset.
        for s in dmags_seds:
            for f in dmags_filters:
                dmags[s][f] = dmags[s][f] - offset[f]
    return dmags

def calcGiColors(mags):
    """Calculate the g-i colors of objects in a set of magnitudes, mags calculated as above."""
    gi = {}
    for s in mags.keys():
        gi[s] = mags[s]['g'] - mags[s]['i']
    return gi


def printDmags(sedlists, dmags, filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Print changes in magnitudes to the screen."""
    print 'Delta mmag:'
    writestring = "object"
    for f in filterlist:
        writestring += '\t %s ' %(f)
    print writestring
    for objtype in sedlists.keys():
        print 'Object type: ', objtype
        for s in sedlists[objtype]:
            writestring = 'dm %s ' %(s)
            for f in filterlist:
                writestring += ' %f ' %(dmags[s][f])
            print writestring
    return


def plotDmags(sedlists, gi, dmags, newfig=True, titletext=None, savefig=False,
              filterlist=('u', 'g', 'r', 'i', 'z', 'y')):
    """Generate a plot of the change in magnitudes. """
    symbs = {'quasar':'o', 'stars':'s', 'sn':'x', 'galaxies':'+', 'photoZ_outliers':'*', 'any':'^'}
    colors = {'quasar':'g', 'stars':'k', 'sn':'b', 'galaxies':'r', 'photoZ_outliers':'m', 'any':'k'}
    if newfig:
        pylab.figure()
    pylab.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
    for f, i in zip(filterlist, range(1, len(filterlist)+1)):
        pylab.subplot(3,2,i)
        for objtype in sedlists.keys():
            for s in sedlists[objtype]:
                pylab.plot(gi[s], dmags[s][f], color=colors[objtype], marker=symbs[objtype])
        pylab.xlabel('g-i')
        pylab.ylabel(r'$\Delta$%s (mmag)' %(f))
    pylab.suptitle(titletext)
    if savefig:
        if titletext != None:
            pylab.savefig('%s.%s' %(titletext, figformat), format=figformat)
        else:
            pylab.savefig('Dmag.%s' %(figformat), format=figformat)
    return
