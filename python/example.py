# This is an example program to show how to use some of the functionality in photTools.py (and a small part of the underlying
# capability of Bandpass and Sed).
# It can be run by entering 'python example.py' at the command line. 

import numpy
import pylab
import photTools as pT

print "This is an example script to show uses of photTools, and also to illustrate some differences in how"
print " we calculate magnitudes. "
print "Instrumental magnitudes include both grey-scale and color-dependent variations as bandpasses and SEDs change."
print "Calibrated natural magnitudes remove the grey-scale variation, and only show color-dependent/wavelength-dependent"
print " variations in measured magnitudes as the shape of the bandpass or the shape of the SED change. "
print "If the blue stars are used as calibration stars for the instrumental magnitudes (i.e. the zeropoint variation is"
print " removed with the blue stars, so the grey-scale variation is removed), the changes in instrumental magnitudes then"
print " show only the color-dependent changes in reported magnitudes. (this can be done here using matchBlue=True in calcDeltaMags)."
print "The reason that the natural magnitudes show a small change in magnitude even for blue objects is that these are being"
print " referenced against a flat Fnu spectrum, instead of the blue star spectrum. If matchBlue is turned on for changes in magnitudes"
print "  calculated using natural magnitudes (calcNatMags), the results are identical to the results for changes in instrumental"
print "  magnitudes calculated using matchBlue."


filterlist = ('u', 'g', 'r', 'i', 'z', 'y')

# Read the SED files.
seds, sedlists = pT.readPhotSeds(sedDir='../')

# And redshift the galaxies, quasar, and SN.
# These are just sample redshifts, not supposed to be 'typical'. 
redshifts = {}
redshifts['galaxies'] = numpy.array([0.5, 1.0], 'float')
redshifts['quasar'] = numpy.array([1.0, 1.5, 2.5], 'float')
redshifts['sn'] = numpy.array([0.3, 0.8, 1.2, 1.5], 'float')
redshifts['photoZ_outliers'] = numpy.array([0, 0.2, 2.0], 'float')
seds, sedlists = pT.makeRedshiftedSeds(seds, sedlists, redshifts)

# Okay, now let's look at two different bandpass comparisons. 
# We will compare the 'normal' base throughputs to a set of throughputs where the bandpasses are shifted by 1%.
# First read the basic components, except for the filters.
componentList_common = ['atmos.dat', 'detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat',
                        'm1_ProtAl_Aged.dat', 'm2_ProtAl_Aged.dat', 'm3_ProtAl_Aged.dat']
lsst_base = pT.buildBandpasses(filterlist=filterlist, componentList_common=componentList_common,
                                componentList_filter = None, filterDir='../throughputs/baseline')
# Now let's read the filters (only) from disk. 
componentList_common = None
componentList_filter = ['filter_u.dat', 'filter_g.dat', 'filter_r.dat', 'filter_i.dat',
                        'filter_z.dat', 'filter_y.dat']
lsst_filters = pT.buildBandpasses(filterlist=filterlist, componentList_common=None,
                                  componentList_filter = componentList_filter, filterDir='../throughputs/baseline')

# Now multiply the base with the filters, for the 'standard' ('normal') set. 
lsst_std = pT.multiplyBandpassDict(lsst_base, lsst_filters)

# And also generate a 'shifted' version of the filters, shifting by 1% times the effective wavelength of each filter.
for f in filterlist:
    # Calculate the effective wavelength of this filter (filter alone). 
    effphi, eff_wavelen = lsst_filters[f].calcEffWavelen()
    # Calculate the shift.
    shift = 0.01 * eff_wavelen
    # Apply the shift to the wavelengths of this filter. 
    lsst_filters[f].wavelen = lsst_filters[f].wavelen + shift
    # Resynchronize and rebin the wavelengths / transmission curve of the filter.
    lsst_filters[f].resampleBandpass()

# Then apply these shifted filters to the same base throughput set, to get our 'shifted' bandpasses.
lsst_shifted = pT.multiplyBandpassDict(lsst_base, lsst_filters)

# Plot the transmission curves. 
pT.plotBandpasses(lsst_std, newfig=True, addlegend=True)
pT.plotBandpasses(lsst_shifted, newfig=False, addlegend=False, titletext='Shifted and non-shifted LSST Bandpasses', linestyle=':')

# Calculate the magnitudes in each of these bandpasses, first including all effects (greyscale and color-dependent).
mags_std = pT.calcInstMags(lsst_std, seds, sedlists)
mags_shifted = pT.calcInstMags(lsst_shifted, seds, sedlists)
# And calculate the differences, matching the blue stars in each bandpass.
dmags = pT.calcDeltaMags(mags_std, mags_shifted, mmags=True, matchBlue=False)

print_to_screen = False
if print_to_screen:
    pT.printDmags(sedlists, dmags)
# And generate a plot.
gi = pT.calcGiColors(mags_std)
pT.plotDmags(sedlists, gi, dmags, titletext ='Instrumental Mags')

# And then let's do the same thing, but calculate the delta mags while matching the magnitudes of the blue stars
# in each bandpass. (this is essentially removing the grey-scale effects, but comparing the color-dependent effects
# compared to a blue stellar SED, since those are the objects which provide the zeropoint scale). 

# And calculate the differences, matching the blue stars in each bandpass.                                                             
dmags = pT.calcDeltaMags(mags_std, mags_shifted, mmags=True, matchBlue=True)

print_to_screen = False
if print_to_screen:
    pT.printDmags(sedlists, dmags)
# And generate a plot.                                                                                                                
pT.plotDmags(sedlists, gi, dmags, titletext='Instrumental Mags, Blue Matched')



# And let's do the same thing, only generate natural magnitudes (removing the grey-scale effect and leaving only
# color-dependent terms ... but note that this is now 'color-dependent' compared to a flat Fnu SED, rather than a blue SED). 

# Calculate the magnitudes in each of these bandpasses, first including all effects (greyscale and color-dependent).                  
mags_std = pT.calcNatMags(lsst_std, seds, sedlists)
mags_shifted = pT.calcNatMags(lsst_shifted, seds, sedlists)
# And calculate the differences, matching the blue stars in each bandpass.                                                         
dmags = pT.calcDeltaMags(mags_std, mags_shifted, mmags=True, matchBlue=False)

print_to_screen = False
if print_to_screen:
    pT.printDmags(sedlists, dmags)
# And generate a plot.            
gi = pT.calcGiColors(mags_std)
pT.plotDmags(sedlists, gi, dmags, titletext='Natural Magnitudes')

pylab.show()

# Note that if we repeated the natural magnitudes plot, after doing a scaling so that the blue stars were the same 
# magnitude in each bandpass, we would see the SAME result as when doing this with the InstMag values. 
