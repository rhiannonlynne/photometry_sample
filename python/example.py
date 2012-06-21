# This is an example program to show how to use some of the functionality in photTools.py (and a small part of the underlying
# capability of Bandpass and Sed).
# It can be run by entering 'python example.py' at the command line. 

import numpy
import pylab
import photTools as pT

filterlist = ('u', 'g', 'r', 'i', 'z', 'y')

# Read the bandpass files.
lsstBP = pT.readBandpasses(filterlist=filterlist)

# Read the SED files.
seds, sedlists = pT.readPhotSeds(sedDir='../')

# And redshift the galaxies, quasar, and SN
redshifts = {}
redshifts['galaxies'] = numpy.array([0.5, 1.0], 'float')
redshifts['quasar'] = numpy.array([1.0, 1.5, 2.5], 'float')
redshifts['sn'] = numpy.array([0.3, 0.8, 1.2, 1.5], 'float')
seds, sedlists = pT.makeRedshiftedSeds(seds, sedlists, redshifts)

# Calculate magnitudes for all the objects in this bandpass set.
mags = pT.calcMags(lsstBP, seds, sedlists)

########

# Okay, now let's do something a little cooler - let's look
#  at a comparison case, where we will replace the filter throughputs with a shifted version.
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

# Now multiply the base with the filters.
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

# Then apply these shifted filters to the same base throughput set.
lsst_shifted = pT.multiplyBandpassDict(lsst_base, lsst_filters)

# Calculate the magnitudes in each of these bandpasses.
mags_std = pT.calcMags(lsst_std, seds, sedlists)
mags_shifted = pT.calcMags(lsst_shifted, seds, sedlists)

# And look at the difference in magnitude between each type of bandpass.
# Note that the magnitudes calculated at this point are 'zero-point calibrated' natural magnitudes:
# that is, they include the zeropoint of the system so absorb any grey-scale variation and
# only leave the color or wavelength-dependent changes in the natural magnitude due to the fact that
# the bandpass itself has changed.

dmags = pT.calcDeltaMags(mags_std, mags_shifted, matchBlue=False)

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



# And generate a plot.

symbs = {'quasar':'o', 'stars':'s', 'sn':'x', 'galaxies':'+'}
colors = {'quasar':'g', 'stars':'k', 'sn':'b', 'galaxies':'r'}

gi = pT.calcGiColors(mags_std)

pylab.figure()
pylab.subplots_adjust(top=0.93, wspace=0.32, hspace=0.32, bottom=0.09, left=0.12, right=0.96)
for f, i in zip(filterlist, range(1, len(filterlist)+1)):
    pylab.subplot(3,2,i)
    for objtype in sedlists.keys():
        for s in sedlists[objtype]:        
            pylab.plot(gi[s], dmags[s][f], color=colors[objtype], marker=symbs[objtype])
    pylab.xlabel('g-i')
    pylab.ylabel(r'$\Delta$%s (mmag)' %(f))
pylab.show()

