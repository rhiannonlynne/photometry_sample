# Example python script for calculating m5 values


# Import necessary modules & classes
import os
import numpy
from Bandpass import Bandpass
from Sed import Sed


# Get directory where 'throughputs' live from environment variable
#  (this env variable is set by phot_ref.csh to be the throughputs/baseline directory in photometry_sample)
throughputsDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')

# Read the individual throughput components for the system/hardware and for the total throughputs, for each filter. 
#   (note that m1.dat = m1_ProtAlIdeal.dat .. if you want the Aged aluminum, change the name in the components lists below).

filterlist = ('u', 'g', 'r', 'i', 'z', 'y')

hardware = {}
multifilter_hardware_components = ['detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat']
total = {}
multifilter_total_components = ['detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat', 'atmos_10.dat']
for f in filterlist:
    hardware[f] = Bandpass()
    hardware_components = multifilter_hardware_components + ['filter_' +f +'.dat',]
    hardware[f].readThroughputList(componentList=hardware_components, rootDir = throughputsDir)
    total[f] = Bandpass()
    total_components = multifilter_total_components + ['filter_' +f+'.dat',]
    total[f].readThroughputList(componentList=total_components, rootDir=throughputsDir)

# Calculate and print Tb (sum of total? throughput / lambda, integrated over lambda) and Sb (similar, but for hardware)
Sb = {}
Tb = {}
for f in filterlist:
    Sb[f] = (hardware[f].sb / hardware[f].wavelen).sum() * hardware[f].wavelen_step
    Tb[f] = (total[f].sb / total[f].wavelen).sum() * total[f].wavelen_step

writestring1 = 'Sb: '
writestring2 = 'Tb: ' 
for f in filterlist:
    writestring1 += '%s %.3f ' %(f, Sb[f])
    writestring2 += '%s %.3f ' %(f, Tb[f])
print writestring1
print writestring2

# Set up to calculate m5. 

# Read in the dark sky SED
darksky = Sed()
darksky.readSED_flambda(os.path.join(throughputsDir, 'darksky.dat'))

# Set up a range of exposure times (per exposure, not per visit)
exptimes = numpy.arange(0.1, 100., 5.)
# Set a range of readnoise values and zero out other noise contributions (I think this is what you want, to isolate readnoise completely)
#  (plus the camera team includes 'othernoise' into their 'instrumental noise' value .. this is another potential source of confusion when talking
#  to them about this .. usually when they say 'readnoise' they actually mean the full instrumental noise, but do not consider dark current)
othernoise = 0
darkcurrent = 0
readnoises = [10., 13.]  # noise per VISIT (so set nexp=1 and then exptime = visit time)
nexp = 1

# Calculate m5 values for each of these exposure times. 
for exptime in exptimes:
    for readnoise in readnoises:
        writestring = 'Exptime %.3f Nexp %d Instnoise %.1f --M5' %(exptime, nexp, readnoise)
        for f in filterlist:
            m5 = total[f].calcM5(darksky, hardware[f], expTime=exptime, nexp=nexp, readnoise=readnoise, 
                                 darkcurrent=darkcurrent, othernoise=othernoise, seeing=None, filter=f, gain=1.)
            writestring += ' %s %.3f' %(f, m5)
        print writestring


# Things you can change for calcM5 -- 
#        calcM5(self, skysed, hardware, expTime=EXPTIME, nexp=NEXP, readnoise=RDNOISE,
#               darkcurrent=DARKCURRENT, othernoise=OTHERNOISE,
#               seeing=None, platescale=PLATESCALE, 
#               gain=GAIN, effarea=EFFAREA, filter='r'):
#  .. if seeing==None, but filter is specified, then uses default value of seeing in that filter 
#     (i.e. SEEING = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63, 'y3':0.63, 'y4':0.63}  # Default seeing values (in ")
#  other default values are:
#EXPTIME = 15                      # Default exposure time. (option for method calls).
#NEXP = 2                          # Default number of exposures. (option for methods).
#EFFAREA = numpy.pi*(6.5*100/2.0)**2   # Default effective area of primary mirror. (option for methods).
#GAIN = 2.3                        # Default gain. (option for method call).
#RDNOISE = 5                       # Default value - readnoise electrons or adu per pixel (per exposure)
#DARKCURRENT = 0.2                 # Default value - dark current electrons or adu per pixel per second
#OTHERNOISE = 4.69                 # Default value - other noise electrons or adu per pixel per exposure
#PLATESCALE = 0.2                  # Default value - "/pixel
# ... EFFAREA needs to be updated, but I can't remember/find the actual value that should be used instead to account for overall 
#      vignetting & spider obscurations

