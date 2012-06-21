# This is another example script, although I did use it to generate all the SED plots in ../plots. 

import pylab
import photTools as pT

figformat = 'png'

def plot_sed(sed, sedname=None, newfig=True, savefig=True):
    """Plot a single SED - f_nu vs. lambda."""
    # These are the quantities that translate directly into magnitudes (when multiplied by phi, the 
    #   normalized system response curve = Sb / lambda). 
    if newfig:
        pylab.figure()
    # Just in case fnu calculation hasn't been triggered yet, let's do it here.
    sed.flambdaTofnu()
    # Plot the information from this SED. 
    if sedname != None:
        pylab.plot(sed.wavelen, sed.fnu, label='%s' %(str(sedname)))
    else:
        pylab.plot(sed.wavelen, sed.fnu)
    pylab.xlim(300, 1150)
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel(r'F$_\nu$')
    if sedname != None:
        pylab.title('%s' %(str(sedname)))
    if savefig:
        if sedname != None:
            pylab.savefig('%s.%s' %(str(sedname), figformat), format=figformat)
        else:
            pylab.savefig('fnu.%s' %(figformat), format=figformat)
    return
    

def plot_all_seds(sedDict):
    """Plot all seds in the sedDict."""
    for s in sedDict.keys():
        plot_sed(sedDict[s], sedname=s)
    return


if __name__ == "__main__":
    
    sedDict, sedlists = pT.readPhotSeds()
    plot_all_seds(sedDict)
    
