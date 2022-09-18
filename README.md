### README ###

This is a subset of SEDs which can be used to calculate magnitudes of some relatively widely
ranging and hopefully useful objects, in combination with [rubin_sim](https://github.com/lsst/rubin_sim).

The python package `rubin_sim` can be installed with conda from conda-forge or with pip by following the directions in this [README](https://github.com/lsst/rubin_sim). 

---

This photometry_reference directory contains 5 subdirectories containing a sample of reference SEDs, chosen 
to represent a range of different types of SEDs suitable for addresing
various calibration questions.  These subdirectories are: 
 - galaxies (composite galaxy SEDs from Andy Connolly)
 - photoZ_outliers (containing galaxy templates which are most likely to
        result in catastrophic failures in determining redshifts 
	using photo-Z techniques. Templates provided by Sam Schmidt, via Tony Tyson). 
 - quasar (containing a single quasar SED which has been generated 
	by combining a composite spectrum from 2204 QSO spectra 
	from SDSS (vanden Berk et al) with a simple model for the
	flux behavior beyond the SDSS wavelengths, provided by Zeljko Ivezic).
 - sn (containing SN SEDs at various points in their time evolution, 
	these come from the SN templates created by Peter Nugent),
 - stars (containing several stars from blue to red -- km10_7250.fits_g45 
  	being the bluest and equivalent to about an F0 type main sequence
	star [7250 indicates its effective temperature .. hotter=bluer],
	down to the m3.0Full.dat SED which represents a red Mdwarf.
	The Kurucz models come from Kurucz 1993, while the Mdwarf is
	a composite empirical spectra from SDSS extended in wavelength
	by other observations.)

In each SED file, there may or may not be a header line starting with '#', 
but the data is always contained in columns of wavelength (in nm) and 
flux - 'Flambda' specified, that is, in units of ergs/cm^s/s/nm, although
the flux is arbitrarily scaled (that is, the SEDs are not 
photometrically pre-scaled to any particular magnitude). 


The directory 'mag_info' contains some brief explanations of magnitudes as
can be calculated by `rubin_sim.photUtils`.
