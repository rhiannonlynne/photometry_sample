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

"""
sed -

Class data:
wavelen (nm)
flambda (ergs/cm^s/s/nm)
fnu (Jansky)
zp  (basically translates to units of fnu = -8.9 (if Janskys) or 48.6 (ergs/cm^2/s/hz))
the name of the sed file

It is important to note the units are NANOMETERS, not ANGSTROMS. It is possible to rig this so you can
use angstroms instead of nm, but you should know what you're doing and understand the wavelength grid
limits applied here and in Bandpass.py.

Methods:
 Because of how these methods will be applied for catalog generation, (taking one base SED and then
  applying various dust extinctions and redshifts), many of the methods will either work on,
  and update self, OR they can be given a set of lambda/flambda arrays and then will return
  new versions of these arrays. In general, the methods will not explicitly set flambda or fnu to
  something you (the user) did not specify - so, for example, when calculating magnitudes (which depend on
  a wavelength/fnu gridded to match the given bandpass) the wavelength and fnu used are temporary copies
  and the object itself is not changed.
 In general, the philosophy of Sed.py is to not define the wavelength grid for the object until necessary
  (so, not until needed for the magnitude calculation or resampleSED is called). At that time the min/max/step
  wavelengths or the bandpass wavelengths are used to define a new wavelength grid for the sed object.
 When considering whether to use the internal wavelen/flambda (self) values, versus input values:
  For consistency, anytime self.wavelen/flambda is used, it will be updated if the values are changed
  (except in the special case of calculating magnitudes), and if self.wavelen/flambda is updated,
  self.fnu will be set to None. This is because many operations are typically chained together
  which alter flambda -- so it is more efficient to wait and recalculate fnu at the end, plus it avoids possible
  de-synchronization errors (flambda reflecting the addition of dust while fnu does not, for example).
  If arrays are passed into a method, they will not be altered and the arrays which are returned will be
  allocated new memory.
 Another general philosophy for Sed.py is use separate methods for items which only need to be generated once
  for several objects (such as the dust A_x, b_x arrays). This allows the user to optimize their code for
  faster operation, depending on what their requirements are (see example_SedBandpass_star.py and
  exampleSedBandpass_galaxy for examples).

Method include:
  setSED / setFlatSED / readSED_flambda / readSED_fnu -- to input information into Sed wavelen/flambda.
  getSED_flambda / getSED_fnu -- to return wavelen / flambda or fnu to the user.
  clearSED -- set everything to 0.
  synchronizeSED -- to calculate wavelen/flambda/fnu on the desired grid and calculate fnu.
  _checkUseSelf/needResample -- not expected to be useful to the user, rather intended for internal use.
  resampleSED -- primarily internal use, but may be useful to user. Resamples SED onto specified grid.
  flambdaTofnu / fnuToflambda -- conversion methods, does not affect wavelen gridding.
  redshiftSED -- redshifts the SED, optionally adding dimmingx
  setupCCMab / addCCMDust -- separated into two components, so that a_x/b_x can be reused between SEDS
if the wavelength range and grid is the same for each SED (calculate a_x/b_x with setupCCMab).
  multiplySED -- multiply two SEDS together.
  calcADU / calcMag / calcFlux -- with a Bandpass, calculate the ADU/magnitude/flux of a SED.
  calcFluxNorm / multiplyFluxNorm -- handle fluxnorm parameters (from UW LSST database) properly.
     These methods are intended to give a user an easy way to scale an SED to match an expected magnitude.
  renormalizeSED  -- intended for rescaling SEDS to a common flambda or fnu level.
  writeSED -- keep a file record of your SED.
  calcSNR_psf / calcSNR_mag -- two methods to calculate the SNR of a SED. (_psf is more accurate, but
requires knowing the sky count backgrounds. _mag assumes you know the m5 already).
  calcMagError / calcAstrometricError -- calculated based on SciBook values, but just estimates.
  setPhiArray -- given a list of bandpasses, sets up the 2-d phiArray (for manyMagCalc) and dlambda value.
  manyMagCalc -- given 2-d phiArray and dlambda, this will return an array of magnitudes (in the same
order as the bandpasses) of this SED in each of those bandpasses.

"""

import warnings
import numpy
import gzip
from PhotometricDefaults import PhotometricDefaults

class Sed(object):
    """Class for holding and utilizing spectral energy distributions (SEDs)"""
    def __init__(self, wavelen=None, flambda=None, fnu=None, badval=numpy.NaN, name=None):
        """
        Initialize sed object by giving filename or lambda/flambda array.

        Note that this does *not* regrid flambda and leaves fnu undefined.
        """
        self.fnu = None
        self.wavelen = None
        self.flambda = None
        #self.zp = -8.9  # default units, Jansky.
        self.zp = -2.5*numpy.log10(3631)
        self.name = name
        self.badval = badval
        # If init was given data to initialize class, use it.
        if (wavelen is not None) and ((flambda is not None) or (fnu is not None)):
            if name is None:
                name = 'FromArray'
            self.setSED(wavelen, flambda=flambda, fnu=fnu, name=name)
        return

    ### Methods for getters and setters.

    def setSED(self, wavelen, flambda=None, fnu=None, name='FromArray'):
        """
        Populate wavelen/flambda fields in sed by giving lambda/flambda or lambda/fnu array.

        If flambda present, this overrides fnu. Method sets fnu=None unless only fnu is given.
        Sets wavelen/flambda or wavelen/flambda/fnu over wavelength array given.
        """
        # Check wavelen array for type matches.
        if isinstance(wavelen, numpy.ndarray)==False:
            raise ValueError("Wavelength must be a numpy array")
        # Wavelen type ok - make new copy of data for self.
        self.wavelen = numpy.copy(wavelen)
        self.flambda=None
        self.fnu=None
        # Check if given flambda or fnu.
        if flambda is not None:
            # Check flambda data type and length.
            if (isinstance(flambda, numpy.ndarray)==False) or (len(flambda)!=len(self.wavelen)):
                raise ValueError("Flambda must be a numpy array of same length as Wavelen.")
            # Flambda ok, make a new copy of data for self.
            self.flambda = numpy.copy(flambda)
        else:
            # Were passed fnu instead : check fnu data type and length.
            if fnu is None:
                raise ValueError("Both fnu and flambda are 'None', cannot set the SED.")
            elif (isinstance(fnu, numpy.ndarray)==False) or (len(fnu)!=len(self.wavelen)):
                raise ValueError("(No Flambda) - Fnu must be numpy array of same length as Wavelen.")
            # Convert fnu to flambda.
            self.wavelen, self.flambda = self.fnuToflambda(wavelen, fnu)
        self.name = name
        return

    def setFlatSED(self, wavelen_min=PhotometricDefaults.minwavelen, wavelen_max=PhotometricDefaults.maxwavelen,
                   wavelen_step=PhotometricDefaults.wavelenstep, name='Flat'):
        """
        Populate the wavelength/flambda/fnu fields in sed according to a flat fnu source.
        """
        self.wavelen = numpy.arange(wavelen_min, wavelen_max+wavelen_step, wavelen_step, dtype='float')
        self.fnu = numpy.ones(len(self.wavelen), dtype='float') * 3631 #jansky
        self.fnuToflambda()
        self.name = name
        return

    def readSED_flambda(self, filename, name=None):
        """
        Read a file containing [lambda Flambda] (lambda in nm) (Flambda erg/cm^2/s/nm).

        Does not resample wavelen/flambda onto grid; leave fnu=None.
        """
        # Try to open data file.
        # ASSUME that if filename ends with '.gz' that the file is gzipped. Otherwise, regular file.
        try:
            if filename.endswith('.gz'):
                f = gzip.open(filename, 'r')
            else:
                f = open(filename, 'r')
        #if the above fails, look for the file with and without the gz
        except IOError:
            try:
                if filename.endswith(".gz"):
                    f = open(filename[:-3], 'r')
                else:
                    f = gzip.open(filename+".gz", 'r')
            except IOError:
                raise IOError("The sed file %s does not exist" %(filename))
        # Read source SED from file - lambda, flambda should be first two columns in the file.
        # lambda should be in nm and flambda should be in ergs/cm2/s/nm
        sourcewavelen = []
        sourceflambda = []
        for line in f:
            if line.startswith("#"):
                continue
            values = line.split()
            sourcewavelen.append(float(values[0]))
            sourceflambda.append(float(values[1]))
        f.close()
        self.wavelen = numpy.array(sourcewavelen)
        self.flambda = numpy.array(sourceflambda)
        self.fnu = None
        if name is None:
            self.name = filename
        else:
            self.name = name
        return

    def readSED_fnu(self, filename, name=None):
        """
        Read a file containing [lambda Fnu] (lambda in nm) (Fnu in Jansky).

        Does not resample wavelen/fnu/flambda onto a grid; leaves fnu set.
        """
        # Try to open the data file.
        try:
            if filename.endswith('.gz'):
                f = gzip.open(filename, 'r')
            else:
                f = open(filename, 'r')
        #if the above fails, look for the file with and without the gz
        except IOError:
            try:
                if filename.endswith(".gz"):
                    f = open(filename[:-3], 'r')
                else:
                    f = gzip.open(filename+".gz", 'r')
            except IOError:
                raise IOError("The throughput file %s does not exist" %(filename))
        # Read source SED from file - lambda, fnu should be first two columns in the file.
        # lambda should be in nm and fnu should be in Jansky.
        sourcewavelen = []
        sourcefnu = []
        for line in f:
            if line.startswith("#"):
                continue
            values = line.split()
            sourcewavelen.append(float(values[0]))
            sourcefnu.append(float(values[1]))
        f.close()
        # Convert to numpy arrays.
        sourcewavelen = numpy.array(sourcewavelen)
        sourcefnu = numpy.array(sourcefnu)
        # Convert fnu to flambda
        self.fnuToflambda(sourcewavelen, sourcefnu)
        if name is None:
            self.name = filename
        else:
            self.name = name
        return

    def getSED_flambda(self):
        """
        Return copy of wavelen/flambda.
        """
        # Get new memory copies of the arrays.
        wavelen = numpy.copy(self.wavelen)
        flambda = numpy.copy(self.flambda)
        return wavelen, flambda

    def getSED_fnu(self):
        """
        Return copy of wavelen/fnu, without altering self.
        """
        wavelen = numpy.copy(self.wavelen)
        # Check if fnu currently set.
        if self.fnu is not None:
            # Get new memory copy of fnu.
            fnu = numpy.copy(self.fnu)
        else:
            # Fnu was not set .. grab copy fnu without changing self.
            wavelen, fnu = self.flambdaTofnu(self.wavelen, self.flambda)
            # Now wavelen/fnu (new mem) are gridded evenly, but self.wavelen/flambda/fnu remain unchanged.
        return wavelen, fnu

    ## Methods that update or change self.

    def clearSED(self):
        """
        Reset all data in sed to None.
        """
        self.wavelen = None
        self.fnu = None
        self.flambda = None
        self.zp = -8.9
        self.name = None
        return

    def synchronizeSED(self, wavelen_min=None, wavelen_max=None, wavelen_step=None):
        """
        Set all wavelen/flambda/fnu values, potentially on min/max/step grid.

        Uses flambda to recalculate fnu. If wavelen min/max/step are given, resamples
        wavelength/flambda/fnu onto an even grid with these values.
        """
        # Grid wavelength/flambda/fnu if desired.
        if ((wavelen_min is not None) and (wavelen_max is not None) and (wavelen_step is not None)):
            self.resampleSED(wavelen_min=wavelen_min, wavelen_max=wavelen_max,
                             wavelen_step=wavelen_step)
        # Reset or set fnu.
        self.flambdaTofnu()
        return

    ## Utilities common to several later methods.

    def _checkUseSelf(self, wavelen, flux):
        """
        Simple utility to check if should be using self's data or passed arrays.

        Also does data integrity check on wavelen/flux if not self.
        """
        update_self = False
        if (wavelen is None) or (flux is None):
            # Then one of the arrays was not passed - check if this is true for both arrays.
            if (wavelen is not None) or (flux is not None):
                # Then one of the arrays was passed - raise exception.
                raise ValueError("Must either pass *both* wavelen/flux pair, or use defaults.")
            update_self = True
        else:
            # Both of the arrays were passed in - check their validity.
            if (isinstance(wavelen, numpy.ndarray)==False) or (isinstance(flux, numpy.ndarray)==False):
                raise ValueError("Must pass wavelen/flux as numpy arrays.")
            if len(wavelen)!=len(flux):
                raise ValueError("Must pass equal length wavelen/flux arrays.")
        return update_self

    def _needResample(self, wavelen_match=None, wavelen=None,
                     wavelen_min=None, wavelen_max=None, wavelen_step=None):
        """
        Check if wavelen or self.wavelen matches wavelen or wavelen_min/max/step grid.
        """
        # Check if should use self or passed wavelen.
        if wavelen is None:
            wavelen = self.wavelen
        # Check if wavelength arrays are equal, if wavelen_match passed.
        if wavelen_match is not None:
            if numpy.shape(wavelen_match) != numpy.shape(wavelen):
                need_regrid=True
            else:
                # check the elements to see if any vary
                need_regrid = numpy.any(abs(wavelen_match-wavelen)>1e-10)
        else:
            need_regrid = True
            # Check if wavelen_min/max/step are set - if ==None, then return (no regridding).
            # It's possible (writeSED) to call this routine, even with no final grid in mind.
            if ((wavelen_min is None) and (wavelen_max is None) and (wavelen_step is None)):
                need_regrid = False
            else:
                # Okay, now look at comparison of wavelen to the grid.
                wavelen_max_in = wavelen[len(wavelen)-1]
                wavelen_min_in = wavelen[0]
                # First check match to minimum/maximum :
                if ((wavelen_min_in == wavelen_min) and (wavelen_max_in == wavelen_max)):
                    # Then check on step size in wavelength array.
                    stepsize = numpy.unique(numpy.diff(wavelen))
                    if (len(stepsize) == 1) and (stepsize[0] == wavelen_step):
                        need_regrid = False
        # At this point, need_grid=True unless it's proven to be False, so return value.
        return need_regrid

    def resampleSED(self, wavelen=None, flux=None, wavelen_match=None,
                    wavelen_min=None, wavelen_max=None, wavelen_step=None, force=False):
        """
        Resample flux onto grid defined by min/max/step OR another wavelength array.

        Give method wavelen/flux OR default to self.wavelen/self.flambda.
        Method either returns wavelen/flambda (if given those arrays) or updates wavelen/flambda in self.
         If updating self, resets fnu to None.
         Method will first check if resampling needs to be done or not, unless 'force' is True.
        """
        # Check if need resampling:
        if (self._needResample(wavelen_match=wavelen_match, wavelen=wavelen, wavelen_min=wavelen_min,
                                wavelen_max=wavelen_max, wavelen_step=wavelen_step)) or (force):
            # Is method acting on self.wavelen/flambda or passed in wavelen/flux arrays?
            update_self=self._checkUseSelf(wavelen, flux)
            if update_self:
                wavelen=self.wavelen
                flux=self.flambda
                self.fnu = None
            # Now, on with the resampling.
            # Set up gridded wavelength or copy of wavelen array to match.
            if wavelen_match is None:
                if ((wavelen_min is None) and (wavelen_max is None) and (wavelen_step is None)):
                    raise ValueError('Must set either wavelen_match or wavelen_min/max/step.')
                wavelen_grid = numpy.arange(wavelen_min, wavelen_max+wavelen_step,
                                            wavelen_step, dtype='float')
            else:
                wavelen_grid = numpy.copy(wavelen_match)
            # Check if the wavelength range desired and the wavelength range of the object overlap.
            # If there is any non-overlap, raise exception.
            if (wavelen.max() < wavelen_grid.max()) or (wavelen.min() > wavelen_grid.min()):
                warnings.warn('There is an area of non-overlap between desired wavelength range (%.2f to %.2f) and sed %s (%.2f to %2.f)' % (wavelen_grid.min(), wavelen_grid.max(), self.name, wavelen.min(), wavelen.max()))
            flux_grid = numpy.empty(len(wavelen), dtype='float')
            # Do the interpolation of wavelen/flux onto grid. (type/len failures will die here).
            flux_grid = numpy.interp(wavelen_grid, wavelen, flux, left=numpy.NaN, right=numpy.NaN)
            # Update self values if necessary.
            if update_self:
                self.wavelen = wavelen_grid
                self.flambda = flux_grid
                return
            return wavelen_grid, flux_grid
        else: # wavelength grids already match.
            update_self = self._checkUseSelf(wavelen, flux)
            if update_self:
                return
            return wavelen, flux

    def flambdaTofnu(self, wavelen=None, flambda=None):
        """
        Convert flambda into fnu.

        This routine assumes that flambda is in ergs/cm^s/s/nm and produces fnu in Jansky.
        Can act on self or user can provide wavelen/flambda and get back wavelen/fnu.
        """
        # Change Flamda to Fnu by multiplying Flambda * lambda^2 = Fv
        # Fv dv = Fl dl .. Fv = Fl dl / dv = Fl dl / (dl*c/l/l) = Fl*l*l/c
        # Check - Is the method acting on self.wavelen/flambda/fnu or passed wavelen/flambda arrays?
        update_self = self._checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen = self.wavelen
            flambda=self.flambda
            self.fnu = None
        # Now on with the calculation.
        # Calculate fnu.
        fnu = flambda * wavelen * wavelen * PhotometricDefaults.nm2m / PhotometricDefaults.lightspeed
        fnu = fnu * PhotometricDefaults.ergsetc2jansky
        # If are using/updating self, then *all* wavelen/flambda/fnu will be gridded.
        # This is so wavelen/fnu AND wavelen/flambda can be kept in sync.
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
            self.fnu = fnu
            return
        # Return wavelen, fnu, unless updating self (then does not return).
        return wavelen, fnu

    def fnuToflambda(self, wavelen=None, fnu=None):
        """
        Convert fnu into flambda.

        Assumes fnu in units of Jansky and flambda in ergs/cm^s/s/nm.
        Can act on self or user can give wavelen/fnu and get wavelen/flambda returned.
        """
        # Fv dv = Fl dl .. Fv = Fl dl / dv = Fl dl / (dl*c/l/l) = Fl*l*l/c
        # Is method acting on self or passed arrays?
        update_self = self._checkUseSelf(wavelen, fnu)
        if update_self:
            wavelen = self.wavelen
            fnu = self.fnu
        # On with the calculation.
        # Calculate flambda.
        flambda = fnu / wavelen / wavelen * PhotometricDefaults.lightspeed / PhotometricDefaults.nm2m
        flambda = flambda / PhotometricDefaults.ergsetc2jansky
        # If updating self, then *all of wavelen/fnu/flambda will be updated.
        # This is so wavelen/fnu AND wavelen/flambda can be kept in sync.
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
            self.fnu = fnu
            return
        # Return wavelen/flambda.
        return wavelen, flambda

    ## methods to alter the sed

    def redshiftSED(self, redshift, dimming=False, wavelen=None, flambda=None):
        """
        Redshift an SED, optionally adding cosmological dimming.

        Pass wavelen/flambda or redshift/update self.wavelen/flambda (unsets fnu).
        """
        # Updating self or passed arrays?
        update_self = self._checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen=self.wavelen
            flambda = self.flambda
            self.fnu = None
        else:
            # Make a copy of input data, because will change its values.
            wavelen = numpy.copy(wavelen)
            flambda = numpy.copy(flambda)
        # Okay, move onto redshifting the wavelen/flambda pair.
        # Or blueshift, as the case may be.
        if redshift<0:
            wavelen = wavelen / (1.0-redshift)
        else:
            wavelen = wavelen * (1.0+redshift)
        # Flambda now just has different wavelength for each value.
        # Add cosmological dimming if required.
        if dimming:
            if redshift<0:
                flambda = flambda * (1.0-redshift)
            else:
                flambda = flambda / (1.0+redshift)
        # Update self, if required - but just flambda (still no grid required).
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
            return
        return wavelen, flambda

    def setupCCMab(self, wavelen=None):
        """
        Calculate a(x) and b(x) for CCM dust model. (x=1/wavelen).

        If wavelen not specified, calculates a and b on the own object's wavelength grid.
        Returns a(x) and b(x) can be common to many seds, wavelen is the same.
        """
        # This extinction law taken from Cardelli, Clayton and Mathis ApJ 1989.
        # The general form is A_l / A(V) = a(x) + b(x)/R_V  (where x=1/lambda in microns),
        # then different values for a(x) and b(x) depending on wavelength regime.
        # Also, the extinction is parametrized as R_v = A_v / E(B-V).
        # Magnitudes of extinction (A_l) translates to flux by a_l = -2.5log(f_red / f_nonred).
        if wavelen is None:
            wavelen = numpy.copy(self.wavelen)
        a_x = numpy.zeros(len(wavelen), dtype='float')
        b_x = numpy.zeros(len(wavelen), dtype='float')
        # Convert wavelength to x (in inverse microns).
        x = numpy.empty(len(wavelen), dtype=float)
        nm_to_micron = 1/1000.0
        x = 1.0 / (wavelen * nm_to_micron)
        # Dust in infrared 0.3 /mu < x < 1.1 /mu (inverse microns).
        condition = (x>=0.3) & (x<=1.1)
        if len(a_x[condition]) > 0 :
            y = x[condition]
            a_x[condition] = 0.574 * y**1.61
            b_x[condition] = -0.527 * y**1.61
        # Dust in optical/NIR 1.1 /mu < x < 3.3 /mu region.
        condition = (x >=1.1) & (x<=3.3)
        if len(a_x[condition])>0:
            y = x[condition] - 1.82
            a_x[condition] = 1 + 0.104*y - 0.609*y**2 + 0.701*y**3 + 1.137*y**4
            a_x[condition] = a_x[condition] - 1.718*y**5 - 0.827*y**6 + 1.647*y**7 - 0.505*y**8
            b_x[condition] = 1.952*y + 2.908*y**2 - 3.989*y**3 - 7.985*y**4
            b_x[condition] = b_x[condition] + 11.102*y**5 + 5.491*y**6 - 10.805*y**7 + 3.347*y**8
        # Dust in ultraviolet and UV (if needed for high-z) 3.3 /mu< x< 8 /mu.
        condition = (x>=3.3) & (x<5.9)
        if len(a_x[condition])>0:
            y = x[condition]
            a_x[condition] = 1.752 - 0.316*y - 0.104/((y-4.67)**2 + 0.341)
            b_x[condition] = -3.090 + 1.825*y + 1.206/((y-4.62)**2 + 0.263)
        condition = (x>5.9) & (x<8)
        if len(a_x[condition])>0:
            y = x[condition]
            Fa_x = numpy.empty(len(a_x[condition]), dtype=float)
            Fb_x = numpy.empty(len(a_x[condition]), dtype=float)
            Fa_x = -0.04473*(y-5.9)**2 - 0.009779*(y-5.9)**3
            Fb_x = 0.2130*(y-5.9)**2 + 0.1207*(y-5.9)**3
            a_x[condition] = 1.752 - 0.316*y - 0.104/((y-4.67)**2 + 0.341) + Fa_x
            b_x[condition] = -3.090 + 1.825*y + 1.206/((y-4.62)**2 + 0.263) + Fb_x
        # Dust in far UV (if needed for high-z) 8 /mu < x < 10 /mu region.
        condition = (x >= 8) & (x<= 11.)
        if len(a_x[condition])>0:
            y = x[condition]-8.0
            a_x[condition] = -1.073 - 0.628*(y) + 0.137*(y)**2 - 0.070*(y)**3
            b_x[condition] = 13.670 + 4.257*(y) - 0.420*(y)**2 + 0.374*(y)**3
        return a_x, b_x

    def addCCMDust(self, a_x, b_x, A_v=None, ebv=None, R_v=3.1, wavelen=None, flambda=None):
        """
        Add CCM dust model extinction to the SED, modifying flambda and fnu.

        Specify any two of A_V, E(B-V) or R_V (=3.1 default).
        """
        # The extinction law taken from Cardelli, Clayton and Mathis ApJ 1989.
        # The general form is A_l / A(V) = a(x) + b(x)/R_V  (where x=1/lambda in microns).
        # Then, different values for a(x) and b(x) depending on wavelength regime.
        # Also, the extinction is parametrized as R_v = A_v / E(B-V).
        # The magnitudes of extinction (A_l) translates to flux by a_l = -2.5log(f_red / f_nonred).
        #
        # Figure out if updating self or passed arrays.
        update_self = self._checkUseSelf(wavelen, flambda)
        if update_self:
            wavelen = self.wavelen
            flambda = self.flambda
            self.fnu = None
        else:
            wavelen = numpy.copy(wavelen)
            flambda = numpy.copy(flambda)
        # Input parameters for reddening can include any of 3 parameters; only 2 are independent.
        # Figure out what parameters were given, and see if self-consistent.
        if R_v == 3.1:
            if A_v is None:
                A_v = R_v * ebv
            elif (A_v is not None) and (ebv is not None):
                # Specified A_v and ebv, so R_v should be nondefault.
                R_v = A_v / ebv
        if (R_v != 3.1):
            if (A_v is not None) and (ebv is not None):
                calcRv = A_v / ebv
                if calcRv != R_v:
                    raise ValueError("CCM parametrization expects R_v = A_v / E(B-V);",
                                     "Please check input values, because values are inconsistent.")
            elif A_v is None:
                A_v = R_v * ebv
        # R_v and A_v values are specified or calculated.
        A_lambda = numpy.empty(len(wavelen), dtype=float)
        dust = numpy.empty(len(wavelen), dtype=float)
        A_lambda = (a_x + b_x / R_v) * A_v
        # dmag_red(dust) = -2.5 log10 (f_red / f_nored) : (f_red / f_nored) = 10**-0.4*dmag_red
        dust = numpy.power(10.0, -0.4*A_lambda)
        flambda = flambda * dust
        # Update self if required.
        if update_self:
            self.flambda = flambda
            return
        return wavelen, flambda


    def multiplySED(self, other_sed, wavelen_step=PhotometricDefaults.wavelenstep):
        """
        Multiply two SEDs together - flambda * flambda - and return a new sed object.

        Unless the two wavelength arrays are equal, returns a SED gridded with stepsize wavelen_step
        over intersecting wavelength region. Does not alter self or other_sed.
        """
        # Check if the wavelength arrays are equal (in which case do not resample)
        if (numpy.all(self.wavelen==other_sed.wavelen)):
            flambda = self.flambda * other_sed.flambda
            new_sed = Sed(self.wavelen, flambda=flambda)
        else:
            # Find overlapping wavelength region.
            wavelen_max = min(self.wavelen.max(), other_sed.wavelen.max())
            wavelen_min = max(self.wavelen.min(), other_sed.wavelen.min())
            if wavelen_max < wavelen_min:
                raise Exception('The two SEDS do not overlap in wavelength space.')
            # Set up wavelen/flambda of first object, on grid.
            wavelen_1, flambda_1 = self.resampleSED(self.wavelen, self.flambda,
                                                    wavelen_min=wavelen_min,
                                                    wavelen_max=wavelen_max,
                                                    wavelen_step=wavelen_step)
            # Set up wavelen/flambda of second object, on grid.
            wavelen_2, flambda_2 = self.resampleSED(wavelen=other_sed.wavelen, flux=other_sed.flambda,
                                                    wavelen_min=wavelen_min, wavelen_max=wavelen_max,
                                                    wavelen_step = wavelen_step)
            # Multiply the two flambda together.
            flambda = flambda_1 * flambda_2
            # Instantiate new sed object. wavelen_1 == wavelen_2 as both are on grid.
            new_sed = Sed(wavelen_1, flambda)
        return new_sed

    ## routines related to magnitudes and fluxes

    def calcADU(self, bandpass, wavelen=None, fnu=None,
                expTime=PhotometricDefaults.exptime, effarea=PhotometricDefaults.effarea,
                gain=PhotometricDefaults.gain):
        """
        Calculate the number of adu from camera, using sb and fnu.

        Given wavelen/fnu arrays or use self. Self or passed wavelen/fnu arrays will be unchanged.
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass;
         (temporary values of these are used).
        """
        use_self = self._checkUseSelf(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if self.fnu is None:
                # If fnu not present, calculate. (does not regrid).
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
        # Make sure wavelen/fnu are on the same wavelength grid as bandpass.
        wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate the number of photons.
        dlambda = wavelen[1] - wavelen[0]
        # Nphoton in units of 10^-23 ergs/cm^s/nm.
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        adu = nphoton * (expTime * effarea/gain) * \
              (1/PhotometricDefaults.ergsetc2jansky) * \
              (1/PhotometricDefaults.planck) * dlambda
        return adu


    def calcMag(self, bandpass, wavelen=None, fnu=None):
        """
        Calculate the AB magnitude of an object, using phi the normalized system response.

        Can pass wavelen/fnu arrays or use self. Self or passed wavelen/fnu arrays will be unchanged.
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass;
         (but only temporary values of these are used).
         """
        # Note - the behavior in this first section might be considered a little odd.
        # However, I felt calculating a magnitude should not (unexpectedly) regrid your
        # wavelen/flambda information if you were using self., as this is not obvious from the "outside".
        # To preserve 'user logic', the wavelen/flambda of self are left untouched. Unfortunately
        # this means, this method can be used inefficiently if calculating many magnitudes with
        # the same sed and same bandpass region - in that case, use self.synchronizeSED() with
        # the wavelen min/max/step set to the bandpass min/max/step first ..
        # then you can calculate multiple magnitudes much more efficiently!
        use_self = self._checkUseSelf(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if self.fnu is None:
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
        # Continue with magnitude calculation.
        # Put bandpass and wavelen/fnu are on the same grid.
        wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate bandpass phi value if required.
        if bandpass.phi is None:
            bandpass.sbTophi()
        # Calculate flux in bandpass and then AB magnitude.
        dlambda = wavelen[1] - wavelen[0]
        flux = (fnu*bandpass.phi).sum() * dlambda
        if flux < 1e-300:
            raise Exception("This SED has no flux within this bandpass.")
        mag = -2.5 * numpy.log10(flux) - self.zp
        return mag

    def calcFlux(self, bandpass, wavelen=None, fnu=None):
        """
        Calculate the F_b (integrated flux of an object, **above the atmosphere**), using phi.

        Passed wavelen/fnu arrays will be unchanged, but if uses self will check if fnu is set.
        Calculating the AB mag requires the wavelen/fnu pair to be on the same grid as bandpass;
           (temporary values of these are used).
        """
        use_self = self._checkUseSelf(wavelen, fnu)
        # Use self values if desired, otherwise use values passed to function.
        if use_self:
            # Calculate fnu if required.
            if self.fnu is None:
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
        # Go on with magnitude calculation.
        wavelen, fnu = self.resampleSED(wavelen, fnu, wavelen_match=bandpass.wavelen)
        # Calculate bandpass phi value if required.
        if bandpass.phi is None:
            bandpass.sbTophi()
        # Calculate flux in bandpass and return this value.
        dlambda = wavelen[1] - wavelen[0]
        flux = (fnu*bandpass.phi).sum() * dlambda
        return flux

    def calcFluxNorm(self, magmatch, bandpass, wavelen=None, fnu=None):
        """
        Calculate the fluxNorm (SED normalization value for a given mag) for a sed.

        Equivalent to adjusting a particular f_nu to Jansky's appropriate for the desired mag.
        Can pass wavelen/fnu or apply to self.
        """
        use_self = self._checkUseSelf(wavelen, fnu)
        if use_self:
            # Check possibility that fnu is not calculated yet.
            if self.fnu is None:
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
        # Fluxnorm gets applied to f_nu (fluxnorm * SED(f_nu) * PHI = mag - 8.9 (AB zeropoint).
        # FluxNorm * SED => correct magnitudes for this object.
        # Calculate fluxnorm.
        curmag = self.calcMag(bandpass, wavelen, fnu)
        if curmag == self.badval:
            return self.badval
        dmag = magmatch - curmag
        fluxnorm = numpy.power(10, (-0.4*dmag))
        return fluxnorm

    def multiplyFluxNorm(self, fluxNorm, wavelen=None, fnu=None):
        """
        Multiply wavelen/fnu (or self.wavelen/fnu) by fluxnorm.

        Returns wavelen/fnu arrays (or updates self).
        Note that multiplyFluxNorm does not regrid self.wavelen/flambda/fnu at all.
        """
        # Note that fluxNorm is intended to be applied to f_nu,
        # so that fluxnorm*fnu*phi = mag (expected magnitude).
        update_self = self._checkUseSelf(wavelen, fnu)
        if update_self:
            # Make sure fnu is defined.
            if self.fnu is None:
                self.flambdaTofnu()
            wavelen = self.wavelen
            fnu = self.fnu
        else:
            # Require new copy of the data for multiply.
            wavelen = numpy.copy(wavelen)
            fnu = numpy.copy(fnu)
        # Apply fluxnorm.
        fnu = fnu * fluxNorm
        # Update self.
        if update_self:
            self.wavelen = wavelen
            self.fnu = fnu
            # Update flambda as well.
            self.fnuToflambda()
            return
        # Else return new wavelen/fnu pairs.
        return wavelen, fnu

    def renormalizeSED(self, wavelen=None, flambda=None, fnu=None,
                       lambdanorm=500, normvalue=1, gap=0, normflux='flambda',
                       wavelen_step=PhotometricDefaults.wavelenstep):
        """
        Renormalize sed in flambda to have normflux=normvalue @ lambdanorm or averaged over gap.

        Can normalized in flambda or fnu values. wavelen_step specifies the wavelength spacing when using 'gap'.
        Either returns wavelen/flambda values or updates self.
        """
        # Normalizes the fnu/flambda SED at one wavelength or average value over small range (gap).
        # This is useful for generating SED catalogs, mostly, to make them match schema.
        # Do not use this for calculating specific magnitudes -- use calcfluxNorm and multiplyFluxNorm.
        # Start normalizing wavelen/flambda.
        if normflux=='flambda':
            update_self = self._checkUseSelf(wavelen, flambda)
            if update_self:
                wavelen = self.wavelen
                flambda = self.flambda
            else:
                # Make a copy of the input data.
                wavelen = numpy.copy(wavelen)
                # Look for either flambda or fnu in input data.
                if flambda is None:
                    if fnu is None:
                        raise Exception("If passing wavelength, must also pass fnu or flambda.")
                    # If not given flambda, must calculate from the given values of fnu.
                    wavelen, flambda = self.fnuToflambda(wavelen, fnu)
                # Make a copy of the input data.
                else:
                    flambda = numpy.copy(flambda)
            # Calculate renormalization values.
            # Check that flambda is defined at the wavelength want to use for renormalization.
            if (lambdanorm > wavelen.max()) or (lambdanorm < wavelen.min()):
                raise Exception("Desired wavelength for renormalization, %f, is outside defined wavelength range." %(lambdanorm))
            # "standard" schema have flambda = 1 at 500 nm.
            if gap==0:
                flambda_atpt = numpy.interp(lambdanorm, wavelen, flambda, left=None, right=None)
                gapval = flambda_atpt
            else:
                lambdapt = numpy.arange(lambdanorm-gap, lambdanorm+gap, wavelen_step, dtype=float)
                flambda_atpt = numpy.zeros(len(lambdapt), dtype='float')
                flambda_atpt = numpy.interp(lambdapt, wavelen, flambda, left=None, right=None)
                gapval = flambda_atpt.sum()/len(lambdapt)
            # Now renormalize fnu and flambda in the case of normalizing flambda.
            if gapval == 0:
                raise Exception("Original flambda is 0 at the desired point of normalization. Cannot renormalize.")
            konst = normvalue/gapval
            flambda = flambda * konst
            wavelen, fnu = self.flambdaTofnu(wavelen, flambda)
        elif normflux=='fnu':
            update_self = self._checkUseSelf(wavelen, fnu)
            if update_self:
                wavelen = self.wavelen
                if self.fnu is None:
                    self.flambdaTofnu()
                fnu = self.fnu
            else:
                # Make a copy of the input data.
                wavelen = numpy.copy(wavelen)
                # Look for either flambda or fnu in input data.
                if fnu is None:
                    if flambda is None:
                        raise Exception("If passing wavelength, must also pass fnu or flambda.")
                    wavelen, fnu = self.flambdaTofnu(wavelen, fnu)
                # Make a copy of the input data.
                else:
                    fnu = numpy.copy(fnu)
            # Calculate renormalization values.
            # Check that flambda is defined at the wavelength want to use for renormalization.
            if (lambdanorm > wavelen.max()) or (lambdanorm < wavelen.min()):
                raise Exception("Desired wavelength for renormalization, %f, is outside defined wavelength range." %(lambdanorm))
            if gap==0:
                fnu_atpt = numpy.interp(lambdanorm, wavelen, flambda, left=None, right=None)
                gapval = fnu_atpt
            else:
                lambdapt = numpy.arange(lambdanorm-gap, lambdanorm+gap, wavelen_step, dtype=float)
                fnu_atpt = numpy.zeros(len(lambdapt), dtype='float')
                fnu_atpt = numpy.interp(lambdapt, wavelen, fnu, left=None, right=None)
                gapval = fnu_atpt.sum()/len(lambdapt)
            # Now renormalize fnu and flambda in the case of normalizing fnu.
            if gapval == 0:
                raise Exception("Original fnu is 0 at the desired point of normalization. Cannot renormalize.")
            konst = normvalue/gapval
            fnu = fnu * konst
            wavelen, flambda = self.fnutoflambda(wavelen,fnu)
        if update_self:
            self.wavelen = wavelen
            self.flambda = flambda
            self.fnu = fnu
            return
        new_sed = Sed(wavelen=wavelen, flambda=flambda)
        return new_sed


    def writeSED(self, filename, print_header=None, print_fnu=False,
                 wavelen_min=None, wavelen_max=None, wavelen_step=None):
        """
        Write SED (wavelen, flambda, optional fnu) out to file.

        Option of adding a header line (such as version info) to output file.
        Does not alter self, regardless of grid or presence/absence of fnu.
        """
        # This can be useful for debugging or recording an SED.
        f = open(filename, 'w')
        wavelen = self.wavelen
        flambda = self.flambda
        wavelen, flambda = self.resampleSED(wavelen, flambda, wavelen_min=wavelen_min,
                                            wavelen_max=wavelen_max,
                                            wavelen_step=wavelen_step)
        # Then just use this gridded wavelen/flambda to calculate fnu.
        # Print header.
        if print_header is not None:
            print >>f, "#", print_header
        # Print standard header info.
        if print_fnu:
            wavelen, fnu = self.flambdaTofnu(wavelen, flambda)
            print >>f, "# Wavelength(nm)  Flambda(ergs/cm^s/s/nm)   Fnu(Jansky)"
        else:
            print >>f, "# Wavelength(nm)  Flambda(ergs/cm^s/s/nm)"
        for i in range(0, len(wavelen), 1):
            if print_fnu:
                fnu = self.flambdaTofnu(wavelen=wavelen, flambda=flambda)
                print >> f, wavelen[i], flambda[i], fnu[i]
            else:
                #print >> f, self.wavelen[i], self.flambda[i]
                print >> f, "%.2f %.7g" %(wavelen[i], flambda[i])
        # Done writing, close file.
        f.close()
        return

    def calcNeff(self, seeing, platescale):
        """
        Calculate the effective number of pixels in a double Gaussian PSF

        @param [in] seeing in arcseconds

        @param [in] platescale in arcseconds per pixel

        @param [out] the effective number of pixels contained by a double
        Gaussian PSF

        see equation 31 of the SNR document
        http://www.astro.washington.edu/users/ivezic/Teaching/Astr511/LSST_SNRdoc.pdf
        """
        return 2.436*(seeing/platescale)**2

    def calcInstrNoiseSq(self, readnoise, darkcurrent, expTime, nexp, othernoise):
        """
        Combine all of the noise due to intrumentation into one value

        @param [in] readnoise

        @param [in] darkcurrent

        @param [in] expTime (length of a single exposure in seconds)

        @param [in] nexp (number of exposures)

        @param [in] othernoise (noise from sources not accounted for above)

        @param [out] The noise due to all of these sources added in quadrature
        """
        return nexp*readnoise**2 + darkcurrent*expTime*nexp + nexp*othernoise**2

    def calcNonSourceNoiseSq(self, skySed, hardwarebandpass, readnoise, darkcurrent,
                           othernoise, seeing, effarea, expTime, nexp, platescale, gain):
        """
        Calculate the noise due to things that are not the source being observed
        (i.e. intrumentation and sky background)

        @param [in] skySed -- an instantiation of the Sed class representing the sky

        @param [in] hardwarebandpass -- an instantiation of the Bandpass class representing
        just the instrumentation throughputs

        @param [in] readnoise

        @param [in] darkcurrent

        @param [in] othernoise

        @param [in] seeing in arcseconds

        @param [in] effarea -- effective area of the primary mirror in square centimeters

        @param [in] expTime -- duration of one exposure in seconds

        @param [in] nexp -- number of exposures being combined

        @param [in] platescale in arcseconds per pixel

        @param [in] gain in electrons per adu

        @param [out] total noise squared (in ADU)

        @param [out] noise squared due just to the instrument

        @param [out] noise squared due to the sky

        @param [out] noise squared due to sky measurement (presently set to zero)

        @param [out] the effective number of pixels in a double Gaussian PSF
        """

        #This method outputs all of the parameters calculated along the way
        #so that the verbose version of calcSNR_psf still works

        #Calculate the effective number of pixels for double-Gaussian PSF
        neff = self.calcNeff(seeing, platescale)

        #Calculate the counts form the sky
        skycounts = skySed.calcADU(hardwarebandpass, expTime=expTime*nexp, effarea=effarea, gain=gain) \
                    * platescale * platescale

        #Calculate the square of the noise due to instrumental effects.
        #Include the readout noise as many times as there are exposures
        noise_instr_sq = self.calcInstrNoiseSq(readnoise, darkcurrent, expTime, nexp, othernoise)

        #Calculate the square of the noise due to sky background poisson noise
        noise_sky_sq = skycounts/gain

        #Discount error in sky measurement for now
        noise_skymeasurement_sq = 0

        total_noise_sq = neff*(noise_sky_sq + noise_instr_sq + noise_skymeasurement_sq)
        return total_noise_sq, noise_instr_sq, noise_sky_sq, noise_skymeasurement_sq, skycounts, neff

    def calcSNR_psf(self, totalbandpass, skysed, hardwarebandpass,
                    readnoise=PhotometricDefaults.rdnoise,
                    darkcurrent=PhotometricDefaults.darkcurrent,
                    othernoise=PhotometricDefaults.othernoise,
                    seeing=PhotometricDefaults.seeing['r'],
                    effarea=PhotometricDefaults.effarea,
                    expTime=PhotometricDefaults.exptime, nexp=PhotometricDefaults.nexp,
                    platescale=PhotometricDefaults.platescale,
                    gain=PhotometricDefaults.gain, verbose=False):
        """
        Calculate the signal to noise ratio for a source, given the bandpass(es) and sky SED.

        For a given source, sky sed, total bandpass and hardware bandpass, as well as
        seeing / expTime, calculates the SNR with optimal PSF extraction
        assuming a double-gaussian PSF. Assumes that all values (readnoise/othernoise
        /darkcurrent) are given in appropriate units for gain. (gain=1 -> values are in electrons,
        while if gain!=1 then values given should be in adu, including readnoise/darkcurrent).
        """
        # Calculate the counts from the source.
        sourcecounts = self.calcADU(totalbandpass, expTime=expTime*nexp, effarea=effarea, gain=gain)

        # Calculate the (square of the) noise due to signal poisson noise.
        noise_source_sq = sourcecounts/gain

        non_source_noise_sq, \
        noise_instr_sq, \
        noise_sky_sq, \
        noise_skymeasurement_sq, \
        skycounts, neff = self.calcNonSourceNoiseSq(skysed, hardwarebandpass, readnoise, darkcurrent, \
                                                    othernoise, seeing, effarea, expTime, nexp, platescale, gain)
        # Calculate total noise
        noise = numpy.sqrt(noise_source_sq + non_source_noise_sq)
        # Calculate the signal to noise ratio.
        snr = sourcecounts / noise
        if verbose:
            print "For Nexp %.1f of time %.1f: " % (nexp, expTime)
            print "Counts from source: %.2f  Counts from sky: %.2f" %(sourcecounts, skycounts)
            print "Seeing: %.2f('')  Neff pixels: %.3f(pix)" %(seeing, neff)
            print "Noise from sky: %.2f Noise from instrument: %.2f" \
                %(numpy.sqrt(noise_sky_sq), numpy.sqrt(noise_instr_sq))
            print "Noise from source: %.2f" %(numpy.sqrt(noise_source_sq))
            print " Total Signal: %.2f   Total Noise: %.2f    SNR: %.2f" %(sourcecounts, noise, snr)
            # Return the signal to noise value.
        return snr

    ## Below here, methods that are appropriate for sed, but don't require SED.

    def calcSNR_mag(self, mag, m5):
        """
        Calculate the signal to noise of an object, given only the 5-sigma limiting mag.
        """
        flux_ratio = numpy.power(10, 0.4*(m5-mag))
        snr_obj = 5 * (flux_ratio)
        return snr_obj

    def calcMagError(self, mag, m5, nvisit=1):
        """
        Calculate the photometric error, for object catalog purposes.

        Returns photometric error for a given SNR.
        """
        # Calculate mag errors - see eqns 4-6 astroph/0805.2366.
        # Can either do single image (nvisit=1) or coadded (nvisit>1).
        # m5 should be the 5-sigma limiting magnitude of one image.
        error_sys = 0.005     # systematic error - 0.005
        if (nvisit>50):
            # Allowing for reduction in systematic error after many visits
            error_sys = 0.0025
        # Rgamma value comes from overview paper.
        rgamma = 0.039
        xval = numpy.power(10, 0.4*(mag-m5))
        error_rand = numpy.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)/numpy.sqrt(nvisit)
        mag_error = numpy.sqrt(error_sys*error_sys + error_rand*error_rand)
        return mag_error

    def calcAstrometricError(self, mag, m5, nvisit=1):
        """
        Calculate the astrometric error, for object catalog purposes.

        Returns astrometric error for a given SNR, in mas.
        """
        # The astrometric error can be applied to parallax or proper motion (for nvisit>1).
        # If applying to proper motion, should also divide by the # of years of the survey.
        # This is also referenced in the astroph/0805.2366 paper.
        # D. Monet suggests sqrt(Nvisit/2) for first 3 years, sqrt(N) for longer, in reduction of error
        # because of the astrometric measurement method, the systematic and random error are both reduced.
        # Zeljko says 'be conservative', so removing this reduction for now.
        rgamma = 0.039
        xval = numpy.power(10, 0.4*(mag-m5))
        # The average seeing is 0.7" (or 700 mas).
        error_rand = 700.0 * numpy.sqrt((0.04-rgamma)*xval + rgamma*xval*xval)
        error_rand = error_rand / numpy.sqrt(nvisit)
        # The systematic error floor in astrometry:
        error_sys = 10.0
        # These next few lines are the code removed due to Zeljko's 'be conservative' requirement.
        #if (nvisit<30):
        #    error_sys = error_sys/numpy.sqrt(nvisit/2.0)
        #if (nvisit>30):
        #    error_sys = error_sys/numpy.sqrt(nvisit)
        astrom_error = numpy.sqrt(error_sys * error_sys + error_rand*error_rand)
        return astrom_error


## Bonus, functions for many-magnitude calculation for many SEDs with a single bandpass

    def setupPhiArray(self, bandpasslist):
        """
        Sets up a 2-d numpy phi array from bandpasslist suitable for input to Sed's manyMagCalc.

        This is intended to be used once, most likely before using Sed's manyMagCalc many times on many SEDs.
        Returns 2-d phi array and the wavelen_step (dlambda) appropriate for that array.
        """
        # Calculate dlambda for phi array.
        wavelen_step = bandpasslist[0].wavelen[1] - bandpasslist[0].wavelen[0]
        wavelen_min = bandpasslist[0].wavelen[0]
        wavelen_max = bandpasslist[0].wavelen[len(bandpasslist[0].wavelen)-1]
        # Set up
        phiarray = numpy.empty((len(bandpasslist), len(bandpasslist[0].wavelen)), dtype='float')
        # Check phis calculated and on same wavelength grid.
        i = 0
        for bp in bandpasslist:
            # Be sure bandpasses on same grid and calculate phi.
            bp.resampleBandpass(wavelen_min=wavelen_min, wavelen_max=wavelen_max, wavelen_step=wavelen_step)
            bp.sbTophi()
            phiarray[i] = bp.phi
            i = i + 1
        return phiarray, wavelen_step


    def manyFluxCalc(self, phiarray, wavelen_step, observedBandPassInd=None):
        """
        Calculate fluxes of a single sed for which fnu has been evaluated in a
        set of bandpasses for which phiarray has been set up to have the same
        wavelength grid as the SED in units of ergs/cm^2/sec. It is assumed
        that `self.fnu` is set before calling this method, and that phiArray
        has the same wavelength grid as the Sed.


        Parameters
        ----------
        phiarray: `np.ndarray`, mandatory
            phiarray corresponding to the list of bandpasses in which the band
            fluxes need to be calculated, in the same wavelength grid as the SED

        wavelen_step: `float`, mandatory
            the uniform grid size of the SED

        observedBandPassInd: list of integers, optional, defaults to None
            list of indices of phiarray corresponding to observed bandpasses,
            if None, the original phiarray is returned


        Returns
        -------
        `np.ndarray` with size equal to number of bandpass filters  band flux
        values in units of ergs/cm^2/sec

        .. note: Sed.manyFluxCalc `assumes` phiArray has the same wavelenghth
        grid as the Sed and that `sed.fnu` has been calculated for the sed,
        perhaps using `sed.flambdaTofnu()`. This requires calling
        `sed.setupPhiArray()` first. These assumptions are to avoid error
        checking within this function (for speed), but could lead to errors if
        method is used incorrectly.
        """

        if observedBandPassInd is not None:
            phiarray = phiarray[observedBandPassInd]
        flux = numpy.empty(len(phiarray), dtype='float')
        flux = numpy.sum(phiarray*self.fnu, axis=1)*wavelen_step
        return flux


    def manyMagCalc(self, phiarray, wavelen_step, observedBandPassInd=None):
        """
        Calculate many magnitudes for many bandpasses using a single sed.

        This method assumes that there will be flux within a particular bandpass
        (could return '-Inf' for a magnitude if there is none).
        Use setupPhiArray first, and note that Sed.manyMagCalc *assumes*
        phiArray has the same wavelength grid as the Sed, and that fnu has
        already been calculated for Sed.
        These assumptions are to avoid error checking within this function (for
        speed), but could lead to errors if method is used incorrectly.
        Parameters
        ----------
        phiarray: `np.ndarray`, mandatory
            phiarray corresponding to the list of bandpasses in which the band
            fluxes need to be calculated, in the same wavelength grid as the SED

        wavelen_step: `float`, mandatory
            the uniform grid size of the SED

        observedBandPassInd: list of integers, optional, defaults to None
            list of indices of phiarray corresponding to observed bandpasses,
            if None, the original phiarray is returned

        """
        fluxes = self.manyFluxCalc(phiarray, wavelen_step, observedBandPassInd)
        mags = -2.5*numpy.log10(fluxes) - self.zp
        return mags
