# This is usually done within the LSST environment by EUPS and setting up various packages (and done better). 
# For now, this is a potential shortcut if you do not have eups installed. Please set 'dir' to be the directory
# the phot_ref package was installed into. 
 
setenv dir $PWD

setenv PHOT_REF_SEDS_DIR $dir
setenv LSST_THROUGHPUTS_ATMOS $dir"/throughputs/atmos/"
setenv LSST_THROUGHPUTS_BASELINE $dir"/throughputs/baseline/"
setenv LSST_THROUGHPUTS_DEFAULT $LSST_THROUGHPUTS_BASELINE
setenv SDSS_THROUGHPUTS $dir"/throughputs/sdss/"
setenv MEGACAM_THROUGHPUTS $dir"/throughputs/megacam/"

setenv PYTHONPATH ${PYTHONPATH}:$dir"/python/"
