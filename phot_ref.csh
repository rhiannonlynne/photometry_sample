# This might usuaully be done within the LSST environment by EUPS and setting up various packages. 
# For now, this is a potential shortcut. 

setenv dir $PWD
setenv LSST_THROUGHPUTS_ATMOS $dir"/throughputs/atmos/"
setenv LSST_THROUGHPUTS_BASELINE $dir"/throughputs/baseline/"
setenv LSST_THROUGHPUTS_DEFAULT $LSST_THROUGHPUTS_BASELINE
setenv SDSS_THROUGHPUTS $dir"/throughputs/sdss/"
setenv MEGACAM_THROUGHPUTS $dir"/throughputs/megacam/"

setenv PYTHONPATH ${PYTHONPATH}:$dir"/python/"
