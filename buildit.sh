make distclean
./configure --FC=mpif90 --with-qipack=/home/feathern/devel/BenM/qi_pack -devel --with-fftw=/custom/software/fftw/3.3.7/intel/18.1
make
make install
