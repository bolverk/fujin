cd $FUJIN_ROOT
mkdir -p external_libraries
cd external_libraries
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz
tar xf ./hdf5-1.10.4.tar.gz
cd hdf5-1.10.4
./configure --enable-cxx --prefix=`cd .. && pwd`
make
make check
make install