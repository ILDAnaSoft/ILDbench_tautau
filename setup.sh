unset MARLIN_DLL
source /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/v02-00-01/init_ilcsoft.sh

# this is stuff for tauola/tauspinner
export LHAPDF_DIR=/home/ilc/jeans/LHAPDF/lhapdf-5.9.1
export LHAPDF_INSTALL_DIR=${LHAPDF_DIR}/install
export LHAPATH=/home/ilc/jeans/LHAPDF/datasets
export TAUOLA_DIR=/home/ilc/jeans/tauola/v116c/TAUOLA

THIS_TOP_DIR=/home/ilc/jeans/myProcessors/ILDbench_tautau

export LD_LIBRARY_PATH=${LHAPDF_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${TAUOLA_DIR}/lib:${LD_LIBRARY_PATH}
export MARLIN_DLL=${THIS_TOP_DIR}/lib/libILDbench_tautau.so:$MARLIN_DLL
