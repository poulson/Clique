set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# The serial XL compilers
set(BGXLC_BASE /soft/apps/ibmcmp-aug2011/vacpp/bg/9.0)
set(BGXLF_BASE /soft/apps/ibmcmp-aug2011/xlf/bg/11.1)
set(CMAKE_C_COMPILER       ${BGXLC_BASE}/bin/bgxlc_r)
set(CMAKE_CXX_COMPILER     ${BGXLC_BASE}/bin/bgxlC_r)
set(CMAKE_Fortran_COMPILER ${BGXLF_BASE}/bin/bgxlf_r)

# The MPI wrappers for the XL C and C++ compilers
set(BGMPI_BASE /bgsys/drivers/ppcfloor/comm)
set(MPI_C_COMPILER   ${BGMPI_BASE}/bin/mpixlc_r)
set(MPI_CXX_COMPILER ${BGMPI_BASE}/bin/mpixlcxx_r)

set(CXX_FLAGS "-g -O3")

set(ESSL_BASE "/soft/apps/ESSL-4.3.1-1")
set(IBMCMP_BASE "/soft/apps/ibmcmp-apr2011")
set(XLF_BASE "${IBMCMP_BASE}/xlf/bg/11.1/bglib")
set(XLSMP_BASE "${IBMCMP_BASE}/xlsmp/bg/1.7/bglib")

set(BGP_LAPACK "-L/soft/apps/LAPACK -llapack_bgp")
set(ESSL "-L${ESSL_BASE}/lib -lesslbg")
set(XLF_LIBS "-L${XLF_BASE} -lxlfmath -lxlf90_r")
set(XLOMP_SER "-L${XLSMP_BASE} -lxlomp_ser")

set(MATH_LIBS "${BGP_LAPACK};${ESSL};${XLF_LIBS};${XLOMP_SER}")

# Make sure we can find the ESSL headers
set(ESSL_INC "${ESSL_BASE}/include")
include_directories(${ESSL_INC})

