set(CMAKE_SYSTEM_NAME BlueGeneP-static)

# The serial XL compilers
set(CMAKE_C_COMPILER   /soft/apps/ibmcmp-apr2011/vacpp/bg/9.0/bin/bgxlc_r)
set(CMAKE_CXX_COMPILER /soft/apps/ibmcmp-apr2011/vacpp/bg/9.0/bin/bgxlC_r)

# The MPI wrappers for the XL C and C++ compilers
set(MPI_C_COMPILER   /bgsys/drivers/ppcfloor/comm/bin/mpixlc_r)
set(MPI_CXX_COMPILER /bgsys/drivers/ppcfloor/comm/bin/mpixlcxx_r)

set(CXX_FLAGS "-g -O2 -DBGP")

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

