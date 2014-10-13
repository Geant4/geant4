#!/bin/sh
#
# Last update: 26-Jul-2013
# 
# Script to use only when the CMake build system is utilized,
# after  "cmake ..."  and  "make" , but before  "make install" .
#
# After the use of  EXCLUDE_FROM_ALL  on  test22/CMakeLists.txt
# (suggeseted by Ben on 26-Jul-2013 to solve a CMake problem in Windows
# platforms), the test22 executables need to be built explicitly,
# by doing "make test22_xxx". 
# This is done by this script for all the 47 tests.
# Notice that if you don't run this script, then you get an error during
# "make install".
# 
# You need to copy this script in the cmake-build-directory/tests/test22
# before you execute it as:
#   ./make_test22_4CMake.sh
#
make test22_ChipsX
make test22_E941
make test22_hA100
make test22_hA120
make test22_HARP
make test22_KmPchan
make test22_KpPchan
make test22_NA35
make test22_NA49
make test22_NA61
make test22_Nb_C_750
make test22_pA200
make test22_PbA_1p22_nE
make test22_PbarA_X
make test22_PbA_Rest
make test22_PbarP100
make test22_PbarP14p8
make test22_PbarP22p4
make test22_PbarP32
make test22_PbarP4p6
make test22_PbarP5p7
make test22_PbarP7p3
make test22_PbarP9p1
make test22_PbarPchan
make test22_PbC_608
make test22_PbDNe_Rest
make test22_PbTa_4
make test22_PbU_608
make test22_PimP100
make test22_PimP16
make test22_PimP40
make test22_PimPchan
make test22_PipP100
make test22_PipP16
make test22_PipP175
make test22_PipP8
make test22_PipPchan
make test22_PP100
make test22_PP12
make test22_PP158
make test22_PP175
make test22_PP205
make test22_PP24
make test22_PP360
make test22_PP400
make test22_PP69
make test22_PPchan
