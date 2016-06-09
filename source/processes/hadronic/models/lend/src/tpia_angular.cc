/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#include <string.h>

#include "tpia_target.h"
#include "tpia_misc.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int tpia_angular_initialize( statusMessageReporting *smr, tpia_angular *angular ) {

    memset( angular, 0, sizeof( tpia_angular ) );
    if( tpia_frame_setFromString( smr, "", "", 0, &(angular->frame) ) ) return( 1 );
    angular->binned.numberOfEs = 0;
    angular->binned.energies = NULL;
    return( 0 );
}
/*
************************************************************
*/
int tpia_angular_release( statusMessageReporting *smr, tpia_angular *angular ) {

    //angular->binned.energies = xData_free( smr, angular->binned.energies );
    angular->binned.energies = (tpia_EqualProbableBinSpectrum*) xData_free( smr, angular->binned.energies );
    angular->binned.numberOfEs = 0;
    return( 0 );
}
/*
************************************************************
*/
int tpia_angular_getFromElement( statusMessageReporting *smr, xData_element *angularElement, tpia_angular *angular ) {

    xData_Int n, nBins;

    if( ( tpia_frame_setFromElement( smr, angularElement, 3, &(angular->frame) ) ) != 0 ) return( 1 );
    if( ( angular->binned.energies = tpia_misc_getEqualProbableBin( smr, angularElement, &n, &nBins ) ) != NULL ) {
        angular->binned.iValue = 0;
        angular->binned.dValue = 0;
        angular->binned.nBins = (int) nBins;
        angular->binned.numberOfEs = (int) n;
        return( 0 );
    }
    return( 1 );
}
/*
************************************************************
*/
int tpia_angular_SampleMu( statusMessageReporting *smr, tpia_angular *angular, tpia_decaySamplingInfo *decaySamplingInfo ) {

    tpia_EqualProbableBinSpectra *binned = &(angular->binned);
    int nBins = binned->nBins, status = 1;
    double mu;
/*
Currently, only equal probable binning is supported.
Need to return frame info for mu also.
*/
    if( binned->numberOfEs > 0 ) {
        tpia_misc_sampleEqualProbableBin( smr, decaySamplingInfo, decaySamplingInfo->e_in, nBins, binned, &mu );
        if( mu < -1. ) mu = -1.;
        if( mu >  1. ) mu =  1.;
        decaySamplingInfo->mu = mu;
        status = 0;
    }

    return( status );
}

#if defined __cplusplus
}
#endif
