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
#include <math.h>
#ifdef WIN32
   #include <float.h>
#endif

#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int _tpia_Legendre_getOrder( statusMessageReporting *smr, xData_element *LegendreOrder, tpia_Legendre *Legendre, tpia_EqualProbableBinSpectra *l );
/*
************************************************************
*/
int tpia_Legendre_initialize( statusMessageReporting *smr, tpia_Legendre *Legendre ) {

    memset( Legendre, 0, sizeof( tpia_Legendre ) );
    if( tpia_frame_setFromString( smr, "", "", 0, &(Legendre->frame) ) ) return( 1 );
    Legendre->binned.numberOfLs = 0;
    Legendre->binned.ls = NULL;
    return( 0 );
}
/*
************************************************************
*/
int tpia_Legendre_release( statusMessageReporting *smr, tpia_Legendre *Legendre ) {

    int i;

    for( i = 0; i < Legendre->binned.numberOfLs; i++ ) xData_free( smr, Legendre->binned.ls[i].energies );
    //Legendre->binned.ls = xData_free( smr, Legendre->binned.ls );
    Legendre->binned.ls = (tpia_EqualProbableBinSpectra*) xData_free( smr, Legendre->binned.ls );
    tpia_Legendre_initialize( smr, Legendre );
    return( 0 );
}
/*
************************************************************
*/
int tpia_Legendre_getFromElement( statusMessageReporting *smr, xData_element *LegendreElement, tpia_Legendre *Legendre ) {

    int i, status = 0;
    xData_elementList *list;

    xData_addToAccessed( smr, LegendreElement, 1 );
    if( ( tpia_frame_setFromElement( smr, LegendreElement, 4, &(Legendre->frame) ) ) != 0 ) return( 1 );
    list = xData_getElementsByTagNameAndSort( smr, LegendreElement, "l", NULL, NULL );
    if( list->n == 0 ) {
        status = 1;
        tpia_misc_setMessageError_Element( smr, NULL, LegendreElement, __FILE__, __LINE__, 1, "Legendre element does not contain any l elements" ); }
    else {
        //if( ( Legendre->binned.ls = xData_malloc2( smr, list->n * sizeof( tpia_EqualProbableBinSpectra ), 1, "ls" ) ) != NULL ) {
        if( ( Legendre->binned.ls = (tpia_EqualProbableBinSpectra*) xData_malloc2( smr, list->n * sizeof( tpia_EqualProbableBinSpectra ), 1, "ls" ) ) != NULL ) {
            Legendre->binned.numberOfLs = 0;
            for( i = 0; i < list->n; i++ ) {
                if( ( status = _tpia_Legendre_getOrder( smr, list->items[i].element, Legendre, &(Legendre->binned.ls[i]) ) ) != 0 ) break;
                Legendre->binned.numberOfLs++;
            }
        }
    }
    xData_freeElementList( smr, list );
    return( status );
}
/*
************************************************************
*/
//static int _tpia_Legendre_getOrder( statusMessageReporting *smr, xData_element *LegendreOrder, tpia_Legendre *Legendre, tpia_EqualProbableBinSpectra *l ) {
static int _tpia_Legendre_getOrder( statusMessageReporting *smr, xData_element *LegendreOrder, tpia_Legendre *, tpia_EqualProbableBinSpectra *l ) {

    int status = 1;
    xData_Int n, nBins, lOrder;

    if( xData_convertAttributeTo_xData_Int( smr, LegendreOrder, "value", &lOrder ) != 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, LegendreOrder, __FILE__, __LINE__, 1, "missing or invalid nBins attribute for angular bins" ); }
    else {
        l->iValue = lOrder;
        l->dValue = 0.;
        if( ( l->energies = tpia_misc_getEqualProbableBin( smr, LegendreOrder, &n, &nBins ) ) != NULL ) {
            status = 0;
            l->nBins = (int) nBins;
            l->numberOfEs = (int) n;
        }
    }
    return( status );
}
/*
************************************************************
*/
int tpia_Legendre_SampleEp( statusMessageReporting *smr, tpia_Legendre *Legendre, int sampleMu, tpia_decaySamplingInfo *decaySamplingInfo ) {

    tpia_EqualProbableBinSpectra *binned_l0 = Legendre->binned.ls;
    double Ep;

/*
Currently, only l = 0, equal probable binning is supported.
Need to return frame info for Ep, mu, also.
*/
    if( Legendre->binned.numberOfLs > 0 ) {
        if( sampleMu ) decaySamplingInfo->mu = 2. * tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState ) - 1.;
#ifndef WIN32
        if( decaySamplingInfo->mu <= -1 ) decaySamplingInfo->mu = nextafter( -1., 0. );
#endif
#ifdef WIN32
        if( decaySamplingInfo->mu <= -1 ) decaySamplingInfo->mu = _nextafter( -1., 0. );
#endif

        tpia_misc_sampleEqualProbableBin( smr, decaySamplingInfo, decaySamplingInfo->e_in, binned_l0->nBins, binned_l0, &Ep );
        /* ??? Need to check that e_in > E_Threshold */
        decaySamplingInfo->Ep = Ep; }
    else {
        return( 1 );
    }

    return( 0 );
}

#if defined __cplusplus
}
#endif
