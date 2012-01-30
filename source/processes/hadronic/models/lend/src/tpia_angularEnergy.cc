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

#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int tpia_angularEnergy_initialize( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy ) {

    memset( angularEnergy, 0, sizeof( tpia_angularEnergy ) );
    if( tpia_frame_setFromString( smr, "", "", 0, &(angularEnergy->frame) ) ) return( 1 );
    angularEnergy->binned.numberOfEs = 0;
    angularEnergy->binned.energies = NULL;
    return( 0 );
}
/*
************************************************************
*/
int tpia_angularEnergy_release( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy ) {

    int i;

    for( i = 0; i < angularEnergy->binned.numberOfEs; i++ ) xData_free( smr, angularEnergy->binned.energies[i].energies );
    //angularEnergy->binned.energies = xData_free( smr, angularEnergy->binned.energies );
    angularEnergy->binned.energies = (tpia_EqualProbableBinSpectra*) xData_free( smr, angularEnergy->binned.energies );
    tpia_angularEnergy_initialize( smr, angularEnergy );
    return( 0 );
}
/*
************************************************************
*/
int tpia_angularEnergy_getFromElement( statusMessageReporting *smr, xData_element *angularEnergyElement, tpia_angularEnergy *angularEnergy ) {

    int i, j, status = 1;
    xData_Int nBins, size, n, index;
    xData_element *epbElement, *element;
    xData_elementList *list;
    tpia_EqualProbableBinSpectra *mus;

    xData_addToAccessed( smr, angularEnergyElement, 1 );
    if( ( tpia_frame_setFromElement( smr, angularEnergyElement, 4, &(angularEnergy->frame) ) ) != 0 ) return( 1 );
    //if( ( epbElement = xData_getOneElementByTagName( smr, angularEnergyElement, "equalProbableBins", 0 ) ) == NULL ) return( 1 );
    if( ( epbElement = xData_getOneElementByTagName( smr, angularEnergyElement, (char*)"equalProbableBins", 0 ) ) == NULL ) return( 1 );
    xData_addToAccessed( smr, epbElement, 1 );
    if( xData_convertAttributeTo_xData_Int( smr, epbElement, "nBins", &nBins ) != 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, epbElement, __FILE__, __LINE__, 1, "missing or invalid nBins attribute" );
        return( 1 );
    }
    list = xData_getElementsByTagNameAndSort( smr, epbElement, "energy", NULL, NULL );
    if( list->n == 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, epbElement, __FILE__, __LINE__, 1, "bins does not contain any energy elements" ); }
    else {
        size = list->n * sizeof( tpia_EqualProbableBinSpectra );
        angularEnergy->binned.nBins = (int) nBins;
        //if( ( angularEnergy->binned.energies = xData_malloc2( smr, size, 1, "energies" ) ) != NULL ) {
        if( ( angularEnergy->binned.energies = (tpia_EqualProbableBinSpectra*) xData_malloc2( smr, size, 1, "energies" ) ) != NULL ) {
            status = 0;
            for( i = 0, mus = angularEnergy->binned.energies; i < list->n; i++, mus++ ) {
                mus->iValue = 0;
                element = list->items[i].element;
                if( xData_convertAttributeTo_xData_Int( smr, element, "index", &index ) != 0 ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "missing or invalid index attribute" );
                    status = 1;
                    break;
                }
                if( index != i ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "index = %lld is not incremental", index );
                    status = 1;
                    break;
                }
                if( ( j = xData_convertAttributeToDouble( smr, element, "value", &(mus->dValue) ) ) != 0 ) {
                    if( j == 1 ) {
                        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "element does not have value attribute" ); }
                    else {
                        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "failed to convert value attribute to double" );
                    }
                    status = 1;
                    break;
                }
                if( ( mus->energies = tpia_misc_getEqualProbableBins( smr, element, "mu", nBins, &n ) ) == NULL ) {
                    status = 1;
                    break; }
                else {
                    mus->numberOfEs = n;
                }
                angularEnergy->binned.numberOfEs++;
            }
        }
    }
    xData_freeElementList( smr, list );
    return( status );
}
/*
************************************************************
*/
int tpia_angularEnergy_SampleEp( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy, tpia_decaySamplingInfo *decaySamplingInfo ) {

    int iE1, iE2;
    tpia_EqualProbableBinSpectra *energies = angularEnergy->binned.energies;
    double Ep, Ep1, Ep2, f, e_in = decaySamplingInfo->e_in;
/*
Currently, only equal probable binning is supported.
Need to return frame info for Ep also.
*/

    if( !smr_isOk( smr ) ) return( 1 );
    if( angularEnergy->binned.numberOfEs == 0 ) return( 1 );
    for( iE2 = 0; iE2 < angularEnergy->binned.numberOfEs; iE2++ ) if( energies[iE2].dValue >= e_in ) break;
    if( iE2 == 0 ) {
        iE1 = iE2; }
    else if( iE2 == angularEnergy->binned.numberOfEs ) {
        iE1 = iE2 = angularEnergy->binned.numberOfEs - 1; }
    else {
        iE1 = iE2 - 1;
    }
    tpia_misc_sampleEqualProbableBin( smr, decaySamplingInfo, decaySamplingInfo->mu, angularEnergy->binned.nBins, &(energies[iE1]), &Ep1 );
    if( iE1 == iE2 ) {
        Ep = Ep1; }
    else {
        tpia_misc_sampleEqualProbableBin( smr, decaySamplingInfo, decaySamplingInfo->mu, angularEnergy->binned.nBins, &(energies[iE2]), &Ep2 );
        f = ( energies[iE2].dValue - e_in ) / ( energies[iE2].dValue - energies[iE1].dValue );
        Ep = f * Ep1 + ( 1 - f ) * Ep2;
    }
    decaySamplingInfo->Ep = Ep;

    return( 0 );
}

#if defined __cplusplus
}
#endif
