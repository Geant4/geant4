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

#include "Randomize.hh"

#include <float.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <tpia_target.h>
#include <tpia_misc.h>

#include <string>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

struct ZSymbol {
    int Z;
    //char *symbol;
    std::string symbol; 
};

static const struct ZSymbol ZSymbols[] = {    {   0, "n"  },  {   1, "H"  },  {   2, "He" },  {   3, "Li" },  {   4, "Be" },  {   5, "B"  },  {   6, "C"  },
        {   7, "N"  },  {   8, "O"  },  {   9, "F"  },  {  10, "Ne" },  {  11, "Na" },  {  12, "Mg" },  {  13, "Al" },  {  14, "Si" },  {  15, "P"  },
        {  16, "S"  },  {  17, "Cl" },  {  18, "Ar" },  {  19, "K"  },  {  20, "Ca" },  {  21, "Sc" },  {  22, "Ti" },  {  23, "V"  },  {  24, "Cr" },
        {  25, "Mn" },  {  26, "Fe" },  {  27, "Co" },  {  28, "Ni" },  {  29, "Cu" },  {  30, "Zn" },  {  31, "Ga" },  {  32, "Ge" },  {  33, "As" },
        {  34, "Se" },  {  35, "Br" },  {  36, "Kr" },  {  37, "Rb" },  {  38, "Sr" },  {  39, "Y"  },  {  40, "Zr" },  {  41, "Nb" },  {  42, "Mo" },
        {  43, "Tc" },  {  44, "Ru" },  {  45, "Rh" },  {  46, "Pd" },  {  47, "Ag" },  {  48, "Cd" },  {  49, "In" },  {  50, "Sn" },  {  51, "Sb" },
        {  52, "Te" },  {  53, "I"  },  {  54, "Xe" },  {  55, "Cs" },  {  56, "Ba" },  {  57, "La" },  {  58, "Ce" },  {  59, "Pr" },  {  60, "Nd" },
        {  61, "Pm" },  {  62, "Sm" },  {  63, "Eu" },  {  64, "Gd" },  {  65, "Tb" },  {  66, "Dy" },  {  67, "Ho" },  {  68, "Er" },  {  69, "Tm" },
        {  70, "Yb" },  {  71, "Lu" },  {  72, "Hf" },  {  73, "Ta" },  {  74, "W"  },  {  75, "Re" },  {  76, "Os" },  {  77, "Ir" },  {  78, "Pt" },
        {  79, "Au" },  {  80, "Hg" },  {  81, "Tl" },  {  82, "Pb" },  {  83, "Bi" },  {  84, "Po" },  {  85, "At" },  {  86, "Rn" },  {  87, "Fr" },
        {  88, "Ra" },  {  89, "Ac" },  {  90, "Th" },  {  91, "Pa" },  {  92, "U"  },  {  93, "Np" },  {  94, "Pu" },  {  95, "Am" },  {  96, "Cm" },
        {  97, "Bk" },  {  98, "Cf" },  {  99, "Es" },  { 100, "Fm" },  { 101, "Md" },  { 102, "No" },  { 103, "Lr" },  { 104, "Rf" },  { 105, "Db" },
        { 106, "Sg" }, { 107, "Bh" },  { 108, "Hs" },  { 109, "Mt" } };

/*
************************************************************
*/
int tpia_misc_NumberOfZSymbols( void ) {

    return( sizeof( ZSymbols ) / sizeof( struct ZSymbol ) );
}
/*
************************************************************
*/
const char *tpia_misc_ZToSymbol( int iZ ) {

    if( ( iZ < 0 ) || ( iZ >= tpia_misc_NumberOfZSymbols( ) ) ) return( NULL );
    //return( ZSymbols[iZ].symbol );
    return( ZSymbols[iZ].symbol.c_str() );
}
/*
************************************************************
*/
int tpia_misc_symbolToZ( const char *Z ) {

    int i, n = tpia_misc_NumberOfZSymbols( );

    for( i = 0; i < n; i++ ) {
        //if( !strcmp( Z, ZSymbols[i].symbol ) ) return( ZSymbols[i].Z );
        if( !strcmp( Z, ZSymbols[i].symbol.c_str() ) ) return( ZSymbols[i].Z );
    }
    return( -1 );
}
/*
************************************************************
*/
int tpia_miscNameToZAm( statusMessageReporting *smr, const char *name, int *Z, int *A, int *m ) {

    int i, n;
    const char *p;
    char s[1024] = "", *q, *e;   /* Note 1) routine will fail when parts of a particle name can be longer than 1024. */

    n = sizeof( s ) - 1;

    *Z = *A = *m = 0;
    if( !strncmp( "FissionProduct", name, 14 ) ) {
        *Z = 99;
        *A = 120;
        return( 0 );
    }
    if( !strcmp( "gamma", name ) ) return( 0 );
    for( p = name, q = s, i = 0; ( *p != '_' ) && ( i != n ); p++, q++, i++ ) *q = *p;
    if( i == n ) {              /* See note 1 above. */
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Failed to find first '_' in particle name %s", name ); }
    else {
        *q = 0;
        if( ( *Z = tpia_misc_symbolToZ( s ) ) < 0 ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Particle %s's symbol = '%s' not found", name, s ); }
        else {                  /* Getting here implies that *p == '_'. */
            for( p++, q = s; ( *p != '_' ) && ( *p != 0 ) && ( i != n ); p++, q++, i++ ) *q = *p;
            if( i == n ) {      /* See note 1 above. */
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Failed to find second '_' in particle name %s", name ); }
            else {
                *q = 0;
                if( strcmp( s, "natural" ) == 0 ) {
                    e = s;
                    while( *e ) e++; }
                else {
                    *A = (int) strtol( s, &e, 10 );
                }
                if( *e != 0 ) {
                    smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Failed to convert A to integer in particle name %s", name ); }
                else {          /* Getting here implies that *p == '_' or 0. */
                    *m = 0;
                    if( *p == '_' ) {
                        p++;
                        if( *p != 'm' ) {
                            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Particle name %s missing meta-stable label 'm'", name ); }
                        else {
                            p++;
                            *m = (int) strtol( p, &e, 10 );
                            if( *e != 0 ) smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Failed to convert m to integer in particle name %s", name );
                        }
                    }
                }
            }
        }
    }

    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
char *tpia_misc_pointerToAttributeIfAllOk( statusMessageReporting *smr, xData_element *element, const char *path, int required, 
        xData_attributionList *attributes, const char *name, const char *file, int line ) {

    char *value;

    if( !smr_isOk( smr ) ) return( NULL );
    if( ( value = xData_getAttributesValue( attributes, name ) ) == NULL ) {
        if( required ) {
            if( element != NULL ) {
                tpia_misc_setMessageError_Element( smr, NULL, element, file, line, 1, "element does not have attribute named %s", name ); }
            else {
                smr_setMessageError( smr, NULL, file, line, 1, "element does not have attribute named %s for file = %d", name, path );
            }
        }
    }
    return( value );
}
/*
************************************************************
*/
int tpia_misc_setMessageError_Element( statusMessageReporting *smr, void *userInterface, xData_element *element, const char *file, int line, int code, 
    const char *fmt, ... ) {

    int status = 0;
    va_list args;
    char *msg;

    va_start( args, fmt );
    msg = smr_vallocateFormatMessage( fmt, &args );
    va_end( args );
    if( msg == NULL ) {
        status = 1;
        va_start( args, fmt );
        smr_vsetMessageError( smr, userInterface, file, line, code, fmt, &args );
        va_end( args ); }
    else {
        status = smr_setMessageError( smr, userInterface, file, line, code, "%s for element %s at line %d column %d", msg, element->fullName, 
            (int) element->docInfo.line, (int) element->docInfo.column );
        free( msg );
    }
    return( status );
}
/*
************************************************************
*/
xData_Int tpia_misc_binarySearch( xData_Int n, double *ds, double d ) {

    xData_Int imin = 0, imid, imax = n - 1;

    if( d < ds[0] ) return( -2 );
    if( d > ds[n-1] ) return( -1 );
    while( 1 ) {
        imid = ( imin + imax ) >> 1;
        if( imid == imin ) break;
        if( d < ds[imid] ) {
            imax = imid; }
        else {
            imin = imid;
        }
    }
    return( imin );
}
/*
************************************************************
*/
double *tpia_misc_get2dx_y_data( statusMessageReporting *smr, xData_element *element, xData_Int *length ) {

    xData_element *xDataElement;
    double *d = NULL;
    xData_Int length_;

    xData_addToAccessed( smr, element, 1 );
    //if( ( xDataElement = xData_getOneElementByTagName( smr, element, "xData", 1 ) ) != NULL ) {
    if( ( xDataElement = xData_getOneElementByTagName( smr, element, (char*)"xData", 1 ) ) != NULL ) {
        xData_addToAccessed( smr, xDataElement, 1 );
        if( xData_is_2d_xy( smr, &(xDataElement->xDataTypeInfo), 1 ) ) {
            d = xData_2d_xy_allocateCopyData( smr, xDataElement, &length_ );
            *length = (xData_Int) length_;
        }
    }
    return( d );
}
/*
************************************************************
*/
double *tpia_misc_get2dxindex_y_data( statusMessageReporting *smr, xData_element *element, xData_Int *start, xData_Int *end, double *xValues ) {

    xData_element *xDataElement;
    double *d = NULL;

    xData_addToAccessed( smr, element, 1 );
    //if( ( xDataElement = xData_getOneElementByTagName( smr, element, "xData", 1 ) ) != NULL ) {
    if( ( xDataElement = xData_getOneElementByTagName( smr, element, (char*) "xData", 1 ) ) != NULL ) {
        xData_addToAccessed( smr, xDataElement, 1 );
        if( xData_is_2d_xindex_y( smr, &(xDataElement->xDataTypeInfo), 1 ) ) {
            if( start != NULL ) *start = xDataElement->xDataTypeInfo.start;
            if( end != NULL ) *end = xDataElement->xDataTypeInfo.end;
            d = xData_2d_xindex_y_toFilledYs( smr, xDataElement, xValues );
        }
    }
    return( d );
}
/*
************************************************************
*/  
double *tpia_misc_get2d_xShared_yHistogram_data( statusMessageReporting *smr, xData_element *element, xData_Int *start, xData_Int *end, xData_Int *length ) {

    xData_Int i;
    xData_element *xDataElement;
    double *d = NULL;

    xData_addToAccessed( smr, element, 1 );
    //if( ( xDataElement = xData_getOneElementByTagName( smr, element, "xData", 1 ) ) != NULL ) {
    if( ( xDataElement = xData_getOneElementByTagName( smr, element, (char*) "xData", 1 ) ) != NULL ) {
        xData_addToAccessed( smr, xDataElement, 1 );
        if( ( d = xData_2d_xShared_yHistogram_copyData( smr, xDataElement, &i ) ) != NULL ) {
            if( start != NULL ) *start = xDataElement->xDataTypeInfo.start;
            if( end != NULL ) *end = xDataElement->xDataTypeInfo.end;
            if( length != NULL ) *length = xDataElement->xDataTypeInfo.length;
        }
    }
    return( d );
}
/*
************************************************************
*/  
int tpia_misc_get2d_xShared_yHistogram_data_Grouped( statusMessageReporting *smr, xData_element *element, tpia_1dData *group ) {

    if( ( group->data = tpia_misc_get2d_xShared_yHistogram_data( smr, element, &(group->start), &(group->end), &(group->length) ) ) == NULL ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
//double tpia_misc_getPointwiseCrossSectionAtE( statusMessageReporting *smr, tpia_1dData *crossSection, double *energyGrid, xData_Int index, double e_in ) {
double tpia_misc_getPointwiseCrossSectionAtE( statusMessageReporting *, tpia_1dData *crossSection, double *energyGrid, xData_Int index, double e_in ) {

    double xsec = 0.0, e1, e2;

    if( ( index >= crossSection->start ) && ( index < crossSection->end ) ) {
        e1 = energyGrid[index];
        e2 = energyGrid[index + 1];
        index -= crossSection->start;
        if( e1 == e2 ) {                            /* Allow for future where step function may be allowed. */
            xsec = 0.5 * ( crossSection->data[index] + crossSection->data[index + 1] ); }
        else {
            xsec = ( crossSection->data[index] * ( e2 - e_in ) + crossSection->data[index + 1] * ( e_in - e1 ) ) / ( e2 - e1 );
        }
    }
    return( xsec );
}
/*
************************************************************
*/
tpia_EqualProbableBinSpectrum *tpia_misc_getEqualProbableBin( statusMessageReporting *smr, xData_element *parent, xData_Int *n, xData_Int *nBins ) {

    xData_element *element;

    xData_addToAccessed( smr, parent, 1 );
    //if( ( element = xData_getOneElementByTagName( smr, parent, "equalProbableBins", 0 ) ) == NULL ) return( NULL );
    if( ( element = xData_getOneElementByTagName( smr, parent, (char*) "equalProbableBins", 0 ) ) == NULL ) return( NULL );
    if( xData_convertAttributeTo_xData_Int( smr, element, "nBins", nBins ) != 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "missing or invalid nBins attribute" );
        return( NULL );
    }
    return( tpia_misc_getEqualProbableBins( smr, element, "energy", *nBins, n ) );
}
/*
************************************************************
*/
tpia_EqualProbableBinSpectrum *tpia_misc_getEqualProbableBins( statusMessageReporting *smr, xData_element *parent, const char *name, xData_Int nBins, 
    xData_Int *n ) {
        

    int i, j;
    xData_Int index, size;
    xData_elementList *list;
    xData_element *element, *xData;
    double *d;
    tpia_EqualProbableBinSpectrum *epbs = NULL, *epb;

    xData_addToAccessed( smr, parent, 1 );
    list = xData_getElementsByTagNameAndSort( smr, parent, name, NULL, NULL );
    if( list->n == 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, parent, __FILE__, __LINE__, 1, "bins does not contain any %s elements", name ); }
    else {
        *n = list->n;
        size = list->n * ( sizeof( tpia_EqualProbableBinSpectrum ) + ( nBins + 1 ) * sizeof( double ) );
        //if( ( epbs = xData_malloc2( smr, size, 0, "energies" ) ) != NULL ) {
        if( ( epbs = (tpia_EqualProbableBinSpectrum*) xData_malloc2( smr, size, 0, "energies" ) ) != NULL ) {
            d = (double *) &(epbs[list->n]);
            for( i = 0, epb = epbs; i < list->n; i++, epb++ ) {    /* Loop to test nBins and index are proper. */
                element = list->items[i].element;
                xData_addToAccessed( smr, element, 1 );
                if( xData_convertAttributeTo_xData_Int( smr, element, "index", &index ) != 0 ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "missing or invalid index attribute" );
                    //epbs = xData_free( smr, epbs );
                    epbs = (tpia_EqualProbableBinSpectrum*) xData_free( smr, epbs );
                    break;
                }
                if( index != i ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "index = %lld is not incremental", index );
                    //epbs = xData_free( smr, epbs );
                    epbs = (tpia_EqualProbableBinSpectrum*) xData_free( smr, epbs );
                    break;
                }
                if( ( j = xData_convertAttributeToDouble( smr, element, "value", &(epb->value) ) ) != 0 ) {
                    if( j == 1 ) {
                        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "element does not have value attribute" ); }
                    else {
                        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "failed to convert value attribute to double" );
                    }
                    //epbs = xData_free( smr, epbs );
                    epbs = (tpia_EqualProbableBinSpectrum*) xData_free( smr, epbs );
                    break;
                }
                if( ( xData = xData_getElements_xDataElement( smr, element ) ) == NULL ) {
                    //epbs = xData_free( smr, epbs );
                    epbs = (tpia_EqualProbableBinSpectrum*) xData_free( smr, epbs );
                    break;
                }
                xData_addToAccessed( smr, xData, 1 );
                epb->index = index;
                epb->nBins = nBins;
                epb->bins = d;
                if( xData_1d_x_copyData( smr, xData, ( nBins + 1 ) * sizeof( double ), d ) != 0 ) {
                    //epbs = xData_free( smr, epbs );
                    epbs = (tpia_EqualProbableBinSpectrum*) xData_free( smr, epbs );
                    break;
                }
                d += nBins + 1;
            }
        }
    }
    xData_freeElementList( smr, list );
    return( epbs );
}
/*
************************************************************
*/
double tpia_misc_drng( double (*rng)( void * ), void *rngState ) {

    double r;

    if( rng != NULL ) {
        r = rng( rngState ); }
    else {
        //r = drand48( );
        r = CLHEP::HepRandom::getTheEngine()->flat();
          
    }
    return( r );
}
/*
************************************************************
*/
//int tpia_misc_sampleEqualProbableBin( statusMessageReporting *smr, tpia_decaySamplingInfo *decaySamplingInfo, double e_in, int nBins, 
int tpia_misc_sampleEqualProbableBin( statusMessageReporting *, tpia_decaySamplingInfo *decaySamplingInfo, double e_in, int nBins, 
        tpia_EqualProbableBinSpectra *binned, double *value_ ) {

    int i, j, index1, index2, method = 0;
    double fE = 1., r, value1, value2, value3, P12, P23, value = -2.;

    if( e_in <= binned->energies[0].value ) {
                index1 = 0;
                index2 = 0; }
    else if( e_in >= binned->energies[binned->numberOfEs-1].value ) {
                index1 = binned->numberOfEs - 1;
                index2 = binned->numberOfEs - 1; }
    else {
            for( i = 0; i < binned->numberOfEs - 1; i++ ) {
                if( e_in <= binned->energies[i].value ) break;
            }
            i--;
            index1 = i;
            index2 = i;
            if( e_in != binned->energies[i].value ) {
                index2++;
                fE = ( e_in - binned->energies[i].value ) / ( binned->energies[i+1].value - binned->energies[i].value );
            }
    }
    r = tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );
    j = (int) (r * nBins);
    if( j >= nBins ) j = nBins - 1;
    if( j < 0 ) j = 0;
    r = tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );          // Do not change r until after Point1 below.
    if( tpia_samplingMethods_isLinear( decaySamplingInfo->samplingMethods->angular_equalProbableBinMethod ) ) {
        method = 1;
        if( ( ( j == 0 ) && ( r <= 0.5 ) ) || ( j == ( nBins - 1 ) && r > 0.5 ) ) method = 0;
    }
    if( method == 0 ) {                 /* Constant interpolaton. */
        value1 = ( 1. - fE ) * binned->energies[index1].bins[j] + fE * binned->energies[index2].bins[j];
        value2 = ( 1. - fE ) * binned->energies[index1].bins[j+1] + fE * binned->energies[index2].bins[j+1];
        value = ( 1. - r ) * value1 + r * value2; }
    else {                              /* Linear interpolation. */
        if( r <= 0.5 ) j--;             /* Point1: Above test insures that j is not 0 (nBins-1) if r <= 0.5 (> 0.5); */
        value1 = ( 1. - fE ) * binned->energies[index1].bins[j] + fE * binned->energies[index2].bins[j];
        value2 = ( 1. - fE ) * binned->energies[index1].bins[j+1] + fE * binned->energies[index2].bins[j+1];
        value3 = ( 1. - fE ) * binned->energies[index1].bins[j+2] + fE * binned->energies[index2].bins[j+2];
//
//TK140602 Modified for protecting divided by 0 BEGIN
        if ( value1 == value2 && value2 == value3 ) {
          value = value1;
        } else {

        if ( value2 != value1 ) 
        P12 = 1. / ( value2 - value1 );
        else
        P12 =FLT_MAX;

        if ( value3 != value2 ) 
        P23 = 1. / ( value3 - value2 );
        else
        P23 =FLT_MAX;

        r = tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );
        if( 0.25 * ( 1.0 + 2.0 * ( value2 - value1 ) / ( value3 - value1 ) ) > r ) {
            P23 = 2. / ( value3 - value1 );
            value3 = value2; }
        else {
            P12 = 2. / ( value3 - value1 );
            value1 = value2;
        }
        r = tpia_misc_drng( decaySamplingInfo->rng, decaySamplingInfo->rngState );
        if( P23 != P12 ) r = ( -P12 + std::sqrt( P12 * P12 * ( 1. - r ) + r * P23 * P23 ) ) / ( P23 - P12 );
        value = 0.5 * ( value1 + value2 + r * ( value3 - value1 ) );
        }
//TK140602 Modified for protecting divided by 0 END
    }
    *value_ = value;
    return( 0 );
}

#if defined __cplusplus
}
#endif
