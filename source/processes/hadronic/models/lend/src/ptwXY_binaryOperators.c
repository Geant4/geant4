/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <math.h>
#include <float.h>

#include "ptwXY.h"

static double ptwXY_mod2( double v, double m, int pythonMod );
static nfu_status ptwXY_mul2_s_ptwXY( statusMessageReporting *smr, ptwXYPoints *div, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double x1, double y1, double x2, double y2, int level );
static nfu_status ptwXY_div_s_ptwXY( statusMessageReporting *smr, ptwXYPoints *div, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double x1, double y1, double x2, double y2, int level, int isNAN1, int isNAN2 );
static ptwXYPoints *ptwXY_div_ptwXY_forFlats( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int safeDivide );
static nfu_status ptwXY_getValueAtX_ignore_XOutsideDomainError( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double x, double *y );
static nfu_status ptwXY_getValueAtX_signal_XOutsideDomainError( statusMessageReporting *smr, int line, 
        char const *function, ptwXYPoints *ptwXY1, double x, double *y );
/*
************************************************************
*/
nfu_status ptwXY_slopeOffset( statusMessageReporting *smr, ptwXYPoints *ptwXY, double slope, double offset ) { 

    int64_t i, nonOverflowLength;
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ( nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ) ) < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = slope * p->y + offset;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = slope * o->point.y + offset; 
    return( ptwXY->status );
}   
/*
************************************************************
*/
nfu_status ptwXY_add_double( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) {
    
    if( ptwXY_slopeOffset( smr, ptwXY, 1., value ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_sub_doubleFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) {

    if( ptwXY_slopeOffset( smr, ptwXY,  1., -value ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_sub_fromDouble( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) {

    if( ptwXY_slopeOffset( smr, ptwXY, -1.,  value ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_mul_double( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) {

    if( ptwXY_slopeOffset( smr, ptwXY, value, 0. ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_div_doubleFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) { 

    if( value == 0. ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." );
        ptwXY->status = nfu_divByZero; }
    else {
        if( ptwXY_slopeOffset( smr, ptwXY, 1. / value, 0. ) != nfu_Okay )
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_div_fromDouble( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value ) {
/*
*   This does not do any infilling and it should?????????
*/

    int64_t i, nonOverflowLength;
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ( nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ) ) < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }
    if( ptwXY->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Other interpolation not allowed." );
        return( nfu_otherInterpolation );
    }

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) if( p->y == 0. ) ptwXY->status = nfu_divByZero;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) if( o->point.y == 0. ) ptwXY->status = nfu_divByZero;
    if( ptwXY->status == nfu_divByZero ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." ); }
    else {
        for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = value / p->y;
        for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = value / o->point.y; 
    }
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_mod( statusMessageReporting *smr, ptwXYPoints *ptwXY, double m, int pythonMod ) { 

    int64_t i, nonOverflowLength;
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ( nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ) ) < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( m == 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." );
        ptwXY->status = nfu_divByZero;
    }

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = ptwXY_mod2( p->y, m, pythonMod );
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = ptwXY_mod2( o->point.y, m, pythonMod );
    return( ptwXY->status );
}
/*
************************************************************
*/
static double ptwXY_mod2( double v, double m, int pythonMod ) {

    double r = fmod( fabs( v ), fabs( m ) );

    if( pythonMod ) {
        if( ( v * m ) < 0. ) r = fabs( m ) - fabs( r );
        if( m < 0. ) r *= -1.; }
    else {
        if( v < 0. ) r *= -1.;
    }

    return( r );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_binary_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double v1, double v2, double v1v2 ) {

    int64_t i;
    int unionOptions = ptwXY_union_fill | ptwXY_union_mergeClosePoints;
    double y;
    ptwXYPoints *ptwXYNew;
    ptwXYPoint *p;

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Source1: Other interpolation not allowed." );
        return( NULL );
    }
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Source2: Other interpolation not allowed." );
        return( NULL );
    }

    if( ptwXY1->interpolation != ptwXY2->interpolation ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, 
                "Source1 interpolation '%s' not same as Source2 interpolation '%s'.",
                    ptwXY1->interpolationString, ptwXY2->interpolationString );
        return( NULL );
    }

    if( ( ptwXY1->interpolation != ptwXY_interpolationLinLin ) && ( ptwXY1->interpolation != ptwXY_interpolationLinLog ) &&
            ( ptwXY1->interpolation != ptwXY_interpolationFlat ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation,
                "Only '%s' and '%s' interpolation supported, not '%s' interpolation.", ptwXY_interpolationToString( ptwXY_interpolationLinLin ),
                ptwXY_interpolationToString( ptwXY_interpolationFlat ), ptwXY1->interpolationString );
        return( NULL );
    }

    if( ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 ) != nfu_Okay ) {
        double domainMin1, domainMax1, domainMin2, domainMax2;

        ptwXY_domainMin( NULL, ptwXY1, &domainMin1 );
        ptwXY_domainMax( NULL, ptwXY1, &domainMax1 );
        ptwXY_domainMin( NULL, ptwXY2, &domainMin2 );
        ptwXY_domainMax( NULL, ptwXY2, &domainMax2 );

        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, "Domains not mutual (%.17e, %.17e) vs (%.17e, %.17e).",
                domainMin1, domainMax1, domainMin2, domainMax2 );
        return( NULL );
    }

    if( ( ptwXYNew = ptwXY_union( smr, ptwXY1, ptwXY2, unionOptions ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." ); }
    else {
        for( i = 0, p = ptwXYNew->points; i < ptwXYNew->length; i++, p++ ) {
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, p->x, &y ) != nfu_Okay ) goto Err;
            p->y = v1 * p->y + v2 * y + v1v2 * y * p->y;
        }
    }
    return( ptwXYNew );
Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( ptwXYNew ) ptwXY_free( ptwXYNew );
    return( NULL );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_add_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    ptwXYPoints *sum;

    if( ptwXY1->length == 0 ) {
        sum = ptwXY_clone( smr, ptwXY2 ); }
    else if( ptwXY2->length == 0 ) {
        sum = ptwXY_clone( smr, ptwXY1 ); }
    else {
        sum = ptwXY_binary_ptwXY( smr, ptwXY1, ptwXY2, 1., 1., 0. );
    }
    if( sum == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( sum );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_sub_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    ptwXYPoints *diff = NULL;

    if( ptwXY1->length == 0 ) {
        diff = ptwXY_clone( smr, ptwXY2 );
        if( diff != NULL ) {
            if( ptwXY_neg( smr, diff ) != nfu_Okay ) diff = ptwXY_free( diff );
        } }
    else if( ptwXY2->length == 0 ) {
        diff = ptwXY_clone( smr, ptwXY1 ); }
    else {
        diff = ptwXY_binary_ptwXY( smr, ptwXY1, ptwXY2, 1., -1., 0. );
    }
    if( diff == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( diff );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_mul_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    ptwXYPoints *mul;

    if( ptwXY1->length == 0 ) {
        mul = ptwXY_clone( smr, ptwXY1 ); }
    else if( ptwXY2->length == 0 ) {
        mul = ptwXY_clone( smr, ptwXY2 ); }
    else {
        mul = ptwXY_binary_ptwXY( smr, ptwXY1, ptwXY2, 0., 0., 1. );
    }
    if( mul == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( mul );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_mul2_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    int64_t i, length;
    ptwXYPoints *mul = NULL;
    int found;
    double x1, y1, x2, y2, u1, u2, v1, v2, xz1 = 0, xz2 = 0, x;

    if( ( mul = ptwXY_mul_ptwXY( smr, ptwXY1, ptwXY2 ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
    if( mul->length == 0 ) return( mul );
    if( ptwXY1->interpolation == ptwXY_interpolationFlat ) return( mul );
    if( ptwXY2->interpolation == ptwXY_interpolationFlat ) return( mul );

    length = mul->length - 1;
    if( length > 0 ) {
        x2 = mul->points[length].x;
        for( i = length - 1; i >= 0; i-- ) {             /* Find and add y zeros not currently in mul's. */
            x1 = mul->points[i].x;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x1, &u1 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x2, &u2 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, x1, &v1 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, x2, &v2 ) != nfu_Okay ) goto Err;
            found = 0;
            if( u1 * u2 < 0 ) {
                xz1 = ( u1 * x2 - u2 * x1 ) / ( u1 - u2 );
                if( ptwXY_setValueAtX( smr, mul, xz1, 0. ) != nfu_Okay ) goto Err;
                found = 1;
            }
            if( v1 * v2 < 0 ) {
                xz2 = ( v1 * x2 - v2 * x1 ) / ( v1 - v2 );
                if( ptwXY_setValueAtX( smr, mul, xz2, 0. ) != nfu_Okay ) goto Err;
                found += 1;
            }
            if( found > 1 ) {
                x = 0.5 * ( xz1 + xz2 );
                if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x, &u1 ) != nfu_Okay ) goto Err;
                if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, x, &v1 ) != nfu_Okay ) goto Err;
                if( ptwXY_setValueAtX( smr, mul, x, u1 * v1 ) != nfu_Okay ) goto Err;
            }
            x2 = x1;
        }

        if( ptwXY_simpleCoalescePoints( smr, mul ) != nfu_Okay ) goto Err;
        length = mul->length;
        x2 = mul->points[mul->length-1].x;
        y2 = mul->points[mul->length-1].y;
        for( i = mul->length - 2; i >= 0; i-- ) {           /* Make interpolation fit accuracy. Work backwards so new */
            x1 = mul->points[i].x;                          /* points will not mess up loop. */
            y1 = mul->points[i].y;
            if( ptwXY_mul2_s_ptwXY( smr, mul, ptwXY1, ptwXY2, x1, y1, x2, y2, 0 ) != nfu_Okay ) goto Err;
            x2 = x1;
            y2 = y1;
        }
        ptwXY_update_biSectionMax( mul, (double) length );
    }
    return( mul );

Err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    if( mul ) ptwXY_free( mul );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_mul2_s_ptwXY( statusMessageReporting *smr, ptwXYPoints *mul, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double x1, double y1, double x2, double y2, int level ) {

    double u1, u2, v1, v2, x, y, yp, dx, a1, a2;
    nfu_status status;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( fabs( x1 ) + fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= mul->biSectionMax ) return( nfu_Okay );
    level++;
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x1, &u1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x2, &u2 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, x1, &v1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, x2, &v2 ) ) != nfu_Okay ) return( status );
    if( ( u1 == u2 ) || ( v1 == v2 ) ) return( nfu_Okay );
    a1 = u1 * v1;
    if( y1 == 0 ) a1 = 0.;                                  /* Fix rounding problem. */
    a2 = u2 * v2;
    if( y2 == 0 ) a2 = 0.;                                  /* Fix rounding problem. */
    if( ( a1 == 0. ) || ( a2 == 0. ) ) {                    /* Handle special case of 0 where accuracy can never be met. */
        x = 0.5 * ( x1 + x2 ); }
    else {
        if( ( a1 * a2 < 0. ) ) return( nfu_Okay );          /* Assume rounding error and no point needed as zero */
        a1 = sqrt( fabs( a1 ) );                            /* crossings are set in ptwXY_mul2_ptwXY. */
        a2 = sqrt( fabs( a2 ) );
        x = ( a2 * x1 + a1 * x2 ) / ( a2 + a1 );
    }
    dx = x2 - x1;
    yp = ( u1 * v1 * ( x2 - x ) + u2 * v2 * ( x - x1 ) ) / dx;
    y = ( u1 * ( x2 - x ) + u2 * ( x - x1 ) ) * ( v1 * ( x2 - x ) + v2 * ( x - x1 ) ) / ( dx * dx );
    if( fabs( y - yp ) < fabs( y * mul->accuracy ) ) return( nfu_Okay );
    if( ptwXY_setValueAtX( smr, mul, x, y ) != nfu_Okay ) return( mul->status );
    if( ptwXY_mul2_s_ptwXY( smr, mul, ptwXY1, ptwXY2, x, y, x2, y2, level ) != nfu_Okay ) return( mul->status );
    ptwXY_mul2_s_ptwXY( smr, mul, ptwXY1, ptwXY2, x1, y1, x, y, level );
    return( mul->status );
}
/*
************************************************************
*/
ptwXYPoints *ptwXY_div_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int safeDivide ) {

    int isNAN1, isNAN2;
    int64_t i, j, k, zeros = 0, length, iYs;
    double x1, x2, y1, y2, u1, u2, v1, v2, y, xz, nan = nfu_getNAN( ), s1, s2;
    ptwXYPoints *div = NULL;
    ptwXYPoint *p;
    nfu_status status = ptwXY_simpleCoalescePoints( smr, ptwXY1 );

    if( status != nfu_Okay ) goto Err;
    if( ptwXY_simpleCoalescePoints( smr, ptwXY2 ) != nfu_Okay ) goto Err;

    if( ptwXY1->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Source1: Other interpolation not allowed." );
        return( NULL );
    }
    if( ptwXY2->interpolation == ptwXY_interpolationOther ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_otherInterpolation, "Source2: Other interpolation not allowed." );
        return( NULL );
    }

    if( ( ptwXY1->interpolation == ptwXY_interpolationFlat ) || ( ptwXY1->interpolation == ptwXY_interpolationFlat ) ) {
        div = ptwXY_div_ptwXY_forFlats( smr, ptwXY1, ptwXY2, safeDivide );
        if( div == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( div );
    }

    if( ptwXY_areDomainsMutual( smr, ptwXY1, ptwXY2 ) != nfu_Okay ) {
        double domainMin1, domainMax1, domainMin2, domainMax2;

        ptwXY_domainMin( NULL, ptwXY1, &domainMin1 );
        ptwXY_domainMax( NULL, ptwXY1, &domainMax1 );
        ptwXY_domainMin( NULL, ptwXY2, &domainMin2 );
        ptwXY_domainMax( NULL, ptwXY2, &domainMax2 );

        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_domainsNotMutual, "Domains not mutual (%.17e, %.17e) vs (%.17e, %.17e).",
                domainMin1, domainMax1, domainMin2, domainMax2 );
        return( NULL );
    }

    if( ( div = ptwXY_union( smr, ptwXY1, ptwXY2, ptwXY_union_fill | ptwXY_union_mergeClosePoints ) ) == NULL ) goto Err;
    for( i = 0, p = div->points; i < div->length; i++, p++ ) {
        if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, p->x, &y ) != nfu_Okay ) goto Err;
        if( y == 0. ) {
            if( p->y == 0. ) {
                iYs = 0;
                y1 = 0.;
                y2 = 0.;
                if( i > 0 ) {
                    if( ( status = ptwXY_getSlopeAtX( smr, ptwXY1, p->x, '-', &s1 ) ) != nfu_Okay ) {
                        if( status != nfu_XOutsideDomain ) goto Err;
                        s1 = 0.;
                    }
                    if( ( status = ptwXY_getSlopeAtX( smr, ptwXY2, p->x, '-', &s2 ) ) != nfu_Okay ) goto Err;
                    if( s2 == 0. ) {
                        y1 = nan; }
                    else {
                        y1 = s1 / s2;
                    } 
                    iYs++;
                }
                if( i < ( div->length - 1 ) ) {
                    if( ( status = ptwXY_getSlopeAtX( smr, ptwXY1, p->x, '+', &s1 ) ) != nfu_Okay ) {
                        if( status != nfu_XOutsideDomain ) goto Err;
                        s1 = 0.;
                    }
                    if( ( status = ptwXY_getSlopeAtX( smr, ptwXY2, p->x, '+', &s2 ) ) != nfu_Okay ) goto Err;
                    if( s2 == 0. ) {
                        y2 = nan; }
                    else {
                        y2 = s1 / s2;
                    }
                    iYs++;
                }
                p->y = ( y1 + y2 ) / iYs; 
                if( nfu_isNAN( p->y ) ) zeros++; }
            else {
                if( !safeDivide ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." );
                    goto Err2;
                }
                zeros++;
                p->y = nan;
            } }
        else {
            p->y /= y;
        }
    }
    length = div->length - 1;
    if( length > 0 ) {
        x2 = div->points[length].x;
        for( i = length - 1; i >= 0; i-- ) {             /* Find and add y zeros and NAN not currently in div's. */
            x1 = div->points[i].x;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x1, &u1 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x2, &u2 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x1, &v1 ) != nfu_Okay ) goto Err;
            if( ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x2, &v2 ) != nfu_Okay ) goto Err;
            if( u1 * u2 < 0 ) {
                xz = ( u1 * x2 - u2 * x1 ) / ( u1 - u2 );
                if( ptwXY_setValueAtX( smr, div, xz, 0. ) != nfu_Okay ) goto Err;
            }
            if( v1 * v2 < 0 ) {
                if( !safeDivide ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." );
                    goto Err2;
                }
                zeros++;
                xz = ( v1 * x2 - v2 * x1 ) / ( v1 - v2 );
                if( ptwXY_setValueAtX( smr, div, xz, nan ) != nfu_Okay ) goto Err;
            }
            x2 = x1;
        }
        if( ptwXY_simpleCoalescePoints( smr, div ) != nfu_Okay ) goto Err;
        length = div->length;
        x2 = div->points[div->length-1].x;
        y2 = div->points[div->length-1].y;
        isNAN2 = nfu_isNAN( y2 );
        for( i = div->length - 2; i >= 0; i-- ) {             /* Make interpolation fit accuracy. Work backwards so new points will not mess up loop. */
            x1 = div->points[i].x;
            y1 = div->points[i].y;
            isNAN1 = nfu_isNAN( y1 );
            if( !isNAN1 || !isNAN2 ) {
                if( ptwXY_div_s_ptwXY( smr, div, ptwXY1, ptwXY2, x1, y1, x2, y2, 0, isNAN1, isNAN2 ) != nfu_Okay ) goto Err;
            }
            x2 = x1;
            y2 = y1;
            isNAN2 = isNAN1;
        }
        ptwXY_update_biSectionMax( div, (double) length );
        if( zeros ) {
            if( ptwXY_simpleCoalescePoints( smr, div ) != nfu_Okay ) goto Err;
            for( i = 0; i < div->length; i++ ) if( !nfu_isNAN( div->points[i].y ) ) break;
            if( nfu_isNAN( div->points[0].y ) ) {                     /* Special case for first point. */
                if( i == div->length ) {                              /* They are all nan's, what now? */
                    zeros = 0;
                    for( i = 0; i < div->length; i++ ) div->points[i].y = 0.; }
                else {
                     div->points[0].y = 2. * div->points[i].y;
                    zeros--;
                }
            }
            for( i = div->length - 1; i > 0; i-- ) if( !nfu_isNAN( div->points[i].y ) ) break;
            if( nfu_isNAN( div->points[div->length - 1].y ) ) {         /* Special case for last point. */
                div->points[div->length - 1].y = 2. * div->points[i].y;
                zeros--;
            }
            if( zeros ) {
                for( i = 0; i < div->length; i++ ) if( nfu_isNAN( div->points[i].y ) ) break;
                for( k = i + 1, j = i; k < div->length; k++ ) {
                if( nfu_isNAN( div->points[k].y ) ) continue;
                    div->points[j] = div->points[k];
                    j++;
                }
                div->length = j;
            }
        }
    }
    return( div );

Err:
    if( status == nfu_XOutsideDomain ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_XOutsideDomain, "x-value outside domain." ); }
    else {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    }
Err2:
    if( div ) ptwXY_free( div );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_div_s_ptwXY( statusMessageReporting *smr, ptwXYPoints *div, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double x1, double y1, double x2, double y2, int level, int isNAN1, int isNAN2 ) {

    nfu_status status;
    double u1, u2, v1, v2, v, x, y, yp, dx, a1, a2;

    if( ( x2 - x1 ) < ClosestAllowXFactor * DBL_EPSILON * ( fabs( x1 ) + fabs( x2 ) ) ) return( nfu_Okay );
    if( level >= div->biSectionMax ) return( nfu_Okay );
    level++;
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x1, &u1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x2, &u2 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x1, &v1 ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x2, &v2 ) ) != nfu_Okay ) return( status );
    if( isNAN1 ) {
        x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x, &u1 ) ) != nfu_Okay ) return( status );
        if( ( status = ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x, &v1 ) ) != nfu_Okay ) return( status );
        y = u1 / v1; }
    else if( isNAN2 ) {
        x = 0.5 * ( x1 + x2 );
        if( ( status = ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY1, x, &u2 ) ) != nfu_Okay ) return( status );
        if( ( status = ptwXY_getValueAtX_signal_XOutsideDomainError( smr, __LINE__, __func__, ptwXY2, x, &v2 ) ) != nfu_Okay ) return( status );
        y = u2 / v2; }
    else {
        if( ( u1 == u2 ) || ( v1 == v2 ) ) return( nfu_Okay );
        if( ( y1 == 0. ) || ( y2 == 0. ) ) {                    /* Handle special case of 0 where accuracy can never be met. */
            x = 0.5 * ( x1 + x2 ); }
        else {
            if( ( u1 * u2 < 0. ) ) return( nfu_Okay );  /* Assume rounding error and no point needed. */
            a1 = sqrt( fabs( u1 ) );
            a2 = sqrt( fabs( u2 ) );
            x = ( a2 * x1 + a1 * x2 ) / ( a2 + a1 );
        }
        dx = x2 - x1;
        v = v1 * ( x2 - x ) + v2 * ( x - x1 );
        if( ( v1 == 0. ) || ( v2 == 0. ) || ( v == 0. ) ) return( nfu_Okay );     /* Probably not correct, but I had to do something. */
        yp = ( u1 / v1 * ( x2 - x ) + u2 / v2 * ( x - x1 ) ) / dx;
        y = ( u1 * ( x2 - x ) + u2 * ( x - x1 ) ) / v;
        if( fabs( y - yp ) < fabs( y * div->accuracy ) ) return( nfu_Okay );
    }
    if( ( status = ptwXY_setValueAtX( smr, div, x, y ) ) != nfu_Okay ) return( status );
    if( ( status = ptwXY_div_s_ptwXY( smr, div, ptwXY1, ptwXY2, x, y, x2, y2, level, 0, isNAN2 ) ) != nfu_Okay ) return( status );
    status = ptwXY_div_s_ptwXY( smr, div, ptwXY1, ptwXY2, x1, y1, x, y, level, isNAN1, 0 );
    return( status );
}
/*
************************************************************
*/
static ptwXYPoints *ptwXY_div_ptwXY_forFlats( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int safeDivide ) {

    int64_t i;
    ptwXYPoints *div = NULL;
    ptwXYPoint *p;
    double y;

    if( ptwXY1->interpolation != ptwXY_interpolationFlat ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Source1 interpolation not 'flat' but '%s'.",
                ptwXY1->interpolationString );
        return( NULL );
    }
    if( ptwXY2->interpolation != ptwXY_interpolationFlat ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_invalidInterpolation, "Source2 interpolation not 'flat' but '%s'.",
                ptwXY2->interpolationString );
        return( NULL );
    }

    if( ( div = ptwXY_union( smr, ptwXY1, ptwXY2, ptwXY_union_fill | ptwXY_union_mergeClosePoints ) ) != NULL ) {
        for( i = 0, p = div->points; i < div->length; i++, p++ ) {
            if( ptwXY_getValueAtX_ignore_XOutsideDomainError( smr, ptwXY2, p->x, &y ) != nfu_Okay ) goto Err;
            if( y == 0. ) {
                if( ( safeDivide ) && ( p->y == 0 ) ) {
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_divByZero, "Divide by 0." );
                    goto Err;
                } }
            else {
                p->y /= y;
            }
        }
    }
    return( div );

Err:
    if( div ) ptwXY_free( div );
    return( NULL );
}
/*
************************************************************
*/
static nfu_status ptwXY_getValueAtX_ignore_XOutsideDomainError( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double x, double *y ) {

    nfu_status status = ptwXY_getValueAtX( smr, ptwXY1, x, y );

    if( status == nfu_XOutsideDomain ) status = nfu_Okay;
    return( status );
}
/*
************************************************************
*/
static nfu_status ptwXY_getValueAtX_signal_XOutsideDomainError( statusMessageReporting *smr, int line, 
        char const *function, ptwXYPoints *ptwXY1, double x, double *y ) {

    nfu_status status = ptwXY_getValueAtX( smr, ptwXY1, x, y );

    if( status == nfu_XOutsideDomain ) {
        double domainMin, domainMax;

        ptwXY_domainMin( NULL, ptwXY1, &domainMin );
        ptwXY_domainMax( NULL, ptwXY1, &domainMax );

        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XOutsideDomain, "x-value = %.17e outside domain (%.17e, %.17e)",
                x, domainMin, domainMin );
        smr_setReportError( smr, NULL, __FILE__, line, function, nfu_SMR_libraryID, nfu_Error, "Via." );
    }
    return( status );
}
