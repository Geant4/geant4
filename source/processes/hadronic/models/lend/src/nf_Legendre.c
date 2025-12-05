/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "nf_Legendre.h"

struct nf_Legendre_from_ptwXY_callback_s {
    int l;
    double mu1, mu2, f1, f2;
};

static nfu_status nf_Legendre_to_ptwXY2( statusMessageReporting *smr, double mu, double *P, void *argList );
static nfu_status nf_Legendre_from_ptwXY_callback( double mu, double *f, void *argList );
/*
************************************************************
*/
nf_Legendre *nf_Legendre_new( statusMessageReporting *smr, int initialSize, int maxOrder, double *Cls ) {

    int l;
    nf_Legendre *Legendre = (nf_Legendre *) smr_malloc2( smr, sizeof( nf_Legendre ), 1, "Legendre" );

    if( Legendre == NULL ) return( NULL );
    if( nf_Legendre_initialize( smr, Legendre, initialSize, maxOrder ) != nfu_Okay ) {
        nfu_free( Legendre );
        return( NULL );
    }
    for( l = 0; l <= Legendre->maxOrder; l++ ) Legendre->Cls[l] = Cls[l];
    return( Legendre );
}
/*
************************************************************
*/
nfu_status nf_Legendre_initialize( statusMessageReporting *smr, nf_Legendre *Legendre, int initialSize, int maxOrder ) {

    nfu_status status;

    memset( Legendre, 0, sizeof( nf_Legendre ) );
    Legendre->status = nfu_Okay;
    if( maxOrder < 0 ) maxOrder = -1;
    if( maxOrder > nf_Legendre_maxMaxOrder ) maxOrder = nf_Legendre_maxMaxOrder;
    Legendre->maxOrder = maxOrder;
    if( initialSize < ( maxOrder + 1 ) ) initialSize = maxOrder + 1;
    if( ( status = nf_Legendre_reallocateCls( smr, Legendre, initialSize, 0 ) ) != nfu_Okay )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( status );
}
/*
************************************************************
*/
nfu_status nf_Legendre_release( statusMessageReporting *smr, nf_Legendre *Legendre ) {

    if( Legendre->allocated > 0 ) nfu_free( Legendre->Cls );
    memset( Legendre, 0, sizeof( nf_Legendre ) );
    return( nfu_Okay );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_free( nf_Legendre *Legendre ) {

    nf_Legendre_release( NULL, Legendre );
    nfu_free( Legendre );
    return( NULL );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_clone( statusMessageReporting *smr, nf_Legendre *nfL ) {

    nf_Legendre *Legendre;

    if( nfL->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( ( Legendre = nf_Legendre_new( smr, 0, nfL->maxOrder, nfL->Cls ) ) == NULL )
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( Legendre );
}
/*
************************************************************
*/
nfu_status nf_Legendre_reallocateCls( statusMessageReporting *smr, nf_Legendre *Legendre, int size, int forceSmallerResize ) {

    int i1;

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( size < nf_Legendre_minMaxOrder ) size = nf_Legendre_minMaxOrder;
    if( size > ( nf_Legendre_maxMaxOrder + 1 ) ) size = nf_Legendre_maxMaxOrder + 1;
    if( size != Legendre->allocated ) {
        if( size > Legendre->allocated ) {
            Legendre->Cls = (double *) smr_realloc2( smr, Legendre->Cls, size * sizeof( double ), "Cls" ); }
        else {
            if( size < ( Legendre->maxOrder + 1 ) ) size = Legendre->maxOrder + 1;
            if( ( Legendre->allocated > 2 * size ) || forceSmallerResize ) {
                    Legendre->Cls = (double *) nfu_realloc( size * sizeof( double ), Legendre->Cls ); } 
            else {
                size = Legendre->allocated;
            }
        }
        if( Legendre->Cls == NULL ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            size = 0;
            Legendre->status = nfu_mallocError;
        }
        Legendre->allocated = size;
    }

    if( Legendre->Cls != NULL ) {
        for( i1 = Legendre->maxOrder + 1; i1 < size; ++i1 ) Legendre->Cls[i1] = 0;
    }

    return( Legendre->status );
}
/*
************************************************************
*/
nfu_status nf_Legendre_maxOrder( statusMessageReporting *smr, nf_Legendre *Legendre, int *maxOrder ) {

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    *maxOrder = Legendre->maxOrder;
    return( Legendre->status );
}
/*
************************************************************
*/
nfu_status nf_Legendre_allocated( statusMessageReporting *smr, nf_Legendre *Legendre, int *allocated ) {

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    *allocated = Legendre->allocated;
    return( Legendre->status );
}
/*
************************************************************
*/
nfu_status nf_Legendre_getCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double *Cl ) {

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    *Cl = 0;
    if( ( l < 0 ) || ( l > Legendre->maxOrder ) ) {
        return( nfu_badIndex ); }
    else {
        *Cl = Legendre->Cls[l];
    }
    return( Legendre->status );
}
/*
************************************************************
*/
nfu_status nf_Legendre_setCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double Cl ) {

    nfu_status status;

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( l > nf_Legendre_maxMaxOrder ) return( nfu_Okay );

    if( l < 0 ) {
        Legendre->status = nfu_badIndex;
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badIndex, "Negative l-order %d is not allowed.", l );
        return( nfu_badIndex );
    }
    if( Legendre->allocated <= l ) {
        if( ( status = nf_Legendre_reallocateCls( smr, Legendre, l + nf_Legendre_sizeIncrement, 0 ) ) != nfu_Okay ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( status );
        }
    }

    Legendre->Cls[l] = Cl;
    if( l > Legendre->maxOrder ) Legendre->maxOrder = l;
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status nf_Legendre_normalize( statusMessageReporting *smr, nf_Legendre *Legendre ) {

    int l;
    double norm;

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    if( Legendre->maxOrder >= 0 ) {
        if( ( norm = Legendre->Cls[0] ) == 0 ) {
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
            return( Legendre->status = nfu_divByZero );
        }
        for( l = 0; l <= Legendre->maxOrder; l++ ) Legendre->Cls[l] /= norm;
    }
    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status nf_Legendre_evauluateAtMu( statusMessageReporting *smr, nf_Legendre *Legendre, double mu, double *P ) {

    int l;

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( nfu_badSelf );
    }

    *P = 0;
    if( ( mu >= -1. ) && ( mu <= 1. ) ) {
        for( l = 0; l <= Legendre->maxOrder; l++ ) *P += ( l + 0.5 ) * Legendre->Cls[l] * nf_Legendre_PofL_atMu( l, mu ); }
    else {
        return( nfu_XOutsideDomain );
    }
    return( Legendre->status );
}
/*
************************************************************
*/
double nf_Legendre_PofL_atMu( int l, double mu ) {

    int l_, twoL_plus1;
    double Pl_minus1, Pl, Pl_plus1;

    if( l == 0 ) {
        return( 1. ); }
    else if( l == 1 ) {
        return( mu ); }
/*
    else if( l <= 9 ) {
        double mu2 = mu * mu;
        if ( l == 2 ) {
            return(                1.5 * mu2 - 0.5 ); }
        else if( l == 3 ) {
            return(                2.5 * mu2 -          1.5 ) * mu; }
        else if( l == 4 ) {
            return(              4.375 * mu2 -         3.75 ) * mu2 +       0.375; }
        else if( l == 5 ) {
            return( (            7.875 * mu2 -         8.75 ) * mu2 +       1.875 ) * mu; }
        else if( l == 6 ) {
            return( (          14.4375 * mu2 -      19.6875 ) * mu2 +      6.5625 ) * mu2 -      0.3125; }
        else if( l == 7 ) {
            return( ( (        26.8125 * mu2 -      43.3125 ) * mu2 +     19.6875 ) * mu2 -      2.1875 ) * mu; }
        else if( l == 8 ) {
            return( ( (     50.2734375 * mu2 -     93.84375 ) * mu2 +   54.140625 ) * mu2 -     9.84375 ) * mu2 + 0.2734375; }
        else {
            return( ( ( (   94.9609375 * mu2 -    201.09375 ) * mu2 +  140.765625 ) * mu2 -    36.09375 ) * mu2 + 2.4609375 ) * mu;
        }
    }
*/

    Pl = 0.;
    Pl_plus1 = 1.;
    for( l_ = 0, twoL_plus1 = 1; l_ < l; l_++, twoL_plus1 += 2 ) {
        Pl_minus1 = Pl;
        Pl = Pl_plus1;
        Pl_plus1 = ( twoL_plus1 * mu * Pl - l_ * Pl_minus1 ) / ( l_ + 1 );
    }
    return( Pl_plus1 );
}
/*
************************************************************
*/
ptwXYPoints *nf_Legendre_to_ptwXY( statusMessageReporting *smr, nf_Legendre *Legendre, double accuracy, 
        int biSectionMax, int checkForRoots ) {

    int i, n = 1;
    double dx, xs[1000];
    void *argList = (void *) Legendre;
    ptwXYPoints *ptwXY;

    if( Legendre->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    xs[0] = -1;
    if( Legendre->maxOrder > 1 ) {
        n = Legendre->maxOrder - 1;
        if( n > 249 ) n = 249;
        n = 4 * n + 1;
        dx = 2. / n;
        for( i = 1; i < n; i++ ) xs[i] = xs[i-1] + dx;
    }
    xs[n] = 1.;
    ptwXY = ptwXY_createFromFunction( smr, n + 1, xs, nf_Legendre_to_ptwXY2, argList, accuracy, checkForRoots, biSectionMax );
    if( ptwXY == NULL ) smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    return( ptwXY );
}
/*
************************************************************
*/
static nfu_status nf_Legendre_to_ptwXY2( statusMessageReporting *smr, double mu, double *P, void *argList ) {

    return( nf_Legendre_evauluateAtMu( smr, (nf_Legendre *) argList, mu, P ) );
}
/*
************************************************************
*/
nf_Legendre *nf_Legendre_from_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY, int maxOrder ) {

    int l, i, n = (int) ptwXY_length( NULL, ptwXY );
    nf_Legendre *Legendre;
    double mu1, mu2, f1, f2, Cl, Cls[1] = { 0 }, integral;
    struct nf_Legendre_from_ptwXY_callback_s argList;

    if( ptwXY->status != nfu_Okay ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_badSelf, "Invalid source." );
        return( NULL );
    }

    if( n == 0 ) {
        if( ( Legendre = nf_Legendre_new( smr, maxOrder + 1, -1, Cls ) ) == NULL )
            smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( Legendre );
    }

    if( n == 1 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_tooFewPoints, "ptwXY only has 1 point." );
        return( NULL );
    }

    ptwXY_getXYPairAtIndex( smr, ptwXY, 0,     &mu1, &f1 );
    ptwXY_getXYPairAtIndex( smr, ptwXY, n - 1, &mu2, &f2 );
    if( ( mu1 < -1 ) || ( mu2 > 1 ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_XOutsideDomain, "bad domain for ptwXY: %.17e %17.e", mu1, mu2 );
        return( NULL );
    }

    if( ( Legendre = nf_Legendre_new( smr, maxOrder + 1, -1, Cls ) ) == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }

    if( maxOrder > nf_Legendre_maxMaxOrder ) maxOrder = nf_Legendre_maxMaxOrder;
    for( l = 0; l <= maxOrder; l++ ) {
        ptwXY_getXYPairAtIndex( smr, ptwXY, 0, &mu1, &f1 );
        argList.l = l;
        for( i = 1, Cl = 0; i < n; i++ ) {
            ptwXY_getXYPairAtIndex( smr, ptwXY, i, &mu2, &f2 );
            argList.mu1 = mu1;
            argList.f1 = f1;
            argList.mu2 = mu2;
            argList.f2 = f2;
            if( nf_Legendre_GaussianQuadrature( l + 1, mu1, mu2, nf_Legendre_from_ptwXY_callback, (void *) &argList, 
                    &integral ) != nfu_Okay ) goto err;
            Cl += integral;
            mu1 = mu2;
            f1 = f2;
        }
        if( nf_Legendre_setCl( smr, Legendre, l, Cl ) != nfu_Okay ) goto err;
    }
    return( Legendre );

err:
    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
    nf_Legendre_free( Legendre );
    return( NULL );    
}
/*
************************************************************
*/
static nfu_status nf_Legendre_from_ptwXY_callback( double mu, double *f, void *argList ) {

    struct nf_Legendre_from_ptwXY_callback_s *args = (struct nf_Legendre_from_ptwXY_callback_s *) argList;

    *f = ( args->f1 * ( args->mu2 - mu ) + args->f2 * ( mu - args->mu1 ) ) / ( args->mu2 - args->mu1 );
    *f *= nf_Legendre_PofL_atMu( args->l, mu );
    return( nfu_Okay );
}
