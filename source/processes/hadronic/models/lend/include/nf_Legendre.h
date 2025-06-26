/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef nf_Legendre_h_included
#define nf_Legendre_h_included

#include <nf_utilities.h>
#include <ptwXY.h>

#if defined __cplusplus
    extern "C" {
#endif

#define nf_Legendre_minMaxOrder 4
#define nf_Legendre_maxMaxOrder 128
#define nf_Legendre_sizeIncrement 8

typedef struct nf_Legendre_s nf_Legendre;

struct nf_Legendre_s {
    nfu_status status;
    int maxOrder;
    int allocated;          /* Will never be less than nf_Legendre_minMaxOrder. */
    double *Cls;
};

typedef nfu_status (*nf_Legendre_GaussianQuadrature_callback)( double x, double *y, void *argList );

/*
* Methods in nf_Legendre.c
*/
nf_Legendre *nf_Legendre_new( statusMessageReporting *smr, int initialSize, int maxOrder, double *Cls );
nfu_status nf_Legendre_initialize( statusMessageReporting *smr, nf_Legendre *nfL, int initialSize, int maxOrder );
nfu_status nf_Legendre_release( statusMessageReporting *smr, nf_Legendre *nfL );
nf_Legendre *nf_Legendre_free( nf_Legendre *nfL );
nf_Legendre *nf_Legendre_clone( statusMessageReporting *smr, nf_Legendre *nfL );
nfu_status nf_Legendre_reallocateCls( statusMessageReporting *smr, nf_Legendre *Legendre, int size, int forceSmallerResize );
nfu_status nf_Legendre_maxOrder( statusMessageReporting *smr, nf_Legendre *Legendre, int *maxOrder );
nfu_status nf_Legendre_allocated( statusMessageReporting *smr, nf_Legendre *Legendre, int *allocated );
nfu_status nf_Legendre_getCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double *Cl );
nfu_status nf_Legendre_setCl( statusMessageReporting *smr, nf_Legendre *Legendre, int l, double Cl );
nfu_status nf_Legendre_normalize( statusMessageReporting *smr, nf_Legendre *Legendre );
nfu_status nf_Legendre_evauluateAtMu( statusMessageReporting *smr, nf_Legendre *nfL, double mu, double *P );
double nf_Legendre_PofL_atMu( int l, double mu );
ptwXYPoints *nf_Legendre_to_ptwXY( statusMessageReporting *smr, nf_Legendre *nfL, double accuracy, int biSectionMax, 
        int checkForRoots );
nf_Legendre *nf_Legendre_from_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY, int maxOrder );

/*
* Methods in nf_Legendre_GaussianQuadrature.c
*/
nfu_status nf_Legendre_GaussianQuadrature( int degree, double x1, double x2, nf_Legendre_GaussianQuadrature_callback func, void *argList, double *integral );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_Legendre_h_included. */
