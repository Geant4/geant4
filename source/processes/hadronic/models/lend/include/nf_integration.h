/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef nf_integration_h_included
#define nf_integration_h_included

#include <nf_utilities.h>
#include <nf_Legendre.h>

#if defined __cplusplus
    extern "C" {
#endif

#define nf_GnG_adaptiveQuadrature_MaxMaxDepth 20

typedef nfu_status (*nf_GnG_adaptiveQuadrature_callback)( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, double x1, 
    double x2, double *integral );

nfu_status nf_GnG_adaptiveQuadrature( nf_GnG_adaptiveQuadrature_callback quadratureFunction, nf_Legendre_GaussianQuadrature_callback integrandFunction, 
    void *argList, double x1, double x2, int maxDepth, double tolerance, double *integral, long *evaluations );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_integration_h_included. */

