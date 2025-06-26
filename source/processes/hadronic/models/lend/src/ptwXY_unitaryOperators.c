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

/*
************************************************************
*/
nfu_status ptwXY_abs( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength;
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ( nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ) ) < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = fabs( p->y );
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = fabs( o->point.y );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_neg( statusMessageReporting *smr, ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength;
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ( nonOverflowLength = ptwXY_getNonOverflowLength( smr, ptwXY ) ) < 0 ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( ptwXY->status );
    }

    if( ( ptwXY->interpolation != ptwXY_interpolationLinLin ) && ( ptwXY->interpolation != ptwXY_interpolationLinLog ) &&
            ( ptwXY_interpolationLinLog != ptwXY_interpolationFlat ) ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_unsupportedInterpolation, 
                "Negation of non-linear y-interpolation not allowed: interpolation = '%s'.",
                ptwXY->interpolationString ); 
        return( nfu_unsupportedInterpolation );
    }

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = -p->y;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = -o->point.y;
    return( ptwXY->status );
}
