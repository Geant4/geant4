/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <cmath>
#include <float.h>

#include "ptwXY.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
nfu_status ptwXY_abs( ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = std::fabs( p->y );
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = std::fabs( o->point.y );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_neg( ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = -p->y;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = -o->point.y;
    return( ptwXY->status );
}

#if defined __cplusplus
}
#endif
