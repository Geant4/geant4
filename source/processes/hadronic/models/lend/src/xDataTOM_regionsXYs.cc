/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>

#include "xDataTOM_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

char const *xDataTOM_regionsXYs_ID = "regionsXYs";

/*
************************************************************
*/
int xDataTOM_regionsXYs_free( xDataTOM_xDataInfo *xDI ) {

    int i;
    xDataTOM_regionsXYs *regionsXYs;

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_regionsXYs_ID, xDI->ID ) != 0 ) return( 1 );
    regionsXYs = (xDataTOM_regionsXYs *) xDI->data;
    for( i = 0; i < regionsXYs->length; i++ ) xDataTOM_XYs_release( &(regionsXYs->XYs[i]) );
    smr_freeMemory( (void **) &(regionsXYs->XYs) );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}

#if defined __cplusplus
}
#endif
