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

char const *xDataTOM_XYs_ID = "XYs";

/*
************************************************************
*/
int xDataTOM_XYs_free( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_XYs_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_XYs_release( (xDataTOM_XYs *) xDI->data );
    smr_freeMemory( (void **) &(xDI->data) );
    
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_XYs_release( xDataTOM_XYs *XYs ) {

    xDataTOM_subAxes_release( &(XYs->subAxes) );
    XYs->length = 0;
    smr_freeMemory( (void **) &(XYs->data) );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_XYs_getData( xDataTOM_XYs *XYs, double **data ) {

    *data = XYs->data;
    return( XYs->length );
}
/*
************************************************************
*/
int xDataTOM_XYs_getDataFromXDataInfo( xDataTOM_xDataInfo *xDI, double **data ) {

    return( xDataTOM_XYs_getData( (xDataTOM_XYs *) xDI->data, data ) );
}

#if defined __cplusplus
}
#endif
