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

char const *xDataTOM_W_XYs_ID = "W_XYs";

/*
************************************************************
*/
xDataTOM_W_XYs *xDataTOM_W_XYs_new( statusMessageReporting *smr, int index, int length, double value, xDataTOM_axes *axes, int axesOffset ) {

    xDataTOM_W_XYs *W_XYs;

    if( ( W_XYs = (xDataTOM_W_XYs *) smr_malloc2( smr, sizeof( xDataTOM_W_XYs ), 0, "W_XYs" ) ) == NULL ) return( NULL );
    if( xDataTOM_W_XYs_initialize( smr, W_XYs, index, length, value, axes, axesOffset ) != 0 ) smr_freeMemory( (void **) &W_XYs );
    return( W_XYs );
}
/*
************************************************************
*/
int xDataTOM_W_XYs_initialize( statusMessageReporting *smr, xDataTOM_W_XYs *W_XYs, int index, int length, double value, xDataTOM_axes *axes,
    int axesOffset ) {

    W_XYs->XYs = NULL;
    W_XYs->index = index;
    W_XYs->length = length;
    W_XYs->value = value;
    if( ( W_XYs->XYs = (xDataTOM_XYs *) smr_malloc2( smr, length * sizeof( xDataTOM_XYs ), 1, "W_XYs->XYs" ) ) == NULL ) return( 1 );
    if( xDataTOM_subAxes_initialize( smr, &(W_XYs->subAxes), xDataTOM_subAxesType_proxy, axesOffset, axes, NULL ) != 0 ) {
        smr_freeMemory( (void **) &(W_XYs->XYs) );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
xDataTOM_W_XYs *xDataTOM_W_XYs_free( xDataTOM_W_XYs *W_XYs ) {

    if( W_XYs == NULL ) return( NULL );
    xDataTOM_W_XYs_release( W_XYs );
    smr_freeMemory( (void **) &W_XYs );
    return( (xDataTOM_W_XYs *) NULL );
}
/*
************************************************************
*/
int xDataTOM_W_XYs_freeFrom_xDataInfo( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_W_XYs_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_W_XYs_free( (xDataTOM_W_XYs *) xDI->data );
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_W_XYs_release( xDataTOM_W_XYs *W_XYs ) {

    int i;

    xDataTOM_subAxes_release( &(W_XYs->subAxes) );
    for( i = 0; i < W_XYs->length; i++ ) xDataTOM_XYs_release( &(W_XYs->XYs[i]) );
    W_XYs->length = 0;
    smr_freeMemory( (void **) &(W_XYs->XYs) );

    return( 0 );
}

#if defined __cplusplus
}
#endif
