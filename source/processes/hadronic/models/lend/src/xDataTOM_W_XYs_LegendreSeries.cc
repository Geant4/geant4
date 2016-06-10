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

char const *xDataTOM_W_XYs_LegendreSeries_ID = "W_XYs_LegendreSeries";

/*
************************************************************
*/
int xDataTOM_W_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, int index, 
        int length, double value, enum xDataTOM_subAxesType subAxesType, xDataTOM_axes *axes, xDataTOM_interpolation *interpolation ) {

    W_XYs_LegendreSeries->LegendreSeries = NULL;
    W_XYs_LegendreSeries->index = index;
    W_XYs_LegendreSeries->length = length;
    W_XYs_LegendreSeries->value = value;
    if( xDataTOM_subAxes_initialize( smr, &(W_XYs_LegendreSeries->subAxes), subAxesType, 0, axes, interpolation ) != 0 ) return( 1 );
    if( ( W_XYs_LegendreSeries->LegendreSeries = (xDataTOM_LegendreSeries *) smr_malloc2( smr, length * sizeof( xDataTOM_LegendreSeries ), 1, "W_XYs_LegendreSeries->LegendreSeries" ) ) == NULL ) return( 1 );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_W_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_W_XYs_LegendreSeries_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_W_XYs_LegendreSeries_release( (xDataTOM_W_XYs_LegendreSeries *) xDI->data );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_W_XYs_LegendreSeries_release( xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries ) {

    int i;

    for( i = 0; i < W_XYs_LegendreSeries->length; i++ ) xDataTOM_LegendreSeries_release( &(W_XYs_LegendreSeries->LegendreSeries[i]) );
    smr_freeMemory( (void **) &(W_XYs_LegendreSeries->LegendreSeries) );
    return( 0 );
}
/*
************************************************************
*/
#if 0
xDataTOM_W_XYs *xDataTOM_W_XYs_LegendreSeries_toW_XYs( statusMessageReporting *smr, xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, int axesOffset ) {

    xDataTOM_W_XYs *W_XYs = NULL;

/*
    if( ( W_XYs = xDataTOM_W_XYs_new( smr, W_XYs_LegendreSeries->index, W_XYs_LegendreSeries->length, W_XYs_LegendreSeries->value, axes, axesOffset ) )
        == NULL ) return( NULL );
*/

    return( W_XYs );

err:
    xDataTOM_W_XYs_free( W_XYs );
    return( NULL );
}
#endif

#if defined __cplusplus
}
#endif
