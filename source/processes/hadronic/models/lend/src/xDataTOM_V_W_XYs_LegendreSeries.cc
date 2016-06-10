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

char const *xDataTOM_V_W_XYs_LegendreSeries_ID = "V_W_XYs_LegendreSeries";

/*
************************************************************
*/
int xDataTOM_V_W_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_V_W_XYs_LegendreSeries *V_W_XYs_LegendreSeries,
        int length, xDataTOM_axes *axes ) {

    V_W_XYs_LegendreSeries->W_XYs_LegendreSeries = NULL;
    V_W_XYs_LegendreSeries->length = length;
    if( xDataTOM_subAxes_initialize( smr, &(V_W_XYs_LegendreSeries->subAxes), xDataTOM_subAxesType_proxy, 0, axes, NULL ) != 0 ) return( 1 );
    if( ( V_W_XYs_LegendreSeries->W_XYs_LegendreSeries = (xDataTOM_W_XYs_LegendreSeries *) smr_malloc2( smr, length * sizeof( xDataTOM_W_XYs_LegendreSeries ), 1, 
        "V_W_XYs_LegendreSeries->W_XYs_LegendreSeries" ) ) == NULL ) return( 1 );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_V_W_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI ) {

    int i;
    xDataTOM_V_W_XYs_LegendreSeries *V_W_XYs_LegendreSeries;

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_V_W_XYs_LegendreSeries_ID, xDI->ID ) != 0 ) return( 1 );
    V_W_XYs_LegendreSeries = (xDataTOM_V_W_XYs_LegendreSeries *) xDI->data;
    for( i = 0; i < V_W_XYs_LegendreSeries->length; i++ ) xDataTOM_W_XYs_LegendreSeries_release( &(V_W_XYs_LegendreSeries->W_XYs_LegendreSeries[i]) );
    smr_freeMemory( (void **) &(V_W_XYs_LegendreSeries->W_XYs_LegendreSeries) );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}

#if defined __cplusplus
}
#endif
