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

char const *xDataTOM_regionsW_XYs_LegendreSeries_ID = "regionsW_XYs_LegendreSeries";

/*
************************************************************
*/
int xDataTOM_regionsW_XYs_LegendreSeries_initialize( statusMessageReporting *smr, xDataTOM_regionsW_XYs_LegendreSeries *regionsW_XYs_LegendreSeries, 
        int length, xDataTOM_axes *axes ) {

    regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries = NULL;
    regionsW_XYs_LegendreSeries->length = length;
    regionsW_XYs_LegendreSeries->axes = axes;
    if( ( regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries = (xDataTOM_W_XYs_LegendreSeries *) smr_malloc2( smr, length * sizeof( xDataTOM_W_XYs_LegendreSeries ), 1, 
        "regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries" ) ) == NULL ) return( 1 );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_regionsW_XYs_LegendreSeries_free( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_regionsW_XYs_LegendreSeries_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_regionsW_XYs_LegendreSeries_release( (xDataTOM_regionsW_XYs_LegendreSeries *) xDI->data );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_regionsW_XYs_LegendreSeries_release( xDataTOM_regionsW_XYs_LegendreSeries *regionsW_XYs_LegendreSeries ) {

    int i;

    for( i = 0; i < regionsW_XYs_LegendreSeries->length; i++ ) xDataTOM_W_XYs_LegendreSeries_release( &(regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries[i]) );
    smr_freeMemory( (void **) &(regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries) );
    return( 0 );
}

#if defined __cplusplus
}
#endif
