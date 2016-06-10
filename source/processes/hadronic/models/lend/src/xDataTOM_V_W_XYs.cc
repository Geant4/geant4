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

char const *xDataTOM_V_W_XYs_ID = "V_W_XYs";

/*
************************************************************
*/
int xDataTOM_V_W_XYs_initialize( statusMessageReporting *smr, xDataTOM_V_W_XYs *V_W_XYs, int length, xDataTOM_axes *axes ) {

    V_W_XYs->W_XYs = NULL;
    V_W_XYs->length = length;
    if( ( V_W_XYs->W_XYs = (xDataTOM_W_XYs *) smr_malloc2( smr, length * sizeof( xDataTOM_W_XYs ), 1, "V_W_XYs->W_XYs" ) ) == NULL ) return( 1 );
    if( xDataTOM_subAxes_initialize( smr, &(V_W_XYs->subAxes), xDataTOM_subAxesType_proxy, 0, axes, NULL ) != 0 ) {
        smr_freeMemory( (void **) &(V_W_XYs->W_XYs) );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_V_W_XYs_free( xDataTOM_xDataInfo *xDI ) {

    int i;
    xDataTOM_V_W_XYs *V_W_XYs;

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_V_W_XYs_ID, xDI->ID ) != 0 ) return( 1 );
    if( ( V_W_XYs = (xDataTOM_V_W_XYs *) xDI->data ) != NULL ) {
        for( i = 0; i < V_W_XYs->length; i++ ) xDataTOM_W_XYs_release( &(V_W_XYs->W_XYs[i]) );
        smr_freeMemory( (void **) &(V_W_XYs->W_XYs) );
        smr_freeMemory( (void **) &(xDI->data) );
    }
    return( 0 );
}

#if defined __cplusplus
}
#endif
