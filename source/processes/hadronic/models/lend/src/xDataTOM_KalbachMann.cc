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

char const *xDataTOM_KalbachMann_ID = "KalbachMann";

/*
************************************************************
*/
int xDataTOM_KalbachMann_initialize( statusMessageReporting *smr, xDataTOM_KalbachMann *KalbachMann, int length, xDataTOM_axes *axes ) {

    KalbachMann->coefficients = NULL;
    KalbachMann->numberOfEnergies = 0;
    if( ( KalbachMann->coefficients = (xDataTOM_KalbachMannCoefficients *) smr_malloc2( smr, length * sizeof( xDataTOM_KalbachMannCoefficients ), 1, "KalbachMann->coefficients" ) ) == NULL ) return( 1 );
    if( xDataTOM_subAxes_initialize( smr, &(KalbachMann->subAxes), xDataTOM_subAxesType_proxy, 0, axes, NULL ) != 0 ) {
        smr_freeMemory( (void **) &(KalbachMann->coefficients) );
        return( 1 );
    }

    return( 0 );
}

/*
************************************************************
*/
int xDataTOM_KalbachMann_free( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_KalbachMann_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_KalbachMann_release( (xDataTOM_KalbachMann *) xDI->data );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_KalbachMann_release( xDataTOM_KalbachMann *KalbachMann ) {

    int i;

    xDataTOM_subAxes_release( &(KalbachMann->subAxes) );
    for( i = 0; i < KalbachMann->numberOfEnergies; i++ ) smr_freeMemory( (void **) &(KalbachMann->coefficients[i].coefficients) );
    KalbachMann->numberOfEnergies = 0;
    smr_freeMemory( (void **) &(KalbachMann->coefficients) );

    return( 0 );
}

#if defined __cplusplus
}
#endif
