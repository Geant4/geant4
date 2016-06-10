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

char const *xDataTOM_polynomial_ID = "polynomial";

/*
************************************************************
*/
int xDataTOM_polynomial_initialize( statusMessageReporting *smr, xDataTOM_polynomial *polynomial, int length, xDataTOM_axes *axes ) {

    polynomial->length = length;
    if( ( polynomial->coefficients = (double *) smr_malloc2( smr, length * sizeof( double ), 1, "polynomial->coefficients" ) ) == NULL ) return( 1 );
    if( xDataTOM_subAxes_initialize( smr, &(polynomial->subAxes), xDataTOM_subAxesType_proxy, 0, axes, NULL ) != 0 ) {
        smr_freeMemory( (void **) &(polynomial->coefficients) );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_polynomial_free( xDataTOM_xDataInfo *xDI ) {

    if( xDI == NULL ) return( 0 );
    if( strcmp( xDataTOM_polynomial_ID, xDI->ID ) != 0 ) return( 1 );
    xDataTOM_polynomial_release( (xDataTOM_polynomial *) xDI->data );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_polynomial_release( xDataTOM_polynomial *polynomial ) {

    xDataTOM_subAxes_release( &(polynomial->subAxes) );
    polynomial->length = 0;
    smr_freeMemory( (void **) &(polynomial->coefficients) );

    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_polynomial_getData( xDataTOM_polynomial *polynomial, double **data ) {

    *data = polynomial->coefficients;
    return( polynomial->length );
}
/*
************************************************************
*/
int xDataTOM_polynomial_getDataFromXDataInfo( xDataTOM_xDataInfo *xDI, double **data ) {

    return( xDataTOM_polynomial_getData( (xDataTOM_polynomial *) xDI->data, data ) );
}

#if defined __cplusplus
}
#endif
