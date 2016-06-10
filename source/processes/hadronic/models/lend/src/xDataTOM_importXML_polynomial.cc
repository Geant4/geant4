/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "xDataTOM_importXML_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
int xDataXML_polynomialToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int length, dataProcessed = 0;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataTOM_polynomial *polynomial = NULL;
    xDataXML_element *XMLChild;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_polynomial_ID, sizeof( xDataTOM_polynomial ) ) ) == NULL ) return( 1 );
    if( xDataTOM_polynomial_initialize( smr, (xDataTOM_polynomial  *) xDI->data, length, &(xDI->axes) ) != 0 ) goto err;
    polynomial = (xDataTOM_polynomial  *) xDI->data;

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( "data", XMLChild->name ) == 0 ) {
            if( dataProcessed ) {
                smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "multiple 'data' elements found" );
                goto err;
            }
            dataProcessed = 1;
            if( xDataXML_stringToDoubles( smr, XE, XMLChild->text.text, length, (double *) polynomial->coefficients ) != 0 ) goto err;
        }
    }
    if( dataProcessed == 0 ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "data element missing" );
        goto err;
    }
    return( 0 );

err:
    if( polynomial != NULL ) xDataTOM_polynomial_release( polynomial );
    smr_freeMemory( (void **) &(xDI->data) );
    return( 1 );
}

#if defined __cplusplus
}
#endif
