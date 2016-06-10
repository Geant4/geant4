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

static int xDataXML_W_XYs_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes,
        int axesOffset );
/*
************************************************************
*/
int xDataXML_W_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataTOM_W_XYs *W_XYs;

/* Need to release W_XYs if an error occurs later. */
    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_W_XYs_ID, sizeof( xDataTOM_W_XYs ) ) ) == NULL ) return( 1 );
    W_XYs = (xDataTOM_W_XYs  *) xDI->data;

    if( xDataXML_W_XYsDataToTOM( smr, XE, W_XYs, 0, 0., &(xDI->axes), 0 ) != 0 ) goto err;
    return( 0 );

err:
    smr_freeMemory( (void **) &(xDI->data) );
    return( 1 );
}
/*
************************************************************
*/
int xDataXML_W_XYsDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_W_XYs *W_XYs, int index, double value, xDataTOM_axes *axes,
    int axesOffset ) {

    int length;
    char const *wLabel;
    xDataXML_element *XMLChild;

    if( ( wLabel = xDataTOM_axes_getLabel( smr, axes, axesOffset ) ) == NULL ) goto err;
    length = xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_W_XYs_initialize( smr, W_XYs, index, length, value, axes, axesOffset ) != 0 ) return( 1 );

    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( wLabel, XMLChild->name ) == 0 ) {
            if( xDataXML_W_XYs_XYsToTOM( smr, XMLChild, &(W_XYs->XYs[index]), axes, axesOffset + 1 ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in xData = 'W_XYs'", XMLChild->name );
            goto err;
        }
    }

    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int xDataXML_W_XYs_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes,
        int axesOffset ) {

    int index, length;
    double accuracy, value;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "accuracy", &accuracy, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "value", &value, 1 ) != 0 ) return( 1 );
    return( xDataXML_XYsDataToTOM( smr, XE, XYs, index, length, value, accuracy, xDataTOM_subAxesType_proxy, axesOffset, axes, NULL ) );
}

#if defined __cplusplus
}
#endif
