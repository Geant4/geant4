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

static int xDataXML_V_W_XYs_LegendreSeries_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, 
        xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, xDataTOM_axes *axes );
/*
************************************************************
*/
int xDataXML_V_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int index, length;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataXML_element *XMLChild;
    char const *wLabel;
    xDataTOM_V_W_XYs_LegendreSeries *V_W_XYs_LegendreSeries;

    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_V_W_XYs_LegendreSeries_ID, sizeof( xDataTOM_V_W_XYs_LegendreSeries ) ) ) == NULL ) 
        return( 1 );
    V_W_XYs_LegendreSeries = (xDataTOM_V_W_XYs_LegendreSeries *) xDI->data;
    if( ( wLabel = xDataTOM_axes_getLabel( smr, &(xDI->axes), 0 ) ) == NULL ) goto err;
    length = (int) xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_V_W_XYs_LegendreSeries_initialize( smr, V_W_XYs_LegendreSeries, length, &(xDI->axes) ) != 0 ) goto err;

    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( wLabel, XMLChild->name ) == 0 ) {
            if( xDataXML_V_W_XYs_LegendreSeries_W_XYs_LegendreSeriesToTOM( smr, XMLChild, &(V_W_XYs_LegendreSeries->W_XYs_LegendreSeries[index]), 
                &(xDI->axes) ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                    "invalid element '%s' in xData = 'V_W_XYs_LegendreSeries'", XMLChild->name );
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
static int xDataXML_V_W_XYs_LegendreSeries_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, 
        xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, xDataTOM_axes *axes ) {

    int index, length;
    double value;
    char const *wLabel;
    xDataXML_element *XMLChild;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "value", &value, 1 ) != 0 ) return( 1 );
    if( ( wLabel = xDataTOM_axes_getLabel( smr, axes, 1 ) ) == NULL ) goto err;
    length = xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_W_XYs_LegendreSeries_initialize( smr, W_XYs_LegendreSeries, index, length, value, xDataTOM_subAxesType_proxy, axes, NULL ) != 0 ) goto err;

    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( wLabel, XMLChild->name ) == 0 ) {
            if( xDataXML_W_XYs_LegendreSeries_LegendreSeriesToTOM( smr, XMLChild, &(W_XYs_LegendreSeries->LegendreSeries[index]) ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in xData = 'W_XYs_LegendreSeries'", XMLChild->name );
            goto err;
        }
    }

    return( 0 );

err:
    return( 1 );
}

#if defined __cplusplus
}
#endif
