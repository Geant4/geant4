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

static int xDataXML_LegendreSeriesDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_LegendreSeries *LegendreSeries, 
        int index, int length, double value );
/*
************************************************************
*/
int xDataXML_W_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int index, length;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataXML_element *XMLChild;
    char const *wLabel;
    xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries = NULL;

    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_W_XYs_LegendreSeries_ID, sizeof( xDataTOM_W_XYs_LegendreSeries ) ) ) == NULL ) 
        return( 1 );
    W_XYs_LegendreSeries = (xDataTOM_W_XYs_LegendreSeries *) xDI->data;
    if( ( wLabel = xDataTOM_axes_getLabel( smr, &(xDI->axes), 0 ) ) == NULL ) goto err;
    length = xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_W_XYs_LegendreSeries_initialize( smr, W_XYs_LegendreSeries, 0, length, 0., xDataTOM_subAxesType_proxy, &(xDI->axes), NULL ) != 0 ) goto err;

    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( wLabel, XMLChild->name ) == 0 ) {
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
/*
************************************************************
*/
int xDataXML_W_XYs_LegendreSeries_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_LegendreSeries *LegendreSeries ) {

    int index, length;
    double value;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "value", &value, 1 ) != 0 ) return( 1 );
    return( xDataXML_LegendreSeriesDataToTOM( smr, XE, LegendreSeries, index, length, value ) );
}
/*
************************************************************
*/
static int xDataXML_LegendreSeriesDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_LegendreSeries *LegendreSeries, 
        int index, int length, double value ) {

    if( xDataTOM_LegendreSeries_initialize( smr, LegendreSeries, index, length, value ) != 0 ) return( 1 );
    if( xDataXML_stringToDoubles( smr, XE, XE->text.text, length, LegendreSeries->LegendreSeries ) == 0 ) return( 0 );
    xDataTOM_LegendreSeries_release( LegendreSeries );
    return( 1 );
}

#if defined __cplusplus
}
#endif
