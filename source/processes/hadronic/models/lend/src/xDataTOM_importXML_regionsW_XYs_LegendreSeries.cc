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

static int xDataXML_regionsW_XYs_LegendreSeries_regionToTOM( statusMessageReporting *smr, xDataXML_element *XE, 
    xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, char const *wLabel, xDataTOM_axes *axes );
/*
************************************************************
*/
int xDataXML_regionsW_XYs_LegendreSeriesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int index, length;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataXML_element *XMLChild;
    char const *wLabel;
    xDataTOM_regionsW_XYs_LegendreSeries *regionsW_XYs_LegendreSeries;

    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_regionsW_XYs_LegendreSeries_ID, sizeof( xDataTOM_regionsW_XYs_LegendreSeries ) ) ) 
        == NULL ) return( 1 );
    regionsW_XYs_LegendreSeries = (xDataTOM_regionsW_XYs_LegendreSeries *) xDI->data;
    length = xDataXML_numberOfElementsByTagName( smr, XE, "region" );
    if( xDataTOM_regionsW_XYs_LegendreSeries_initialize( smr, regionsW_XYs_LegendreSeries, length, &(xDI->axes) ) != 0 ) goto err;
    if( ( wLabel = xDataTOM_axes_getLabel( smr, &(xDI->axes), 0 ) ) == NULL ) goto err;

    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( "region", XMLChild->name ) == 0 ) {
            if( xDataXML_regionsW_XYs_LegendreSeries_regionToTOM( smr, XMLChild, &(regionsW_XYs_LegendreSeries->W_XYs_LegendreSeries[index]), 
                wLabel, regionsW_XYs_LegendreSeries->axes ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in xData 'regionsW_XYs_LegendreSeries'", XMLChild->name );
            goto err;
        }
    }

    return( 0 );

err:
/* Need to free things here?????????.*/
    return( 1 );
}
/*
************************************************************
*/
static int xDataXML_regionsW_XYs_LegendreSeries_regionToTOM( statusMessageReporting *smr, xDataXML_element *XE, 
        xDataTOM_W_XYs_LegendreSeries *W_XYs_LegendreSeries, char const *wLabel, xDataTOM_axes *axes ) {

    int index, length;
    xDataXML_element *XMLChild, *interpolationAxesElement = NULL;
    xDataTOM_interpolation interpolation;
    char const *sInterpolation;

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "interpolationAxes", XMLChild->name ) == 0 ) {
            if( interpolationAxesElement != NULL ) {
                smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                    "multiple %s elements in element 'region'", XMLChild->name );
                goto err;
            }
            interpolationAxesElement = XMLChild;
        }
    }
    if( interpolationAxesElement == NULL ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "missing element 'interpolationAxes'" );
        goto err; 
    }
    if( ( sInterpolation = xDataXML_getAttributesValueInElement( interpolationAxesElement, "interpolation" ) ) == NULL ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( interpolationAxesElement ), xDataTOM_smrLibraryID, -1,
                "missing attribute 'interpolation'" );
        goto err; 
    }
    if( xDataTOM_interpolation_setFromString( smr, &interpolation, sInterpolation ) != 0 ) goto err;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    length = xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_W_XYs_LegendreSeries_initialize( smr, W_XYs_LegendreSeries, index, length, 0., xDataTOM_subAxesType_intepolationAxes, axes,
        &interpolation ) != 0 ) goto err;
    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "interpolationAxes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( wLabel, XMLChild->name ) == 0 ) {
            if( xDataXML_W_XYs_LegendreSeries_LegendreSeriesToTOM( smr, XMLChild, &(W_XYs_LegendreSeries->LegendreSeries[index]) ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in element 'region'", XMLChild->name );
            goto err;
        }
    }
    return( 0 );

err:
/* Need to free things here?????????.*/

    return( 1 );
}

#if defined __cplusplus
}
#endif
