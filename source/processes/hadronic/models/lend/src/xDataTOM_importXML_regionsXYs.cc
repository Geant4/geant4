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

static int xDataXML_regionsXYs_regionToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes );
static int xDataXML_regionsXYs_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes, 
        xDataTOM_interpolation *interpolation, int index, int length, double accuracy );
/*
************************************************************
*/
int xDataXML_regionsXYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int index;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataXML_element *XMLChild;
    xDataTOM_regionsXYs *regionsXYs;

    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_regionsXYs_ID, sizeof( xDataTOM_regionsXYs ) ) ) == NULL ) return( 1 );
    regionsXYs = (xDataTOM_regionsXYs *) xDI->data;
    regionsXYs->axes = &(xDI->axes);
    regionsXYs->length = xDataXML_numberOfElementsByTagName( smr, XE, "region" );
    if( ( regionsXYs->XYs = (xDataTOM_XYs *) smr_malloc2( smr, regionsXYs->length * sizeof( xDataTOM_XYs ), 1, "regionsXYs->XYs" ) ) == NULL ) goto err;


    for( XMLChild = xDataXML_getFirstElement( XE ), index = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( "region", XMLChild->name ) == 0 ) {
            if( xDataXML_regionsXYs_regionToTOM( smr, XMLChild, &(regionsXYs->XYs[index]), regionsXYs->axes ) != 0 ) goto err;
            index++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in xData 'regionsXYs'", XMLChild->name );
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
static int xDataXML_regionsXYs_regionToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes ) {

    int index, length;
    double accuracy;
    xDataXML_element *XMLChild, *interpolationAxesElement = NULL, *dataElement = NULL;
    xDataTOM_interpolation interpolation;
    char const *sInterpolation;

    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "accuracy", &accuracy, 1 ) != 0 ) return( 1 );

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "interpolationAxes", XMLChild->name ) == 0 ) {
            if( interpolationAxesElement != NULL ) {
                smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                    "multiple %s elements in element 'region'", XMLChild->name );
                goto err;
            }
            interpolationAxesElement = XMLChild; }
        else if( strcmp( "data", XMLChild->name ) == 0 ) {
            if( dataElement != NULL ) {
                smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                    "multiple %s elements in element 'region'", XMLChild->name );
                goto err;
            }
            dataElement = XMLChild; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in element 'region'", XMLChild->name );
            goto err;
        }
    }
    if( interpolationAxesElement == NULL ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                "missing 'interpolationAxes' element in element 'region'" );
        goto err;
    }
    if( ( sInterpolation = xDataXML_getAttributesValueInElement( interpolationAxesElement, "interpolation" ) ) == NULL ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( interpolationAxesElement ), xDataTOM_smrLibraryID, -1,
                "missing attribute 'interpolation'" );
        goto err; 
    }
    if( xDataTOM_interpolation_setFromString( smr, &interpolation, sInterpolation ) != 0 ) goto err;
    if( dataElement == NULL ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                "missing 'data' element in element 'region'" );
        goto err;
    }
    xDataXML_regionsXYs_XYsToTOM( smr, dataElement, XYs, axes, &interpolation, index, length, accuracy );
    return( 0 );

err:
/* Need to free things here?????????.*/

    return( 1 );
}
/*
************************************************************
*/
static int xDataXML_regionsXYs_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, xDataTOM_axes *axes, 
        xDataTOM_interpolation *interpolation, int index, int length, double accuracy ) {

    return( xDataXML_XYsDataToTOM( smr, XE, XYs, index, length, 0., accuracy, xDataTOM_subAxesType_intepolationAxes, 0, axes, interpolation ) );
}

#if defined __cplusplus
}
#endif
