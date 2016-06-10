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

static int xDataXML_XYsDataToTOM2( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_xDataInfo *xDI, int index, int length, double value, 
    double accuracy );
/*
************************************************************
*/
int xDataXML_XYsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int dataProcessed = 0, length;
    double accuracy;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataXML_element *XMLChild;

    xDI->element = TE;
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "accuracy", &accuracy, 1 ) != 0 ) return( 1 );
    if( xDataXML_axesElememtToTOM( smr, XE, &(xDI->axes) ) != 0 ) return( 1 );
    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( "data", XMLChild->name ) == 0 ) {
            if( dataProcessed ) {
                smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "multiple 'data' elements found" );
                goto err;
            }
            dataProcessed = 1;
            if( xDataXML_XYsDataToTOM2( smr, XMLChild, xDI, -1, length, 0., accuracy ) != 0 ) goto err;
        }
    }
    if( dataProcessed == 0 ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "data element missing" );
        goto err;
    }
    return( 0 );

err:
    return( 1 );
}
/*
************************************************************
*/
static int xDataXML_XYsDataToTOM2( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_xDataInfo *xDI, int index, int length, double value, 
        double accuracy ) {

    xDataTOM_XYs *XYs;

    xDI->ID = xDataTOM_XYs_ID;
    if( ( xDI->data = (xDataTOM_XYs *) smr_malloc2( smr, sizeof( xDataTOM_XYs ), 1, "xDI->data" ) ) == NULL ) goto err;
    XYs = (xDataTOM_XYs *) xDI->data;

    if( xDataXML_XYsDataToTOM( smr, XE, XYs, index, length, value, accuracy, xDataTOM_subAxesType_proxy, 0, &(xDI->axes), NULL ) != 0 ) goto err;
    return( 0 );

err:
    smr_freeMemory( (void **) &(xDI->data) );
    return( 1 );   
}
/*
************************************************************
*/
int xDataXML_XYsDataToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_XYs *XYs, int index, int length, double value, double accuracy, 
        enum xDataTOM_subAxesType subAxesType, int axesOffest, xDataTOM_axes *axes, xDataTOM_interpolation *interpolation ) {

    XYs->index = index;
    XYs->length = length;
    XYs->value = value;
    XYs->accuracy = accuracy;
    if( xDataTOM_subAxes_initialize( smr, &(XYs->subAxes), subAxesType, axesOffest, axes, interpolation ) != 0 ) return( 1 );
    if( ( XYs->data = (double *) smr_malloc2( smr, 2 * length * sizeof( double ), 0, "XYs->data" ) ) == NULL ) goto err;

    if( xDataXML_stringToDoubles( smr, XE, XE->text.text, 2 * length, (double *) XYs->data ) != 0 ) goto err;
    return( 0 );

err:
    smr_freeMemory( (void **) &(XYs->data) );
    return( 1 );
}
/*
************************************************************
*/
int xDataXML_stringToDoubles( statusMessageReporting *smr,  xDataXML_element *XE, char const *s1, int length, double *d1 ) {

    char *e1 = (char *) s1;
    int i1;

    for( i1 = 0; i1 < length; i1++, d1++, s1 = e1 ) {
        if( xDataXML_stringTo_double( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), s1, d1, " \n", &e1 ) ) return( 1 );
    }
    while( isspace( *e1 ) ) e1++;     /* There should be nothing but white spaces left in the string. */ // Loop checking, 11.06.2015, T. Koi
    if( *e1 != 0 ) {
        smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "text contains extra data = %s", e1 );
        return( 1 );
    }
    return( 0 );
}

#if defined __cplusplus
}
#endif
