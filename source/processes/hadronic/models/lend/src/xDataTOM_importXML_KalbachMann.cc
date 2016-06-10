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

static int xDataXML_KalbachMannCoefficientsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_KalbachMannCoefficients *coefficients );
/*
************************************************************
*/
int xDataXML_KalbachMannToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE ) {

    int length;
    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);
    xDataTOM_KalbachMann *KalbachMann;
    char const *wLabel, *form;
    xDataXML_element *XMLChild;
    xDataTOM_axes *axes = &(xDI->axes);

/* Need to release KalbachMann if an error occurs later. */
    if( ( xDI->data = xDataXML_initializeData( smr, XE, TE, xDataTOM_KalbachMann_ID, sizeof( xDataTOM_KalbachMann ) ) ) == NULL ) return( 1 );
    KalbachMann = (xDataTOM_KalbachMann *) xDI->data;

    if( ( form = xDataXML_getAttributesValueInElement( XE, "form" ) ) == NULL ) goto err;
    if( strcmp( form, "fr" ) == 0 ) {
        KalbachMann->type = xDataTOM_KalbachMannType_fr; }
    else if( strcmp( form, "fra" ) == 0 ) {
        KalbachMann->type = xDataTOM_KalbachMannType_fra; }
    else {
        smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid KalbachMann type - '%s'", form );
        goto err;
    }
    if( ( wLabel = xDataTOM_axes_getLabel( smr, axes, 0 ) ) == NULL ) goto err;
    length = xDataXML_numberOfElementsByTagName( smr, XE, wLabel );
    if( xDataTOM_KalbachMann_initialize( smr, KalbachMann, length, axes ) != 0 ) return( 1 );

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            continue; }
        else if( strcmp( wLabel, XMLChild->name ) == 0 ) {
            if( xDataXML_KalbachMannCoefficientsToTOM( smr, XMLChild, &(KalbachMann->coefficients[KalbachMann->numberOfEnergies]) ) != 0 ) goto err;
            KalbachMann->numberOfEnergies++; }
        else {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1,
                "invalid element '%s' in xData = 'KalbachMann'", XMLChild->name );
            goto err;
        }
    }

    return( 0 );

err:
    smr_freeMemory( (void **) &(xDI->data) );
    return( 1 );
}
/*
************************************************************
*/
static int xDataXML_KalbachMannCoefficientsToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_KalbachMannCoefficients *coefficients ) {

    int index, length;
    double value;

    coefficients->coefficients = NULL;
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "index", &index, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XE, "length", &length, 1 ) != 0 ) return( 1 );
    if( xDataXML_convertAttributeToDouble( smr, XE, "value", &value, 1 ) != 0 ) return( 1 );
    coefficients->index = index;
    coefficients->length = length;
    coefficients->value = value;
    if( ( coefficients->coefficients = (double *) smr_malloc2( smr, length * sizeof( double ), 0, "coefficients->coefficients" ) ) == NULL ) goto err;
    if( xDataXML_stringToDoubles( smr, XE, XE->text.text, length, (double *) coefficients->coefficients ) != 0 ) goto err;
    return( 0 );

err:
    if( coefficients->coefficients != NULL ) smr_freeMemory( (void **) &(coefficients->coefficients) );
    return( 1 );
}

#if defined __cplusplus
}
#endif
