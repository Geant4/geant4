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
int xDataXML_axesElememtToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_axes *axes ) {

    int axesProcessed = 0;
    xDataXML_element *XMLChild;

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axes", XMLChild->name ) == 0 ) {
            if( axesProcessed ) {
                smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "multiple 'axes' elements found" );
                return( 1 );
            }
            axesProcessed = 1;
            if( xDataXML_axesToTOM( smr, XMLChild, axes ) != 0 ) return( 1 );
        }
    }
    if( axesProcessed == 0 ) {
        smr_setReportError3p( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "axes element missing" );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
int xDataXML_axesToTOM( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_axes *axes ) {

    int i = 0, n = 0, index;
    xDataXML_element *XMLChild;
    char const *label, *unit, *sInterpolation, *attribute;
    xDataTOM_interpolation interpolation;

    for( XMLChild = xDataXML_getFirstElement( XE ); XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ) ) {
        if( strcmp( "axis", XMLChild->name ) != 0 ) {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                "non axis element found: name = %s", XMLChild->name );
            return( 1 );
        }
        n++;
    }
    if( xDataTOM_axes_initialize( smr, axes, n ) != 0 ) return( 1 );

    for( XMLChild = xDataXML_getFirstElement( XE ), i = 0; XMLChild != NULL; XMLChild = xDataXML_getNextElement( XMLChild ), i++ ) {
        attribute = "index";
        if( xDataXML_convertAttributeTo_xDataTOM_Int( smr, XMLChild, attribute, &index, 1 ) != 0 ) goto errA;
        attribute = "label";
        if( ( label = xDataXML_getAttributesValueInElement( XMLChild, attribute ) ) == NULL ) goto errA;
        attribute = "unit";
        if( ( unit  = xDataXML_getAttributesValueInElement( XMLChild, attribute  ) ) == NULL ) goto errA;
        if( i < ( n - 1 ) ) {
            attribute = "interpolation";
            if( ( sInterpolation = xDataXML_getAttributesValueInElement( XMLChild, attribute ) ) == NULL ) goto errA;
            if( xDataTOM_interpolation_setFromString( smr, &interpolation, sInterpolation ) != 0 ) goto err; }
        else {
            sInterpolation = "";
            if( xDataTOM_interpolation_set( smr, &interpolation, xDataTOM_interpolationFlag_linear, xDataTOM_interpolationFlag_linear, 
                xDataTOM_interpolationQualifier_dependent ) != 0 ) goto err;
        }
        xDataTOM_axis_initialize( smr, &(axes->axis[i]), index, label, unit, &interpolation );
    }
    return( 0 );

errA:
    smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, "axis missing attribute '%s'", attribute );
err:
    n = i;
    for( i = 0; i < n; i++ ) xDataTOM_axis_release( smr, &(axes->axis[i]) );
    smr_freeMemory( (void **) &(axes->axis) );
    return( 1 );
}

#if defined __cplusplus
}
#endif
