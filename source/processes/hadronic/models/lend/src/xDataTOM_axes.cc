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

static char const *xDataTOM_frame_labString = "lab";
static char const *xDataTOM_frame_centerOfMassString = "centerOfMass";
static char const *xDataTOM_frame_invalidString = "invalid";
/*
************************************************************
*/
int xDataTOM_axes_initialize( statusMessageReporting *smr, xDataTOM_axes *axes, int numberOfAxes ) {

    axes->numberOfAxes = 0;
    if( ( axes->axis = (xDataTOM_axis *) smr_malloc2( smr, numberOfAxes * sizeof( xDataTOM_axis ), 1, "axes->axis" ) ) == NULL ) return( 1 );
    axes->numberOfAxes = numberOfAxes;
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_axes_release( xDataTOM_axes *axes ) {

    int i;

    for( i = 0; i < axes->numberOfAxes; i++ ) {
        xDataTOM_axis_release( NULL, &(axes->axis[i]) );
    }
    smr_freeMemory( (void **) &(axes->axis) );
    return( 0 );
}
/*
************************************************************
*/
char const *xDataTOM_axes_getLabel( statusMessageReporting *smr, xDataTOM_axes *axes, int index ) {

    if( ( index < 0 ) || ( index >= axes->numberOfAxes ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid axes index = %d", index );
        return( NULL );
    }
    return( axes->axis[index].label );
}
/*
************************************************************
*/
char const *xDataTOM_axes_getUnit( statusMessageReporting *smr, xDataTOM_axes *axes, int index ) {

    if( ( index < 0 ) || ( index >= axes->numberOfAxes ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid axes index = %d", index );
        return( NULL );
    }
    return( axes->axis[index].unit );
}
/*
************************************************************
*/
int xDataTOM_axes_getInterpolation( statusMessageReporting *smr, xDataTOM_axes *axes, int index, enum xDataTOM_interpolationFlag *independent, 
        enum xDataTOM_interpolationFlag *dependent, enum xDataTOM_interpolationQualifier *qualifier ) {

    xDataTOM_interpolation *interpolation;

    if( ( index < 0 ) || ( index >= axes->numberOfAxes ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid axes index = %d", index );
        return( 1 );
    }
    interpolation = &(axes->axis[index].interpolation);
    *independent = interpolation->independent;
    *dependent = interpolation->dependent;
    *qualifier = interpolation->qualifier;

    return( 0 );
}

/*
c   subAxes functions.
*/
/*
************************************************************
*/
int xDataTOM_subAxes_initialize( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, enum xDataTOM_subAxesType type, int offset, 
        xDataTOM_axes *axes, xDataTOM_interpolation *interpolation ) {

    subAxes->type = type;
    if( axes == NULL ) {
        smr_setReportError2p( smr, xDataTOM_smrLibraryID, -1, "Axes must not be NULL" );
        return( 1 );
    }
    subAxes->offset = offset;
    if( ( offset < 0 ) || ( offset >= axes->numberOfAxes ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "offset = %d < 0 or >= axes->numberOfAxes = %d", offset, axes->numberOfAxes );
        return( 1 );
    }
    if( type == xDataTOM_subAxesType_intepolationAxes ) {
        if( interpolation == NULL ) {
            smr_setReportError2p( smr, xDataTOM_smrLibraryID, -1, "Interpolation must not be NULL for intepolationAxes" );
            return( 1 );
        }
        if( xDataTOM_interpolation_copy( smr, &(subAxes->interpolation), interpolation ) ) return( 1 ); }
    else {      /* Not used but fill in anyway. */
        xDataTOM_interpolation_set( smr, &(subAxes->interpolation), xDataTOM_interpolationFlag_linear, xDataTOM_interpolationFlag_linear,
            xDataTOM_interpolationQualifier_none );
    }
    subAxes->axes = axes;
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_subAxes_release( xDataTOM_subAxes *subAxes ) {

    subAxes->axes = NULL;
    return( 0 );
}
/*
************************************************************
*/
char const *xDataTOM_subAxes_getLabel( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, int index ) {

    return( xDataTOM_axes_getLabel( smr, subAxes->axes, index + subAxes->offset ) );
}
/*
************************************************************
*/
char const *xDataTOM_subAxes_getUnit( statusMessageReporting *smr, xDataTOM_subAxes *subAxes, int index ) {

    return( xDataTOM_axes_getUnit( smr, subAxes->axes, index + subAxes->offset ) );
}

/*
c   Axis functions.
*/
/*
************************************************************
*/
xDataTOM_axis *xDataTOM_axis_new( statusMessageReporting *smr, int index, char const *label, char const *unit, xDataTOM_interpolation *interpolation ) {

    xDataTOM_axis *axis = NULL;

    if( ( axis = (xDataTOM_axis *) smr_malloc2( smr, sizeof( xDataTOM_axis ), 0, "axis" ) ) == NULL ) return( NULL );
    if( xDataTOM_axis_initialize( smr, axis, index, label, unit, interpolation ) != 0 ) smr_freeMemory( (void **) &axis );
    return( axis );
}
/*
************************************************************
*/
int xDataTOM_axis_initialize( statusMessageReporting *smr, xDataTOM_axis *axis, int index, char const *label, char const *unit, xDataTOM_interpolation *interpolation ) {

    axis->index = index;
    if( ( axis->label = smr_allocateCopyString2( smr, label, "label" ) ) == NULL ) goto err;
    if( ( axis->unit = smr_allocateCopyString2( smr, unit, "unit" ) ) == NULL ) goto err;
    if( xDataTOM_interpolation_copy( smr, &(axis->interpolation), interpolation ) != 0 ) goto err;

    return( 0 );

err:
    smr_freeMemory( (void **) &(axis->label) );
    smr_freeMemory( (void **) &(axis->unit) );
    return( 1 );
}
/*
************************************************************
*/
xDataTOM_axis *xDataTOM_axis_release( statusMessageReporting * /*smr*/, xDataTOM_axis *axis ) {

    axis->index = -1;
    smr_freeMemory( (void **) &(axis->label) );
    smr_freeMemory( (void **) &(axis->unit) );
    return( NULL );
}
/*
************************************************************
*/
enum xDataTOM_frame xDataTOM_axis_stringToFrame( statusMessageReporting * /*smr*/, char const *frame ) {

    if( strcmp( "lab", frame ) == 0 ) return( xDataTOM_frame_lab );
    if( strcmp( "centerOfMass", frame ) == 0 ) return( xDataTOM_frame_centerOfMass );
    return( xDataTOM_frame_invalid );
}
/*
************************************************************
*/
char const *xDataTOM_axis_frameToString( statusMessageReporting * /*smr*/, enum xDataTOM_frame frame ) {

    switch( frame ) {
    case xDataTOM_frame_lab : return( xDataTOM_frame_labString );
    case xDataTOM_frame_centerOfMass : return( xDataTOM_frame_centerOfMassString );
    default :
        break;
    }
    return( xDataTOM_frame_invalidString );
}

#if defined __cplusplus
}
#endif
