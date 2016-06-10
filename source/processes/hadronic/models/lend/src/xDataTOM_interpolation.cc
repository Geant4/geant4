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

#define dependentAxis 1
#define allowByRegion 2

static enum xDataTOM_interpolationFlag xDataTOM_interpolation_getFromString( statusMessageReporting *smr, char const *s, char const **e,
    char const *str, int flag );
/*
************************************************************
*/
int xDataTOM_interpolation_set( statusMessageReporting *smr, xDataTOM_interpolation *interpolation, enum xDataTOM_interpolationFlag independent, 
    enum xDataTOM_interpolationFlag dependent,  enum xDataTOM_interpolationQualifier qualifier ) {

    if( ( independent < xDataTOM_interpolationFlag_linear ) || ( independent > xDataTOM_interpolationFlag_byRegion ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid independent interpolation = %d", independent );
        return( 1 );
    }
    if( ( dependent < xDataTOM_interpolationFlag_linear ) || ( dependent > xDataTOM_interpolationFlag_flat ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid dependent interpolation = %d", dependent );
        return( 1 );
    }
    if( ( qualifier <= xDataTOM_interpolationQualifier_invalid ) || ( qualifier > xDataTOM_interpolationQualifier_correspondingPoints ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid interpolation qualifier = %d", qualifier );
        return( 1 );
    }

    interpolation->independent = independent;
    interpolation->dependent = dependent;
    interpolation->qualifier = qualifier;
    return( 0 );
}
/*
************************************************************
*/
int xDataTOM_interpolation_setFromString( statusMessageReporting *smr, xDataTOM_interpolation *interpolation, char const *str ) {

    int flag = 0;
    char const *c, *e;
    enum xDataTOM_interpolationQualifier qualifier = xDataTOM_interpolationQualifier_none;
    enum xDataTOM_interpolationFlag independent, dependent;

    if( ( c = strchr( str, ':' ) ) != NULL ) {
        if( strncmp( "unitBase:", str, 9 ) == 0 ) {
            qualifier = xDataTOM_interpolationQualifier_unitBase; }
        else if( strncmp( "correspondingPoints:", str, 20 ) == 0 ) {
            qualifier = xDataTOM_interpolationQualifier_correspondingPoints; }
        else {
            smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid interpolation string qualifier '%s'", str );
            return( 1 );
        }
        c++; }
    else {
        c = str;
    }
    if( ( independent = xDataTOM_interpolation_getFromString( smr, c, &e, str, flag ) ) == xDataTOM_interpolationFlag_invalid ) return( 1 );
    if( *e != ',' ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "missing ',' separator in interpolation string'%s'", str );
        return( 1 );
    }
    c = ++e;
    flag |= dependentAxis;
    if( ( dependent   = xDataTOM_interpolation_getFromString( smr, c, &e, str, flag ) ) == xDataTOM_interpolationFlag_invalid ) return( 1 );
    xDataTOM_interpolation_set( smr, interpolation, independent, dependent, qualifier );
    return( 0 );
}
/*
************************************************************
*/
static enum xDataTOM_interpolationFlag xDataTOM_interpolation_getFromString( statusMessageReporting *smr, char const *s, char const **e,
    char const *str, int flag ) {

    if( strncmp( "linear", s, 6 ) == 0 ) { *e = &(s[6]); return( xDataTOM_interpolationFlag_linear ); }
    if( strncmp( "log", s, 3 ) == 0 ) { *e = &(s[3]); return( xDataTOM_interpolationFlag_log ); }
    if( flag | allowByRegion ) {
        if( strncmp( "byRegion", s, 8 ) == 0 ) { *e = &(s[8]); return( xDataTOM_interpolationFlag_byRegion ); }
    }
    if( flag | dependentAxis ) {
        if( strncmp( "flat", s, 4 ) == 0 ) { *e = &(s[4]); return( xDataTOM_interpolationFlag_flat ); }
    }
    smr_setReportError2( smr, xDataTOM_smrLibraryID, -1, "invalid interpolation component '%s' in string '%s'", s, str );
    return( xDataTOM_interpolationFlag_invalid );

/*  Currently not supported.
    otherToken = 'other'
    chargedParticleToken = 'charged-particle'
*/
}
/*
************************************************************
*/
int xDataTOM_interpolation_copy( statusMessageReporting *smr, xDataTOM_interpolation *desc, xDataTOM_interpolation *src ) {

    return( xDataTOM_interpolation_set( smr, desc, src->independent, src->dependent, src->qualifier ) );
}

#if defined __cplusplus
}
#endif
