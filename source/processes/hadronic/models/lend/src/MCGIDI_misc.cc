/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <ctype.h>

#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <ptwXY.h>
#include <xDataTOM_importXML_private.h>

#include "MCGIDI.h"
#include "MCGIDI_misc.h"
#include "MCGIDI_fromTOM.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

struct ZSymbol {
    int Z;
    char const *symbol;
};

static struct ZSymbol ZSymbols[] = {    {   0, "n"  },  {   1, "H"  },  {   2, "He" },  {   3, "Li" },  {   4, "Be" },  {   5, "B"  },  {   6, "C"  },
        {   7, "N"  },  {   8, "O"  },  {   9, "F"  },  {  10, "Ne" },  {  11, "Na" },  {  12, "Mg" },  {  13, "Al" },  {  14, "Si" },  {  15, "P"  },
        {  16, "S"  },  {  17, "Cl" },  {  18, "Ar" },  {  19, "K"  },  {  20, "Ca" },  {  21, "Sc" },  {  22, "Ti" },  {  23, "V"  },  {  24, "Cr" },
        {  25, "Mn" },  {  26, "Fe" },  {  27, "Co" },  {  28, "Ni" },  {  29, "Cu" },  {  30, "Zn" },  {  31, "Ga" },  {  32, "Ge" },  {  33, "As" },
        {  34, "Se" },  {  35, "Br" },  {  36, "Kr" },  {  37, "Rb" },  {  38, "Sr" },  {  39, "Y"  },  {  40, "Zr" },  {  41, "Nb" },  {  42, "Mo" },
        {  43, "Tc" },  {  44, "Ru" },  {  45, "Rh" },  {  46, "Pd" },  {  47, "Ag" },  {  48, "Cd" },  {  49, "In" },  {  50, "Sn" },  {  51, "Sb" },
        {  52, "Te" },  {  53, "I"  },  {  54, "Xe" },  {  55, "Cs" },  {  56, "Ba" },  {  57, "La" },  {  58, "Ce" },  {  59, "Pr" },  {  60, "Nd" },
        {  61, "Pm" },  {  62, "Sm" },  {  63, "Eu" },  {  64, "Gd" },  {  65, "Tb" },  {  66, "Dy" },  {  67, "Ho" },  {  68, "Er" },  {  69, "Tm" },
        {  70, "Yb" },  {  71, "Lu" },  {  72, "Hf" },  {  73, "Ta" },  {  74, "W"  },  {  75, "Re" },  {  76, "Os" },  {  77, "Ir" },  {  78, "Pt" },
        {  79, "Au" },  {  80, "Hg" },  {  81, "Tl" },  {  82, "Pb" },  {  83, "Bi" },  {  84, "Po" },  {  85, "At" },  {  86, "Rn" },  {  87, "Fr" },
        {  88, "Ra" },  {  89, "Ac" },  {  90, "Th" },  {  91, "Pa" },  {  92, "U"  },  {  93, "Np" },  {  94, "Pu" },  {  95, "Am" },  {  96, "Cm" },
        {  97, "Bk" },  {  98, "Cf" },  {  99, "Es" },  { 100, "Fm" },  { 101, "Md" },  { 102, "No" },  { 103, "Lr" },  { 104, "Rf" },  { 105, "Db" },
        { 106, "Sg" }, { 107, "Bh" },  { 108, "Hs" },  { 109, "Mt" } };

static int MCGIDI_miscNameToZAm_getLevel( statusMessageReporting *smr, const char *name, const char *p );
static ptwXYPoints *MCGIDI_misc_Data2ptwXYPointsInUnitsOf( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
    int length, double *data, char const *fromUnits[2], char const *toUnits[2] );
/*
************************************************************
*/
int MCGIDI_misc_NumberOfZSymbols( void ) {

    return( sizeof( ZSymbols ) / sizeof( struct ZSymbol ) );
}
/*
************************************************************
*/
const char *MCGIDI_misc_ZToSymbol( int iZ ) {

    if( ( iZ < 0 ) || ( iZ >= MCGIDI_misc_NumberOfZSymbols( ) ) ) return( NULL );
    return( ZSymbols[iZ].symbol );
}
/*
************************************************************
*/
int MCGIDI_misc_symbolToZ( const char *Z ) {

    int i, n = MCGIDI_misc_NumberOfZSymbols( );

    for( i = 0; i < n; i++ ) {
        if( strcmp( Z, ZSymbols[i].symbol ) == 0 ) return( ZSymbols[i].Z );
    }
    return( -1 );
}
/*
************************************************************
*/
int MCGIDI_miscNameToZAm( statusMessageReporting *smr, const char *name, int *Z, int *A, int *m, int *level ) {

    const char *p;
    char s[1024] = "", *q, *e;   /* Note 1) routine will fail when parts of a particle name can be longer than 1024. */

    if( strlen( name ) >= ( sizeof( s ) - 1 ) ) {
        smr_setReportError2( smr, smr_unknownID, 0, "particle name too long: '%s'", name );
        return( 1 );
    }

    *Z = *A = *m = *level = 0;
    if( ( !strncmp( "FissionProduct", name, 14 ) ) || !strncmp( "99120", name, 5 ) ) {
        *Z = 99;
        *A = 120;
        return( 0 );
    }
    if( strcmp( "gamma", name ) == 0 ) return( 0 );
    if( strcmp( "n", name ) == 0 ) { *A = 1; return( 0 ); }

    for( p = name, q = s; ( *p != 0 ) && !isdigit( *p ) && ( *p != '_' ); p++, q++ ) *q = *p;   /* '_' only for "natural". */
    if( *p == 0 ) {
        smr_setReportError2( smr, smr_unknownID, 0, "unsupport particle name = '%s'", name );
        return( 1 );
    }
    *q = 0;
    if( ( *Z = MCGIDI_misc_symbolToZ( s ) ) < 0 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "Particle %s's symbol = '%s' not found", name, s ); }
    else {                  /* Getting here implies that *p is a digit. */
        if( *p == '_' ) {
            if( strncmp( p, "_natural", 8 ) == 0 ) {
                p += 8;
                if( *p ) *level = MCGIDI_miscNameToZAm_getLevel( smr, name, p ); }
            else {
                smr_setReportError2( smr, smr_unknownID, 0, "expecting 'natural': %s", name );
            } }
        else {
            for( q = s; isdigit( *p ); p++, q++ ) *q = *p;
            *q = 0;
            if( strcmp( s, "natural" ) == 0 ) {
                e = s;
                while( *e ) e++; /* Loop checking, 11.06.2015, T. Koi*/ }
            else {
                *A = (int) strtol( s, &e, 10 );
            }
            if( *e != 0 ) {
                smr_setReportError2( smr, smr_unknownID, 1, "Failed to convert A to integer in particle name %s", name ); }
            else {          /* Getting here implies that *p == '_' or 0. */
                if( *p ) *level = MCGIDI_miscNameToZAm_getLevel( smr, name, p );
            }
        }
    }

    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static int MCGIDI_miscNameToZAm_getLevel( statusMessageReporting *smr, const char *name, const char *p ) {

    int level = 0;
    char *e;

    if( *p == '_' ) {
        p++;
        switch( *p ) {
        case 'e' :
            p++;
            level = (int) strtol( p, &e, 10 );
            if( *e != 0 ) smr_setReportError2( smr, smr_unknownID, 1, "Failed to convert level to integer in particle name %s", name );
            break;
        case 'c' :
            level = MCGIDI_particleLevel_continuum;
            break;
        case 's' :
            level = MCGIDI_particleLevel_sum;
            break;
        default :
            smr_setReportError2( smr, smr_unknownID, 0, "invalid 'natural': %s", name );
        } }
    else {
        smr_setReportError2( smr, smr_unknownID, 0, "invalid level specifier: %s", name );
    }
    return( level );
}
/*
************************************************************
*/
char const *MCGIDI_misc_pointerToTOMAttributeIfAllOk( statusMessageReporting *smr, const char *path, int required, 
        xDataTOM_attributionList *attributes, const char *name, const char *file, int line ) {

    char const *value;

    if( !smr_isOk( smr ) ) return( NULL );
    if( ( value = xDataTOMAL_getAttributesValue( attributes, name ) ) == NULL ) {
        if( required ) {
            smr_setReportError( smr, NULL, file, line, __func__, smr_unknownID, 1, "element does not have attribute named %s for file = %d", name, path );
        }
    }
    return( value );
}
/*
************************************************************
*/
char const *MCGIDI_misc_pointerToAttributeIfAllOk( statusMessageReporting *smr, xDataXML_element *element, const char *path, int required, 
        xDataTOM_attributionList *attributes, const char *name, const char *file, int line ) {

    char const *value;

    if( !smr_isOk( smr ) ) return( NULL );
    if( ( value = xDataTOMAL_getAttributesValue( attributes, name ) ) == NULL ) {
        if( required ) {
            if( element != NULL ) {
                MCGIDI_misc_setMessageError_Element( smr, NULL, element, file, line, 1, "element does not have attribute named %s", name ); }
            else {
                smr_setReportError( smr, NULL, file, line, __func__, smr_unknownID, 1, "element does not have attribute named %s for file = %d", name, path );
            }
        }
    }
    return( value );
}
/*
************************************************************
*/
int MCGIDI_misc_setMessageError_Element( statusMessageReporting *smr, void *userInterface, xDataXML_element *element, const char *file, int line, int code, 
    const char *fmt, ... ) {

    int status = 0;
    va_list args;
    char *msg;

    va_start( args, fmt );
    msg = smr_vallocateFormatMessage( fmt, &args );
    va_end( args );
    if( msg == NULL ) {
        status = 1;
        va_start( args, fmt );
        smr_vsetReportError( smr, userInterface, file, line, __func__, smr_unknownID, code, fmt, &args );
        va_end( args ); }
    else {
        status = smr_setReportError( smr, userInterface, file, line, __func__, smr_unknownID, code, "%s for element %s", msg, element->name );
        smr_freeMemory( (void **) &msg );
    }
    return( status );
}
/*
************************************************************
*/
xDataTOM_Int MCGIDI_misc_binarySearch( xDataTOM_Int n, double *ds, double d ) {
/*
*   Returns -2 is d < first point of ds, -1 if > last point of ds and the lower index of ds otherwise.
*/
    xDataTOM_Int imin = 0, imid, imax = n - 1;

    if( d < ds[0] ) return( -2 );
    if( d > ds[n-1] ) return( -1 );
    while( 1 ) { // Loop checking, 11.06.2015, T. Koi
        imid = ( imin + imax ) >> 1;
        if( imid == imin ) break;
        if( d < ds[imid] ) {
            imax = imid; }
        else {
            imin = imid;
        }
    }
    return( imin );
}
/*
************************************************************
*/
char *MCGIDI_misc_getAbsPath( statusMessageReporting *smr, const char *fileName ) {
/*
*   User must free returned string.
*/
    int n = (int) strlen( fileName ) + 1, nCwd = 0;
    char *absPath, cwd[4 * 1024] = "", *p, *needle;

    if( fileName[0] != '/' ) {
        //if( getcwd( cwd, sizeof( cwd ) + 1 ) == NULL ) {
        //TK modified above line for compiler(gcc.4.8) warning message
        if( getcwd( cwd, sizeof( cwd ) ) == NULL ) {
            smr_setReportError2p( smr, smr_unknownID, -1, "hardwired cwd too small" );
            return( NULL );
        }
        nCwd = (int) strlen( cwd );
        n += nCwd + 1;                                  /* cwd + '/'. */
    }
    if( ( absPath = (char *) smr_malloc2( smr, n, 0, "absPath" ) ) == NULL ) return( NULL );
    if( fileName[0] != '/' ) {
        strcpy( absPath, cwd );
        strcat( absPath, "/" );
        strcat( absPath, fileName ); }
    else {
        strcpy( absPath, fileName );
    }

    while( 1 ) {                                        /* Remove all ./ from path. */ // Loop checking, 11.06.2015, T. Koi
        if( ( needle = strstr( absPath, "/./" ) ) == NULL ) break;
        p = needle;
        for( needle += 2; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    }
    while( 1 ) {                                        /* Remove all ../ from path. */ // Loop checking, 11.06.2015, T. Koi
        if( ( needle = strstr( absPath, "/../" ) ) == NULL ) break;
        p = needle - 1;
        while( ( p > absPath ) && ( *p != '/' ) ) p--; // Loop checking, 11.06.2015, T. Koi
        if( *p != '/' ) break;                           /* This should not happen if path is legit, I think, and I do not know what to do so will leave it. */
        if( p == absPath ) break;                       /* Ditto. */
        for( needle += 3; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    }
    return( absPath );
}
/*
************************************************************
*/
int MCGIDI_misc_copyXMLAttributesToTOM( statusMessageReporting *smr, xDataTOM_attributionList *TOM, xDataXML_attributionList *XML ) {

    int i;
    xDataXML_attribute *attribute;

    xDataTOMAL_initial( smr, TOM );
    for( i = 0; ; i++ ) {
        if( ( attribute = xDataXML_attributeByIndex( XML, i ) ) == NULL ) break;
        if( xDataTOMAL_addAttribute( smr, TOM, attribute->name, attribute->value ) != 0 ) goto err;
    }
    return( 0 );

err:
    xDataTOMAL_release( TOM );
    return( 1 );
}
/*
************************************************************
*/
enum xDataTOM_frame MCGIDI_misc_getProductFrame( statusMessageReporting *smr, xDataTOM_element *frameElement ) {

    char const *frameString;
    enum xDataTOM_frame frame = xDataTOM_frame_invalid;

    if( ( frameString = xDataTOM_getAttributesValueInElement( frameElement, MCGIDI_token_productFrame ) ) != NULL ) {
        if( ( frame = xDataTOM_axis_stringToFrame( smr, frameString ) ) == xDataTOM_frame_invalid ) {
            smr_setReportError2( smr, smr_unknownID, 1, "Invalid frame = '%s'", frameString );
        }
    }
    return( frame );
}
/*
************************************************************
*/
int MCGIDI_misc_PQUStringToDouble( statusMessageReporting *smr, char const *str, char const *unit, double conversion, double *value ) {
/*
*   Currently, white spaces are not allowed after the unit.
*
*   Examples of allowed strings are: "2.39e6 eV", " 2.39e6eV" and " 2.39e6 eV".
*/
    char const *s = str;
    char *e;


    while( isspace( *s ) ) s++; // Loop checking, 11.06.2015, T. Koi
    *value = strtod( s, &e ) * conversion;
    if( e == s ) {
        smr_setReportError2( smr, smr_unknownID, 1, "no number at start of string = <%s>", str );
        return( 1 );
    }
    while( isspace( *e ) ) e++; // Loop checking, 11.06.2015, T. Koi
    if( strcmp( e, unit ) != 0 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "unit = '%s' not '%s' in '%s'", e, unit, str );
        return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_misc_PQUStringToDoubleInUnitOf( statusMessageReporting *smr, char const *str, char const *toUnit, double *value ) {
/*
*   Currently, white spaces are not allowed after the unit.
*
*   Examples of allowed strings are: "2.39e6 eV", " 2.39e6eV" and " 2.39e6 eV".
*/
    char const *s1 = str;
    char *e1;
    double factor;

    while( isspace( *s1 ) ) s1++; // Loop checking, 11.06.2015, T. Koi
    *value = strtod( s1, &e1 );
    if( e1 == s1 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "no number at start of string = <%s>", str );
        return( 1 );
    }
    while( isspace( *e1 ) ) e1++; // Loop checking, 11.06.2015, T. Koi

    factor = MCGIDI_misc_getUnitConversionFactor( smr, e1, toUnit );
    *value *= factor;
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
double MCGIDI_misc_getUnitConversionFactor( statusMessageReporting *smr, char const *fromUnit, char const *toUnit ) {
/*
*   This is a kludge until units are better supported.
*/
    if( strcmp( fromUnit, toUnit ) == 0 ) return( 1.0 );

    if( strcmp( fromUnit, "eV" ) == 0 ) {
        if( strcmp( toUnit, "MeV" ) == 0 ) return( 1e-6 ); }
    else if( strcmp( fromUnit, "MeV" ) == 0 ) {
        if( strcmp( toUnit, "eV" ) == 0 ) return( 1e+6 ); }
    else if( strcmp( fromUnit, "1/eV" ) == 0 ) {
        if( strcmp( toUnit, "1/MeV" ) == 0 ) return( 1e+6 ); }
    else if( strcmp( fromUnit, "1/MeV" ) == 0 ) {
        if( strcmp( toUnit, "1/eV" ) == 0 ) return( 1e-6 ); }
    else if( strcmp( fromUnit, "K" ) == 0 ) {
        if( strcmp( toUnit, "MeV/k" ) == 0 ) return( 8.617343183775137e-11 );
    }

    smr_setReportError2( smr, smr_unknownID, 1, "Cannot convert unit '%s' to unit '%s'", fromUnit, toUnit );
    return( 1.0 );
}
/*
************************************************************
*/
ptwXYPoints *MCGIDI_misc_dataFromXYs2ptwXYPointsInUnitsOf( statusMessageReporting *smr, xDataTOM_XYs *XYs, 
        ptwXY_interpolation interpolation, char const *toUnits[2] ) {

    int length;
    double *data;
    char const *fromUnits[2];

    fromUnits[0] = xDataTOM_subAxes_getUnit( smr, &(XYs->subAxes), 0 );
    if( !smr_isOk( smr ) ) return( NULL );
    fromUnits[1] = xDataTOM_subAxes_getUnit( smr, &(XYs->subAxes), 1 );
    if( !smr_isOk( smr ) ) return( NULL );

    length = xDataTOM_XYs_getData( XYs, &data );

    return( MCGIDI_misc_Data2ptwXYPointsInUnitsOf( smr, interpolation, length, data, fromUnits, toUnits ) );
}
/*
************************************************************
*/
ptwXYPoints *MCGIDI_misc_dataFromElement2ptwXYPointsInUnitsOf( statusMessageReporting *smr, xDataTOM_element *linear, char const *toUnits[2] ) {

    int length;
    double *data;
    xDataTOM_axes *axes = &(linear->xDataInfo.axes);
    char const *fromUnits[2];
    ptwXY_interpolation interpolation;

    if( axes->numberOfAxes != 2 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "axes must have 2 axis, it has %d", axes->numberOfAxes );
        return( NULL );
    }

    if( MCGIDI_fromTOM_interpolation( smr, linear, 0, &interpolation ) != 0 ) return( NULL );
    fromUnits[0] = axes->axis[0].unit;
    fromUnits[1] = axes->axis[1].unit;

    length = xDataTOM_XYs_getDataFromXDataInfo( (xDataTOM_xDataInfo *) &(linear->xDataInfo), &data );
    return( MCGIDI_misc_Data2ptwXYPointsInUnitsOf( smr, interpolation, length, data, fromUnits, toUnits ) );
}
/*
************************************************************
*/
static ptwXYPoints *MCGIDI_misc_Data2ptwXYPointsInUnitsOf( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
        int length, double *data, char const *fromUnits[2], char const *toUnits[2] ) {

    double xFactor, yFactor;
    ptwXYPoints *ptwXY = NULL;
    nfu_status status;

    xFactor = MCGIDI_misc_getUnitConversionFactor( smr, fromUnits[0], toUnits[0] );
    if( !smr_isOk( smr ) ) goto err;
    yFactor = MCGIDI_misc_getUnitConversionFactor( smr, fromUnits[1], toUnits[1] );
    if( !smr_isOk( smr ) ) goto err;


    ptwXY = ptwXY_create( interpolation, NULL, 2., 1e-3, length, 10, length, data, &status, 0 );
    if( status != nfu_Okay ) {
        smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_create err = %d: %s\n", status, nfu_statusMessage( status ) );
        goto err;
    }

    if( ( xFactor != 1. ) || ( yFactor != 1. ) ) {
        if( ( status = ptwXY_scaleOffsetXAndY( ptwXY, xFactor, 0., yFactor, 0. ) ) != nfu_Okay ) {
            smr_setReportError2( smr, smr_unknownID, 1, "ptwXY_scaleOffsetXAndY err = %d: %s\n", status, nfu_statusMessage( status ) );
            goto err;
        }
    }

    return( ptwXY );

err:
    if( ptwXY != NULL ) ptwXY_free( ptwXY );
    return( NULL );
}
/*
************************************************************
*/
void MCGIDI_misc_updateTransportabilitiesMap( transportabilitiesMap *transportabilities, int PoPID, enum MCGIDI_transportability transportability ) {

    transportabilitiesMap::iterator iter = transportabilities->find( PoPID );

    if( iter != transportabilities->end( ) ) {
        switch ( iter->second ) {
        case MCGIDI_transportability_unknown :
            break;
        case MCGIDI_transportability_none :
            switch( transportability ) {
                case MCGIDI_transportability_unknown :
                case MCGIDI_transportability_none :
                    transportability = MCGIDI_transportability_none;
                    break;
                case MCGIDI_transportability_partial :
                    break;
                case MCGIDI_transportability_full :
                    transportability = MCGIDI_transportability_partial;
                    break;
            }
            break;
        case MCGIDI_transportability_partial :
            transportability = MCGIDI_transportability_partial;
            break;
        case MCGIDI_transportability_full :
            switch( transportability ) {
                case MCGIDI_transportability_none :
                case MCGIDI_transportability_partial :
                    transportability = MCGIDI_transportability_partial;
                    break;
                case MCGIDI_transportability_unknown :
                case MCGIDI_transportability_full :
                    break;
            }
            break;
        }
    }
    (*transportabilities)[PoPID] = transportability;
}
/*
************************************************************
*/
void MCGIDI_misc_updateTransportabilitiesMap2( transportabilitiesMap *transportabilities, int PoPID, int transportable ) {

    MCGIDI_misc_updateTransportabilitiesMap( transportabilities, PoPID, ( transportable ? MCGIDI_transportability_full : MCGIDI_transportability_none ) );
}

#if defined __cplusplus
}
#endif

