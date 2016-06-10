/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include "xDataTOM_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

#define nameValueLength 1024

static xDataTOM_element *xDataTOM_getLinksElement2( statusMessageReporting *smr, xDataTOM_element *element, char const *link, char const *fullLink );
static int xDataTOM_getLinksElement3( statusMessageReporting *smr, char const *nameValue, char *name, char *value, char const *fullLink );
/*
************************************************************
*/
char *xDataTOMMisc_getAbsPath( statusMessageReporting *smr, const char *fileName ) {
/*
*   User must free returned string.
*/
    int n = (int) strlen( fileName ) + 1, nCwd = 0;
    char *absPath, cwd[4 * 1024] = "", *p, *needle;

    if( fileName[0] != '/' ) {
        //if( getcwd( cwd, sizeof( cwd ) + 1 ) == NULL ) {
        //TK modified above line for compiler(gcc.4.8) warning message
        if( getcwd( cwd, sizeof( cwd ) ) == NULL ) {
            smr_setReportError2p( smr, xDataTOM_smrLibraryID, -1, "hardwired cwd too small" );
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

    while( 1 ) {                                        /* Remove all ./ from path. */
        if( ( needle = strstr( absPath, "/./" ) ) == NULL ) break; 
        p = needle;
        for( needle += 2; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    } // Loop checking, 11.06.2015, T. Koi

    while( 1 ) {                                        /* Remove all ../ from path. */
        if( ( needle = strstr( absPath, "/../" ) ) == NULL ) break;
        p = needle - 1;
        while( ( p > absPath ) && ( *p != '/' ) ) p--; // Loop checking, 11.06.2015, T. Koi
        if( *p != '/' ) break;                           /* This should not happen if path is legit, I think, and I do not know what to do so will leave it. */
        if( p == absPath ) break;                       /* Ditto. */
        for( needle += 3; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    } // Loop checking, 11.06.2015, T. Koi
    return( absPath );
}
/*
************************************************************
*/
int xDataTOM_setMessageError_ReturnInt( int value, statusMessageReporting *smr, void *userInterface, const char *packageName, int lineNumber, int code, 
    const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_setReportError( smr, userInterface, packageName, lineNumber, __func__, xDataTOM_smrLibraryID, code, fmt, args );
    va_end( args );
    return( value );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOM_getLinksElement( statusMessageReporting *smr, xDataTOM_element *element, char const *link ) {

    xDataTOM_element *linkedElement = NULL;

    if( link[0] == '/' ) {
        for( linkedElement = element; linkedElement->parent != NULL;  ) linkedElement = linkedElement->parent;
        linkedElement = xDataTOM_getLinksElement2( smr, linkedElement, &(link[1]), link ); }
    else {
        smr_setReportError2( smr, smr_unknownID, 1, "Only absolute link currently supported: requested link = '%s'", link );
    }
    return( linkedElement );
}
/*
************************************************************
*/
static xDataTOM_element *xDataTOM_getLinksElement2( statusMessageReporting *smr, xDataTOM_element *element, char const *link, char const *fullLink ) {

    int n = (int) strlen( link );
    char const *slash = strchr( link, '/' ), *bracket = strchr( link, '[' ), *attributesValue;
    char name[nameValueLength], value[nameValueLength];
    xDataTOM_element *child;

    if( bracket != NULL ) n = (int) ( bracket - link );
    if( slash != NULL ) {
        if( (int) ( slash - link ) < n ) {
            n = (int) ( slash - link );
            bracket = NULL;
        }
    }
    for( child = element->children; child != NULL; child = child->next ) {
        if( strncmp( link, child->name, n ) == 0 ) {
            if( bracket != NULL ) {
                if( bracket[1] != '@' ) {
                    smr_setReportError2( smr, smr_unknownID, 1, "bad link info at '%s' of '%s'", bracket, fullLink );
                    return( NULL );
                }
                if( xDataTOM_getLinksElement3( smr, &(bracket[2]), name, value, fullLink ) ) return( NULL );
                if( ( attributesValue = xDataTOM_getAttributesValueInElement( child, name ) ) == NULL ) continue;
                if( strcmp( value, attributesValue ) ) continue;
            }
            if( slash == NULL ) return( child );
            return( xDataTOM_getLinksElement2( smr, child, &(slash[1]), fullLink ) );
        }
    }
    return( NULL );
}
/*
************************************************************
*/
static int xDataTOM_getLinksElement3( statusMessageReporting *smr, char const *nameValue, char *name, char *value, char const *fullLink ) {

    int n;
    char const *equal = strchr( nameValue, '=' ), *p;
    char quote = '\'';

    if( equal == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "link qualifier missing '=' character at '%s' of '%s'", nameValue, fullLink );
        return( 1 );
    }
    n = (int) ( equal - nameValue );
    if( n >= ( nameValueLength - 1 ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "link's name qualifier too long at '%s' of '%s'", nameValue, fullLink );
        return( 1 );
    }
    strncpy( name, nameValue, n );
    name[n] = 0;

    equal++;
    if( *equal != quote ) quote = '"';
    if( *equal != quote ) {
        smr_setReportError2( smr, smr_unknownID, 1, "link's name qualifier missing quote at '%s' of '%s'", nameValue, fullLink );
        return( 1 );
    }

    equal++;
    p = strchr( equal, quote );
    if( p == NULL ) {
        smr_setReportError2( smr, smr_unknownID, 1, "link's name qualifier missing end quote at '%s' of '%s'", nameValue, fullLink );
        return( 1 );
    }

    n = (int) ( p - equal );
    if( n >= ( nameValueLength - 1 ) ) {
        smr_setReportError2( smr, smr_unknownID, 1, "link's value qualifier too long at '%s' of '%s'", nameValue, fullLink );
        return( 1 );
    }
    strncpy( value, equal, n );
    value[n] = 0;

    return( 0 );
}

#if defined __cplusplus
}
#endif
