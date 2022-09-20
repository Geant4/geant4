/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#if defined(WIN32) || defined(__MINGW32__)
#include <windows.h>
#define realpath( a, b ) GetFullPathName( a, PATH_MAX, b, NULL )
#define strtoll _strtoi64
#else
#include <unistd.h>
#endif

#include "xDataTOM_importXML_private.h"
#include "xDataTOM_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

int xDataTOM_smrLibraryID = smr_unknownID;

static int xDataTOME_initializeElement( statusMessageReporting *smr, xDataTOM_element *element, xDataTOM_element *parent, int ordinal, int index, 
    char const *name );
static void xDataTOME_displayTree( statusMessageReporting *smr, xDataTOM_element *element, int printAttributes, int level );

static int xDataTOM_initialize_xData( statusMessageReporting *smr, xDataTOM_xDataInfo *xDataInfo );
/*
************************************************************
*/
xDataTOM_TOM *xDataTOM_importFile( statusMessageReporting *smr, const char *fileName ) {
/*
*   Returns NULL is an error occurred.
*/

    return( xDataXML_importFile( smr, fileName ) );
}
/*
************************************************************
*/
xDataTOM_TOM *xDataTOM_mallocTOM( statusMessageReporting *smr ) {
/*
*   Returns NULL is an error occurred.
*/
    xDataTOM_TOM *TOM;

    if( ( TOM = (xDataTOM_TOM *) smr_malloc2( smr, sizeof( xDataTOM_TOM ), 1, "xDataTOM_TOM" ) ) != NULL ) {
        if( xDataTOM_initializeTOM( smr, TOM ) ) smr_freeMemory( (void **) &TOM );
    }
    return( TOM );
}
/*
************************************************************
*/
int xDataTOM_initializeTOM( statusMessageReporting *smr, xDataTOM_TOM *TOM ) {

    TOM->fileName = NULL;
    TOM->realFileName = NULL;
    xDataTOME_initializeElement( smr, &(TOM->root), NULL, 0, 0, "" );
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
void *xDataTOM_freeTOM( statusMessageReporting * /*smr*/, xDataTOM_TOM **TOM ) {

    xDataTOM_TOM *TOMp;

    if( TOM == NULL ) return( NULL );
    TOMp = *TOM;
    if( TOMp == NULL ) return( NULL );
    xDataTOM_releaseElement( &(TOMp->root) );
    smr_freeMemory( (void **) &(TOMp->fileName) );
    smr_freeMemory( (void **) &(TOMp->realFileName) );
    smr_freeMemory( (void **) TOM );
    return( NULL );
}
/*
************************************************************
*/
int xDataTOM_setFileNameTOM( statusMessageReporting *smr, xDataTOM_TOM *TOM, const char *fileName ) {
/*
*   Returns not zero value if error occurred.
*/

    char realPath[PATH_MAX+1];

    smr_freeMemory( (void **) &(TOM->fileName) );
    smr_freeMemory( (void **) &(TOM->realFileName) );
        if( fileName != NULL ) {
        if( ( TOM->fileName = smr_allocateCopyString2( smr, fileName, "fileName" ) ) == NULL ) return( 1 );
        if( realpath( fileName, realPath ) != NULL ) {
            if( ( TOM->realFileName = smr_allocateCopyString2( smr, realPath, "realFileName" ) ) == NULL ) return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
void xDataTOM_displayTree( statusMessageReporting *smr, xDataTOM_TOM *TOM, int printAttributes ) {

    if( TOM->root.children != NULL ) xDataTOME_displayTree( smr, TOM->root.children, printAttributes, 0 );
}

/****************************************
*       Element functions.
****************************************/
/*
************************************************************
*/
xDataTOM_element *xDataTOM_mallocElement( statusMessageReporting *smr, xDataTOM_element *parent, int ordinal, int index, char const *name ) {
/*
*   Returns NULL is an error occurred.
*/
    xDataTOM_element *element;

    if( ( element = (xDataTOM_element *) smr_malloc2( smr, sizeof( xDataTOM_element ), 1, "xDataTOM_elelument" ) ) != NULL ) {
        if( xDataTOME_initializeElement( smr, element, parent, ordinal, index, name ) ) smr_freeMemory( (void **) &element );
    }
    return( element );
}
/*
************************************************************
*/
void xDataTOM_freeElement( xDataTOM_element **element ) {

    if( element == NULL ) return;
    xDataTOM_releaseElement( *element );
    smr_freeMemory( (void **) element );
}
/*
************************************************************
*/
void xDataTOM_releaseElement( xDataTOM_element *element ) {

    xDataTOM_element *child, *nextChild;

    if( element == NULL ) return;
    xDataTOMAL_release( &(element->attributes) );
    for( child = element->children; child != NULL; child = nextChild ) {
        nextChild = child->next;
        xDataTOM_freeElement( &child );
    }
    if( element->xDataInfo.ID != NULL ) {
        xDataTOM_axes_release( &(element->xDataInfo.axes) );
        if( strcmp( element->xDataInfo.ID, xDataTOM_XYs_ID ) == 0 ) {
            xDataTOM_XYs_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_regionsXYs_ID ) == 0 ) {
            xDataTOM_regionsXYs_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_W_XYs_ID ) == 0 ) {
            xDataTOM_W_XYs_freeFrom_xDataInfo( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_V_W_XYs_ID ) == 0 ) {
            xDataTOM_V_W_XYs_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_W_XYs_LegendreSeries_ID ) == 0 ) {
            xDataTOM_W_XYs_LegendreSeries_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_regionsW_XYs_LegendreSeries_ID ) == 0 ) {
            xDataTOM_regionsW_XYs_LegendreSeries_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_V_W_XYs_LegendreSeries_ID ) == 0 ) {
            xDataTOM_V_W_XYs_LegendreSeries_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_KalbachMann_ID ) == 0 ) {
            xDataTOM_KalbachMann_free( &(element->xDataInfo ) ); }
        else if( strcmp( element->xDataInfo.ID, xDataTOM_polynomial_ID ) == 0 ) {
            xDataTOM_polynomial_free( &(element->xDataInfo ) ); }
        else {
            printf( "not freed for %s\n", element->xDataInfo.ID );
        }
    }
    element->parent = NULL;
    smr_freeMemory( (void **) &(element->name) );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOM_addElementInElement( statusMessageReporting *smr, xDataTOM_element *parent, int index, char const *name ) {

    xDataTOM_element *element;

    if( ( element = xDataTOM_mallocElement( smr, parent, parent->numberOfChildren, index, name ) ) == NULL ) return( NULL ); 
    if( parent->children == NULL ) {
        parent->children = element; }
    else {
        xDataTOM_element *last;

        for( last = parent->children; last->next != NULL; last = last->next ) ;
        last->next = element;
    }
    (parent->numberOfChildren)++;
    return( element );
}
/*
************************************************************
*/
static int xDataTOME_initializeElement( statusMessageReporting *smr, xDataTOM_element *element, xDataTOM_element *parent, int ordinal, int index, 
        char const *name ) {

    element->ordinal = ordinal;
    element->index = index;
    element->parent = parent;
    element->next = NULL;
    element->name = smr_allocateCopyString2( smr, name, "element->name" );
    xDataTOMAL_initial( smr, &(element->attributes) );
    element->numberOfChildren = 0;
    element->children = NULL;
    return( ( xDataTOM_initialize_xData( smr, &(element->xDataInfo) ) || ( element->name == NULL ) ) ? 1 : 0 );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOM_getDocumentsElement( xDataTOM_TOM *TOM ) {

    return( TOM->root.children );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOME_getFirstElement( xDataTOM_element *element ) {

    if( element != NULL ) element = element->children;
    return( element );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOME_getNextElement( xDataTOM_element *element ) {

    if( element != NULL ) element = element->next;
    return( element );
}
/*
************************************************************
*/
xDataTOM_element *xDataTOME_getOneElementByName( statusMessageReporting *smr, xDataTOM_element *element, char const *name, int required ) {

    int n = 0;
    xDataTOM_element *child, *desired = NULL;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        if( strcmp( child->name, name ) == 0 ) {
            if( n == 0 ) desired = child;
            n++;
        }
    }
    if( n == 0 ) {
        if( required ) smr_setReportError2( smr, smr_unknownID, 1, "elements '%s' not found in element '%s'", name, element->name ); }
    else if( n > 1 ) {
        smr_setReportError2( smr, smr_unknownID, 1, "multiple (= %d) elements '%s' found in element '%s'", name, element->name );
        desired = NULL;
    }
    return( desired );
}
/*
************************************************************
*/
int xDataTOM_numberOfElementsByName( statusMessageReporting * /*smr*/, xDataTOM_element *element, char const *name ) {

    int n = 0;
    xDataTOM_element *child;

    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) if( !strcmp( child->name, name ) ) n++;
    return( n );
}
/*
************************************************************
*/
int xDataTOME_addAttribute( statusMessageReporting *smr, xDataTOM_element *element, char const *name, char const *value ) {

    return( xDataTOMAL_addAttribute( smr, &(element->attributes), name, value ) );
}
/*
************************************************************
*/
char const *xDataTOM_getAttributesValueInElement( xDataTOM_element *element, char const *name ) {

    return( xDataTOMAL_getAttributesValue( &(element->attributes), name ) );
}
/*
************************************************************
*/
int xDataTOME_copyAttributionList( statusMessageReporting *smr, xDataTOM_attributionList *desc, xDataTOM_element *element ) {

    return( xDataTOMAL_copyAttributionList( smr, desc, &(element->attributes) ) );
}
/*
************************************************************
*/
int xDataTOME_convertAttributeToInteger( statusMessageReporting *smr, xDataTOM_element *element, char const *name, int *n ) {

    return( xDataTOMAL_convertAttributeToInteger( smr, &(element->attributes), name, n ) );
}
/*
************************************************************
*/
int xDataTOME_convertAttributeToDouble( statusMessageReporting *smr, xDataTOM_element *element, char const *name, double *d ) {

    return( xDataTOMAL_convertAttributeToDouble( smr, &(element->attributes), name, d ) );
}
/*
************************************************************
*/
int xDataTOME_getInterpolation( statusMessageReporting *smr, xDataTOM_element *element, int index, enum xDataTOM_interpolationFlag *independent, 
    enum xDataTOM_interpolationFlag *dependent, enum xDataTOM_interpolationQualifier *qualifier ) {

    xDataTOM_xDataInfo *xDI = &(element->xDataInfo);

    if( xDI->ID == NULL ) return( 1 );

    return( xDataTOM_axes_getInterpolation( smr, &(xDI->axes), index, independent, dependent, qualifier ) );
}
/*
************************************************************
*/
static void xDataTOME_displayTree( statusMessageReporting *smr, xDataTOM_element *element, int printAttributes, int level ) {

    int i;
    xDataTOM_element *child;

    for( i = 0; i < level; i++ ) printf( "    " );
    printf( "/%s", element->name );
    if( element->index >= 0 ) printf( " (%d)", element->index );
    if( printAttributes ) {
        xDataTOM_attribute *attribute;

        for( attribute = element->attributes.attributes; attribute != NULL; attribute = attribute->next ) {
            printf( " (%s, \"%s\")", attribute->name, attribute->value );
        }
    }
    printf( "\n" );
    for( child = xDataTOME_getFirstElement( element ); child != NULL; child = xDataTOME_getNextElement( child ) ) {
        xDataTOME_displayTree( smr, child, printAttributes, level + 1 );
    }
}

/****************************************
*       Attribute functions.
****************************************/
/*
************************************************************
*/
void xDataTOMAL_initial( statusMessageReporting * /*smr*/, xDataTOM_attributionList *attributes ) {

    attributes->number = 0;
    attributes->attributes = NULL;
}
/*
************************************************************
*/
void xDataTOMAL_release( xDataTOM_attributionList *attributes ) {

    xDataTOM_attribute *attribute, *next;

    for( attribute = attributes->attributes; attribute != NULL; attribute = next ) {
        next = attribute->next;
        smr_freeMemory( (void **) &(attribute->name) );
        smr_freeMemory( (void **) &(attribute->value) );
        smr_freeMemory( (void **) &(attribute) );
    }
    xDataTOMAL_initial( NULL, attributes );
}
/*
************************************************************
*/
int xDataTOMAL_addAttribute( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, char const *value ) {

    xDataTOM_attribute *attribute;

    if( ( attribute = (xDataTOM_attribute *) smr_malloc2( smr, sizeof( xDataTOM_attribute ), 1, "xDataTOM_attribute" ) ) == NULL ) return( 1 );
    if( ( attribute->name = smr_allocateCopyString2( smr, name, "name" ) ) == NULL ) goto err;
    if( ( attribute->value = smr_allocateCopyString2( smr, value, "value" ) ) == NULL ) goto err;
    if( attributes->attributes == NULL ) {
        attributes->attributes = attribute; }
    else {
        xDataTOM_attribute *last;

        for( last = attributes->attributes; last->next != NULL; last = last->next ) ;
        last->next = attribute;
    }
    attributes->number++;
    return( 0 );

err:
    smr_freeMemory( (void **) &(attribute->name) );
    smr_freeMemory( (void **) &(attribute->value) );
    smr_freeMemory( (void **) &(attribute) );
    return( 1 );
}
/*
************************************************************
*/
char const *xDataTOMAL_getAttributesValue( xDataTOM_attributionList *attributes, char const *name ) {

    xDataTOM_attribute *attribute;

    for( attribute = attributes->attributes; attribute != NULL; attribute = attribute->next ) {
        if( !strcmp( attribute->name, name ) ) return( attribute->value );
    }
    return( NULL );
}
/*
************************************************************
*/
int xDataTOMAL_copyAttributionList( statusMessageReporting *smr, xDataTOM_attributionList *desc, xDataTOM_attributionList *src ) {

    xDataTOM_attribute *attribute;

    xDataTOMAL_initial( smr, desc );
    for( attribute = src->attributes; attribute != NULL; attribute = attribute->next ) {
        if( xDataTOMAL_addAttribute( smr, desc, attribute->name, attribute->value ) != 0 ) goto err;

    }
    return( 0 );

err:
    xDataTOMAL_release( desc );
    return( 1 );
}
/*
************************************************************
*/
int xDataTOMAL_convertAttributeToInteger( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, int *n ) {

    char const *value = xDataTOMAL_getAttributesValue( attributes, name );
    char *e;

    if( value != NULL ) {
        *n = (int) strtoll( value, &e, 10 );
        if( *e == 0 ) return( 0 );
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "could not convert attribute %s's value = '%s' to an integer", name, value ); }
    else {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "no attribute named '%s'", name );
    }
    return( 1 );
}
/*
************************************************************
*/
int xDataTOMAL_convertAttributeToDouble( statusMessageReporting *smr, xDataTOM_attributionList *attributes, char const *name, double *d ) {

    char const *value = xDataTOMAL_getAttributesValue( attributes, name );
    char *e;

    if( value != NULL ) {
        *d = strtod( value, &e );
        if( *e == 0 ) return( 0 );
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "could not convert attribute %s's values = '%s' to a double", name, value ); }
    else {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "no attribute named '%s'", name );
    }
    return( 1 );
}


/****************************************
*       xData functions.
****************************************/
/*
************************************************************
*/
static int xDataTOM_initialize_xData( statusMessageReporting * /*smr*/, xDataTOM_xDataInfo * /*xDataInfo*/ ) {

    return( 0 );
}
/*
************************************************************
*/
void *xData_initializeData( statusMessageReporting *smr, xDataTOM_element *TE, char const *ID, size_t size ) {

    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);

    xDI->data = NULL;
    xDI->ID = ID;
    xDI->element = TE;
    return( xDI->data = (void *) smr_malloc2( smr, size, 1, "xDI->data" ) );
}
/*
************************************************************
*/
int xDataTOM_isXDataID( xDataTOM_element *TE, char const *ID ) {

    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);

    if( xDI->ID != NULL ) {
        return( !strcmp( xDI->ID, ID ) );
    }

    return( 0 );
}
/*
************************************************************
*/
xDataTOM_xDataInfo *xDataTOME_getXData( xDataTOM_element *TE ) {

    if( TE->xDataInfo.ID == NULL ) return( NULL );
    return( &(TE->xDataInfo) );
}
/*
************************************************************
*/
void *xDataTOME_getXDataIfID( statusMessageReporting *smr, xDataTOM_element *TE, char const *ID ) {

    xDataTOM_xDataInfo *xDI = xDataTOME_getXData( TE );

    if( xDI == NULL ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "element '%s' does not have xData", TE->name );
        return( NULL );
    }
    if( strcmp( ID, xDI->ID ) ) {
        smr_setReportError2( smr, xDataTOM_smrLibraryID, 1, "xData has ID = '%s' not '%s' for element %s", xDI->ID, ID, TE->name );
        return( NULL );
    }
    return( xDI->data );

}

#if defined __cplusplus
}
#endif
