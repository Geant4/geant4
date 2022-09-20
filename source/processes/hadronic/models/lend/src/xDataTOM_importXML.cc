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
#include <BaseTsd.h>
#include <io.h>
#include <windows.h>
#define realpath( a, b ) GetFullPathName( a, PATH_MAX, b, NULL )
#define strtoll _strtoi64
typedef SSIZE_T ssize_t;
#else
#include <unistd.h>
#endif

#include "xDataTOM_importXML_private.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

static xDataXML_document *xDataXML_mallocDoc( statusMessageReporting *smr );
static int xDataXML_initializeDoc( statusMessageReporting *smr, xDataXML_document *doc );
static int xDataXML_endXMLParsing( statusMessageReporting *smr, xDataXML_document *doc );
static void *xDataXML_freeElement( statusMessageReporting *smr, xDataXML_element *element );
static void xDataXML_freeElementItems( statusMessageReporting *smr, xDataXML_element *element );
static int xDataXML_parse( xDataXML_document *doc, char const *s );
static void XMLCALL xDataXML_parseStartElement( void *userData, char const *name, char const **attris );
static void XMLCALL xDataXML_parseEndElement( void *userData, char const *name );
static void XMLCALL xDataXML_parseCharacterData( void *userData, XML_Char const *s, int len );
static void xDataXML_initializeRootElement( xDataXML_document *doc, xDataXML_rootElement *re, xDataXML_element *parentElement, int depth );
static int xDataXML_parseInitializeText( xDataXML_document *doc, xDataXML_text *text );
static int xDataXML_addElementToRoot( statusMessageReporting *smr, xDataXML_rootElement *parentRoot, char const *name, char const **attris );
static enum xDataXML_errorCodes xDataXML_parseGetCurrentPosition( xDataXML_document *doc, xDataXML_docInfo *docInfo );
static int xDataXML_init_xDataTypeNone( xDataXMLType *xDT, xDataXML_element *element );
static char *xDataXML_getTraceback( statusMessageReporting *smr, xDataXML_element *element );
static char *xDataXML_getTraceback2( statusMessageReporting *smr, xDataXML_rootElement *parentRoot, int n );
static int xDataXML_setFileName( statusMessageReporting *smr, xDataXML_document *doc, char const *fileName );

static int xDataXML_smrUserInterfaceInitialize( xDataXML_document *doc );
static int xDataXML_smrUserInterfaceFree( xDataXML_document *doc );
static char *xDataXML_smrUserInterface( void *userData );
static char const *xDataXML_shortStringForMessage( size_t size, char *Out, char const *In );

static int xDataXML_constructTOM( statusMessageReporting *smr, xDataTOM_element *TE, xDataXML_element *element );
/*
************************************************************
*/
xDataTOM_TOM *xDataXML_importFile( statusMessageReporting *smr, char const *fileName ) {
/*
*   Returns NULL is any error occurred. If an error occurs in an expat routine, xDataXML_endXMLParsing will set smr appropriately.
*/
    xDataTOM_TOM *TOM = NULL;
    xDataXML_document *XML = NULL;
    xDataXML_element *element;

    if( ( XML = xDataXML_importFile2( smr, fileName ) ) == NULL ) return( NULL );

    if( ( TOM = xDataTOM_mallocTOM( smr ) ) == NULL ) goto Err;
    if( xDataTOM_setFileNameTOM( smr, TOM, fileName ) != 0 ) goto Err;

    element = xDataXML_getDocumentsElement( XML );
    if( xDataXML_constructTOM( smr, (&TOM->root), element ) != 0 ) goto Err;

    xDataXML_freeDoc( smr, XML );
    return( TOM );

Err:
    if( XML != NULL ) xDataXML_freeDoc( smr, XML );
    if( TOM != NULL ) xDataTOM_freeTOM( smr, &TOM );
    return( NULL );
}
/*
************************************************************
*/
xDataXML_document *xDataXML_importFile2( statusMessageReporting *smr, char const *fileName ) {
/*
*   Returns NULL is any error occurred. If an error occurs in an expat routine, xDataXML_endXMLParsing will set smr appropriately.
*/
    int f;
    char buffer[10 * 1000];
    ssize_t count, n = sizeof( buffer ) - 1;
    xDataXML_document *doc;

    if( ( doc = xDataXML_mallocDoc( smr ) ) == NULL ) return( NULL );
    if( xDataXML_setFileName( smr, doc, fileName ) == 0 ) {
        f = open( fileName, O_RDONLY );
        if( f == -1 ) {
                xDataXML_endXMLParsing( smr, doc );
                smr_setReportError2( smr, xDataTOM_smrLibraryID, xDataXML_errFileError, "could not open XML file %s", fileName ); }
        else {
                while( ( count = read( f, buffer, n ) ) > 0 ) {
                    buffer[count] = 0;
                    if( xDataXML_parse( doc, buffer ) ) break;
                    if( !smr_isOk( doc->smr ) ) break;
                }  // Loop checking, 11.06.2015, T. Koi
                close( f );
                xDataXML_endXMLParsing( smr, doc );
                if( count < 0 ) smr_setReportError2( smr, xDataTOM_smrLibraryID, xDataXML_errFileError, "read failed with errno = %d for XML %s", 
                    errno, fileName );
        }
    }
    if( doc != NULL ) {
        if( !smr_isOk( smr ) ) {
            xDataXML_freeDoc( smr, doc );
            doc = NULL;
        }
    }
    return( doc );
}
/*
************************************************************
*/
static xDataXML_document *xDataXML_mallocDoc( statusMessageReporting *smr ) {

    xDataXML_document *doc;

    if( ( doc = (xDataXML_document *) smr_malloc2( smr, sizeof( xDataXML_document ), 1, "xDataXML_document" ) ) != NULL ) {
        if( xDataXML_initializeDoc( smr, doc ) ) doc = (xDataXML_document *) xDataXML_freeDoc( smr, doc );
    }
    return( doc );
}
/*
************************************************************
*/
static int xDataXML_initializeDoc( statusMessageReporting *smr, xDataXML_document *doc ) {

    doc->status = xDataXML_statusParsing;
    doc->error = xDataXML_errNone;
    doc->err = XML_ERROR_NONE;
    doc->err_line = 0;
    doc->err_column = 0;
    doc->fileName = NULL;
    doc->realFileName = NULL;
    xDataXML_smrUserInterfaceInitialize( doc );
    doc->smr= smr;
    if( ( doc->xmlParser = XML_ParserCreate( NULL ) ) == NULL ) {
        smr_setReportError2p( smr, xDataTOM_smrLibraryID, xDataXML_errXML_ParserCreate, "XML_ParserCreate failed" ); }
    else {
        XML_SetUserData( doc->xmlParser, doc  );
        xDataXML_initializeRootElement( doc, &(doc->root), NULL, 0 );
        doc->currentRoot = &(doc->root);
        XML_SetElementHandler( doc->xmlParser, xDataXML_parseStartElement, xDataXML_parseEndElement );
        XML_SetCharacterDataHandler( doc->xmlParser, xDataXML_parseCharacterData );
    }
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
static int xDataXML_endXMLParsing( statusMessageReporting *smr, xDataXML_document *doc ) {

    if( doc->xmlParser ) {
        doc->err = XML_GetErrorCode( doc->xmlParser );
        doc->err_line = XML_GetCurrentLineNumber( doc->xmlParser );
        doc->err_column = XML_GetCurrentColumnNumber( doc->xmlParser );
        if( smr_isOk( smr ) && ( XML_Parse( doc->xmlParser, NULL, 0, 1 ) == XML_STATUS_ERROR ) ) {
            doc->status = xDataXML_statusError;
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromDocument( doc ), xDataTOM_smrLibraryID, xDataXML_errXMLParser, 
                "status = %d\nXML_Error code = %d\nXML_ErrorString = %s\nerror line, column = %d, %d", xDataXML_errXMLParser, 
                doc->err, XML_ErrorString( doc->err ), doc->err_line, doc->err_column );
        }
        XML_ParserFree( doc->xmlParser );
        doc->xmlParser = NULL;
        if( doc->status != xDataXML_statusError ) doc->status = xDataXML_statusCompleted;
    }
    return( 0 );
}
/*
************************************************************
*/
void *xDataXML_freeDoc( statusMessageReporting *smr, xDataXML_document *doc ) {

    xDataXML_endXMLParsing( smr, doc );
    doc->root.children = (xDataXML_element *) xDataXML_freeElement( smr, doc->root.children );
    smr_freeMemory( (void **) &(doc->fileName) );
    smr_freeMemory( (void **) &(doc->realFileName) );
    xDataXML_smrUserInterfaceFree( doc );
    smr_freeMemory( (void **) &doc );
    return( NULL );
}
/*
************************************************************
*/
static void *xDataXML_freeElement( statusMessageReporting *smr, xDataXML_element *element ) {
    
    xDataXML_element *next;

    for( ; element != NULL; element = next ) {
        next = element->next;
        xDataXML_freeElementItems( smr, element );
        smr_freeMemory( (void **) &element );
    }
    return( NULL );
}
/*
************************************************************
*/
static void xDataXML_freeElementItems( statusMessageReporting *smr, xDataXML_element *element ) {

    element->childrenRoot.children = (xDataXML_element *) xDataXML_freeElement( smr, element->childrenRoot.children );
/* BRB, The next line needs work */
    if( ( !strcmp( element->name, "xData" ) ) && ( element->xDataTypeInfo.release != NULL ) ) element->xDataTypeInfo.release( smr, &(element->xDataTypeInfo) );
    smr_freeMemory( (void **) &(element->name) );
    smr_freeMemory( (void **) &(element->fullName) );
    if( element->attributes.attributes ) smr_freeMemory( (void **) &(element->attributes.attributes) );
    if( element->text.text ) smr_freeMemory( (void **) &(element->text.text) );
}
/*
************************************************************
*/
static int xDataXML_parse( xDataXML_document *doc, char const *s ) {

    if( doc->status != xDataXML_statusParsing ) return( doc->status );
    if( XML_Parse( doc->xmlParser, s, (int) strlen( s ), 0 ) == XML_STATUS_ERROR ) return( -1 );
    return( 0 );
}
/*
************************************************************
*/
static void XMLCALL xDataXML_parseStartElement( void *userData, char const *name, char const **attris ) {

    xDataXML_document *doc = (xDataXML_document *) userData;

    if( !smr_isOk( doc->smr ) ) return;
    xDataXML_addElementToRoot( doc->smr, doc->currentRoot, name, attris );
}
/*
************************************************************
*/
static void XMLCALL xDataXML_parseEndElement( void *userData, char const * /*name*/ ) {

    xDataXML_document *doc = (xDataXML_document *) userData;

    doc->currentRoot->currentChild = NULL;
    doc->currentRoot = doc->currentRoot->parentRoot;
}
/*
************************************************************
*/
static void XMLCALL xDataXML_parseCharacterData( void *userData, XML_Char const *s, int len ) {
/*
*   Always terminates text with a 0.
*/
    xDataXML_document *doc = (xDataXML_document *) userData;
    xDataXML_text *text = &(doc->currentRoot->parentRoot->currentChild->text);
    size_t needSize = text->length + len + 1, l; 
    char *p;

    if( !smr_isOk( doc->smr ) ) return;
    if( needSize < 8  ) needSize = 8;
    if( needSize > text->allocated ) {
        if( text->allocated != 0 ) {
            l = ( 20 * text->allocated ) / 100;
            if( l < 100 ) l = 100;
            if( needSize < ( text->allocated + l ) ) needSize = text->allocated + l;
        }
        text->allocated = needSize;
        text->text = (char *) smr_realloc2( doc->smr, text->text, text->allocated, "text" );
        if( !smr_isOk( doc->smr ) ) return;
    }
    p = &(text->text[text->length]);
    strncpy( p, s, len );
    text->length += len;
    p[len] = 0;
}
/*
************************************************************
*/
static void xDataXML_initializeRootElement( xDataXML_document *doc, xDataXML_rootElement *re, xDataXML_element *parentElement, int depth ) {

    re->xData_doc = doc;
    re->parentElement = parentElement;
    re->parentRoot = NULL;
    if( parentElement != NULL ) re->parentRoot = parentElement->parentRoot;
    re->depth = depth;
    re->numberOfElements = 0;
    re->children = NULL;
    re->currentChild = NULL;
}
/*
************************************************************
*/
static int xDataXML_parseInitializeText( xDataXML_document *doc, xDataXML_text *text ) {

    xDataXML_parseGetCurrentPosition( doc, &(text->docInfo) );
    text->allocated = 0;
    text->length = 0;
    text->text = NULL;
    return( 0 );
}
/*
************************************************************
*/
static int xDataXML_addElementToRoot( statusMessageReporting *smr, xDataXML_rootElement *parentRoot, char const *name, char const **attris ) {

    xDataXML_document *doc = parentRoot->xData_doc;
    xDataXML_element *element;
    int i, n, status = 1;
    size_t lens;
    char *p, *e;
    char const **pAttris;
    xDataXML_attribute *a;
    void *smrUser;

    element = (xDataXML_element *) smr_malloc2( doc->smr, sizeof( xDataXML_element ), 1, "xDataXML_element" );
    if( element == NULL ) return( 1 );
    xDataXML_parseGetCurrentPosition( doc, &(element->docInfo) );
    element->ordinal = parentRoot->numberOfElements;
    element->index = -1;
    element->accessed = 0;
    element->parentRoot = parentRoot;
    xDataXML_initializeRootElement( doc, &(element->childrenRoot), element, parentRoot->depth + 1 );
    element->next = NULL;
    if( ( element->name = (char *) smr_malloc2( doc->smr, strlen( name ) + 1, 0, "name" ) ) == NULL ) {
        smr_freeMemory( (void **) &element );
        return( 1 );
    }
    strcpy( element->name, name );
    if( ( element->fullName = xDataXML_getTraceback( smr, element ) ) == NULL ) {
        smr_freeMemory( (void **) &(element->name) );
        smr_freeMemory( (void **) &element );
        return( 1 );
    }
    for( i = 0, lens = 0, pAttris = attris; *pAttris; i++, pAttris++ ) lens += strlen( *pAttris ) + 1;
    n = i / 2;
    element->attributes.size = n * sizeof( xDataXML_attribute ) + lens;
    element->attributes.number = n;
    element->attributes.attributes = NULL;
    smrUser = xDataXML_get_smrUserInterfaceFromElement( element );
    if( element->attributes.size  ) {
        if( ( element->attributes.attributes = (xDataXML_attribute *) smr_malloc2( doc->smr, element->attributes.size, 0, "attributes") ) == NULL ) {
            status = 0; }
        else {
            a = element->attributes.attributes;
            p = (char *) &(element->attributes.attributes[n]);
            for( i = 0, pAttris = attris; ( i < n ) && status; i++, a++, pAttris++ ) {
                lens = strlen( *pAttris ) + 1;
                a->name = p;
                strcpy( p, *pAttris );
                p += lens;
                pAttris++;
                lens = strlen( *pAttris ) + 1;
                a->value= p;
                strcpy( p, *pAttris );
                p += lens;
                if( !strcmp( "index", a->name ) ) {
                    element->index = (int) strtoll( a->value, &e, 10 );
                    if( *e != 0 ) {
                        status = 0;
                        smr_setReportError3( doc->smr, smrUser, xDataTOM_smrLibraryID, -1, "could not convert index attribute = %s to integer", a->value );
                    }
                }
            }
        }
    }
    if( !status ) {
        smr_freeMemory( (void **) &(element->attributes.attributes) );
        smr_freeMemory( (void **) &(element->name) );
        smr_freeMemory( (void **) &(element->fullName) );
        smr_freeMemory( (void **) &element );
        return( 1 );
    }
    xDataXML_init_xDataTypeNone( &(element->xDataTypeInfo), element );
    element->textOffset = 0;
    xDataXML_parseInitializeText( doc, &(element->text) );
    if( parentRoot->parentElement != NULL ) element->textOffset = parentRoot->parentElement->text.length;
    if( parentRoot->currentChild == NULL ) {
        parentRoot->children = element; }
    else {
        parentRoot->currentChild->next = element;
    }
    parentRoot->numberOfElements++;
    parentRoot->currentChild = element;
    doc->currentRoot = &(element->childrenRoot);
    return( 0 );
}
/*
************************************************************
*/
static enum xDataXML_errorCodes xDataXML_parseGetCurrentPosition( xDataXML_document *doc, xDataXML_docInfo *docInfo ) {

    docInfo->column = XML_GetCurrentColumnNumber( doc->xmlParser );
    docInfo->line = XML_GetCurrentLineNumber( doc->xmlParser );
    return( xDataXML_errNone );
}
/*
************************************************************
*/
int xDataXML_parseIsError( xDataXML_document *doc ) {

    return( doc->status == xDataXML_statusError );
}
/*
************************************************************
*/
xDataXML_element *xDataXML_getDocumentsElement( xDataXML_document *doc ) { return( doc->root.children ); }
xDataXML_element *xDataXML_getFirstElement( xDataXML_element *element ) { return( element->childrenRoot.children ); }
xDataXML_element *xDataXML_getNextElement( xDataXML_element *element ) { return( element->next ); }
/*
************************************************************
*/
enum xDataXML_itemMode xDataXML_getFirstItem( xDataXML_element *element, xDataXML_item *item ) {

    item->parentElement = element;
    item->element = xDataXML_getFirstElement( element );
    if( item->element == NULL ) {
        item->mode = xDataXML_itemModeText;
        if( element->text.length == 0 ) item->mode = xDataXML_itemModeEnd; }
    else {
        item->mode = xDataXML_itemModeElement;
        if( 0 < item->element->textOffset ) item->mode = xDataXML_itemModeText;
    }
    item->textOffset = 0;
    item->textLength = element->text.length;
    if( item->element != NULL ) item->textLength = item->element->textOffset;
    item->text = element->text.text;
    return( item->mode );
}
/*
************************************************************
*/
enum xDataXML_itemMode xDataXML_getNextItem( xDataXML_item *item ) {

    if( item->mode != xDataXML_itemModeEnd ) {
        if( item->mode == xDataXML_itemModeText ) {
            item->mode = xDataXML_itemModeElement;
            if( item->element == NULL ) item->mode = xDataXML_itemModeEnd;
            item->textOffset += item->textLength;
            item->textLength = 0;
            item->text = &(item->parentElement->text.text[item->textOffset]); }
        else {
            item->element = item->element->next;
            item->mode = xDataXML_itemModeText;
            if( item->element == NULL ) {
                if( item->textOffset < item->parentElement->text.length ) {
                    item->textLength = item->parentElement->text.length - item->textOffset; }
                else {
                    item->mode = xDataXML_itemModeEnd;
                } }
            else {
                item->textLength = item->element->textOffset - item->textOffset;
            }
        }
    }
    return( item->mode );
}
/*
************************************************************
*/
int xDataXML_isAttributeInList( xDataXML_attributionList *attributes, char const *name ) {

    int i;

    for( i = 0; i < attributes->number; i++ ) {
        if( !strcmp( attributes->attributes[i].name, name ) ) return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int xDataXML_isAttributeInElement( xDataXML_element *element, char const *name ) {

    return( xDataXML_isAttributeInList( &(element->attributes), name ) );
}
/*
************************************************************
*/
char *xDataXML_getAttributesValue( xDataXML_attributionList *attributes, char const *name ) {

    int i;

    for( i = 0; i < attributes->number; i++ ) {
        if( !strcmp( attributes->attributes[i].name, name ) ) return( attributes->attributes[i].value );
    }
    return( NULL );
}
/*
************************************************************
*/
char const *xDataXML_getAttributesValueInElement( xDataXML_element *element, char const *name ) {

    return( (char const *) xDataXML_getAttributesValue( &(element->attributes), name ) );
}
/*
************************************************************
*/
int xDataXML_attributeListLength( xDataXML_attributionList *attributes ) {

    return( attributes->number );
}
/*
************************************************************
*/
xDataXML_attribute *xDataXML_attributeByIndex( xDataXML_attributionList *attributes, int index ) {

    if( index >= attributes->number ) return( NULL );
    return( &(attributes->attributes[index]) );
}
/*
************************************************************
*/
static int xDataXML_init_xDataTypeNone( xDataXMLType *xDT, xDataXML_element *element ) {

    xDT->status = xDataXML_xDataType_Ok;
    xDT->ID = NULL;
    xDT->element = element;
    xDT->toData = NULL;
    xDT->toString = NULL;
    xDT->release = NULL;
    xDT->indexPresent = 1;                  /* The following describes the meaning of present variables. */
    xDT->startPresent = 1;                  /* If < 0, an error occured in converting value to an integer. */
    xDT->endPresent = 1;                    /* If > 0, not present as an attribute. */
    xDT->lengthPresent = 1;                 /* Else, if 0, present and converted without an error. */
    xDT->index = -1;
    xDT->start = -1;
    xDT->end = -1;
    xDT->length = -1;
    xDT->data = NULL;
    return( 0 );
}
/*
************************************************************
*/
int xDataXML_getCommonData( statusMessageReporting *smr, xDataXML_element *element, xDataTOM_Int *index, xDataTOM_Int *start, xDataTOM_Int *end,
        xDataTOM_Int *length ) {

    if( element->xDataTypeInfo.ID == NULL ) {
        smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 1, 
            "element %s is not xData", element->fullName );
        return( 1 );
    }
    *index = element->xDataTypeInfo.index;
    *start = element->xDataTypeInfo.start;
    *end = element->xDataTypeInfo.end;
    *length = element->xDataTypeInfo.length;
    return( 0 );
}
/*
************************************************************
*/
int xDataXML_xDataTypeConvertAttributes( statusMessageReporting *smr, xDataXML_element *element ) {

    xDataXMLType *xDT = &(element->xDataTypeInfo);
    void *smrUser = xDataXML_get_smrUserInterfaceFromElement( element );

    xDT->index = -1;
    xDT->start = -1;
    xDT->end = -1;
    xDT->length = -1;
    if( ( xDT->indexPresent = xDataXML_convertAttributeTo_xDataTOM_Int( smr, element, "index", &(xDT->index), 0 ) ) < 0 ) return( 1 );
    if( ( xDT->startPresent = xDataXML_convertAttributeTo_xDataTOM_Int( smr, element, "start", &(xDT->start), 0 ) ) < 0 ) return( 1 );
    if( ( xDT->endPresent = xDataXML_convertAttributeTo_xDataTOM_Int( smr, element, "end", &(xDT->end), 0 ) ) < 0 ) return( 1 );
    if( ( xDT->lengthPresent = xDataXML_convertAttributeTo_xDataTOM_Int( smr, element, "length", &(xDT->length), 0 ) ) < 0 ) return( 1 );
    if( ( xDT->endPresent > 0 ) ) {
        if( xDT->lengthPresent > 0 ) {
            smr_setReportError3p( smr, smrUser, xDataTOM_smrLibraryID, 1, "missing length (or end) in xData" );
            return( 1 );
        }
        xDT->end = xDT->length; }
    else {
        if( xDT->lengthPresent > 0 ) xDT->length = xDT->end;
    }

    if( xDT->startPresent > 0 ) xDT->start = 0;
    if( xDT->start < 0 ) {
        smr_setReportError3( smr, smrUser, xDataTOM_smrLibraryID, 1, "start = %d < 0", xDT->start );
        return( 1 );
    }
    if( xDT->end < xDT->start ) {
        smr_setReportError3( smr, smrUser, xDataTOM_smrLibraryID, 1, "start = %d >= end = %d", xDT->start, xDT->end );
        return( 1 );
    }
    if( xDT->length < 0 ) {
        smr_setReportError3( smr, smrUser, xDataTOM_smrLibraryID, 1, "length = %d < 0", xDT->length );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
xDataTOM_Int xDataXML_convertAttributeTo_xDataTOM_Int( statusMessageReporting *smr, xDataXML_element *element, char const *name, xDataTOM_Int *n, int required ) {
/*
*   Returns 1 if no such attribute, -1 if error converting to xDataTOM_Int and 0 if successful.
*/
    char const *value;
    char *e;

    if( ( value = xDataXML_getAttributesValueInElement( element, name ) ) == NULL ) {
        if( required ) smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 1,
            "missing required attribute '%s'", name );
        return( 1 );
    }
    *n = (xDataTOM_Int) strtoll( value, &e, 10 );
    if( *e != 0 ) {
        smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 1, 
            "could not convert attribute %s's value = %s to an integer", name, value );
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int xDataXML_convertAttributeToDouble( statusMessageReporting *smr, xDataXML_element *element, char const *name, double *d, int required ) {
/*
*   Returns 1 if no such attribute, -1 if error converting to double and 0 if successful.
*/
    char const *value;
    char *e;

    if( ( value = xDataXML_getAttributesValueInElement( element, name ) ) == NULL ) {
        if( required ) smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 1,
            "missing required attribute '%s'", name );
        return( 1 );
    }
    *d = strtod( value, &e );
    if( *e != 0 ) {
        smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element) , xDataTOM_smrLibraryID, 1, 
            "could not convert attribute %s's values = %s to a double", name, value );
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int xDataXML_numberOfElementsByTagName( statusMessageReporting * /*smr*/, xDataXML_element *element, char const *tagName ) {

    int n = 0;
    xDataXML_element *child;

    for( child = xDataXML_getFirstElement( element ); child != NULL; child = xDataXML_getNextElement( child ) ) if( !strcmp( child->name, tagName ) ) n++;
    return( n );
}
/*
************************************************************
*/
xDataXML_elementList *xDataXML_getElementsByTagName( statusMessageReporting *smr, xDataXML_element *element, char const *tagName ) {

    int n = xDataXML_numberOfElementsByTagName( smr, element, tagName );
    size_t size;
    xDataXML_element *child;
    xDataXML_elementListItem *p;
    xDataXML_elementList *list = NULL;


    size = sizeof( xDataXML_elementList ) + n * sizeof( xDataXML_elementListItem );
    if( ( list = (xDataXML_elementList *) smr_malloc2( smr, size, 0, "list" ) ) != NULL ) {
        list->n = n;
        p = list->items = (xDataXML_elementListItem *) &(list[1]);
        for( child = xDataXML_getFirstElement( element ); child != NULL; child = xDataXML_getNextElement( child ) ) {
            if( !strcmp( child->name, tagName ) ) {
                p->element = child;
                p->sortString = NULL;
                p++;
            }
        }
    }
    return( list );
}
/*
************************************************************
*/
xDataXML_element *xDataXML_getOneElementByTagName( statusMessageReporting *smr, xDataXML_element *element, char *name, int required ) {

    xDataXML_elementList *list;
    xDataXML_element *xData = NULL;

    if( ( list = xDataXML_getElementsByTagName( smr, element, name ) ) != NULL ) {
        if( list->n == 0 ) {
            if( required ) smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 
                1, "element %s does not have sub-element named %s", element->fullName, name ); }
        else if( list->n > 1 ) {
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( element ), xDataTOM_smrLibraryID, 1, 
                "element %s contains more than one sub-element named %s", element->fullName, name ); }
        else {
            xData = list->items[0].element;
        }
        xDataXML_freeElementList( smr, list );
    }
    return( xData );
}
/*
************************************************************
*/
void xDataXML_freeElementList( statusMessageReporting * /*smr*/, xDataXML_elementList *list ) {

    smr_freeMemory( (void **) &list );
}
/*
************************************************************
*/
static char *xDataXML_getTraceback( statusMessageReporting *smr, xDataXML_element *element ) {
/*
*   Returned string must be freed by calling routine.
*/
    int size;
    char *s, *name;

    name = element->name;
    size = (int) strlen( name ) + 1;
    if( ( s = xDataXML_getTraceback2( smr, element->parentRoot, size ) ) != NULL ) {
        strcat( s, "/" );
        strcat( s, name );
    }
    return( s );
}
/*
************************************************************
*/
static char *xDataXML_getTraceback2( statusMessageReporting *smr, xDataXML_rootElement *parentRoot, int n ) {

    int size;
    char *s, *name;

    if( parentRoot->parentRoot == NULL ) {
        s = (char *) smr_malloc2( smr, n + 1, 0, "traceback string" );
        *s = 0; }
    else {
        name = parentRoot->parentElement->name;
        size = (int) strlen( name ) + 1;
        n += size;
        if( ( s = xDataXML_getTraceback2( smr, parentRoot->parentRoot, n ) ) != NULL ) {
            strcat( s, "/" );
            strcat( s, name );
        }
    }
    return( s );
}
/*
************************************************************
*/
int xDataXML_is_xDataType( statusMessageReporting *smr, xDataXMLType *xDT, char const * const ID, int setMsg ) {

    if( xDT->ID == NULL ) {
        if( setMsg ) smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( xDT->element ), xDataTOM_smrLibraryID, 1, 
            "element %s not xData object", xDT->element->fullName ); }
    else if( xDT->ID != ID ) {
        if( setMsg ) smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( xDT->element ), xDataTOM_smrLibraryID, 1, 
            "Element %s is not xData object of ID %s but %s", xDT->element->fullName, ID, xDT->ID );
    }
    return( xDT->ID == ID );
}
/*
************************************************************
*/
char const *xDataXML_getFileName( xDataXML_document *doc ) {

    return( doc->fileName );
}
/*
************************************************************
*/
char const *xDataXML_getRealFileName( xDataXML_document *doc ) {

    return( doc->realFileName );
}
/*
************************************************************
*/
static int xDataXML_setFileName( statusMessageReporting *smr, xDataXML_document *doc, char const *fileName ) {

    char realPath[PATH_MAX+1];

    smr_freeMemory( (void **) &(doc->fileName) );
    smr_freeMemory( (void **) &(doc->realFileName) );
    if( fileName != NULL ) {
        if( ( doc->fileName = smr_allocateCopyString2( smr, fileName, "fileName" ) ) == NULL ) return( 1 );
        if( realpath( fileName, realPath ) != NULL ) {
            if( ( doc->realFileName = smr_allocateCopyString2( smr, realPath, "realFileName" ) ) == NULL ) return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
xDataXML_document *xDataXML_getElementsDocument( xDataXML_element *element ) {

    xDataXML_rootElement* root = element->parentRoot;

    while( root->parentRoot != NULL ) root = root->parentRoot; // Loop checking, 11.06.2015, T. Koi
    return( root->xData_doc );
}
/*
************************************************************
*/
void *xDataXML_get_smrUserInterfaceFromDocument( xDataXML_document *doc ) {

    if( doc == NULL ) return( NULL );
    return( &(doc->smrUserInterface ) );
}
/*
************************************************************
*/
void *xDataXML_get_smrUserInterfaceFromElement( xDataXML_element *element ) {

    return( xDataXML_get_smrUserInterfaceFromDocument( xDataXML_getElementsDocument( element ) ) );
}
/*
************************************************************
*/
static int xDataXML_smrUserInterfaceInitialize( xDataXML_document *doc ) {

    doc->smrUserInterface.smrUserInterface = xDataXML_smrUserInterface;
    doc->smrUserInterface.doc = doc;
    return( 0 );
}
/*
************************************************************
*/
static int xDataXML_smrUserInterfaceFree( xDataXML_document *doc ) {

    doc->smrUserInterface.smrUserInterface = NULL;
    doc->smrUserInterface.doc = NULL;
    return( 0 );
}
/*
************************************************************
*/
static char *xDataXML_smrUserInterface( void *userData ) {

    xDataXML_smr *smrUserInterface = (xDataXML_smr *) userData;
    xDataXML_rootElement *currentRoot = smrUserInterface->doc->currentRoot;

    if( currentRoot->parentElement != NULL ) {
        return( smr_allocateFormatMessage( "\nat line %d and column %d of file %s\nin element %s", currentRoot->parentElement->docInfo.line, 
            currentRoot->parentElement->docInfo.column, smrUserInterface->doc->fileName, currentRoot->parentElement->fullName ) ); }
    else if( smrUserInterface->doc->fileName != NULL ) {
        return( smr_allocateFormatMessage( "\nof file %s", smrUserInterface->doc->fileName ) );
    }
    return( smr_allocateFormatMessage( "\nat line %d and column %d\nin element %s", currentRoot->parentElement->docInfo.line,
        currentRoot->parentElement->docInfo.column, currentRoot->parentElement->fullName ) );
}
/*
************************************************************
*/
int xDataXML_stringTo_xDataTOM_Int( statusMessageReporting *smr, void *smrUserInterface, char const *c, xDataTOM_Int *value, char const *endings, char **e ) {

    char const *s;
    char tmp[64];
    int status = 1, n = sizeof( tmp );

    for( s = c; *s != 0; s++ ) if( !isspace( *s ) ) break;
    *value = (xDataTOM_Int) strtoll( s, e, 10 );
    if( *e == s ) {
        smr_setReportError3(smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "could not convert \"%s\" to an integer", xDataXML_shortStringForMessage( n, tmp, c ));}
    else {
        if( *endings == 0 ) while( isspace( **e ) ) (*e)++; // Loop checking, 11.06.2015, T. Koi
        if( **e == 0 ) {
            status = 0; }
        else {
            if( *endings == 0 ) {
                smr_setReportError3( smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "integer string \"%s\" does not end with a '\\0'", 
                    xDataXML_shortStringForMessage( n, tmp, c ) ); }
            else {
                if( strchr( endings, **e ) == NULL ) {
                    smr_setReportError3( smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "integer string \"%s\" does not end with a white space or a '\\0\'", 
                        xDataXML_shortStringForMessage( n, tmp, c ) ); }
                else {
                    status = 0;
                }
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
int xDataXML_stringTo_double( statusMessageReporting *smr, void *smrUserInterface, char const *c, double *value, char const *endings, char **e ) {

    char const *s;
    char tmp[64];
    int status = 1, n = sizeof( tmp );

    for( s = c; *s != 0; s++ ) if( !isspace( *s ) ) break;
    *value = strtod( s, e );
    if( *e == s ) {
        smr_setReportError3( smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "could not convert \"%s\" to an double", 
            xDataXML_shortStringForMessage( n, tmp, c ));}
    else {
        if( *endings == 0 ) while( isspace( **e ) ) (*e)++; // Loop checking, 11.06.2015, T. Koi
        if( **e == 0 ) {
            status = 0; }
        else {
            if( *endings == 0 ) {
                smr_setReportError3( smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "double string \"%s\" does not end with a '\\0'", 
                    xDataXML_shortStringForMessage( n, tmp, c ) ); }
            else {
                if( strchr( endings, **e ) == NULL ) {
                    smr_setReportError3( smr, smrUserInterface, xDataTOM_smrLibraryID, 1, "double string \"%s\" does not end with a white space or a '\\0\'", 
                        xDataXML_shortStringForMessage( n, tmp, c ) ); }
                else {
                    status = 0;
                }
            }
        }
    }
    return( status );
}
/*
************************************************************
*/
int xDataXML_addToAccessed( statusMessageReporting * /*smr*/, xDataXML_element *element, int increment ) {

    element->accessed += increment;
    return( element->accessed );
}
/*
************************************************************
*/
int xDataXML_getAccessed( statusMessageReporting * /*smr*/, xDataXML_element *element ) {

    return( element->accessed );
}
/*
************************************************************
*/
static char const *xDataXML_shortStringForMessage( size_t size, char *Out, char const *In ) {

    if( strlen( In ) > size ) {
        strncpy( Out, In, size - 5 );
        Out[size-5] = 0;
        strcat( Out, " ..." );
        return( Out );
    }
    return( In );
}
/*
************************************************************
*/
static int xDataXML_constructTOM( statusMessageReporting *smr, xDataTOM_element *TE, xDataXML_element *XE ) {

    int i, status = 0;
    xDataTOM_element *TOMChild;
    xDataXML_element *XMLChild;
    xDataXML_attribute *attribute;
    char const *xDataValue = xDataXML_getAttributesValueInElement( XE, "xData" );

    if( !smr_isOk( smr ) ) return( 1 );
    if( ( TOMChild = xDataTOM_addElementInElement( smr, TE, XE->index, XE->name ) ) == NULL ) return( 1 );
    for( i = 0; 1; i++ ) {
        if( ( attribute = xDataXML_attributeByIndex( &(XE->attributes), i ) ) == NULL ) break;
        if( xDataTOME_addAttribute( smr, TOMChild, attribute->name, attribute->value ) != 0 ) return( 1 );
    }

    if( !strcmp( XE->name, xDataTOM_KalbachMann_ID ) ) {
        xDataValue = xDataTOM_KalbachMann_ID;
    }

    if( xDataValue == NULL ) {
        for( XMLChild = xDataXML_getFirstElement( XE ); ( status == 0 ) && ( XMLChild != NULL ); XMLChild = xDataXML_getNextElement( XMLChild ) ) {
            status = xDataXML_constructTOM( smr, TOMChild, XMLChild );
        } }
    else {
        if( strcmp( xDataValue, xDataTOM_XYs_ID ) == 0 ) {
            status = xDataXML_XYsToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_regionsXYs_ID ) == 0 ) {
            status = xDataXML_regionsXYsToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_W_XYs_ID ) == 0 ) {
            status = xDataXML_W_XYsToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_V_W_XYs_ID ) == 0 ) {
            status = xDataXML_V_W_XYsToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_W_XYs_LegendreSeries_ID ) == 0 ) {
            status = xDataXML_W_XYs_LegendreSeriesToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_regionsW_XYs_LegendreSeries_ID ) == 0 ) {
            status = xDataXML_regionsW_XYs_LegendreSeriesToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_V_W_XYs_LegendreSeries_ID ) == 0 ) {
            status = xDataXML_V_W_XYs_LegendreSeriesToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_KalbachMann_ID ) == 0 ) {
            status = xDataXML_KalbachMannToTOM( smr, XE, TOMChild ); }
        else if( strcmp( xDataValue, xDataTOM_polynomial_ID ) == 0 ) {
            status = xDataXML_polynomialToTOM( smr, XE, TOMChild ); }
        else {
            printf( "Unsupported xData type '%s' in element '%s'\n", xDataValue, XE->name );
#if 0
            smr_setReportError3( smr, xDataXML_get_smrUserInterfaceFromElement( XE ), xDataTOM_smrLibraryID, -1, 
                "Unsupported xData type = \"%s\"", xDataValue );
            status = 1;
#endif
        }
    }
    return( status );
}
/*
************************************************************
*/
void *xDataXML_initializeData( statusMessageReporting *smr, xDataXML_element *XE, xDataTOM_element *TE, char const *ID, size_t size ) {

    xDataTOM_xDataInfo *xDI = &(TE->xDataInfo);

    if( xData_initializeData( smr, TE, ID, size ) == NULL ) return( NULL );
    if( xDataXML_axesElememtToTOM( smr, XE, &(xDI->axes) ) != 0 ) smr_freeMemory( (void **) &(xDI->data) );
    return( xDI->data );
}

#if defined __cplusplus
}
#endif
