/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifndef WIN32
   #include <unistd.h>
#endif
#ifdef WIN32
   #include <io.h>
   #ifndef _SSIZE_T_DEFINED
      typedef int ssize_t;
      #define _SSIZE_T_DEFINED
   #endif
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "xData.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static void *xData_parseFreeElement( statusMessageReporting *smr, xData_element *element );
static void xData_parseFreeElementItems( statusMessageReporting *smr, xData_element *element );
static void XMLCALL xData_parseStartElement( void *userData, const char *name, const char **attris );
static void XMLCALL xData_parseEndElement( void *userData, const char *name );
static void XMLCALL xData_parseCharacterData( void *userData, const XML_Char *s, int len );
static void xData_parseInitializeRootElement( xData_document *doc, xData_rootElement *re, xData_element *parentElement, 
    int depth );
static int xData_parseInitializeText( xData_document *doc, xData_text *text );
static int xData_parseAddElementToRoot( statusMessageReporting *smr, xData_rootElement *parentRoot, const char *name, const char **attris );
static enum xData_errorCodes xData_parseGetCurrentPosition( xData_document *doc, xData_docInfo *docInfo );
static int xData_init_xDataTypeNone( xDataType *xDT, xData_element *element );
static char *xData_getTraceback( statusMessageReporting *smr, xData_element *element );
static char *xData_getTraceback2( statusMessageReporting *smr, xData_rootElement *parentRoot, int n );
static int xData_elementList_defaultSorter( void const *p1, void const *p2 );
static int xData_elementList_indexSorter( void const *p1, void const *p2 );
static int xData_smrUserInterfaceInitialize( xData_document *doc );
static int xData_smrUserInterfaceFree( xData_document *doc );
static int xData_smrUserInterface( void *userData, char **smr );
static char const *xData_shortStringForMessage( size_t size, char *Out, char const *In );
/*
************************************************************
*/
xData_document *xData_parseReadFile( statusMessageReporting *smr, const char *fileName, xData_xDataTypeOk func, void *userData ) {
/*
*   Returns NULL is any error occurred. If an error occurs in an expat routine, xData_parseEndOfXML will set smr appropriately.
*/
    int f;
    char buffer[10 * 1000];
    ssize_t count, n = sizeof( buffer ) - 1;
    ssize_t s = 0;
    xData_document *doc = NULL;

    if( ( doc = xData_parseMalloc( smr, func, userData ) ) != NULL ) {
        if( xData_setFileName( smr, doc, fileName ) == 0 ) {
            f = open( fileName, O_RDONLY );
            if( f == -1 ) {
                xData_parseEndOfXML( smr, doc );
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, xData_errFileError, "could not open file %s", fileName ); }
            else {
                while( ( count = read( f, buffer, n ) ) > 0 ) {
                    s += count;
                    buffer[count] = 0;
                    if( xData_parse( doc, buffer ) ) break;
                    if( !smr_isOk( doc->smr ) ) break;
                }
                close( f );
                xData_parseEndOfXML( smr, doc );
                if( count < 0 ) smr_setMessageError( smr, NULL, __FILE__, __LINE__, xData_errFileError, 
                    "read failed with errno = %d for %s", errno, fileName );
            }
        }
        if( ( doc != NULL ) && ( !smr_isOk( smr ) ) ) {
            xData_parseFree( smr, doc );
            doc = NULL;
        }
    }
    return( doc );
}
/*
************************************************************
*/
xData_document *xData_parseString( statusMessageReporting *smr, const char *str, xData_xDataTypeOk func, void *userData ) {
/*
*   Returns NULL is any error occurred. If an error occurs in an expat routine, xData_parseEndOfXML will set smr appropriately.
*/
    xData_document *doc = NULL;
    if( ( doc = xData_parseMalloc( smr, func, userData ) ) != NULL ) {
        xData_parse( doc, str );
        xData_parseEndOfXML( smr, doc );
        if( !smr_isOk( smr ) ) {
            xData_parseFree( smr, doc );
            doc = NULL;
        }
    }
    return( doc );
}
/*
************************************************************
*/
xData_document *xData_parseMalloc( statusMessageReporting *smr, xData_xDataTypeOk func, void *userData ) {
/*
*   Returns NULL is any error occurred.
*/
    xData_document *doc;

    //if( ( doc = xData_malloc2( smr, sizeof( xData_document ), 1, "xData_document" ) ) != NULL ) {
    if( ( doc = (xData_document*) xData_malloc2( smr, sizeof( xData_document ), 1, "xData_document" ) ) != NULL ) {
        //if( xData_parseInitialize( smr, doc, func, userData ) ) doc = xData_parseFree( smr, doc );
        if( xData_parseInitialize( smr, doc, func, userData ) ) doc = (xData_document*) xData_parseFree( smr, doc );
    }
    return( doc );
}
/*
************************************************************
*/
int xData_parseInitialize( statusMessageReporting *smr, xData_document *doc, xData_xDataTypeOk func, void *userData ) {

    XML_Parser xmlParser;

    doc->status = xData_statusParsing;
    doc->error = xData_errNone;
    //doc->err = 0;
    doc->err = (XML_Error) 0;
    doc->err_line = 0;
    doc->err_column = 0;
    doc->fileName = NULL;
    doc->xDataTypeOk_userFunction = func;
    doc->xDataTypeOk_userData = userData;
    xData_smrUserInterfaceInitialize( doc );
    doc->smr= smr;
    doc->xmlParser = xmlParser = XML_ParserCreate( NULL );
    if( xmlParser == NULL ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, xData_errXML_ParserCreate, "XML_ParserCreate failed" ); }
    else {
        XML_SetUserData( doc->xmlParser, doc  );
        xData_parseInitializeRootElement( doc, &(doc->root), NULL, 0 );
        doc->currentRoot = &(doc->root);
        XML_SetElementHandler( xmlParser, xData_parseStartElement, xData_parseEndElement );
        XML_SetCharacterDataHandler( xmlParser, xData_parseCharacterData );
    }
    return( !smr_isOk( smr ) );
}
/*
************************************************************
*/
int xData_parseEndOfXML( statusMessageReporting *smr, xData_document *doc ) {

    if( doc->xmlParser ) {
        doc->err = XML_GetErrorCode( doc->xmlParser );
        doc->err_line = XML_GetCurrentLineNumber( doc->xmlParser );
        doc->err_column = XML_GetCurrentColumnNumber( doc->xmlParser );
        if( smr_isOk( smr ) && ( XML_Parse( doc->xmlParser, NULL, 0, 1 ) == XML_STATUS_ERROR ) ) {
            doc->status = xData_statusError;
            smr_setMessageError( smr, xData_get_smrUserInterfaceFromDocument( doc ), __FILE__, __LINE__, xData_errXMLParser, 
                "status = %d\nXML_Error code = %d\nXML_ErrorString = %s\nerror line, column = %d, %d", xData_errXMLParser, 
                doc->err, XML_ErrorString( doc->err ), doc->err_line, doc->err_column );
        }
        XML_ParserFree( doc->xmlParser );
        doc->xmlParser = NULL;
        if( doc->status != xData_statusError ) doc->status = xData_statusCompleted;
    }
    return( 0 );
}
/*
************************************************************
*/
void *xData_parseFree( statusMessageReporting *smr, xData_document *doc ) {

    xData_parseEndOfXML( smr, doc );
    //doc->root.children = xData_parseFreeElement( smr, doc->root.children );
    doc->root.children = (xData_element*) xData_parseFreeElement( smr, doc->root.children );
    //doc->fileName = xData_free( smr, doc->fileName );
    doc->fileName = (char*) xData_free( smr, doc->fileName );
    xData_smrUserInterfaceFree( doc );
    xData_free( smr, doc );
    return( NULL );
}
/*
************************************************************
*/
static void *xData_parseFreeElement( statusMessageReporting *smr, xData_element *element ) {
    
    xData_element *next;

    if( element == NULL ) return( NULL );
    for( ; element != NULL; element = next ) {
        next = element->next;
        xData_parseFreeElementItems( smr, element );
        xData_free( smr, element );
    }
    return( NULL );
}
/*
************************************************************
*/
static void xData_parseFreeElementItems( statusMessageReporting *smr, xData_element *element ) {

    //element->childrenRoot.children = xData_parseFreeElement( smr, element->childrenRoot.children );
    element->childrenRoot.children = (xData_element*) xData_parseFreeElement( smr, element->childrenRoot.children );
    if( ( !strcmp( element->name, "xData" ) ) && ( element->xDataTypeInfo.release != NULL ) ) element->xDataTypeInfo.release( smr, 
        &(element->xDataTypeInfo) );
    xData_free( smr, element->name );
    xData_free( smr, element->fullName );
    if( element->attributes.attributes ) xData_free( smr, element->attributes.attributes );
    if( element->text.text ) xData_free( smr, element->text.text );
}
/*
************************************************************
*/
int xData_parse( xData_document *doc, const char *s ) {

    if( doc->status != xData_statusParsing ) return( doc->status );
    if( XML_Parse( doc->xmlParser, s, strlen( s ), 0 ) == XML_STATUS_ERROR ) return( -1 );
    return( 0 );
}
/*
************************************************************
*/
static void XMLCALL xData_parseStartElement( void *userData, const char *name, const char **attris ) {

    xData_document *doc = (xData_document *) userData;

    if( !smr_isOk( doc->smr ) ) return;
    xData_parseAddElementToRoot( doc->smr, doc->currentRoot, name, attris );
}
/*
************************************************************
*/
static void XMLCALL xData_parseEndElement( void *userData, const char *name ) {

    int status = 0;
    xData_document *doc = (xData_document *) userData;

    if( !smr_isOk( doc->smr ) ) return;
    if( !strcmp( name, "xData" ) ) {
        xData_element *element = doc->currentRoot->parentRoot->currentChild;
        xData_text *text = &(element->text);
        const char *value = xData_getAttributesValueInElement( element, "type" );
        if( !strcmp( value, xData_oned_x_ID ) ) {
            xData_init_1d_x( doc->smr, element ); }
        else if( !strcmp( value, xData_twod_xy_ID ) ) {
            xData_init_2d_xy( doc->smr, element ); }
        else if( !strcmp( value, xData_twod_xindex_y_ID ) ) {
            xData_init_2d_xindex_y( doc->smr, element ); }
        else if( !strcmp( value, xData_twod_xShared_yHistogram_ID ) ) {
            xData_init_2d_xShared_yHistogram( doc->smr, element ); }
        else if( !strcmp( value, xData_matrix_ID ) ) {
            xData_init_matrix( doc->smr, element ); }
        else {
            status = -1;
            if( doc->xDataTypeOk_userFunction != NULL ) {
                if( doc->xDataTypeOk_userFunction( value, doc, doc->xDataTypeOk_userData ) ) status = 1;
            }
            if( status < 0 ) smr_setMessageError( doc->smr, xData_get_smrUserInterfaceFromElement( element ), __FILE__, __LINE__, -1, 
                "Unsupported xData type = %s", value );
        }
        if( ( status == 0 ) && ( smr_isOk( doc->smr ) ) ) status = element->xDataTypeInfo.toData( doc->smr, &(element->xDataTypeInfo), 
            &(element->attributes), text->text );
    }
    doc->currentRoot->currentChild = NULL;
    doc->currentRoot = doc->currentRoot->parentRoot;
}
/*
************************************************************
*/
static void XMLCALL xData_parseCharacterData( void *userData, const XML_Char *s, int len ) {
/*
*   Always terminates text with a 0.
*/
    //xData_document *doc = userData;
    xData_document *doc = (xData_document*) userData;
    xData_text *text = &(doc->currentRoot->parentRoot->currentChild->text);
    int needSize = text->length + len + 1, l; 
    char *p;

    if( !smr_isOk( doc->smr ) ) return;
    if( needSize < 8  ) needSize = 8;
    //if( needSize > text->allocated ) {
    if( needSize > (int) text->allocated ) {
        if( text->allocated != 0 ) {
            l = ( 20 * text->allocated ) / 100;
            if( l < 100 ) l = 100;
            //if( needSize < ( text->allocated + l ) ) needSize = text->allocated + l;
            if( needSize < ( (int) text->allocated + l ) ) needSize = text->allocated + l;
        }
        text->allocated = needSize;
        //text->text = xData_realloc2( doc->smr, text->text, text->allocated, "text" );
        text->text = (char*) xData_realloc2( doc->smr, text->text, text->allocated, "text" );
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
static void xData_parseInitializeRootElement( xData_document *doc, xData_rootElement *re, 
        xData_element *parentElement, int depth ) {

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
static int xData_parseInitializeText( xData_document *doc, xData_text *text ) {

    xData_parseGetCurrentPosition( doc, &(text->docInfo) );
    text->allocated = 0;
    text->length = 0;
    text->text = NULL;
    return( 0 );
}
/*
************************************************************
*/
static int xData_parseAddElementToRoot( statusMessageReporting *smr, xData_rootElement *parentRoot, const char *name, const char **attris ) {

    xData_document *doc = parentRoot->xData_doc;
    xData_element *element;
    int i, n, status = 1;
    size_t lens;
    char *p, *e;
    const char **pAttris;
    xData_attribute *a;
    void *smrUser;

    //element = xData_malloc2( doc->smr, sizeof( xData_element ), 1, "xData_element" );
    element = (xData_element*) xData_malloc2( doc->smr, sizeof( xData_element ), 1, "xData_element" );
    if( element == NULL ) return( 1 );
    xData_parseGetCurrentPosition( doc, &(element->docInfo) );
    element->ordinal = parentRoot->numberOfElements;
    element->index = -1;
    element->accessed = 0;
    element->parentRoot = parentRoot;
    xData_parseInitializeRootElement( doc, &(element->childrenRoot), element, parentRoot->depth + 1 );
    element->next = NULL;
    //if( ( element->name = xData_malloc2( doc->smr, strlen( name ) + 1, 0, "name" ) ) == NULL ) {
    if( ( element->name = (char*) xData_malloc2( doc->smr, strlen( name ) + 1, 0, "name" ) ) == NULL ) {
        xData_free( smr, element );
        return( 1 );
    }
    strcpy( element->name, name );
    if( ( element->fullName = xData_getTraceback( smr, element ) ) == NULL ) {
        xData_free( smr, element->name );
        xData_free( smr, element );
        return( 1 );
    }
    for( i = 0, lens = 0, pAttris = attris; *pAttris; i++, pAttris++ ) lens += strlen( *pAttris ) + 1;
    n = i / 2;
    element->attributes.size = n * sizeof( xData_attribute ) + lens;
    element->attributes.number = n;
    element->attributes.attributes = NULL;
    smrUser = xData_get_smrUserInterfaceFromElement( element );
    if( element->attributes.size  ) {
        //if( ( element->attributes.attributes = xData_malloc2( doc->smr, element->attributes.size, 0, "attributes") ) == NULL ) {
        if( ( element->attributes.attributes = (xData_attribute*) xData_malloc2( doc->smr, element->attributes.size, 0, "attributes") ) == NULL ) {
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
#ifndef WIN32
                    element->index = strtoll( a->value, &e, 10 );
#endif
#ifdef WIN32
                    element->index = strtol( a->value, &e, 10 );
#endif
                    if( *e != 0 ) {
                        status = 0;
                        smr_setMessageError( doc->smr, smrUser, __FILE__, __LINE__, -1, "could not convert index attribute = %s to integer", a->value );
                    }
                }
            }
        }
    }
    if( !status ) {
        xData_free( smr, element->attributes.attributes );
        xData_free( smr, element->name );
        xData_free( smr, element->fullName );
        xData_free( smr, element );
        return( 1 );
    }
    xData_init_xDataTypeNone( &(element->xDataTypeInfo), element );
    element->textOffset = 0;
    xData_parseInitializeText( doc, &(element->text) );
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
static enum xData_errorCodes xData_parseGetCurrentPosition( xData_document *doc, xData_docInfo *docInfo ) {

    docInfo->column = XML_GetCurrentColumnNumber( doc->xmlParser );
    docInfo->line = XML_GetCurrentLineNumber( doc->xmlParser );
    return( xData_errNone );
}
/*
************************************************************
*/
int xData_parseIsError( xData_document *doc ) {

    return( doc->status == xData_statusError );
}
/*
************************************************************
*/
xData_element *xData_getDocumentsElement( xData_document *doc ) { return( doc->root.children ); }
xData_element *xData_getFirstElement( xData_element *element ) { return( element->childrenRoot.children ); }
xData_element *xData_getNextElement( xData_element *element ) { return( element->next ); }
/*
************************************************************
*/
enum xData_itemMode xData_getFirstItem( xData_element *element, xData_item *item ) {

    item->parentElement = element;
    item->element = xData_getFirstElement( element );
    if( item->element == NULL ) {
        item->mode = xData_itemModeText;
        if( element->text.length == 0 ) item->mode = xData_itemModeEnd; }
    else {
        item->mode = xData_itemModeElement;
        if( 0 < item->element->textOffset ) item->mode = xData_itemModeText;
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
enum xData_itemMode xData_getNextItem( xData_item *item ) {

    if( item->mode != xData_itemModeEnd ) {
        if( item->mode == xData_itemModeText ) {
            item->mode = xData_itemModeElement;
            if( item->element == NULL ) item->mode = xData_itemModeEnd;
            item->textOffset += item->textLength;
            item->textLength = 0;
            item->text = &(item->parentElement->text.text[item->textOffset]); }
        else {
            item->element = item->element->next;
            item->mode = xData_itemModeText;
            if( item->element == NULL ) {
                if( item->textOffset < item->parentElement->text.length ) {
                    item->textLength = item->parentElement->text.length - item->textOffset; }
                else {
                    item->mode = xData_itemModeEnd;
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
char *xData_getAttributesValue( xData_attributionList *attributes, const char *name ) {

    int i;
    char *value = NULL;

    for( i = 0; i < attributes->number; i++ ) {
        if( !strcmp( attributes->attributes[i].name, name ) ) {
            value = attributes->attributes[i].value;
            break;
        }
    }
    return( value );
}
/*
************************************************************
*/
const char *xData_getAttributesValueInElement( xData_element *element, const char *name ) {

    return( (const char *) xData_getAttributesValue( &(element->attributes), name ) );
}
/*
************************************************************
*/
//int xData_initializeAttributionList( statusMessageReporting *smr, xData_attributionList *attributes ) {
int xData_initializeAttributionList( statusMessageReporting *, xData_attributionList *attributes ) {

    attributes->number = 0;
    attributes->size = 0;
    attributes->attributes = NULL;
    return( 0 );
}
/*
************************************************************
*/
int xData_copyAttributionList( statusMessageReporting *smr, xData_attributionList *dest, xData_attributionList *src ) {
/*
*   The dest must not be active, else a memory leak will occur, as this routine does no free-ing.
*/
    int i;
    size_t lens;
    xData_attribute *d, *s;
    char *p;

    //if( ( dest->attributes = xData_malloc2( smr, src->size, 0, "attributes" ) ) == NULL ) return( 1 );
    if( ( dest->attributes = (xData_attribute*) xData_malloc2( smr, src->size, 0, "attributes" ) ) == NULL ) return( 1 );
    dest->number = src->number;
    dest->size = src->size;
    d = dest->attributes;
    s = src->attributes;
    p = (char *) &(dest->attributes[src->number]);
    for( i = 0; i < src->number; i++, s++, d++ ) {
        lens = strlen( s->name ) + 1;
        d->name = p;
        strcpy( p, s->name );
        p += lens;
        lens = strlen( s->value ) + 1;
        d->value= p;
        strcpy( p, s->value );
        p += lens;
    }

    return( 0 );
}
/*
************************************************************
*/
int xData_attributeListLength( xData_attributionList *attributes ) {

    return( attributes->number );
}
/*
************************************************************
*/
xData_attribute *xData_attributeByIndex( xData_attributionList *attributes, int index ) {

    if( index >= attributes->number ) return( NULL );
    return( &(attributes->attributes[index]) );
}
/*
************************************************************
*/
int xData_releaseAttributionList( statusMessageReporting *smr, xData_attributionList *attributes ) {

    attributes->number = 0;
    attributes->size = 0;
    //attributes->attributes = xData_free( smr, attributes->attributes );
    attributes->attributes = (xData_attribute*) xData_free( smr, attributes->attributes );
    return( 0 );
}
/*
************************************************************
*/
xData_element *xData_getElements_xDataElement( statusMessageReporting *smr, xData_element *element ) {

    //return( xData_getOneElementByTagName( smr, element, "xData", 1 ) );
    return( xData_getOneElementByTagName( smr, element, (char*) "xData", 1 ) );
}
/*
************************************************************
*/
static int xData_init_xDataTypeNone( xDataType *xDT, xData_element *element ) {

    xDT->status = xData_xDataType_Ok;
    xDT->typeString = NULL;
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
int xData_getCommonData( statusMessageReporting *smr, xData_element *element, xData_Int *index, xData_Int *start, xData_Int *end,
        xData_Int *length ) {

    if( element->xDataTypeInfo.typeString == NULL ) {
        smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( element ), __FILE__, __LINE__, 1, "element %s is not xData", element->fullName );
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
int xData_xDataTypeConvertAttributes( statusMessageReporting *smr, xData_element *element ) {

    xDataType *xDT = &(element->xDataTypeInfo);
    void *smrUser = xData_get_smrUserInterfaceFromElement( element );

    xDT->index = -1;
    xDT->start = -1;
    xDT->end = -1;
    xDT->length = -1;
    if( ( xDT->indexPresent = xData_convertAttributeTo_xData_Int( smr, element, "index", &(xDT->index) ) ) < 0 ) return( 1 );
    if( ( xDT->startPresent = xData_convertAttributeTo_xData_Int( smr, element, "start", &(xDT->start) ) ) < 0 ) return( 1 );
    if( ( xDT->endPresent = xData_convertAttributeTo_xData_Int( smr, element, "end", &(xDT->end) ) ) < 0 ) return( 1 );
    if( ( xDT->lengthPresent = xData_convertAttributeTo_xData_Int( smr, element, "length", &(xDT->length) ) ) < 0 ) return( 1 );
    if( ( xDT->endPresent > 0 ) ) {
        if( xDT->lengthPresent > 0 ) {
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "missing length (or end) in xData" );
            return( 1 );
        }
        xDT->end = xDT->length; }
    else {
        if( xDT->lengthPresent > 0 ) xDT->length = xDT->end;
    }

    if( xDT->startPresent > 0 ) xDT->start = 0;
    if( xDT->start < 0 ) {
        smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "start = %d < 0", xDT->start );
        return( 1 );
    }
    if( xDT->end < xDT->start ) {
        smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "start = %d >= end = %d", xDT->start, xDT->end );
        return( 1 );
    }
    if( xDT->length < 0 ) {
        smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "length = %d < 0", xDT->length );
        return( 1 );
    }

    return( 0 );
}
/*
************************************************************
*/
xData_Int xData_convertAttributeTo_xData_Int( statusMessageReporting *smr, xData_element *element, const char *name, xData_Int *n ) {
/*
*   Returns 1 if no such attribute, -1 if error converting to xData_Int and 0 if successful.
*/
    const char *value;
    char *e;

    if( ( value = xData_getAttributesValueInElement( element, name ) ) == NULL ) return( 1 );
    //*n = strtoll( value, &e, 10 );
    *n = strtol( value, &e, 10 );
    if( *e != 0 ) {
        smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( element ), __FILE__, __LINE__, 1, 
            "could not convert attribute %s's value = %s to an integer", name, value );
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int xData_convertAttributeToDouble( statusMessageReporting *smr, xData_element *element, const char *name, double *d ) {
/*
*   Returns 1 if no such attribute, -1 if error converting to double and 0 if successful.
*/
    const char *value;
    char *e;

    if( ( value = xData_getAttributesValueInElement( element, name ) ) == NULL ) return( 1 );
    *d = strtod( value, &e );
    if( *e != 0 ) {
        smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( element) , __FILE__, __LINE__, 1, 
            "could not convert attribute %s's values = %s to a double", name, value );
        return( -1 );
    }
    return( 0 );
}
/*
************************************************************
*/
//int xData_numberOfElementsByTagName( statusMessageReporting *smr, xData_element *element, const char *tagName ) {
int xData_numberOfElementsByTagName( statusMessageReporting *, xData_element *element, const char *tagName ) {

    int n = 0;
    xData_element *child;

    for( child = xData_getFirstElement( element ); child != NULL; child = xData_getNextElement( child ) ) if( !strcmp( child->name, tagName ) ) n++;
    return( n );
}
/*
************************************************************
*/
xData_elementList *xData_getElementsByTagName( statusMessageReporting *smr, xData_element *element, const char *tagName ) {

    int n = xData_numberOfElementsByTagName( smr, element, tagName );
    size_t size;
    xData_element *child;
    xData_elementListItem *p;
    xData_elementList *list = NULL;


    size = sizeof( xData_elementList ) + n * sizeof( xData_elementListItem );
    //if( ( list = xData_malloc2( smr, size, 0, "list" ) ) != NULL ) {
    if( ( list = (xData_elementList*) xData_malloc2( smr, size, 0, "list" ) ) != NULL ) {
        list->n = n;
        p = list->items = (xData_elementListItem *) &(list[1]);
        for( child = xData_getFirstElement( element ); child != NULL; child = xData_getNextElement( child ) ) {
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
xData_elementList *xData_getElementsByTagNameAndSort( statusMessageReporting *smr, xData_element *element, const char *tagName, 
    const char *sortAttributeName, xData_sortElementFunc sortFunction ) {

    int i;
    xData_elementList *list = xData_getElementsByTagName( smr, element, tagName );
    xData_elementListItem *p;

    if( list != NULL ) {
        if( sortFunction == NULL ) {
            sortFunction = (xData_sortElementFunc) xData_elementList_defaultSorter;
            if( sortAttributeName == NULL ) sortFunction = (xData_sortElementFunc) xData_elementList_indexSorter;
        }
        if( sortAttributeName == NULL ) sortAttributeName = "index";
        for( i = 0, p = list->items; i < list->n; i++, p++ ) p->sortString = xData_getAttributesValueInElement( p->element, sortAttributeName );
        qsort( list->items, list->n, sizeof( xData_elementListItem ), sortFunction );
    }

    return( list );
}
/*
************************************************************
*/
xData_element *xData_getOneElementByTagName( statusMessageReporting *smr, xData_element *element, char *name, int required ) {

    xData_elementList *list;
    xData_element *xData = NULL;

    if( ( list = xData_getElementsByTagName( smr, element, name ) ) != NULL ) {
        if( list->n == 0 ) {
            if( required ) smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( element ), __FILE__, __LINE__, 
                1, "element %s does not have sub-element named %s", element->fullName, name ); }
        else if( list->n > 1 ) {
            smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( element ), __FILE__, __LINE__, 1, 
                "element %s contains more than one sub-element named %s", element->fullName, name ); }
        else {
            xData = list->items[0].element;
        }
        xData_freeElementList( smr, list );
    }
    return( xData );
}
/*
************************************************************
*/
void xData_freeElementList( statusMessageReporting *smr, xData_elementList *list ) {

    xData_free( smr, list );
}
/*
************************************************************
*/
static char *xData_getTraceback( statusMessageReporting *smr, xData_element *element ) {
/*
*   Returned string must be freed by calling routine.
*/
    int size;
    char *s, *name;

    name = element->name;
    size = strlen( name ) + 1;
    if( ( s = xData_getTraceback2( smr, element->parentRoot, size ) ) != NULL ) {
        strcat( s, "/" );
        strcat( s, name );
    }
    return( s );
}
/*
************************************************************
*/
static char *xData_getTraceback2( statusMessageReporting *smr, xData_rootElement *parentRoot, int n ) {

    int size;
    char *s, *name;

    if( parentRoot->parentRoot == NULL ) {
        //s = xData_malloc2( smr, n + 1, 0, "traceback string" );
        s = (char*) xData_malloc2( smr, n + 1, 0, "traceback string" );
        *s = 0; }
    else {
        name = parentRoot->parentElement->name;
        size = strlen( name ) + 1;
        n += size;
        if( ( s = xData_getTraceback2( smr, parentRoot->parentRoot, n ) ) != NULL ) {
            strcat( s, "/" );
            strcat( s, name );
        }
    }
    return( s );
}
/*
************************************************************
*/
static int xData_elementList_defaultSorter( void const *p1, void const *p2 ) {

    const char *s1 = ((xData_elementListItem *) p1)->sortString, *s2 = ((xData_elementListItem *) p2)->sortString;

    if( s2 == NULL ) return( -1 );
    if( s1 == NULL ) return( 1 );
    return( strcmp( s1, s2 ) );
}
/*
************************************************************
*/
static int xData_elementList_indexSorter( void const *p1, void const *p2 ) {

    xData_element *e1 = ((xData_elementListItem *) p1)->element, *e2 = ((xData_elementListItem *) p2)->element;

    return( e1->index - e2->index );
}
/*
************************************************************
*/
int xData_is_xDataType( statusMessageReporting *smr, xDataType *xDT, char const * const type, int setMsg ) {

    if( xDT->typeString == NULL ) {
        if( setMsg ) smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( xDT->element ), __FILE__, __LINE__, 1, 
            "element %s not xData object", xDT->element->fullName ); }
    else if( xDT->typeString != type ) {
        if( setMsg ) smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( xDT->element ), __FILE__, __LINE__, 1, 
            "Element %s is not xData object of type %s", type );
    }
    return( xDT->typeString == type );
}
/*
************************************************************
*/
char const *xData_getFileName( xData_document *doc ) {

    return( doc->fileName );
}
/*
************************************************************
*/
int xData_setFileName( statusMessageReporting *smr, xData_document *doc, char const *fileName ) {

    doc->fileName = (char*) xData_free( smr, doc->fileName );
    if( fileName != NULL ) {
        //if( ( doc->fileName = xData_malloc2( smr, strlen( fileName ) + 1, 0, "doc->fileName" ) ) == NULL ) return( 1 );
        if( ( doc->fileName = (char*) xData_malloc2( smr, strlen( fileName ) + 1, 0, "doc->fileName" ) ) == NULL ) return( 1 );
        strcpy( doc->fileName, fileName );
    }
    return( 0 );
}
/*
************************************************************
*/
xData_document *xData_getElementsDocument( xData_element *element ) {

    xData_rootElement* root = element->parentRoot;

    while( root->parentRoot != NULL ) root = root->parentRoot;
    return( root->xData_doc );
}
/*
************************************************************
*/
void *xData_get_smrUserInterfaceFromDocument( xData_document *doc ) {

    if( doc == NULL ) return( NULL );
    return( &(doc->smrUserInterface ) );
}
/*
************************************************************
*/
void *xData_get_smrUserInterfaceFromElement( xData_element *element ) {

    return( xData_get_smrUserInterfaceFromDocument( xData_getElementsDocument( element ) ) );
}
/*
************************************************************
*/
static int xData_smrUserInterfaceInitialize( xData_document *doc ) {

    doc->smrUserInterface.smrUserInterface = xData_smrUserInterface;
    doc->smrUserInterface.doc = doc;
    return( 0 );
}
/*
************************************************************
*/
static int xData_smrUserInterfaceFree( xData_document *doc ) {

    doc->smrUserInterface.smrUserInterface = NULL;
    doc->smrUserInterface.doc = NULL;
    return( 0 );
}
/*
************************************************************
*/
static int xData_smrUserInterface( void *userData, char **str ) {

    int size, fileNameSize = 0, elementSize = 0;
    xData_smr *smrUserInterface = (xData_smr *) userData;
    static const char lcl[] = "\nat line %d and column %d", el[] = "\nin element ", fl[] = "\nof file ";
    char str_lcl[sizeof( lcl ) + 40];
    xData_rootElement *currentRoot = smrUserInterface->doc->currentRoot;

    if( smrUserInterface->doc->fileName != NULL ) fileNameSize = strlen( smrUserInterface->doc->fileName ) + strlen( fl );
    if( currentRoot != NULL ) {
        if( currentRoot->parentElement != NULL ) {
            sprintf( str_lcl, lcl, (int)currentRoot->parentElement->docInfo.line, (int)currentRoot->parentElement->docInfo.column );
            elementSize = strlen( str_lcl ) + strlen( currentRoot->parentElement->fullName ) + strlen( el );
        }
    }
    size = fileNameSize + elementSize;
    if( ( fileNameSize != 0 ) && ( elementSize != 0 ) ) size++;
    if( ( size > 0 ) && ( str != NULL ) ) {
        //if( ( *str = malloc( size + 1 ) ) == NULL ) return( -1 );
        if( ( *str = (char*) malloc( size + 1 ) ) == NULL ) return( -1 );
        if( ( size != 0 ) && ( elementSize != 0 ) ) {
            sprintf( *str, "%s%s%s%s%s", str_lcl, el, currentRoot->parentElement->fullName, fl, smrUserInterface->doc->fileName ); }
        else if( size != 0 ) {
            sprintf( *str, "%s%s", fl, smrUserInterface->doc->fileName ); }
        else {
            sprintf( *str, "%s%s%s", str_lcl, el, currentRoot->parentElement->fullName );
        }
    }
    return( size );
}
/*
************************************************************
*/
int xData_stringTo_xData_Int( statusMessageReporting *smr, void *smrUserInterface, char const *c, xData_Int *value, char const *endings, char **e ) {

    char const *s;
    char tmp[64];
    int status = 1, n = sizeof( tmp );

    for( s = c; *s != 0; s++ ) if( !isspace( *s ) ) break;
    //*value = strtoll( s, e, 10 );
    *value = strtol( s, e, 10 );
    if( *e == s ) {
        smr_setMessageError(smr, smrUserInterface, __FILE__, __LINE__, 1, "could not convert \"%s\" to an integer", xData_shortStringForMessage( n, tmp, c ));}
    else {
        if( *endings == 0 ) while( isspace( **e ) ) (*e)++;
        if( **e == 0 ) {
            status = 0; }
        else {
            if( *endings == 0 ) {
                smr_setMessageError( smr, smrUserInterface, __FILE__, __LINE__, 1, "integer string \"%s\" does not end with a '\\0'", 
                    xData_shortStringForMessage( n, tmp, c ) ); }
            else {
                if( strchr( endings, **e ) == NULL ) {
                    smr_setMessageError( smr, smrUserInterface, __FILE__, __LINE__, 1, "integer string \"%s\" does not end with a white space or a '\\0\'", 
                        xData_shortStringForMessage( n, tmp, c ) ); }
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
int xData_stringTo_double( statusMessageReporting *smr, void *smrUserInterface, char const *c, double *value, char const *endings, char **e ) {

    char const *s;
    char tmp[64];
    int status = 1, n = sizeof( tmp );

    for( s = c; *s != 0; s++ ) if( !isspace( *s ) ) break;
    *value = strtod( s, e );
    if( *e == s ) {
        smr_setMessageError(smr, smrUserInterface, __FILE__, __LINE__, 1, "could not convert \"%s\" to an double", xData_shortStringForMessage( n, tmp, c ));}
    else {
        if( *endings == 0 ) while( isspace( **e ) ) (*e)++;
        if( **e == 0 ) {
            status = 0; }
        else {
            if( *endings == 0 ) {
                smr_setMessageError( smr, smrUserInterface, __FILE__, __LINE__, 1, "double string \"%s\" does not end with a '\\0'", 
                    xData_shortStringForMessage( n, tmp, c ) ); }
            else {
                if( strchr( endings, **e ) == NULL ) {
                    smr_setMessageError( smr, smrUserInterface, __FILE__, __LINE__, 1, "double string \"%s\" does not end with a white space or a '\\0\'", 
                        xData_shortStringForMessage( n, tmp, c ) ); }
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
//int xData_addToAccessed( statusMessageReporting *smr, xData_element *element, int increment ) {
int xData_addToAccessed( statusMessageReporting *, xData_element *element, int increment ) {

    element->accessed += increment;
    return( element->accessed );
}
/*
************************************************************
*/
//int xData_getAccessed( statusMessageReporting *smr, xData_element *element ) {
int xData_getAccessed( statusMessageReporting *, xData_element *element ) {

    return( element->accessed );
}
/*
************************************************************
*/
static char const *xData_shortStringForMessage( size_t size, char *Out, char const *In ) {

    if( strlen( In ) > size ) {
        strncpy( Out, In, size - 5 );
        Out[size-5] = 0;
        strcat( Out, " ..." );
        return( Out );
    }
    return( In );
}

#if defined __cplusplus
}
#endif
