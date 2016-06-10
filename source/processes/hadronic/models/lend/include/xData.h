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
#ifndef xData_h_included
#define xData_h_included

#if defined __cplusplus
    extern "C" {
#endif

#include <expat.h>
#include "statusMessageReporting.h"

#if defined __cplusplus
    namespace GIDI {
#endif

typedef int xData_Int;

#ifndef XMLCALL
#define XMLCALL
#endif

#ifndef XML_Char
#define XML_Char char
#endif

#ifndef XML_Size
#define XML_Size long
#endif

//TK move follwoing reference from source 
//extern char const * const xData_oned_x_ID;
//extern char const * const xData_twod_xy_ID;
//extern char const * const xData_twod_xindex_y_ID;
//extern char const * const xData_twod_xShared_yHistogram_ID;
//extern char const * const xData_matrix_ID;

const char * const xData_oned_x_ID = "1d.x";
const char * const xData_twod_xy_ID = "2d.xy";
const char * const xData_twod_xindex_y_ID = "2d.xindex_y";
const char * const xData_twod_xShared_yHistogram_ID = "2d_xShared_yHistogram";
const char * const xData_matrix_ID = "matrix";

enum xData_status { xData_statusParsing = 1, xData_statusCompleted, xData_statusError };
enum xData_errorCodes { xData_errNone, xData_errXML_ParserCreate, xData_errFileError,
    xData_errXMLParser };
enum xData_itemMode { xData_itemModeEnd = 0, xData_itemModeElement, xData_itemModeText };
enum xData_xDataType { xData_xDataType_Ok, xData_xDataType_ConvertingToData, xData_xDataType_ConvertingToString };

typedef struct xData_docInfo_s xData_docInfo;
typedef struct xData_attribute_s xData_attribute;
typedef struct xData_attributionList_s xData_attributionList;
typedef struct xData_rootElement_s xData_rootElement;
typedef struct xData_element_s xData_element;
typedef struct xData_text_s xData_text;
typedef struct xDataType_s xDataType;
typedef struct xData_smr_s xData_smr;
typedef struct xData_document_s xData_document;
typedef struct xData_item_s xData_item;
typedef struct xData_2d_xindex_y_s xData_2d_xindex_y;
typedef struct xData_matrix_s xData_matrix;
typedef struct xData_matrix_rowStartEnd_s xData_matrix_rowStartEnd;
typedef struct xData_elementListItem_s xData_elementListItem;
typedef struct xData_elementList_s xData_elementList;

typedef int (*xData_sortElementFunc)( const void *, const void * );

struct xData_2d_xindex_y_s {
    xData_Int index;
    double value;
};

struct xData_matrix_rowStartEnd_s {
    xData_Int row, start, end;
};

struct xData_matrix_s {
    xData_Int rows, columns;
    xData_matrix_rowStartEnd *rowStartEnds;
    double *values;
};

struct xData_docInfo_s {
    size_t line;
    size_t column;
};

struct xData_attribute_s {
    char *name;
    char *value;
};

struct xData_attributionList_s {
    int number;
    size_t size;
    xData_attribute *attributes;
};

struct xData_text_s {
    xData_docInfo docInfo;
    size_t allocated;
    size_t length;
    char *text;
};

typedef int (*xData_xDataTypeOk)( char const *name, xData_document *doc, void *userData );
typedef int (*xDT_toDataFunction)( statusMessageReporting *smr, xDataType *, xData_attributionList *, const char * );
typedef char *(*xDT_toStringFunction)( statusMessageReporting *smr, xDataType * );
typedef int (*xDT_releaseFunction)( statusMessageReporting *smr, xDataType * );

struct xData_elementListItem_s {
    xData_element *element;
    const char *sortString;
};

struct xData_elementList_s {
    int n;
    xData_elementListItem *items;
};

struct xDataType_s {
    enum xData_xDataType status;
    const char *typeString;
    xData_element *element;
    xDT_toDataFunction toData;
    xDT_toStringFunction toString;
    xDT_releaseFunction release;
    signed char indexPresent, startPresent, endPresent, lengthPresent;
    xData_Int index, start, end, length;
    void *data;
};

struct xData_rootElement_s {
    xData_document *xData_doc;
    xData_element *parentElement;
    xData_rootElement *parentRoot;
    int depth;
    int numberOfElements;
    xData_element *children;
    xData_element *currentChild;
};

struct xData_element_s {
    xData_docInfo docInfo;
    int ordinal;                                    /* Counting from 0. */
    int index;                                      /* Value from "index" attribute if present or -1 */
    int accessed;                                   /* For the convenience of the users, not used internally. */
    xData_rootElement *parentRoot;
    xData_rootElement childrenRoot;
    xData_element *next;
    char *name;                                     /* Allocated in xData_parseAddElementToRoot. */
    char *fullName;                                 /* Allocated in xData_parseAddElementToRoot. */
    xData_attributionList attributes;               /* attributes->abbributes is allocated in xData_parseAddElementToRoot. */
    xDataType xDataTypeInfo;
    size_t textOffset;
    xData_text text;
};

struct xData_smr_s {
    smr_userInterface smrUserInterface;
    xData_document *doc;
};

struct xData_document_s {
    enum xData_status status;
    enum xData_errorCodes error;
    enum XML_Error err;
    XML_Size err_line, err_column;
    char *fileName;
    xData_xDataTypeOk xDataTypeOk_userFunction;
    void *xDataTypeOk_userData;
    xData_smr smrUserInterface;
    statusMessageReporting *smr;
    XML_Parser xmlParser;
    xData_rootElement root;
    xData_rootElement *currentRoot;
};

struct xData_item_s {
    xData_element *parentElement;
    xData_element *element;
    enum xData_itemMode mode;
    size_t textOffset;
    size_t textLength;
    char *text;
};

xData_document *xData_parseReadFile( statusMessageReporting *smr, const char *fileName, xData_xDataTypeOk func, void *userData );
xData_document *xData_parseString( statusMessageReporting *smr, const char *str, xData_xDataTypeOk func, void *userData );
xData_document *xData_parseMalloc( statusMessageReporting *smr, xData_xDataTypeOk func, void *userData );
int xData_parseInitialize( statusMessageReporting *smr, xData_document *xData_doc, xData_xDataTypeOk func, void *userData );
int xData_parseEndOfXML( statusMessageReporting *smr, xData_document *xData_doc );
void *xData_parseFree( statusMessageReporting *smr, xData_document *xData_doc );
int xData_parse( xData_document *xData_doc, const char *s );
int xData_parseIsError( xData_document *xData_doc );
xData_element *xData_getDocumentsElement( xData_document *xData_doc );
xData_element *xData_getFirstElement( xData_element *element );
xData_element *xData_getNextElement( xData_element *element );
enum xData_itemMode xData_getFirstItem( xData_element *element, xData_item *item );
enum xData_itemMode xData_getNextItem( xData_item *item );
char *xData_getAttributesValue( xData_attributionList *attributes, const char *name );
const char *xData_getAttributesValueInElement( xData_element *element, const char *name );
int xData_initializeAttributionList( statusMessageReporting *smr, xData_attributionList *attributes );
int xData_copyAttributionList( statusMessageReporting *smr, xData_attributionList *dest, xData_attributionList *src );
int xData_releaseAttributionList( statusMessageReporting *smr, xData_attributionList *attributes );
int xData_attributeListLength( xData_attributionList *attributes );
xData_attribute *xData_attributeByIndex( xData_attributionList *attributes, int index );
xData_element *xData_getElements_xDataElement( statusMessageReporting *smr, xData_element *element );
int xData_getCommonData( statusMessageReporting *smr, xData_element *element, xData_Int *index, xData_Int *start, xData_Int *end,
        xData_Int *length );
int xData_xDataTypeConvertAttributes( statusMessageReporting *smr, xData_element *element );
xData_Int xData_convertAttributeTo_xData_Int( statusMessageReporting *smr, xData_element *element, const char *name, xData_Int *n );
int xData_convertAttributeToDouble( statusMessageReporting *smr, xData_element *element, const char *name, double *d );
int xData_numberOfElementsByTagName( statusMessageReporting *smr, xData_element *element, const char *tagName );
xData_elementList *xData_getElementsByTagName( statusMessageReporting *smr, xData_element *element, const char *tagName );
xData_elementList *xData_getElementsByTagNameAndSort( statusMessageReporting *smr, xData_element *element, const char *tagName, 
    const char *sortAttributeName, xData_sortElementFunc sortFunction );
xData_element *xData_getOneElementByTagName( statusMessageReporting *smr, xData_element *element, char *name, int required );
void xData_freeElementList( statusMessageReporting *smr, xData_elementList *list );
int xData_addToAccessed( statusMessageReporting *smr, xData_element *element, int increment );
int xData_getAccessed( statusMessageReporting *smr, xData_element *element );

int xData_init_1d_x( statusMessageReporting *smr, xData_element *element );
int xData_is_1d_x( statusMessageReporting *smr, xDataType *xDT, int setMsg );
int xData_isElement_1d_x( statusMessageReporting *smr, xData_element *element, int setMsg );
int xData_1d_x_copyData( statusMessageReporting *smr, xData_element *element, xData_Int nAllocatedBytes, double *d );
double *xData_1d_x_allocateCopyData( statusMessageReporting *smr, xData_element *element );
int xData_1d_x_free_copyData( statusMessageReporting *smr, void *data );

int xData_init_2d_xindex_y( statusMessageReporting *smr, xData_element *element );
int xData_is_2d_xindex_y( statusMessageReporting *smr, xDataType *xDT, int setMsg );
int xData_isElement_2d_xindex_y( statusMessageReporting *smr, xData_element *element, int setMsg );
xData_Int *xData_2d_xindex_y_rawIndices( statusMessageReporting *smr, xData_element *element );
int xData_2d_xindex_y_free_rawIndices( statusMessageReporting *smr, void *data );
double *xData_2d_xindex_y_toXYs( statusMessageReporting *smr, xData_element *element, double *Xs );
double *xData_2d_xindex_y_toFilledYs( statusMessageReporting *smr, xData_element *element, double *Xs );
int xData_2d_xindex_y_free_toFilledYs( statusMessageReporting *smr, void *data );
double *xData_2d_xindex_y_toFilledXYs( statusMessageReporting *smr, xData_element *element, double *Xs );

int xData_init_2d_xy( statusMessageReporting *smr, xData_element *element );
int xData_is_2d_xy( statusMessageReporting *smr, xDataType *xDT, int setMsg );
int xData_isElement_2d_xy( statusMessageReporting *smr, xData_element *element, int setMsg );
double *xData_2d_xy_allocateCopyData( statusMessageReporting *smr, xData_element *element, xData_Int *length );
int xData_2d_xy_free_copyData( statusMessageReporting *smr, void *data );

int xData_init_2d_xShared_yHistogram( statusMessageReporting *smr, xData_element *element );
int xData_is_2d_xShared_yHistogram( statusMessageReporting *smr, xDataType *xDT, int setMsg );
int xData_isElement_2d_xShared_yHistogram( statusMessageReporting *smr, xData_element *element, int setMsg );
double *xData_2d_xShared_yHistogram_copyData( statusMessageReporting *smr, xData_element *element, xData_Int *n );
int xData_2d_xShared_yHistogram_free_copyData( statusMessageReporting *smr, void *data );
double *xData_2d_xShared_yHistogram_toFilledXYs( xDataType *xDT, xData_Int nXs, double *Xs );

int xData_init_matrix( statusMessageReporting *smr, xData_element *element );
int xData_is_matrix( statusMessageReporting *smr, xDataType *xDT, int setMsg );
int xData_isElement_matrix( statusMessageReporting *smr, xData_element *element, int setMsg );
xData_matrix *xData_matrix_copyData( statusMessageReporting *smr, xData_element *element );
int xData_matrix_free_copyData( statusMessageReporting *smr, void *data );
int getRowStartEndAtIndex( statusMessageReporting *smr, xDataType *xDT, xData_Int index, xData_Int *row, xData_Int *start, xData_Int *end );

int xData_is_xDataType( statusMessageReporting *smr, xDataType *xDT, char const * const type, int setMsg );
char const *xData_getFileName( xData_document *doc );
int xData_setFileName( statusMessageReporting *smr, xData_document *doc, char const *fileName );
xData_document *xData_getElementsDocument( xData_element *element );
void *xData_get_smrUserInterfaceFromDocument( xData_document *doc );
void *xData_get_smrUserInterfaceFromElement( xData_element *element );

int xData_stringTo_xData_Int( statusMessageReporting *smr, void *smrUserInterface, char const *c, xData_Int *value, char const *endings, char **e );
int xData_stringTo_double( statusMessageReporting *smr, void *smrUserInterface, char const *c, double *value, char const *endings, char **e );

/*
* Stuff in xDataMisc.c
*/
void *xData_malloc( statusMessageReporting *smr, size_t size, int zero, const char *forItem, const char *file, int line );
void *xData_realloc( statusMessageReporting *smr, void *pOld, size_t size, const char *forItem, const char *routine, int line );
void *xData_free( statusMessageReporting *smr, void *p );
char *xDataMisc_allocateCopyString( statusMessageReporting *smr, const char *s, const char *forItem, const char *routine, int line );
char *xDataMisc_getAbsPath( statusMessageReporting *smr, const char *fileName );
int xData_setMessageError_ReturnInt( int value, statusMessageReporting *smr, void *userData, const char *file, int line, int code, const char *fmt, ... );

#define xData_malloc2( smr, size, zero, forItem ) xData_malloc( smr, size, zero, forItem, __FILE__, __LINE__ )
#define xData_realloc2( smr, old, size, forItem ) xData_realloc( smr, old, size, forItem, __FILE__, __LINE__ )
#define xDataMisc_allocateCopyString2( smr, s, forItem ) xDataMisc_allocateCopyString( smr, s, forItem, __FILE__, __LINE__ )

#if defined __cplusplus
    }
    }
#endif

#endif              /* End of xData_h_included. */
