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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xData.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

typedef struct xmlTextStruct_s {
    long allocated;
    long length;
    char *text;
} xmlTextStruct;

static int xData_parseReconstructXML2( xData_rootElement *root, xmlTextStruct *XML, char *text, long textLength );
static int addStringToXML( xmlTextStruct *XML, char *s, long len );
static void xData_parseOutlinePrintRoot( FILE *f, xData_rootElement *root, int outputText );
/*
************************************************************
*/
char *xData_parseReconstructXML( xData_document *xData_doc ) {

    int err;
    xmlTextStruct XML = { 0, 0, NULL };

    err = xData_parseReconstructXML2( &(xData_doc->root), &XML, NULL, 0 );
    //if( err == 0 ) addStringToXML( &XML, "\n", -1 );
    if( err == 0 ) addStringToXML( &XML, (char*) "\n", -1 );
    return( XML.text );
}
/*
************************************************************
*/
static int xData_parseReconstructXML2( xData_rootElement *root, xmlTextStruct *XML, char *text, long textLength ) {

    long i, textOffset = 0;
    xData_element *child;

    for( child = root->children; child != NULL; child = child->next ) {
        //if( textOffset < child->textOffset ) {
        if( textOffset < (int) child->textOffset ) {
            if( addStringToXML( XML, &(text[textOffset]), child->textOffset - textOffset ) ) return( 1 );
            textOffset = child->textOffset;
        }
        //if( addStringToXML( XML, "<", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, (char*) "<", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, child->name, -1 ) != 0 ) return( 1 );
        for( i = 0; i < child->attributes.number; i++ ) {
            //if( addStringToXML( XML, " ", -1 ) != 0 ) return( 1 );
            if( addStringToXML( XML, (char*) " ", -1 ) != 0 ) return( 1 );
            if( addStringToXML( XML, child->attributes.attributes[i].name, -1 ) != 0 ) return( 1 );
            //if( addStringToXML( XML, "=\"", -1 ) != 0 ) return( 1 );
            if( addStringToXML( XML, (char*) "=\"", -1 ) != 0 ) return( 1 );
            if( addStringToXML( XML, child->attributes.attributes[i].value, -1 ) != 0 ) return( 1 );
            //if( addStringToXML( XML, "\"", -1 ) != 0 ) return( 1 );
            if( addStringToXML( XML, (char*) "\"", -1 ) != 0 ) return( 1 );
        }
        //if( addStringToXML( XML, ">", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, (char*) ">", -1 ) != 0 ) return( 1 );
        if( xData_parseReconstructXML2( &(child->childrenRoot), XML, child->text.text, child->text.length ) != 0 ) return( 1 );
        //if( addStringToXML( XML, "</", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, (char*) "</", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, child->name, -1 ) != 0 ) return( 1 );
        //if( addStringToXML( XML, ">", -1 ) != 0 ) return( 1 );
        if( addStringToXML( XML, (char*) ">", -1 ) != 0 ) return( 1 );
    }
    if( textOffset < textLength ) if( addStringToXML( XML, &(text[textOffset]), textLength - textOffset ) ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
static int addStringToXML( xmlTextStruct *XML, char *s, long len ) {

    long lenS, length, inc;
    char *p;

    if( len >= 0 ) {
        lenS  = len; }
    else  {
        lenS = strlen( s );
    }
    length = XML->length + lenS + 1;
    if( XML->allocated < length ) {
        inc = ( 140 * XML->allocated ) / 100;
        if( inc < 10000 ) inc = 10000;
        if( length < inc ) length = inc;
        XML->text = (char *) realloc( XML->text, length );
        XML->allocated  = length;
        if( XML->text == NULL ) return( 1 );
    }
    p = &(XML->text[XML->length]);
    if( len >= 0 ) {
        strncpy( p, s, len ); }
    else {
        strcpy( p, s );
    }
    XML->length += lenS;
    return( 0 );
}
/*
************************************************************
*/
int xData_parseOutline( FILE *f, xData_document *xData_doc, int outputText ) {

    xData_parseOutlinePrintRoot( f, &(xData_doc->root), outputText );
    return( 0 );
}
/*
************************************************************
*/
static void xData_parseOutlinePrintRoot( FILE *f, xData_rootElement *root, int outputText ) {

    int i, depth = root->depth;
    xData_element *child;

    for( child = root->children; child != NULL; child = child->next ) {
        for( i = 0; i < depth; i++ ) fprintf( f, "  " );
        fprintf( f, "name = %s at line = %ld column = %ld textOffset = %ld; attris:", 
            //child->name, child->docInfo.line, child->docInfo.column, child->textOffset );
            child->name, (long int)child->docInfo.line, (long int)child->docInfo.column, (long int)child->textOffset );
        for( i = 0; i < child->attributes.number; i++ ) {
            fprintf( f, " %s = \"%s\"", child->attributes.attributes[i].name, child->attributes.attributes[i].value );
        }
        if( outputText && child->text.text != NULL ) fprintf( f, "; text: %s", child->text.text );
        fprintf( f, "\n" );
        xData_parseOutlinePrintRoot( f, &(child->childrenRoot), outputText );
    }
}

#if defined __cplusplus
}
#endif
