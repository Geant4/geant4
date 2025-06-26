/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <algorithm>

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

static std::size_t startIndexAttribute( HAPI::Node const &a_node );

/* *********************************************************************************************************//**
 * This function searchs the list of ascending values *a_Xs* for the two values that bound *a_x* using a bi-section search.
 * If *a_x* is less than the first value, -2 is returned. If *a_x* is greater than the last value, -1 is returned.
 * Otherwise, the returned index will be such that *a_Xs*[index] <= *a_x* < *a_Xs*[index+1].
 *
 * @param a_x       [in]    The value to search.
 * @param a_Xs      [in]    The list of ascending values to 
 *
 * @return          [in]    The index within the *a_Xs* list that bounds *a_x*.
 ***********************************************************************************************************/

long binarySearchVector( double a_x, std::vector<double> const &a_Xs ) {
/*
*   Returns -2 is a_x < first point of a_Xs, -1 if > last point of a_Xs, and the lower index of a_Xs otherwise.
*/
    long size = a_Xs.size( );
    long imin = 0, imid, imax = size - 1;

    if( a_x < a_Xs[0] ) return( -2 );
    if( a_x > a_Xs[size-1] ) return( -1 );
    while( 1 ) {
        imid = ( imin + imax ) >> 1;
        if( imid == imin ) break;
        if( a_x < a_Xs[imid] ) {
            imax = imid; }
        else {
            imin = imid;
        }
    }
    return( imin );
}

/* *********************************************************************************************************//**
 * Adds the list of integers to the list of XML lines in *a_writeInfo*.
 *
 * @param a_writeInfo           [in/out]    Instance containing incremental indentation, values per line and other information and stores the appended lines.
 * @param a_indent              [in]        The amount to indent *this* node.
 * @param a_values              [in]        The list of integers to convert to strings and add to *a_writeInfo*.
 * @param a_attributes          [in]        String representation of the attributes for the GNDS **values** node.
 ***********************************************************************************************************/

void intsToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, std::vector<int> a_values, std::string const &a_attributes ) {

    a_writeInfo.addNodeStarter( a_indent, GIDI_valuesChars, a_attributes );

    std::string intString;
    std::string sep( "" );

    for( std::size_t i1 = 0; i1 < a_values.size( ); ++i1 ) {
        intString += sep + intToString( a_values[i1] );
        if( i1 == 0 ) sep = a_writeInfo.m_sep;
    }

    a_writeInfo.m_lines.back( ) += intString;
    a_writeInfo.addNodeEnder( GIDI_valuesChars );
}

/* *********************************************************************************************************//**
 * This function converts the text of a **HAPI::Node** into a list of doubles.
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** node whose text is to be converted into a list of doubles.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_values              [in]    The list to fill with the converted values.
 ***********************************************************************************************************/

void parseValuesOfDoubles( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
            nf_Buffer<double> &a_values) {

    parseValuesOfDoubles( a_node, a_setupInfo, a_values, a_construction.useSystem_strtod( ) );
}

/* *********************************************************************************************************//**
 * This function converts the text of a **HAPI::Node** into a list of doubles.
 *
 * @param a_node                [in]    The **HAPI::Node** node whose text is to be converted into a list of doubles.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_values              [in]    The list to fill with the converted values.
 * @param a_useSystem_strtod    [in]    Flag passed to the function nfu_stringToListOfDoubles.
 ***********************************************************************************************************/

void parseValuesOfDoubles( HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<double> &a_values, LUPI_maybeUnused int a_useSystem_strtod ) {

    std::string href = a_node.attribute_as_string( GIDI_hrefChars );

    if( href != "" ) {
        std::size_t startIndex = startIndexAttribute( a_node );
        std::size_t count = a_node.attribute_as_long( GIDI_countChars );
        if( a_setupInfo.m_protare->dataManager( ) == nullptr )
            throw Exception( "parseValuesOfDoubles: Cannot read from HDF5 file as GIDI+ was compiled without HDF5 support." );
        a_setupInfo.m_protare->dataManager( )->getDoubles( a_values, startIndex, startIndex + count ); }
    else {
        HAPI::Data data = a_node.data( );
        data.getDoubles( a_values ); // FIXME overload getDoubles() to take std::vector argument, avoid extra copy?
    }
/*
    int64_t numberConverted = p1.size( );
    a_values.resize( numberConverted );
    for( int64_t i1 = 0; i1 < numberConverted; ++i1 ) a_values[i1] = p1[i1];
*/
}

/* *********************************************************************************************************//**
 * This function converts the text of a **HAPI::Node** into a list of ints.
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The **HAPI::Node** node whose text is to be converted into a list of ints.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_values              [in]    The list to fill with the converted values.
 ***********************************************************************************************************/

void parseValuesOfInts( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<int> &a_values) {

    parseValuesOfInts( a_node, a_setupInfo, a_values );
}

/* *********************************************************************************************************//**
 * This function converts the text of a **HAPI::Node** into a list of ints.
 *
 * @param a_node                [in]    The **HAPI::Node** node whoses text is to be converted into a list of ints.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_values              [in]    The list to fill with the converted values.
 ***********************************************************************************************************/

void parseValuesOfInts( HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<int> &a_values ) {

    std::string href = a_node.attribute_as_string( GIDI_hrefChars );

    if( href != "" ) {
        std::size_t startIndex = startIndexAttribute( a_node );
        std::size_t count = a_node.attribute_as_long( GIDI_countChars );
        if( a_setupInfo.m_protare->dataManager( ) == nullptr )
            	throw Exception( "parseValuesOfInts: Cannot read from HDF5 file as GIDI+ was compiled without HDF5 support." );
        a_setupInfo.m_protare->dataManager( )->getInts( a_values, startIndex, startIndex + count ); }
    else {
        HAPI::Data data = a_node.data( );
        data.getInts( a_values ); // FIXME overload getDoubles() to take std::vector argument, avoid extra copy?
    }

/*
    int64_t numberConverted = p1.size( );
    a_values.resize( numberConverted );
    for( int64_t i1 = 0; i1 < numberConverted; ++i1 ) a_values[i1] = p1[i1];
*/
//  a_values.swap( p1 );
}

/* *********************************************************************************************************//**
 * Adds the list of doubles to the list of XML lines in *a_writeInfo*.
 *
 * @param a_writeInfo           [in/out]    Instance containing incremental indentation, values per line and other information and stores the appended lines.
 * @param a_indent              [in]        The amount to indent *this* node.
 * @param a_values              [in]        The list of doubles to convert to strings and add to *a_writeInfo*.
 * @param a_start               [in]        The value for the *start* attribute.
 * @param a_newLine             [in]        If *false*, the first *a_writeInfo.m_valuesPerLine* values are added to the last line with no indentation; otherwise, they are put on a new line with indentation.
 * @param a_valueType           [in]        The value for the *valueType* attribute.
 ***********************************************************************************************************/

void doublesToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, std::vector<double> a_values, std::size_t a_start, bool a_newLine, std::string const &a_valueType ) {

    int valuesPerLine( a_writeInfo.m_valuesPerLine );
    std::string indent( a_indent );
    std::string attributes;
    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string XMLLine;
    std::string sep = "";

    if( !a_newLine ) indent = "";
    if( a_valueType != "" ) attributes += a_writeInfo.addAttribute( GIDI_valueTypeChars, a_valueType );
    if( a_start != 0 ) attributes += a_writeInfo.addAttribute( GIDI_startChars, size_t_ToString( a_start ) );
    XMLLine = a_writeInfo.nodeStarter( indent, GIDI_valuesChars, attributes );

    if( valuesPerLine < 1 ) valuesPerLine = 1;
    int numberOfValuesInLine = 0;
    for( std::size_t i1 = 0; i1 < a_values.size( ); ++i1 ) {
        XMLLine += sep + LUPI::Misc::doubleToShortestString( a_values[i1] );
        sep = a_writeInfo.m_sep;
        ++numberOfValuesInLine;

        if( numberOfValuesInLine == valuesPerLine ) {
            if( a_newLine ) {
                a_writeInfo.push_back( XMLLine ); }
            else {
                a_writeInfo.m_lines.back( ) += XMLLine;
            }
            numberOfValuesInLine = 0;
            XMLLine.clear( );
            XMLLine = indent2;
            a_newLine = true;
            sep = "";
        }
    }
    if( numberOfValuesInLine > 0 ) {
        if( a_newLine ) {
            a_writeInfo.push_back( XMLLine ); }
        else {
            a_writeInfo.m_lines.back( ) += XMLLine;
        } }
    else if( a_values.size( ) == 0 ) {
        a_writeInfo.push_back( XMLLine );
    }

    a_writeInfo.addNodeEnder( GIDI_valuesChars );
}

/* *********************************************************************************************************//**
 * This function returns an frame enum representing a **HAPI::Node**'s attribute with name *a_name*.
 *
 * @param a_node        [in]    The **HAPI::Node** node whoses attribute named *a_node* is to be parsed to determine the frame.
 * @param a_setupInfo   [in]    Information create my the Protare constructor to help in parsing.
 * @param a_name        [in]    The name of the attribute to parse.
 *
 * @return                      The *frame* enum representing the node's frame.
 ***********************************************************************************************************/

Frame parseFrame( HAPI::Node const &a_node, LUPI_maybeUnused SetupInfo &a_setupInfo, std::string const &a_name ) {

    Frame frame = Frame::lab;
    if( strcmp( a_node.attribute_as_string( a_name.c_str( ) ).c_str( ), GIDI_centerOfMassChars ) == 0 ) frame = Frame::centerOfMass;
    return( frame );
}

/* *********************************************************************************************************//**
 * This function converts the y-values from the Gridded1d into a Ys1d instance.
 *
 * @param a_function1d  [in]    The Gridded1d whoses y-values are converted into a Ys1d instance.
 *
 * @return                      A Ys1d instance of the y-values.
  ***********************************************************************************************************/

Functions::Ys1d gridded1d2GIDI_Ys1d( Functions::Function1dForm const &a_function1d ) {

    std::vector<double> ys;
    Functions::Ys1d ys1d( a_function1d.axes( ), ptwXY_interpolationFlat, 0, ys );

    switch( a_function1d.type( ) ) {
    case FormType::gridded1d :
        {
            Functions::Gridded1d const &gridded1d = static_cast<Functions::Gridded1d const &>( a_function1d );
            Vector const &data = gridded1d.data( );
            std::size_t start = 0;

            for( ; start < data.size( ); ++start ) {
                if( data[start] != 0 ) break;
            }
            ys1d.setStart( start );

            for( std::size_t i1 = start; i1 < data.size( ); ++i1 ) ys1d.push_back( data[i1] );
        }
        break;
    default :
        throw Exception( "gridded1d2GIDI_Ys1d: unsupported 1d function type " + a_function1d.label( ) );
    }

    return( ys1d );
}

/* *********************************************************************************************************//**
 * This function converts the values of a Vector into a Ys1d instance.
 *
 * @param a_axes        [in]    The Axes for the returned Ys1d instance.
 * @param a_vector      [in]    The Vector whoses values are converted into a Ys1d instance.
 *
 * @return                      A Ys1d instance of the values.
  ***********************************************************************************************************/

Functions::Ys1d vector2GIDI_Ys1d( Axes const &a_axes, Vector const &a_vector ) {

    std::size_t start = 0;
    for( ; start < a_vector.size( ); ++start ) {
        if( a_vector[start] != 0 ) break;
    }

    std::vector<double> ys;
    Functions::Ys1d ys1d( a_axes, ptwXY_interpolationLinLin, start, ys );

    for( std::size_t i1 = start; i1 < a_vector.size( ); ++i1 ) ys1d.push_back( a_vector[i1] );

    return( ys1d );
}

/* *********************************************************************************************************//**
 * This function converts an integer gid value (i.e., group id) into the LLNL legacy bdfls label.
 *
 * @param a_gid         [in]    The integer gid used to construct the LLNL legacy bdfls label.
 *
 * @return                      The LLNL legacy bdfls label.
  ***********************************************************************************************************/

std::string LLNL_gidToLabel( int a_gid ) {

    return( LUPI::Misc::argumentsToString( "LLNL_gid_%d", a_gid ) );
}

/* *********************************************************************************************************//**
 * This function converts an integer fid value (i.e., flux id) into the LLNL legacy bdfls label.
 *
 * @param a_fid         [in]    The integer fid used to construct the LLNL legacy bdfls label.
 *
 * @return                      The LLNL legacy bdfls label.
  ***********************************************************************************************************/

std::string LLNL_fidToLabel( int a_fid ) {

    return( LUPI::Misc::argumentsToString( "LLNL_fid_%d", a_fid ) );
}

/* *********************************************************************************************************//**
 * This function returns an instance of *std::vector<std::string>* with only a_string as an item.
 * 
 * @param a_string      [in]    The string to add to the returned *std::vector<std::string>* instance.
 *
 * @return                      A *std::vector<std::string>* instance.
  ***********************************************************************************************************/

std::vector<std::string> vectorOfStrings( std::string const &a_string ) {
    std::vector<std::string> vectorOfStrings1;

    vectorOfStrings1.push_back( a_string );
    return( vectorOfStrings1 );
}

/* *********************************************************************************************************//**
 * This function returns a sorted instance of the strings in *a_strings*.
 *
 * @param a_strings             [in]    The string to add to the returned *std::vector<std::string>* instance.
 * @param a_orderIsAscending    [in]    If *true* the strings are sorted in ascending order; otherwise, descending order.
 *
 * @return                      A *std::vector<std::string>* instance.
  ***********************************************************************************************************/

std::vector<std::string> sortedListOfStrings( std::vector<std::string> const &a_strings, bool a_orderIsAscending ) {

    std::vector<std::string> keys( a_strings );

    std::sort( keys.begin( ), keys.end( ) );

    if( a_orderIsAscending ) return( keys );

    std::vector<std::string> keys2;

    for( std::vector<std::string>::reverse_iterator iter = keys.rbegin( ); iter != keys.rend( ); ++iter ) keys2.push_back( *iter );
    return( keys2 );
}

/* *********************************************************************************************************//**
 * This function returns a std::string representation of a *frame*.
 *
 * @param a_frame       [in]    The frame to convert to a string.
 *
 * @return                      A *std::string* instance.
  ***********************************************************************************************************/

std::string frameToString( Frame a_frame ) {

    if( a_frame == Frame::lab ) return( GIDI_labChars );
    return( GIDI_centerOfMassChars );
}

/* *********************************************************************************************************//**
 * Create the XML list for xs, pdf or cdf for an Xs_pdf_cdf1d instance.
 *
 * @param a_writeInfo       [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param a_nodeName        [in]        The name of the node (e.g., "xs" );
 * @param a_values          [in]        The list of doubles to wrap.
 *
 * @return                      A *std::string* instance.
  ***********************************************************************************************************/

std::string nodeWithValuesToDoubles( GUPI::WriteInfo &a_writeInfo, std::string const &a_nodeName, std::vector<double> const &a_values ) {

    std::string xml = a_writeInfo.nodeStarter( "", a_nodeName );
    std::string sep( "" );

    xml += a_writeInfo.nodeStarter( "", GIDI_valuesChars );
    for( std::size_t i1 = 0; i1 < a_values.size( ); ++i1 ) {
        xml += sep + LUPI::Misc::doubleToShortestString( a_values[i1] );
        if( i1 == 0 ) sep = " ";
    }
    xml += a_writeInfo.nodeEnder( GIDI_valuesChars );
    xml += a_writeInfo.nodeEnder( a_nodeName );

    return( xml );
}

/* *********************************************************************************************************//**
 * Returns a string representation of int *a_value*.
 *
 * @param a_value               [in]        The int value to convert to a string.
 *
 * @return                      A *std::string* instance.
  ***********************************************************************************************************/

std::string intToString( int a_value ) {

    return( LUPI::Misc::argumentsToString( "%d", a_value ) );
}

/* *********************************************************************************************************//**
 * Returns a string representation of std::size_t *a_value*.
 *
 * @param a_value               [in]        The std::size value to convert to a string.
 *
 * @return                      A *std::string* instance.
  ***********************************************************************************************************/

std::string size_t_ToString( std::size_t a_value ) {

    return( LUPI::Misc::argumentsToString( "%zu", a_value ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_moniker           [in]        The moniker for the energy type.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_function          [in]        The energy function whose information is converted to XML.
 ***********************************************************************************************************/

void energy2dToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_moniker, std::string const &a_indent, Functions::Function1dForm *a_function ) {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    if( a_function == nullptr ) return;

    a_writeInfo.addNodeStarter( a_indent, a_moniker, "" );
    a_function->toXMLList( a_writeInfo, indent2 );
    a_writeInfo.addNodeEnder( a_moniker );
}

/* *********************************************************************************************************//**
 * Returns a new **ExcludeReactionsSet** with 
 *
 * @param a_protare     [in]            The **Protare** instance used to determine the number of reactions to adjust the new indices by.
 *
 * @return                              returns the startIndex attribute of *a_node*.
 ***********************************************************************************************************/

void excludeReactionsSetAdjust( ExcludeReactionsSet a_excludeReactionsSet, Protare const &a_protare ) {

    ExcludeReactionsSet excludeReactionsSet;

    for( auto iter = a_excludeReactionsSet.begin( ); iter != a_excludeReactionsSet.end( ); ++iter ) {
        int index = (*iter) - a_protare.numberOfReactions( );
        if( index > -1 ) excludeReactionsSet.insert( index );
    }

    a_excludeReactionsSet = excludeReactionsSet;
}

/* *********************************************************************************************************//**
 * For internal use only.
 *
 * @param a_node                [in]    The **HAPI::Node** node whose text is to be converted into a list of doubles.
 *
 * @return                              returns the startIndex attribute of *a_node*.
 ***********************************************************************************************************/

static std::size_t startIndexAttribute( HAPI::Node const &a_node ) {

    std::size_t startIndex = 0;

    std::string attribute = a_node.attribute_as_string( GIDI_startIndexChars );
    if( attribute != "" ) {
        startIndex = a_node.attribute_as_long( GIDI_startIndexChars ); }
    else {
        attribute = a_node.attribute_as_string( GIDI_offsetChars );
        if( attribute != "" ) startIndex = a_node.attribute_as_long( GIDI_offsetChars );
    }

    return( startIndex );
}

}               // End namespace GIDI.
