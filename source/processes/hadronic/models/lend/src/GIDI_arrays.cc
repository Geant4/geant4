/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <algorithm>

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

namespace Array {

/* *********************************************************************************************************//**
 * Helper function to parse the *shape* attribute for an Array instance.
 * 
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the shape. This must be a reference to the array node.
 *
 * @return                              Returns a *std::vector<int>* of the shape.
 ***********************************************************************************************************/

static std::vector<std::size_t> parseArrayShape( HAPI::Node const &a_node ) {

    std::string shape = a_node.attribute_as_string( GIDI_shapeChars );
    std::vector<std::size_t> shapeInts;

    auto shapeItems = LUPI::Misc::splitString( shape, ',', true );
    for( auto iter = shapeItems.begin( ); iter != shapeItems.end( ); ++iter ) {
        shapeInts.push_back( static_cast<std::size_t>( atoi( (*iter).c_str( ) ) ) );
    }

    return( shapeInts );
}

/*! \class FullArray
 * The class for storing any **GNDS** array as a flattened full array. It is the array returned by the 
 * Array.constructArray method. Note, this class is never used to store a **GNDS** array node, but to
 * provide a convenient way to access data in any of the ways an array can be stored compactly in **GNDS**.
 */

/* *********************************************************************************************************//**
 * Constructs a FullArray using the arguments *a_shape*. The values for *m_flattenedValues* are set to 0.0.
 *
 * @param a_shape               [in]    The shape of the full array.
 ***********************************************************************************************************/

FullArray::FullArray( std::vector<std::size_t> const &a_shape ) :
        m_shape( a_shape ) {

    std::size_t size = 1;
    for( auto iter = a_shape.begin( ); iter != a_shape.end( ); ++iter ) size *= *iter;
    m_flattenedValues.resize( size, 0.0 );
}

/* *********************************************************************************************************//**
 * Constructs a FullArray using the arguments *a_shape* and *a_flattenedValues*.
 *
 * @param a_shape               [in]    The shape of the full array.
 * @param a_flattenedValues     [in]    The values of the array. Its size must be the same has that specified by *a_shape*.
 ***********************************************************************************************************/

FullArray::FullArray( std::vector<std::size_t> const &a_shape, std::vector<double> const &a_flattenedValues ) :
        m_shape( a_shape ) {

    std::size_t size = 1;
    for( auto iter = a_shape.begin( ); iter != a_shape.end( ); ++iter ) size *= *iter;

    if( size != static_cast<std::size_t>( a_flattenedValues.size( ) ) ) throw Exception( "FullArray::FullArray: a_shape and a_flattenedValues are not the same size." );

    m_flattenedValues.resize( size );
    for( std::size_t index = 0; index < size; ++index ) m_flattenedValues[index] = a_flattenedValues[index];
}

/*! \class Array
 * The class for storing any **GNDS** array.
 */

/* *********************************************************************************************************//**
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Array2d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_useSystem_strtod    [in]    Flag passed to the function nfu_stringToListOfDoubles.
 ***********************************************************************************************************/

Array::Array( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod ) :
        Form( a_node, a_setupInfo, FormType::array ),
        m_shape( parseArrayShape( a_node ) ),
        m_compression( a_node.attribute_as_string( GIDI_compressionChars ) ),
        m_symmetry( a_node.attribute_as_string( GIDI_symmetryChars ) ),
        m_permutation( a_node.attribute_as_string( GIDI_permutationChars ) ),
        m_storageOrder( a_node.attribute_as_string( GIDI_storageOrderChars ) ) {

    if( m_compression == "" ) m_compression = GIDI_noneChars;

    if( ( m_compression == GIDI_noneChars ) || ( m_compression == GIDI_diagonalChars ) || ( m_compression == GIDI_flattenedChars ) ) {
        for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
            std::string name( child.name( ) );
            if( name != GIDI_valuesChars ) throw Exception( "For '" + moniker( ) + "' array, child node not 'values' but '" + name + "'." );

            std::string label1( child.attribute_as_string( GIDI_labelChars ) );
            if( label1 == "" ) {
                parseValuesOfDoubles( child, a_setupInfo, m_values, a_useSystem_strtod ); }
            else if( label1 == GIDI_startingIndices ) {
                if( m_compression != GIDI_diagonalChars ) throw Exception( "Invalid values' label '" + label1 + "'." );
                parseValuesOfInts( child, a_setupInfo, m_starts ); }
            else if( label1 == GIDI_startsChars ) {
                if( m_compression != GIDI_flattenedChars ) throw Exception( "Invalid values' label " + label1 + "." );
                parseValuesOfInts( child, a_setupInfo, m_starts ); }
            else if( label1 == GIDI_lengthsChars ) {
                if( m_compression != GIDI_flattenedChars ) throw Exception( "Invalid values' label '" + label1 + "'." );
                parseValuesOfInts( child, a_setupInfo, m_length ); }
            else {
                throw Exception( "Invalid values' label '" + label1 + "'." );
            }
        } }
    else {                                              // m_compression == GIDI_embeddedChars
        for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
            std::string name( child.name( ) );
            if( name != GIDI_arrayChars ) throw Exception( "For 'embedded' array, child node not 'array' but '" + name + "'." );

            m_array.push_back( new Array( child, a_setupInfo, a_useSystem_strtod ) );
        }
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Array::~Array( ) {

    for( auto iter = m_array.begin( ); iter != m_array.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Calculates the size of *this* which is the product of all the *m_shape* values.
 *
 * @return              The product of the *m_shape* values.
 ***********************************************************************************************************/

std::size_t Array::size( ) const {

    std::size_t size1 = 1;
    for( auto iter = m_shape.begin( ); iter != m_shape.end( ); ++iter ) size1 *= *iter;

    return( size1 );
}

/* *********************************************************************************************************//**
 * Returns an array that is a FullArray instance of *this*. Note, this does not currently handle symmetry, permutation
 * or storageOrder.
 *
 * @return                              Return a FullArray of *this*.
 ***********************************************************************************************************/

FullArray Array::constructArray( ) const {

    std::vector<double> values( size( ) );

    if( m_compression == GIDI_noneChars ) {
        values = m_values.vector( ); }
    else if( m_compression == GIDI_diagonalChars ) {
        throw Exception( "Array::constructArray: compression '" + m_compression + "' not supported." ); }
    else if( m_compression == GIDI_flattenedChars ) {
        std::size_t index = 0;
        for( std::size_t index1 = 0; index1 < m_starts.size( ); ++index1 ) {
            std::size_t length = static_cast<std::size_t>( m_length[index1] );
            std::size_t start = static_cast<std::size_t>( m_starts[index1] );

            for( std::size_t index2 = 0; index2 < length; ++index2, ++index ) values[start+index2] = m_values[index];
        } }
    else {
        throw Exception( "Array::constructArray: compression '" + m_compression + "' not supported." );
    }

    FullArray fullArray( m_shape, values );
    return( fullArray );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 * Currently does nothing.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Array::toXMLList( LUPI_maybeUnused GUPI::WriteInfo &a_writeInfo, LUPI_maybeUnused std::string const &a_indent ) const {

}

}           // End namespace Array.

/* *********************************************************************************************************//**
 * Function to parse a one-d flattened array.
 *
 * @param a_construction        [in]    Used to pass user options for parsing.
 * @param a_node                [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_data                [out]   An empty GIDI::Vector that is filled with the data.
 *
 * @return                              0 if successfull and 1 otherwise.
 ***********************************************************************************************************/

int parseFlattened1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Vector &a_data ) {

    FlattenedArrayData arrayData( a_node, a_setupInfo, 1, a_construction.useSystem_strtod( ) );

    std::size_t size = (std::size_t) arrayData.m_shape[0];
    a_data.resize( size );

    std::size_t n1 = 0, n2 = size;
    for( std::size_t i1 = 0; i1 < arrayData.m_numberOfStarts; ++i1 ) {
        std::size_t offset = (std::size_t) arrayData.m_starts[i1];
        for( int32_t i2 = 0; i2 < arrayData.m_lengths[i1]; ++i2, ++n1, ++offset ) {
            if( n1 >= arrayData.m_dValues.size( ) ) throw Exception( "Too many values in flattened array." );
            if( offset >= n2 ) throw Exception( "requested size is too small." );
            a_data[offset] = arrayData.m_dValues[n1];
        }
    }
    return( 0 );
}

/* *********************************************************************************************************//**
 * Function to parse a flattened array of dimension **a_dimensions**.
 *
 * @param a_node                [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_dimensions          [in]    The dimension of the flattened array to be parsed.
 * @param a_useSystem_strtod    [in]    Flag passed to the function nfu_stringToListOfDoubles.
 ***********************************************************************************************************/

FlattenedArrayData::FlattenedArrayData( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_dimensions, int a_useSystem_strtod ) :
        Form( a_node, a_setupInfo, FormType::flattenedArrayData ),
        m_numberOfStarts( 0 ), 
        m_numberOfLengths( 0 ),
        m_starts( ),
        m_lengths( ) {

    bool m_dValuesPresent( false );

    std::string shape( a_node.attribute_as_string( GIDI_shapeChars ) );
    auto shapeItems = LUPI::Misc::splitString( shape, ',', true );
    for( auto iter = shapeItems.begin( ); iter != shapeItems.end( ); ++iter ) {
        m_shape.push_back( static_cast<std::size_t>( atoi( (*iter).c_str( ) ) ) );
    }

    if( a_dimensions != (int) m_shape.size( ) ) throw Exception( "a_dimensions != m_shape.size( )" );

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );

        if( name == GIDI_valuesChars ) {
            std::string label( child.attribute_as_string( GIDI_labelChars ) );
            if( label == GIDI_startsChars ) {
                parseValuesOfInts( child, a_setupInfo, m_starts );
                m_numberOfStarts = (std::size_t) m_starts.size(); }
            else if( label == GIDI_lengthsChars ) {
              parseValuesOfInts( child, a_setupInfo, m_lengths );
              m_numberOfLengths = (std::size_t) m_lengths.size(); }
            else if( label == "" ) {
                m_dValuesPresent = true;
                parseValuesOfDoubles( child, a_setupInfo, m_dValues, a_useSystem_strtod ); }
            else {
                throw Exception( "unknown label for flatteded array sub-element" );
            } }
        else {
            throw Exception( "unknown flattened array sub-element" );
        }
    }
    if( m_starts.data() == nullptr ) throw Exception( "array missing starts element" );
    if( m_lengths.data() == nullptr ) throw Exception( "array missing lengths element" );
    if( !m_dValuesPresent ) throw Exception( "array missing dValues element" );
    if( m_numberOfStarts != m_numberOfLengths ) throw Exception( "m_numberOfStarts != m_numberOfLengths for array" );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

FlattenedArrayData::~FlattenedArrayData( ) {

//    smr_freeMemory2( m_starts );
//    smr_freeMemory2( m_lengths );
}

/* *********************************************************************************************************//**
 * Sets all elements in the range [*a_start*,*a_end*) to *a_value*.
 *
 * @param a_start       [in]    The starting flat-cell index of *this* to fill with *a_value*.
 * @param a_end         [in]    One after the last flat-cell index of *this* to fill with *a_value*.
 * @param a_value       [in]    The value to set each double in the range to.
 ***********************************************************************************************************/

void FlattenedArrayData::setToValueInFlatRange( LUPI_maybeUnused std::size_t a_start, std::size_t a_end, LUPI_maybeUnused double a_value ) {

    std::size_t size = 1;
    for( auto iter = m_shape.begin( ); iter != m_shape.end( ); ++iter ) size *= *iter;
    a_end = std::min( a_end, size );

    long numberOfValuesToSet = 0;
    for( std::size_t startIndex = 0; startIndex < m_numberOfStarts; ++startIndex ) {
        std::size_t start = static_cast<std::size_t>( m_starts[startIndex] );
        std::size_t length = static_cast<std::size_t>( m_lengths[startIndex] );

        if( ( start + length ) > a_end ) {
            if( a_end < start ) break;
            length = a_end - start;
        }
        numberOfValuesToSet += length;
    }

    for( long index = 0; index < numberOfValuesToSet; ++index ) m_dValues[index] = 0.0;
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void FlattenedArrayData::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::vector<int> ints;

    std::string shapeString;
    std::string sep = "";
    for( std::size_t i1 = 0; i1 < m_shape.size( ); ++i1 ) {
        shapeString += sep + intToString( static_cast<int>( m_shape[i1] ) );
        sep = ",";
    }

    std::string attributes = a_writeInfo.addAttribute( GIDI_shapeChars, shapeString ) + a_writeInfo.addAttribute( "compression", "flattened" );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );

    for( std::size_t i1 = 0; i1 < m_numberOfStarts; ++i1 ) ints.push_back( m_starts[i1] );
    intsToXMLList( a_writeInfo, indent2, ints, " valueType=\"Integer32\" label=\"starts\"" );

    ints.clear( );
    for( std::size_t i1 = 0; i1 < m_numberOfLengths; ++i1 ) ints.push_back( m_lengths[i1] );
    intsToXMLList( a_writeInfo, indent2, ints, " valueType=\"Integer32\" label=\"lengths\"" );

    doublesToXMLList( a_writeInfo, indent2, m_dValues.vector() );
    a_writeInfo.addNodeEnder( moniker( ) );
}

}           // End namespace GIDI.
