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

/*! \class Array3d
 * Class to store a 3d array.
 */

/* *********************************************************************************************************//**
 *
 * @param a_node                [in]    The **HAPI::Node** to be parsed and used to construct the Array3d.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_useSystem_strtod    [in]    Flag passed to the function nfu_stringToListOfDoubles.
 ***********************************************************************************************************/

Array3d::Array3d( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod ) :
        Form( a_node, a_setupInfo, FormType::array3d ),
        m_array( a_node, a_setupInfo, 3, a_useSystem_strtod ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Array3d::~Array3d( ) {

}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the ProtareSingle that
 * contains the TNSL data covers the lower energy multi-group data.
 *
 * @param a_maxTNSL_index           [in]    All elements up to "row" *a_maxTNSL_index* exclusive are zero-ed.
 ***********************************************************************************************************/

void Array3d::modifiedMultiGroupElasticForTNSL( int a_maxTNSL_index ) {

    std::vector<int> const &m_shape = m_array.shape( );
    int maxFlatIndex = a_maxTNSL_index * m_shape[1] * m_shape[2];

    m_array.setToValueInFlatRange( 0, maxFlatIndex, 0.0 );
}


/* *********************************************************************************************************//**
 * Returns the matrix that represents the specified 3rd dimension. That is the matrix M[i][j] for all i, j of A2d[i][j][*a_index*].
 * This is mainly used for multi-group, Legendre expanded transfer matrices where a specific Legendre order is requested. This is,
 * the matrix represent the *energy_in* as rows and the *energy_outp* as columns for a specific Legendre order.
 *
 * @param a_index           [in]     The requested *index* for the 3rd dimension.
 ***********************************************************************************************************/

Matrix Array3d::matrix( std::size_t a_index ) const {

    if( size( ) <= a_index ) {
        Matrix matrix( 0, 0 );
        return( matrix );
    }

    std::size_t numberOfOrders = m_array.m_shape[2], rows = m_array.m_shape[0], columns = m_array.m_shape[1];
    Matrix matrix( rows, columns );

    std::size_t lengthSum = 0;
    for( std::size_t i1 = 0; i1 < m_array.m_numberOfStarts; ++i1 ) {
        std::size_t start = m_array.m_starts[i1];
        std::size_t length = m_array.m_lengths[i1];

        std::size_t energyInIndex = start / ( numberOfOrders * columns );
        std::size_t energyOutIndex = start % ( numberOfOrders * columns );
        std::size_t orderIndex = energyOutIndex % numberOfOrders;
        energyOutIndex /= numberOfOrders;

        std::size_t step = a_index - orderIndex;
        if( orderIndex > a_index ) {
            ++energyOutIndex;
            if( energyOutIndex >= columns ) {
                energyOutIndex = 0;
                ++energyInIndex;
            }
            step += numberOfOrders;
        }
        std::size_t dataIndex = lengthSum + step;
        for( ; step < length; step += numberOfOrders ) {
            matrix.set( energyInIndex, energyOutIndex, m_array.m_dValues[dataIndex] );
            ++energyOutIndex;
            if( energyOutIndex >= columns ) {
                energyOutIndex = 0;
                ++energyInIndex;
            }
            dataIndex += numberOfOrders;
        }
        lengthSum += length;
    }

    return( matrix );
}

}
