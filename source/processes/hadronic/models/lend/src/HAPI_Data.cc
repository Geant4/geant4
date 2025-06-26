/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"

namespace HAPI {

/*
=========================================================
 *
 * @return
 */
Data::Data() :
        m_data(NULL) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
Data::Data( Data_internal *a_data ) :
        m_data(a_data) {

}
/*
=========================================================
*/
Data::~Data( ) {

    delete m_data;

}

int Data::length( ) const {

    return m_data->length();

}

void Data::getDoubles(nf_Buffer<double> &buffer)
{

    m_data->getDoubles(buffer);

}

void Data::getInts(nf_Buffer<int> &buffer)
{

    m_data->getInts(buffer);

}

}
