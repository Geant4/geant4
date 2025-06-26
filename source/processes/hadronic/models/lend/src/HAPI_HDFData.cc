/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"

#ifdef HAPI_USE_HDF5
namespace HAPI {

/*
=========================================================
 *
 * @return
 */
HDFData::HDFData() :
        m_node_id(-1),
        m_dataspace_id(-1),
        m_length(-1) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
HDFData::HDFData( hid_t node_id ) :
        m_node_id(node_id) {

    m_dataspace_id = H5Dget_space(m_node_id);
    m_length = H5Sget_simple_extent_npoints(m_dataspace_id);

}
/*
=========================================================
*/
HDFData::~HDFData( ) {

}

int HDFData::length( ) const {

    return m_length;
}

void HDFData::getDoubles(nf_Buffer<double> &buffer)
{
    buffer.resize(m_length);
    H5Dread(m_node_id, H5T_NATIVE_DOUBLE, H5S_ALL, m_dataspace_id, H5P_DEFAULT, buffer.data());
}

void HDFData::getInts(nf_Buffer<int> &buffer)
{
    buffer.resize(m_length);
    H5Dread(m_node_id, H5T_NATIVE_INT, H5S_ALL, m_dataspace_id, H5P_DEFAULT, buffer.data());
}

}
#endif
