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
HDFFile::HDFFile() :
        m_name( "" ),
        m_doc( 0 ),
        m_doc_as_node(nullptr){

}
/*
=========================================================
 *
 * @param filename
 * @return
 */
HDFFile::HDFFile(char const *filename) :
        m_name( filename ) {

    m_doc = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
    m_doc_as_node = new HDFNode( m_doc );

}
/*
=========================================================
*/
HDFFile::~HDFFile( ) {

    H5Fclose(m_doc);
    delete m_doc_as_node;

}
/*
============================================================
===================== get child element ====================
============================================================
 *
 * @return
 */
Node HDFFile::child(char const *a_name) {

    return Node( m_doc_as_node->child(a_name) );

}
/*
============================================================
===================== get first child node =================
============================================================
*/
Node HDFFile::first_child() {

    return Node( m_doc_as_node->first_child() );

}
/*
=========================================================
*/
std::string HDFFile::name() const {

    return m_name;

}

}
#endif
