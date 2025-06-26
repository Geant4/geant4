/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <HAPI.hpp>
#include <LUPI.hpp>

namespace HAPI {

/*
=========================================================
 *
 * @return
 */
Node::Node() :
        m_node(NULL) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
Node::Node( Node_internal *a_node ) :
        m_node(a_node) {

}

/* *********************************************************************************************************//**
 * Copy constructor.
 *
 * @param a_node        [In]   The node to copy.
 ***********************************************************************************************************/

Node::Node( Node const &a_node ) :
        m_node( nullptr ) {

#ifdef HAPI_USE_PUGIXML
    if( a_node.m_node->type( ) == NodeInteralType::pugiXML ) {
        m_node = new PugiXMLNode( static_cast<PugiXMLNode const &>( *a_node.m_node ) );
    }
#endif
#ifdef HAPI_USE_HDF5
    if( a_node.m_node->type( ) == NodeInteralType::HDF5 ) {
        m_node = new HDFNode( static_cast<HDFNode const &>( *a_node.m_node ) );
    }
#endif
    if( m_node == nullptr ) throw LUPI::Exception( "Unsupported m_node type." );
}

/*
=========================================================
*/
Node::~Node( ) {

    delete m_node;

}
/*
============================================================
===================== get child element ====================
============================================================
 *
 * @return
 */
Node Node::child(char const *a_name) const {

    if (NULL == m_node)
        return Node();
    return Node( m_node->child( a_name ) );

}
/*
=========================================================
*/
Node Node::first_child() const {

    if (NULL == m_node)
        return Node();
    return Node( m_node->first_child( ) );

}
/*
============================================================
===================== get sibling element ==================
============================================================
 *
 * @return
 */
Node Node::next_sibling() const {

    if (NULL == m_node)
        return Node();
    Node_internal *sibling = m_node->next_sibling( );
    delete m_node;
    return Node( sibling );

}
/*
============================================================
============ update self to point to next sibling ==========
============================================================
 *
 * @return
 */
void Node::to_next_sibling() const {

    m_node->to_next_sibling( );

}
/*
============================================================
===================== assignment operator ==================
============================================================
 */
Node& Node::operator=(const Node &other) {

    if (NULL != other.m_node)
        this->m_node = other.m_node->copy();
    return *this;

}
/*
============================================================
===================== get tag name =========================
============================================================
 *
 * @return
 */
std::string Node::name() const {

    if (NULL == m_node)
        return std::string("");
    return m_node->name();

}
/*
============================================================
================== test for empty node =====================
============================================================
 *
 * @return
 */
bool Node::empty() const {

    if (NULL == m_node)
        return true;
    return m_node->empty();

}
/*
============================================================
======================= text data ==========================
============================================================
 *
 * @return
 */
Text Node::text() const {

    if (NULL == m_node)
        return Text();
    return m_node->text();

}
/*
============================================================
===================== numeric data =========================
============================================================
 *
 * @return
 */
Data Node::data() const {

    if (NULL == m_node)
        return Data();
    return Data( m_node->data() );

}

}
