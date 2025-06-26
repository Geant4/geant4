/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "HAPI.hpp"

#ifdef HAPI_USE_PUGIXML
namespace HAPI {

/*
=========================================================
 *
 * @return
 */
PugiXMLNode::PugiXMLNode() :
        Node_internal( NodeInteralType::pugiXML ),
        m_node( pugi::xml_node() ) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
PugiXMLNode::PugiXMLNode( pugi::xml_node a_node ) :
        Node_internal( NodeInteralType::pugiXML ),
        m_node( a_node ) {

}
/*
============================================================
====================== copy constructor ====================
============================================================
 */
PugiXMLNode::PugiXMLNode(const PugiXMLNode &other) :
        Node_internal( other ),
        m_node( other.m_node ) {

}
/*
=========================================================
*/
PugiXMLNode::~PugiXMLNode( ) {

}
/*
============================================================
=================== get attribute by name ==================
============================================================
 *
 * @param a_name
 * @return
 */
std::string PugiXMLNode::attribute(char const *a_name) {

    pugi::xml_attribute attr = m_node.attribute( a_name );

    return std::string(attr.value( ));

}


int PugiXMLNode::attribute_as_int(const char* a_name){
  pugi::xml_attribute attr = m_node.attribute( a_name );

  return atoi(attr.value( ));
}
long PugiXMLNode::attribute_as_long(const char* a_name){
  pugi::xml_attribute attr = m_node.attribute( a_name );

  return atol(attr.value( ));
}
double PugiXMLNode::attribute_as_double(const char* a_name){
  pugi::xml_attribute attr = m_node.attribute( a_name );

  return atof(attr.value( ));
}

/*
============================================================
===================== get child element ====================
============================================================
 *
 * @return
 */
Node_internal *PugiXMLNode::child(char const *a_name) {

    return new PugiXMLNode( m_node.child( a_name ) );

}
/*
=========================================================
*/
Node_internal *PugiXMLNode::first_child() {

    return new PugiXMLNode( m_node.first_child( ) );

}
/*
============================================================
===================== get sibling element ==================
============================================================
 *
 * @return
 */
Node_internal *PugiXMLNode::next_sibling() {

    return new PugiXMLNode( m_node.next_sibling( ) );

}
/*
============================================================
============= update self to point to next sibling =========
============================================================
 *
 * @return
 */
void PugiXMLNode::to_next_sibling() {

    m_node = m_node.next_sibling( );

}
/*
============================================================
======================== make a copy =======================
============================================================
 *
 * @return
 */
Node_internal *PugiXMLNode::copy() {

    return new PugiXMLNode( m_node );

}
/*
============================================================
===================== assignment operator ==================
============================================================
 */
Node_internal &PugiXMLNode::operator=(const PugiXMLNode &other) {

    this->m_node = other.m_node;
    return *this;

}
/*
============================================================
===================== get tag name =========================
============================================================
 *
 * @return
 */
std::string PugiXMLNode::name() const {

    return std::string(m_node.name());

}
/*
============================================================
================== test for empty node =====================
============================================================
 *
 * @return
 */
bool PugiXMLNode::empty() const {

    return m_node.empty();

}
/*
============================================================
======================= text data ==========================
============================================================
 *
 * @return
 */
Text PugiXMLNode::text() const {

    return Text( std::string(m_node.text().get()) );

}
/*
============================================================
===================== numeric data =========================
============================================================
 *
 * @return
 */
Data_internal *PugiXMLNode::data() const {

    return new PugiXMLData( m_node );

}

}
#endif
