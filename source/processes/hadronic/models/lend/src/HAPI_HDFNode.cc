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

hid_t getNodeId(hid_t loc_id, std::string name);
herr_t map_children(hid_t loc_id, char const *name, const H5L_info_t *info, void *opdata);

namespace HAPI {

/*
=========================================================
 *
 * @return
 */
HDFNode::HDFNode() :
        Node_internal( NodeInteralType::HDF5 ),
        m_node_id( 0 ),
        m_parent_id( 0 ),
        m_index( 0 ) {

}
/*
=========================================================
 *
 * @param a_node
 * @return
 */
HDFNode::HDFNode( hid_t a_node_id, hid_t a_parent_id, size_t a_index, std::vector<childInfo> a_siblings ) :
        Node_internal( NodeInteralType::HDF5 ),
        m_node_id( a_node_id ),
        m_parent_id( a_parent_id ),
        m_index( a_index ),
        m_siblings( a_siblings ) {

    H5I_type_t id_type = H5Iget_type( m_node_id );
    if (H5I_GROUP == id_type) {
        hsize_t start_idx = 0;
        H5Literate(m_node_id, H5_INDEX_NAME, H5_ITER_NATIVE, &start_idx, map_children, &m_children);
    }
    else if (H5I_DATASET == id_type) {
        return; // has no child elements
    }
    else {
        throw "Unknown object type encountered in HDF file";
    }
}
/*
=========================================================
 *
 * Alternate constructor to allow treating H5File as a Node
 * 
 * @param a_file
 * @return
 */
HDFNode::HDFNode( hid_t a_file_id) :
        Node_internal( NodeInteralType::HDF5 ),
        m_node_id( a_file_id ),
        m_parent_id( 0 ),
        m_index( 0 ) {

  hsize_t start_idx = 0;
  H5Literate(m_node_id, H5_INDEX_NAME, H5_ITER_NATIVE, &start_idx, map_children, &m_children);

}
/*
============================================================
====================== copy constructor ====================
============================================================
 */
HDFNode::HDFNode(const HDFNode &other) :
        Node_internal( other ),
        m_node_id( other.m_node_id ),
        m_parent_id( other.m_parent_id ),
        m_index( other.m_index ),
        m_siblings( other.m_siblings ),
        m_children( other.m_children ) {

}
/*
=========================================================
*/
HDFNode::~HDFNode( ) {

}
/*
============================================================
=================== get attribute by name ==================
============================================================
 *
 * @param a_name
 * @return
 */
std::string HDFNode::attribute(char const *a_name) {

    if ((m_node_id == 0) || (!H5Aexists(m_node_id, a_name)))
        return "";

    hid_t attr_id = H5Aopen(m_node_id, a_name, H5P_DEFAULT);
    size_t attr_len = H5Aget_storage_size(attr_id);
    hid_t atype = H5Aget_type(attr_id);

    std::vector<char> buffer(attr_len+1, 0);
    H5Aread(attr_id, atype, buffer.data());

    std::string attrValue(buffer.data());

    return attrValue;
}


int HDFNode::attribute_as_int(const char* a_name){
  hid_t attr_id = H5Aopen(m_node_id, a_name, H5P_DEFAULT);
  size_t attr_len = H5Aget_storage_size(attr_id);
  hid_t atype = H5Aget_type(attr_id);

  std::vector<char> buffer(attr_len+1, 0);
  H5Aread(attr_id, atype, buffer.data());

  return atoi(buffer.data());

}
long HDFNode::attribute_as_long(const char* a_name){
  hid_t attr_id = H5Aopen(m_node_id, a_name, H5P_DEFAULT);
  size_t attr_len = H5Aget_storage_size(attr_id);

  std::vector<char> buffer(attr_len+1, 0);
  H5Aread(attr_id, H5T_NATIVE_CHAR, buffer.data());

  return atol(buffer.data());

}
double HDFNode::attribute_as_double(const char* a_name){
  hid_t attr_id = H5Aopen(m_node_id, a_name, H5P_DEFAULT);
  size_t attr_len = H5Aget_storage_size(attr_id);
  hid_t atype = H5Aget_type(attr_id);

  std::vector<char> buffer(attr_len+1, 0);
  H5Aread(attr_id, atype, buffer.data());

  return atof(buffer.data());
}

/*
============================================================
===================== get child element ====================
============================================================
 *
 * @return
 */
Node_internal *HDFNode::child(char const *a_name) {

    for (size_t idx=0; idx<m_children.size(); idx++)
    {
        if ( m_children[idx].xmlName == a_name ) {
            hid_t child_id = getNodeId(m_node_id, m_children[idx].name.c_str());
            HDFNode *child = new HDFNode(child_id, m_node_id, idx, m_children);
            return child;
        }
    }
    return new HDFNode();

}
/*
============================================================
===================== get first child node =================
============================================================
*/
Node_internal *HDFNode::first_child() {

    if (m_children.size() == 0)
        return new HDFNode();
    hid_t child_id = getNodeId(m_node_id, m_children[0].name);
    HDFNode *child = new HDFNode(child_id, m_node_id, 0, m_children);
    return child;
}
/*
============================================================
===================== get sibling element ==================
============================================================
 *
 * @return
 */
Node_internal *HDFNode::next_sibling() {

    size_t nextIndex = m_index + 1;
    if (nextIndex >= this->m_siblings.size())
        return new HDFNode();
    hid_t sibling_id = this->m_siblings[nextIndex].node_id;
    HDFNode *sibling = new HDFNode(sibling_id, m_parent_id, nextIndex, m_siblings);
    return sibling;

}
/*
============================================================
============ update self to point to next sibling ==========
============================================================
 *
 * @return
 */
void HDFNode::to_next_sibling() {

    size_t nextIndex = m_index + 1;
    size_t nsibs = this->m_siblings.size();
    if (nextIndex >= nsibs)
    {
        m_node_id = 0;
    }
    else
    {
        m_node_id = this->m_siblings[nextIndex].node_id;
        m_index = nextIndex;
    }
}
/*
============================================================
======================== make a copy =======================
============================================================
 *
 * @return
 */
Node_internal *HDFNode::copy() {

    if (m_node_id == 0)
      return new HDFNode();

    HDFNode *copy = new HDFNode( m_node_id, m_parent_id, m_index, m_siblings );
    copy->m_children = m_children;
    return copy;

}
/*
============================================================
===================== assignment operator ==================
============================================================
 */
Node_internal &HDFNode::operator=(const HDFNode &other) {

    this->m_node_id = other.m_node_id;
    this->m_parent_id = other.m_parent_id;
    this->m_index = other.m_index;
    this->m_siblings = other.m_siblings;
    this->m_children = other.m_children;
    return *this;

}
/*
============================================================
===================== get tag name =========================
============================================================
 *
 * @return
 */
std::string HDFNode::name() const {

    if (m_node_id == 0)
        return "";

    hid_t attr_id = H5Aopen(m_node_id, "_xmltag", H5P_DEFAULT);
    size_t attr_len = H5Aget_storage_size(attr_id);
    hid_t atype = H5Aget_type(attr_id);

    std::vector<char> buffer(attr_len+1, 0);
    H5Aread(attr_id, atype, buffer.data());

    return std::string(buffer.data());

}
/*
============================================================
================== test for empty node =====================
============================================================
 *
 * @return
 */
bool HDFNode::empty() const {

    return (m_node_id == 0);

}
/*
============================================================
======================= text data ==========================
============================================================
 *
 * @return
 */
Text HDFNode::text() const {

#if H5_VERSION_GE(1,12,0)
    H5O_info2_t infobuf;
    herr_t status = H5Oget_info3(m_node_id, &infobuf, H5O_INFO_NUM_ATTRS);
#else
    H5O_info_t infobuf;
    herr_t status = H5Oget_info(m_node_id, &infobuf);
#endif
    if (status != 0) throw "unable to extract text from HDF";

    switch (infobuf.type) {
    case H5O_TYPE_GROUP:
        return Text( );
        break;
    case H5O_TYPE_DATASET: {
        size_t len = H5Dget_storage_size(m_node_id);
        hid_t dtype = H5Dget_type(m_node_id);
        std::vector<char> buffer(len+1,0);
        H5Dread(m_node_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());

        return Text(std::string(buffer.data())); }
        break;
    default:
        throw "encountered unexpected type in HDF::text()!";
    }
}
/*
============================================================
===================== numeric data =========================
============================================================
 *
 * @return
 */
Data_internal *HDFNode::data() const {

  return new HDFData( m_node_id );

}

}


/*
============================================================
=============== get child HDF5 id by name ==================
============================================================
 */
hid_t getNodeId(hid_t loc_id, std::string name)
{
#if H5_VERSION_GE(1,12,0)
    H5O_info2_t infobuf;
    herr_t status = H5Oget_info_by_name3( loc_id, name.c_str( ), &infobuf, H5O_INFO_NUM_ATTRS, H5P_DEFAULT );
#else
    H5O_info_t infobuf;
    herr_t status = H5Oget_info_by_name(loc_id, name.c_str(), &infobuf, H5P_DEFAULT);
#endif
    if (status != 0) {
        throw "requested child node not found in getNodeId";
    }

    hid_t node_id;

    switch (infobuf.type) {
        case H5O_TYPE_GROUP:
            node_id = H5Gopen(loc_id, name.c_str(), H5P_DEFAULT);
            break;
        case H5O_TYPE_DATASET:
            node_id = H5Dopen(loc_id, name.c_str(), H5P_DEFAULT);
            break;
        default:
            throw "encountered unexpected type in getNodeId!";
    }

    return node_id;
}

/*
============================================================
========== called by iterateElems in constructor ===========
============================================================
 */
herr_t map_children(hid_t loc_id, char const *name, LUPI_maybeUnused const H5L_info_t *info, void *opdata)
{
  std::vector<HAPI::childInfo> *children =
          static_cast< std::vector<HAPI::childInfo>* >(opdata);

  hid_t node_id = H5Oopen(loc_id, name, H5P_DEFAULT);

  std::string xmlTag;
  {
    hid_t attr_id = H5Aopen(node_id, "_xmltag", H5P_DEFAULT);
    hid_t atype = H5Aget_type(attr_id);
    size_t attr_len = H5Aget_storage_size(attr_id);

    std::vector<char> buffer(attr_len+1, 0);
    H5Aread(attr_id, atype, buffer.data());

    xmlTag = std::string(buffer.data());
  }


  uint index;
  {
    hid_t attr_id = H5Aopen(node_id, "_xmlindex", H5P_DEFAULT);

    H5Aread(attr_id, H5T_NATIVE_UINT16_g, &index);
  }

  HAPI::childInfo infos = {
          std::string(name),
          xmlTag,
          index,
          node_id
  };

  if (index >= children->size())
  {
      children->resize(index+1);
  }
  (*children)[index] = infos;

  //printf("ITER: id=%d, index=%d, name=%s\n", (int)node_id, (int)index, name);

  //H5Oclose(node_id);

  return 0;
 }
#endif
