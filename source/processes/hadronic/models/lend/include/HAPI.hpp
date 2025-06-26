/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef HAPI_hpp_included
#define HAPI_hpp_included 1

#include <string>
#include <stdlib.h>
#include <vector>
#include <stdexcept>

#include <nf_buffer.h>
#include <nf_utilities.h>

#define HAPI_USE_PUGIXML 1

#ifdef HAPI_USE_PUGIXML
#include <pugixml.hpp>
#endif


#ifdef HAPI_USE_HDF5
#include <hdf5.h>
#endif

#include <LUPI.hpp>

namespace HAPI {

enum class NodeInteralType { pugiXML, HDF5 };

// container classes for reading in from various data sources:

/*
============================================================
========================= Attribute ========================
============================================================
 */

class Node_internal;
class Attribute {

    private:
        Node_internal *m_node;
        std::string m_name;
        //std::string m_value;

    public:
        inline Attribute() : m_node(nullptr), m_name() {}

        inline Attribute(Node_internal *a_node, std::string const a_name) :
          m_node(a_node),
          m_name(a_name)
        {
        }
        ~Attribute() = default;
        //std::string const &name() const { return( m_name ); }
        inline std::string const value() const;
        inline int as_int() const;
        inline long as_long() const;
        inline double as_double() const;
};

/*
============================================================
=========================== Text ===========================
============================================================
 */
class Text {

    private:
        std::string m_text;

    public:
        Text();
        Text(std::string const a_text);
        ~Text();
        std::string const &get() const { return( m_text ); }
};

/*
============================================================
================== Data_internal (base class) ==============
============================================================
 */
class Data_internal {

    public:
        Data_internal() { };
        virtual ~Data_internal() = 0;
        //std::string const getDataType();
        //virtual template <typename T> T read() = 0;
        virtual void getDoubles(nf_Buffer<double> &buffer) = 0;
        virtual void getInts(nf_Buffer<int> &buffer) = 0;
        virtual int length() const = 0;
};

/*
============================================================
=================== Node_internal (base class) =============
============================================================
 */
class Node_internal {

    private:
        NodeInteralType m_type;

    public:
        Node_internal( NodeInteralType a_type );
        Node_internal( Node_internal const &a_node );
        virtual ~Node_internal() = 0;

        NodeInteralType type( ) const { return( m_type ); }

        virtual std::string attribute(const char* name) = 0;
        virtual int attribute_as_int(const char* name) = 0;
        virtual long attribute_as_long(const char* name) = 0;
        virtual double attribute_as_double(const char* name) = 0;
        virtual Node_internal *child(const char* name) = 0;
        virtual Node_internal *first_child() = 0;
        virtual Node_internal *next_sibling() = 0;
        virtual void to_next_sibling() = 0;
        virtual Node_internal *copy() = 0;
        virtual std::string name() const = 0;
        virtual bool empty() const = 0;
        virtual Text text() const = 0;
        virtual Data_internal *data() const = 0;
};


inline std::string const Attribute::value() const { return( m_node->attribute( m_name.c_str()) ); }
inline int Attribute::as_int() const { return( m_node->attribute_as_int(m_name.c_str()) ); }
inline long Attribute::as_long() const { return( m_node->attribute_as_long(m_name.c_str()) ); }
inline double Attribute::as_double() const { return( m_node->attribute_as_double(m_name.c_str()) ); }

/*
============================================================
=========================== Data ===========================
============================================================
 */
class Data {

    private:
        Data_internal *m_data;

    public:
        Data();
        Data( Data_internal *a_data );
        ~Data();
        void getDoubles(nf_Buffer<double> &buffer);
        void getInts(nf_Buffer<int> &buffer);
        int length() const;
};

/*
============================================================
============================ Node ==========================
============================================================
 */
class Node {

    private:
        Node_internal *m_node;

    public:
        Node();
        Node( Node_internal *a_node );
        Node( Node const &a_node );
        ~Node();
        inline Attribute attribute(const char* a_name) const{
            return Attribute(m_node, a_name);
        }
        inline std::string attribute_as_string(const char* a_name) const{
            if(m_node == nullptr){
              return "";
            }
            return m_node->attribute(a_name);
        }
        inline int attribute_as_int(const char* a_name) const{
          if(m_node == nullptr){
            return 0;
          }
          return m_node->attribute_as_int(a_name);
        }
        inline long attribute_as_long(const char* a_name) const{
          if(m_node == nullptr){
            return 0;
          }
          return m_node->attribute_as_long(a_name);
        }
        inline double attribute_as_double(const char* a_name) const{
          if(m_node == nullptr){
            return 0.0;
          }
          return m_node->attribute_as_double(a_name);
        }
        Node child(const char* name) const;
        Node first_child() const;
        Node next_sibling() const;
        void to_next_sibling() const;
        Node &operator=(const Node &other);
        std::string name() const;
        bool empty() const;
        Text text() const;
        Data data() const;
};


/*
============================================================
======================= File (base class) ==================
============================================================
 */
class File {

    public:
        File() { };
        virtual ~File() = 0;
        virtual Node child(const char* name) = 0;
        virtual Node first_child() = 0;
        virtual std::string name() const = 0;
};


/*
============================================================
=============== Data Manager (for hybrid files) ============
============================================================
 */
class DataManager {

    public:
        DataManager() {};
        virtual ~DataManager() {};
        static DataManager* m_instance;
    
    public:
        virtual void getDoubles(nf_Buffer<double> &result, size_t startIndex, size_t endIndex) = 0;
        virtual void getInts(nf_Buffer<int> &result, size_t startIndex, size_t endIndex) = 0;
};


#ifdef HAPI_USE_PUGIXML
/*
============================================================
===================== XML using Pugi =======================
============================================================
 */
class PugiXMLNode : public Node_internal {

    private:
        pugi::xml_node m_node;

    public:
        PugiXMLNode();
        PugiXMLNode(pugi::xml_node a_node);
        PugiXMLNode(const PugiXMLNode &other);
        virtual ~PugiXMLNode();

        std::string attribute(const char* name);
        int attribute_as_int(const char* name);
        long attribute_as_long(const char* name);
        double attribute_as_double(const char* name);
        Node_internal *child(char const *name);
        Node_internal *first_child();
        Node_internal *next_sibling();
        void to_next_sibling();
        Node_internal *copy();
        Node_internal &operator=(const PugiXMLNode &other);
        std::string name() const;
        bool empty() const;
        Text text() const;
        Data_internal *data() const;
};

class PugiXMLData : public Data_internal {

    private:
        pugi::xml_node m_node;
        int m_length;

    public:
        PugiXMLData();
        PugiXMLData(pugi::xml_node a_node);
        virtual ~PugiXMLData();
        void getDoubles(nf_Buffer<double> &buffer);
        void getInts(nf_Buffer<int> &buffer);
        int length() const;
};

class PugiXMLFile : public File {

    private:
        std::string m_name;
        pugi::xml_document m_doc;

    public:
        PugiXMLFile();
        PugiXMLFile(char const *filename, std::string const &a_callingFunctionName);
        virtual ~PugiXMLFile();
        Node child(char const *name);
        Node first_child();
        std::string name() const;
};
#endif


#ifdef HAPI_USE_HDF5
/*
============================================================
=========================== HDF ============================
============================================================
 */
typedef struct {
    std::string name;
    std::string xmlName;  // HDF sometimes mangles names, need original name here
    size_t index;
    hid_t node_id;
} childInfo;

class HDFNode : public Node_internal {

    private:
        hid_t m_node_id;
        hid_t m_parent_id;
        size_t m_index;
        std::vector<childInfo> m_siblings;
        std::vector<childInfo> m_children;

    public:
        HDFNode();
        HDFNode(hid_t a_node_id, hid_t a_parent_id, size_t a_index, std::vector<childInfo> a_siblings);
        explicit HDFNode(hid_t a_file_id);
        HDFNode(const HDFNode &other);
        virtual ~HDFNode();

        //Attribute attribute(char const *name);
        std::string attribute(const char* name);
        int attribute_as_int(const char* name);
        long attribute_as_long(const char* name);
        double attribute_as_double(const char* name);
        Node_internal *child(char const *name);
        Node_internal *first_child();
        Node_internal *next_sibling();
        void to_next_sibling();
        Node_internal *copy();
        Node_internal &operator=(const HDFNode &other);
        std::string name() const;
        bool empty() const;
        Text text() const;
        Data_internal *data() const;

};

class HDFData : public Data_internal {

    private:
        hid_t m_node_id;
        hid_t m_dataspace_id;
        int m_length;

    public:
        HDFData();
        explicit HDFData(hid_t node_id);
        virtual ~HDFData();
        void getDoubles(nf_Buffer<double> &buffer);
        void getInts(nf_Buffer<int> &buffer);
        int length() const;
};

class HDFFile : public File {

    private:
        std::string m_name;
        hid_t m_doc;
        HDFNode *m_doc_as_node;

    public:
        HDFFile();
        explicit HDFFile(char const *filename);
        virtual ~HDFFile();
        Node child(char const *name);
        Node first_child();
        std::string name() const;
};

class HDFDataManager : public DataManager{

    private:
        std::string m_filename;
        bool m_iDataPresent;
        bool m_dDataPresent;
        hid_t m_file_id;
        hid_t m_dataset_ints, m_dataset_doubles;
        hid_t m_dataspace_ints, m_dataspace_doubles;

        hsize_t m_stride[1], m_block[1];

        size_t m_num_double_reads;
        size_t m_num_double_elem;
        size_t m_num_int_reads;
        size_t m_num_int_elem;

    public:
        HDFDataManager(std::string const &filename);
        virtual ~HDFDataManager();
        virtual void getDoubles(nf_Buffer<double> &result, size_t startIndex, size_t endIndex);
        virtual void getInts(nf_Buffer<int> &result, size_t startIndex, size_t endIndex);
};
#endif

}               // end of namespace 'HAPI'

#endif		// End of HAPI_hpp_included
