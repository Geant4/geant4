/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef GUPI_hpp_included
#define GUPI_hpp_included 1

#include <list>
#include <map>

#include <LUPI_dataBuffer.hpp>
#include <LUPI.hpp>
#include <HAPI.hpp>

namespace GUPI {

class Entry;
class Suite;

typedef Entry *(*GUPI_parseSuite)( Suite *a_parent, HAPI::Node const &a_node );

#define GUPI_documentationChars "documentation"
#define GUPI_titleChars "title"
#define GUPI_abstractChars "abstract"
#define GUPI_bodyChars "body"
#define GUPI_endfCompatibleChars "endfCompatible"
#define GUPI_doiChars "doi"
#define GUPI_publicationDateChars "publicationDate"
#define GUPI_versionChars "version"

/*
============================================================
======================== WriteInfo =========================
============================================================
*/

class WriteInfo {

    public:
        std::list<std::string> m_lines;
        std::string m_incrementalIndent;
        int m_valuesPerLine;
        std::string m_sep;

        WriteInfo( std::string const &a_incrementalIndent = "  ", int a_valuesPerLine = 100, std::string const &a_sep = " " );

        std::string incrementalIndent( std::string const &indent ) { return( indent + m_incrementalIndent ); }
        void push_back( std::string const &a_line ) { m_lines.push_back( a_line ); }

        void addNodeStarter( std::string const &indent, std::string const &a_moniker, std::string const &a_attributes = "" ) {
                m_lines.push_back( indent + "<" + a_moniker + a_attributes + ">" ); }
        void addNodeStarterEnder( std::string const &indent, std::string const &a_moniker, std::string const &a_attributes = "" ) {
                m_lines.push_back( indent + "<" + a_moniker + a_attributes + "/>" ); }
        void addNodeEnder( std::string const &a_moniker ) { m_lines.back( ) += "</" + a_moniker + ">"; }
        std::string addAttribute( std::string const &a_name, std::string const &a_value ) const { return( " " + a_name + "=\"" + a_value + "\"" ); }

        std::string nodeStarter( std::string const &indent, std::string const &a_moniker, std::string const &a_attributes = "" ) 
                { return( indent + "<" + a_moniker + a_attributes + ">" ); }
        std::string nodeEnder( std::string const &a_moniker ) { return( "</" + a_moniker + ">" ); }

        void print( );
        void clear( ) { m_lines.clear( ); }      /**< Clears the contents of *m_lines*. */
};

/*
============================================================
========================= Ancestry =========================
============================================================
*/
class Ancestry {

    public:
            /* *********************************************************************************************************//**
             * Constructs and returns the key name/value for the *this* node.
             *
             * @return          The constructed key name/value.
             ***********************************************************************************************************/
        static std::string buildXLinkItemKey( std::string const &a_name, std::string const &a_key ) {

            if( a_key.size( ) == 0 ) return( "" );
            return( "[@" + a_name + "='" + a_key + "']" );
        }

    private:
        std::string m_moniker;                                  /**< The node's name (i.e., moniker). */
        Ancestry *m_ancestor;                                   /**< The parent node of *this*. */
        std::string m_attribute;                                /**< The name of the attribute in the node that uniquely identifies the node when the parent node containing other child nodes with the same moniker. */

        Ancestry *findInAncestry2( std::size_t a_index, std::vector<std::string> const &a_segments );
        Ancestry const *findInAncestry2( std::size_t a_index, std::vector<std::string> const &a_segments ) const ;

    public:
        Ancestry( std::string const &a_moniker, std::string const &a_attribute = "" );
        virtual ~Ancestry( );
        Ancestry &operator=( Ancestry const &a_ancestry );

        std::string const &moniker( ) const { return( m_moniker ); }                               /**< Returns the value of the *m_moniker* member. */
        void setMoniker( std::string const &a_moniker ) { m_moniker = a_moniker; }          /**< Set the value of the *m_moniker* member to *a_moniker*. */
        Ancestry *ancestor( ) { return( m_ancestor ); }                                     /**< Returns the value of the *m_ancestor* member. */
        Ancestry const *ancestor( ) const { return( m_ancestor ); }                                     /**< Returns the value of the *m_ancestor* member. */
        void setAncestor( Ancestry *a_ancestor ) { m_ancestor = a_ancestor; }               /**< Sets the *m_ancestor* member to *a_ancestor*. */
        std::string attribute( ) const { return( m_attribute ); }                           /**< Returns the value of the *m_attribute* member. */

        Ancestry *root( );
        Ancestry const *root( ) const ;
        bool isChild( Ancestry *a_instance ) { return( this == a_instance->m_ancestor ); }  /**< Returns true if *a_instance* is a child of *this*. */
        bool isParent( Ancestry *a_parent ) { return( this->m_ancestor == a_parent ); }     /**< Returns true if *a_instance* is the parent of *this*. */
        bool isRoot( ) const { return( this->m_ancestor == nullptr ); }                     /**< Returns true if *this* is the root ancestor. */

        Ancestry *findInAncestry( std::string const &a_href );
        Ancestry const *findInAncestry( std::string const &a_href ) const ;

            /* *********************************************************************************************************//**
             * Used to tranverse **GNDS** nodes. This method returns a pointer to a derived class' *a_item* member or nullptr if none exists.
             *
             * @param a_item    [in]    The name of the class member whose pointer is to be return.
             * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
             ***********************************************************************************************************/
        virtual Ancestry *findInAncestry3( std::string const &a_item ) = 0;
        virtual Ancestry const *findInAncestry3( std::string const &a_item ) const = 0;

        virtual LUPI_HOST void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        virtual std::string xlinkItemKey( ) const { return( "" ); }                         /**< Returns the value of *this*'s key. */
        std::string toXLink( ) const ;

        virtual void toXMLList( WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
        void printXML( ) const ;
};

/*
============================================================
========================== Entry ===========================
============================================================
*/
class Entry : public Ancestry {

    private:
        std::string m_keyName;                                              /**< The name of the key used by the parent suite to reference *this* entry. */
        std::string m_keyValue;                                             /**< The key used by the parent suite to reference *this* entry. */

    public:
        Entry( std::string const &a_moniker, std::string const &a_keyName, std::string const &a_keyValue );
        Entry( HAPI::Node const &a_node, std::string const &a_keyName );
        ~Entry( );

        std::string const &keyName( ) const { return( m_keyName ); }        /**< Returns a const reference to the *m_keyName* member. */
        std::string const &keyValue( ) const { return( m_keyValue ); }      /**< Returns a const reference to the *m_keyValue* member. */

        Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }
        LUPI_HOST void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        std::string xlinkItemKey( ) const {

            if( m_keyValue == "" ) return( "" );
            return( buildXLinkItemKey( m_keyName, m_keyValue ) );
        } 
};

/*
============================================================
========================== Text ============================
============================================================
*/

class Text : public Ancestry {

    public:
        enum class Encoding {
            utf8,
            ascii
        };

        enum class Markup {
            none,
            xml,
            html,
            latex
        };
    
    private:
        std::string m_body;
        Encoding m_encoding;
        Markup m_markup;
        std::string m_label;

    public:
        Text( HAPI::Node const &a_node );
        ~Text( );

        std::string const &body( ) const { return m_body; }
        Encoding encoding( ) const { return m_encoding; }
        Markup markup( ) const { return m_markup; }
        std::string const &label( ) const { return m_label; }

        Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }

};

/*
============================================================
===================== Documentation ========================
============================================================
*/
class Documentation : public Ancestry {

    private:
        std::string m_doi;                                              /**< The name of the key used by the parent suite to reference *this* entry. */
        std::string m_publicationDate;                                             /**< The key used by the parent suite to reference *this* entry. */
        std::string m_version;

        Text m_title;
        Text m_abstract;
        Text m_body;     

    public:
        // Documentation(std::string const &a_moniker, Text const &a_doi, std::string const &a_publicationDate, Text const &a_version);
        Documentation(HAPI::Node const &a_node);
        ~Documentation( );

        std::string const &doi( ) const { return m_doi; }
        std::string const &publicationDate( ) const { return m_publicationDate; }
        std::string const &version( ) const { return m_version; }

        Text const &title( ) const { return m_title; }
        Text const &abstract( ) const { return m_abstract; }
        Text const &body( ) const { return m_body; }

        Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }

};

/*
============================================================
=========================== Suite ==========================
============================================================
*/
class Suite : public Ancestry {

    public:
        typedef std::vector<Entry *> Entries;                           /**< The typedef the the *m_entries* member. */

    private:
        std::string m_keyName;                                          /**< The name of the key used to look up items in the suite. */
        mutable Entries m_entries;                                      /**< The list of nodes stored within *this*. */
        std::map<std::string,std::size_t> m_map;                        /**< A map of *this* node labels to their index in *m_entries*. */

        Suite( Suite const *a_suite );                  // FIXME, should we make public or private copy constructor? Making private for now.

    public:
        Suite( std::string const &a_keyName );
        Suite( std::string const &a_moniker, std::string const &a_keyName );
        Suite( HAPI::Node const &a_node, std::string const &a_keyName, GUPI_parseSuite a_parseSuite );
        ~Suite( );

        std::string const &keyName( ) const { return( m_keyName ); }                        /**< Returns a const reference to the *m_keyName* member. */
        std::size_t size( ) const { return( m_entries.size( ) ); }                            /**< Returns the number of node contained by *this*. */

        std::size_t operator[]( std::string const &a_label ) const ;
        typedef Entries::iterator iterator;
        typedef Entries::const_iterator const_iterator;
        iterator begin( ) { return m_entries.begin( ); }                                      /**< The C++ begin iterator for *this*. */
        const_iterator begin( ) const { return m_entries.begin( ); }                          /**< The C++ const begin iterator for *this*. */
        iterator end( ) { return m_entries.end( ); }                                          /**< The C++ end iterator for *this*. */
        const_iterator end( ) const { return m_entries.end( ); }                              /**< The C++ const end iterator for *this*. */

        template<typename T> T       *get( std::size_t a_Index );
        template<typename T> T const *get( std::size_t a_Index ) const ;
        template<typename T> T       *get( std::string const &a_label );
        template<typename T> T const *get( std::string const &a_label ) const ;

        void parse( HAPI::Node const &a_node, GUPI_parseSuite a_parseSuite );
        void add( Entry *a_entry );
        iterator find( std::string const &a_label );
        const_iterator find( std::string const &a_label ) const ;
        bool has( std::string const &a_label ) const { return( find( a_label ) != m_entries.end( ) ); }

        Ancestry *findInAncestry3( std::string const &a_item );
        Ancestry const *findInAncestry3( std::string const &a_item ) const ;
        std::vector<iterator> findAllOfMoniker( std::string const &a_moniker ) ;
        std::vector<const_iterator> findAllOfMoniker( std::string const &a_moniker ) const ;

        void toXMLList( WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
        void printEntryLabels( std::string const &a_header ) const ;
};

/* *********************************************************************************************************//**
 * Returns the node at index *a_index*.
 *
 * @param a_index               [in]    The index of the node to return.
 *
 * @return                              The node at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T *Suite::get( std::size_t a_index ) {

    Entry *entry = m_entries[a_index];
    T *object = dynamic_cast<T *>( entry );

    if( object == nullptr ) throw LUPI::Exception( "GIDI::Suite::get( std::size_t ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node at index *a_index*.
 *
 * @param a_index               [in]    The index of the node to return.
 *
 * @return                              The node at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T const *Suite::get( std::size_t a_index ) const {

    Entry *entry = m_entries[a_index];
    T *object = dynamic_cast<T *>( entry );

    if( object == nullptr ) throw LUPI::Exception( "GIDI::Suite::get( std::size_t ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node with label *a_label*.
 *
 * @param a_label               [in]    The label of the node to return.
 *
 * @return                              The node with label *a_label*.
 ***********************************************************************************************************/

template<typename T> T *Suite::get( std::string const &a_label ) {

    auto index = (*this)[a_label];
    Entry *entry = m_entries[index];
    T *object = dynamic_cast<T *>( entry );

    if( object == nullptr ) throw LUPI::Exception( "GIDI::Suite::get( std::string const & ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node with label *a_label*.
 *
 * @param a_label               [in]    The label of the node to return.
 *
 * @return                              The node with label *a_label*.
 ***********************************************************************************************************/

template<typename T> T const *Suite::get( std::string const &a_label ) const {

    auto index = (*this)[a_label];
    Entry *entry = m_entries[index];
    T *object = dynamic_cast<T *>( entry );

    if( object == nullptr ) throw LUPI::Exception( "GIDI::Suite::get( std::string const & ): invalid cast" );

    return( object );
}

}               // End of namespace GUPI.

#endif          // GUPI_hpp_included
