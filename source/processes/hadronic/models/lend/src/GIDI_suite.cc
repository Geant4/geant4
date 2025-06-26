/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/*! \class Suite
 * This class is used to store a list (i.e., suite) of similar type **GNDS** nodes.
*/

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::Suite( std::string const &a_keyName ) :
        GUPI::Ancestry( "" ),
        m_keyName( a_keyName ),
        m_styles( nullptr ),
        m_allowsLazyParsing( false ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** moniker for the Suite instance.
 * @param a_keyName             [in]    The name of the key for elements of *this*.
 ***********************************************************************************************************/

Suite::Suite( std::string const &a_moniker, std::string const &a_keyName ) :
        GUPI::Ancestry( a_moniker ),
        m_keyName( a_keyName ),
        m_styles( nullptr ),
        m_allowsLazyParsing( false ) {

}

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_moniker             [in]    The **GNDS** moniker for the Suite instance.
 * @param a_node                [in]    The HAPI::Node to be parsed and used to construct the Suite.
 * @param a_keyName             [in]    The name of the key for referencing up child nodes.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs        [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                      This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parseSuite          [in]    This function to call to parse each sub-node.
 * @param a_styles              [in]    The <**styles**> node under the <**reactionSuite**> node.
 * @param a_allowsLazyParsing   [in]    Boolean stating if the suite allows lazy parsing.
 ***********************************************************************************************************/

Suite::Suite( Construction::Settings const &a_construction, std::string const &a_moniker, std::string const &a_keyName, HAPI::Node const &a_node, 
                SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, parseSuite a_parseSuite,
                Styles::Suite const *a_styles, bool a_allowsLazyParsing ) :
        GUPI::Ancestry( a_moniker ),
        m_keyName( a_keyName ),
        m_styles( a_styles ),
        m_allowsLazyParsing( a_allowsLazyParsing ),
        m_href( "" ) {

    HAPI::Node const node = a_node.child( a_moniker.c_str( ) );
    m_href = node.attribute_as_string( GIDI_hrefChars );

    if( !node.empty( ) ) parse( a_construction, node, a_setupInfo, a_pops, a_internalPoPs, a_parseSuite, a_styles );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Suite::~Suite( ) {

    for( std::vector<Form *>::const_iterator iter = m_forms.begin( ); iter < m_forms.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * This methods parses all the child nodes of *a_node*.
 *
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_node                [in]    The HAPI::Node to be parsed and used to construct the Product.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs        [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                      This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parseSuite          [in]    This function to call to parse each sub-node.
 * @param a_styles              [in]    The <**styles**> node under the <**reactionSuite**> node.
 ***********************************************************************************************************/

void Suite::parse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
		PoPI::Database const &a_internalPoPs, parseSuite a_parseSuite, GIDI::Styles::Suite const *a_styles ) {

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );

        Form *form = nullptr;

        if( m_allowsLazyParsing && a_construction.lazyParsing( ) ) {
            form = new LazyParsingHelperForm( a_construction, this, child, a_setupInfo, a_pops, a_internalPoPs, name, a_styles, a_parseSuite ); }
        else {
            form = a_parseSuite( a_construction, this, child, a_setupInfo, a_pops, a_internalPoPs, name, a_styles );
        }
        if( form != nullptr ) add( form );
    }
}

/* *********************************************************************************************************//**
 * Returns the index of the node in *this* that has keyValue *a_keyValue*.

 * @return                      [in]    The index of the node with keyValue *a_keyValue* in *this*.
 ***********************************************************************************************************/

int Suite::operator[]( std::string const &a_keyValue ) const {

    std::map<std::string, int>::const_iterator iter = m_map.find( a_keyValue );
    if( iter == m_map.end( ) ) {
        throw Exception( "form '" + a_keyValue + "' not in suite " + toXLink( ) + "." );
    }

    return( iter->second );
}

/* *********************************************************************************************************//**
 * Adds the node *a_form* to *this*.
 *
 * @param a_form                [in]    The form to add.
 ***********************************************************************************************************/

void Suite::add( Form *a_form ) {

    int i1 = 0;

    for( Suite::iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter, ++i1 ) {
        if( (*iter)->keyValue( ) == a_form->keyValue( ) ) {
            m_forms[i1] = a_form;
            a_form->setAncestor( this );
            return;
        }
    }
    m_map[a_form->keyValue( )] = (int) m_forms.size( );
    m_forms.push_back( a_form );
    a_form->setAncestor( this );
}

/* *********************************************************************************************************//**
 * Check to see if the form is a **LazyParsingHelperForm**, it so, parses the form, loads it before returning the requested form.
 *
 * @param a_index               [in]    The index of the child to return.
 *
 * @return                              The form at index *a_index*.
 ***********************************************************************************************************/

Form *Suite::checkLazyParsingHelperForm( std::size_t a_index ) {

    Form *form = m_forms[a_index];

    if( form->type( ) == FormType::lazyParsingHelperForm ) {
        LazyParsingHelperForm *lazyParsingHelperForm = static_cast<LazyParsingHelperForm *>( form );
        form = lazyParsingHelperForm->parse( );
        if( form == nullptr ) {             // Happens because several forms (e.g., CoulombPlusNuclearElastic) are not needed by transport codes and are not parsed.
            form = lazyParsingHelperForm; }
        else {
            form->setAncestor( this );
            m_forms[a_index] = form;
            delete lazyParsingHelperForm;
        }
    }

    return( form );
}

/* *********************************************************************************************************//**
 * Check to see if the form is a **LazyParsingHelperForm**, it so, parses the form, loads it before returning the requested form.
 *
 * @param a_index               [in]    The index of the child to return.
 *
 * @return                              The form at index *a_index*.
 ***********************************************************************************************************/

Form *Suite::checkLazyParsingHelperForm( std::size_t a_index ) const {

    Form *form = m_forms[a_index];

    if( form->type( ) == FormType::lazyParsingHelperForm ) {
        LazyParsingHelperForm *lazyParsingHelperForm = static_cast<LazyParsingHelperForm *>( form );
        form = lazyParsingHelperForm->parse( );
        if( form != nullptr ) {
            form->setAncestor( const_cast<Suite *>( this ) );
            m_forms[a_index] = form;
            delete lazyParsingHelperForm;
        }
    }

    return( form );
}
/* *********************************************************************************************************//**
 * Check to see if the form is a **LazyParsingHelperForm**, it so, parses the form, loads it before returning the requested form.
 *
 * @param a_iter                [in]    Iterator to the **Form** to check.
 *
 * @return                              The iterator to the old or converted form.
 ***********************************************************************************************************/

Suite::iterator Suite::checkLazyParsingHelperFormIterator( Suite::iterator a_iter ) {

    if( a_iter == end( ) ) return( a_iter );

    std::size_t index = (*this)[(*a_iter)->keyValue()];
    ++a_iter;
    checkLazyParsingHelperForm( index );

    return( --a_iter );
}

/* *********************************************************************************************************//**
 * Check to see if the form is a **LazyParsingHelperForm**, it so, parses the form, loads it before returning the requested form.
 *
 * @param a_iter                [in]    Iterator to the **Form** to check.
 *
 * @return                              The form at index *a_index*.
 ***********************************************************************************************************/

Suite::const_iterator Suite::checkLazyParsingHelperFormIterator( Suite::const_iterator a_iter ) const {
    
    if( a_iter == end( ) ) return( a_iter );
    
    std::size_t index = (*this)[(*a_iter)->keyValue()];
    a_iter++;
    checkLazyParsingHelperForm( index );

    return( --a_iter );
}

/* *********************************************************************************************************//**
 * Returns the iterator to the node with keyValue *a_keyValue*.
 *
 * @param a_keyValue                        [in]    The keyValue of the node to find.
 * @param a_convertLazyParsingHelperForm    [in]    If true and requested form is a LazyParsingHelperForm instance, that instance is replaced with the parsed form.
 *
 * @return                              The iterator to the node with keyValue *a_keyValue*.
 ***********************************************************************************************************/

Suite::iterator Suite::find( std::string const &a_keyValue, bool a_convertLazyParsingHelperForm ) {

    for( Suite::iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) {
        if( (*iter)->keyValue( ) == a_keyValue ) {
            if( a_convertLazyParsingHelperForm ) return( checkLazyParsingHelperFormIterator( iter ) );
            return( iter );
        }
    }
    return( m_forms.end( ) );
}

/* *********************************************************************************************************//**
 * Returns the iterator to the node with keyValue *a_keyValue*.
 *
 * @param a_keyValue                        [in]    The keyValue of the node to find.
 * @param a_convertLazyParsingHelperForm    [in]    If true and requested form is a LazyParsingHelperForm instance, that instance is replaced with the parsed form.
 *
 * @return                              The iterator to the node with keyValue *a_keyValue*.
 ***********************************************************************************************************/

Suite::const_iterator Suite::find( std::string const &a_keyValue, bool a_convertLazyParsingHelperForm ) const {

    for( Suite::const_iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) {
        if( (*iter)->keyValue( ) == a_keyValue ) {
            if( a_convertLazyParsingHelperForm ) return( checkLazyParsingHelperFormIterator( iter ) );
            return( iter );
        }
    }
    return( m_forms.end( ) );
}

/* *********************************************************************************************************//**
 * Returns a list of iterators to the nodes in *this* that have **GNDS** moniker *a_moniker*.
 *
 * @param a_moniker             [in]    The moniker to search for.
 *
 * @return                              List of iterators to the nodes in *this* that have moniker *a_moniker*.
 ***********************************************************************************************************/

std::vector<Suite::iterator> Suite::findAllOfMoniker( std::string const &a_moniker ) {

    std::vector<Suite::iterator> iters;

    for( Suite::iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) {
        if( (*iter)->moniker( ) == a_moniker ) iters.push_back( iter );
    }

    return( iters );
}

/* *********************************************************************************************************//**
 * Returns a list of iterators to the nodes in *this* that have **GNDS** moniker *a_moniker*.
 *
 * @param a_moniker             [in]    The moniker to search for.
 *
 * @return                              List of iterators to the nodes in *this* that have moniker *a_moniker*.
 ***********************************************************************************************************/

std::vector<Suite::const_iterator> Suite::findAllOfMoniker( std::string const &a_moniker ) const {

    std::vector<Suite::const_iterator> iters;

    for( Suite::const_iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) {
        if( (*iter)->moniker( ) == a_moniker ) iters.push_back( iter );
    }

    return( iters );
}

/* *********************************************************************************************************//**
 * This method finds the nearest form of instance Functions::XYs1d in *this* that is prior to the form with label *a_label*.
 *
 * @param a_label               [in]    The label of the form to start from when looking backwards.
 * @param a_formType            [in]    The type of form to return.
 *
 * @return                              Pointer to an Functions::XYs1d instance of nullptr if one not found.
 ***********************************************************************************************************/

Form const *Suite::findInstanceOfTypeInLineage( std::string const &a_label, std::string const &a_moniker ) const {

    Form const *form1 = nullptr;
    auto formIter = m_forms.end( );

    for( auto iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) {
        if( (*iter)->label( ) == a_label ) break;
        if( (*iter)->actualMoniker( ) == a_moniker ) formIter  = iter;
    }

    if( formIter != m_forms.end( ) ) {
        form1 = *checkLazyParsingHelperFormIterator( formIter );
    }

    return( form1 );
}

/* *********************************************************************************************************//**
 * This method finds the nearest form of instance Functions::XYs1d in *this* that is prior to the form with label *a_label*.
 *
 * @param a_styles              [in]    The styles suite for the protare.
 * @param a_label               [in]    The label of the form to start from when looking backwards.
 * @param a_formType            [in]    The type of form to return.
 *
 * @return                              Pointer to an Functions::XYs1d instance of nullptr if one not found.
 ***********************************************************************************************************/

Form *Suite::findInstanceOfTypeInLineage( Styles::Suite const &a_styles, std::string const &a_label, std::string const &a_moniker ) {

    auto stylesIter = a_styles.find( a_label );
    if( stylesIter != a_styles.end( ) ) {
        auto suiteIter = find( a_label );
        if( suiteIter != end( ) ) {
            if( (*suiteIter)->actualMoniker( ) == a_moniker ) return( *checkLazyParsingHelperFormIterator( suiteIter ) );
        }
        Styles::Base const *style = static_cast<Styles::Base const *>( *stylesIter );
        return( findInstanceOfTypeInLineage( a_styles, style->getDerivedStyle( )->keyValue( ), a_moniker ) );
    }

    return( nullptr );
}

/* *********************************************************************************************************//**
 * Only for internal use. Called by ProtareTNSL instance to zero the lower energy multi-group data covered by the ProtareSingle that
 * contains the TNSL data covers the lower energy multi-group data.
 *
 * @param a_maximumTNSL_MultiGroupIndex     [in]    A map that contains labels for heated multi-group data and the last valid group boundary
 *                                                  for the TNSL data for that boundary.
 ***********************************************************************************************************/

void Suite::modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex ) {

    for( auto iter = a_maximumTNSL_MultiGroupIndex.begin( ); iter != a_maximumTNSL_MultiGroupIndex.end( ); ++iter ) {
        auto formIter = find( iter->first, true );

        if( formIter == m_forms.end( ) ) continue;

        if( (*formIter)->type( ) == FormType::gridded1d ) {
            reinterpret_cast<Functions::Gridded1d *>( (*formIter) )->modifiedMultiGroupElasticForTNSL( iter->second ); }
        else if( (*formIter)->type( ) == FormType::gridded3d ) {
            reinterpret_cast<Functions::Gridded3d *>( (*formIter) )->modifiedMultiGroupElasticForTNSL( iter->second );
        }
    }
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry *Suite::findInAncestry3( std::string const &a_item ) {

    std::size_t index( a_item.find( '=' ) ), lastQuote = a_item.size( ) - 2;

    if( index == std::string::npos ) return( nullptr );
    ++index;
    if( index > lastQuote ) throw Exception( "Suite::findInAncestry3: invalide xlink" );
    if( a_item[index] != '\'' ) throw Exception( "Suite::findInAncestry3: invalid xlink, missing '." );
    ++index;
    if( a_item[lastQuote]  != '\'' ) throw Exception( "Suite::findInAncestry3: invalid xlink, missing endl '." );

    std::string keyValue( a_item.substr( index, lastQuote - index ) );

    return( get<GUPI::Ancestry>( keyValue ) );
}

/* *********************************************************************************************************//**
 * Used by GUPI::Ancestry to tranverse GNDS nodes. This method returns a pointer to a derived class' a_item member or nullptr if none exists.
 *
 * @param a_item    [in]    The name of the class member whose pointer is to be return.
 * @return                  The pointer to the class member or nullptr if class does not have a member named a_item.
 ***********************************************************************************************************/

GUPI::Ancestry const *Suite::findInAncestry3( std::string const &a_item ) const {

    std::size_t index( a_item.find( '=' ) ), lastQuote = a_item.size( ) - 2;

    if( index == std::string::npos ) return( nullptr );
    ++index;
    if( index > lastQuote ) throw Exception( "Suite::findInAncestry3: invalide xlink" );
    if( a_item[index] != '\'' ) throw Exception( "Suite::findInAncestry3: invalid xlink, missing '." );
    ++index;
    if( a_item[lastQuote]  != '\'' ) throw Exception( "Suite::findInAncestry3: invalid xlink, missing endl '." );

    std::string keyValue( a_item.substr( index, lastQuote - index ) );

    return( get<GUPI::Ancestry>( keyValue ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 ***********************************************************************************************************/

void Suite::toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );

    if( size( ) == 0 ) return;

    std::string XMLLine( a_indent + "<" + moniker( ) + ">" );
    a_writeInfo.push_back( XMLLine );

    for( Suite::const_iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) (*iter)->toXMLList( a_writeInfo, indent2 );

    a_writeInfo.addNodeEnder( moniker( ) );
}
 
/* *********************************************************************************************************//**
 * Prints the list of node keyValues to std::cout.
 *
 * @param a_header              [in]    A string printed before the list of keyValues is printed.
 ***********************************************************************************************************/

void Suite::printFormLabels( std::string const &a_header ) const {

    std::cout << a_header << ": size = " << size( ) << std::endl;
    
    for( Suite::const_iterator iter = m_forms.begin( ); iter != m_forms.end( ); ++iter ) 
            std::cout << "    " << (*iter)->keyValue( ) << std::endl;
}

/*! \class Component
 * This class is used to store a list (i.e., suite) of similar type **GNDS** form nodes.
*/

/* *********************************************************************************************************//**
 * @param a_construction        [in]    Used to pass user options to the constructor.
 * @param a_moniker             [in]    The **GNDS** moniker for the Suite instance.
 * @param a_keyName             [in]    The key name for elements of *this*.
 * @param a_node                [in]    The HAPI::Node to be parsed and used to construct the Product.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                [in]    The *external* PoPI::Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs        [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                      This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_parseSuite          [in]    This function to call to parse each sub-node.
 * @param a_styles              [in]    The <**styles**> node under the <**reactionSuite**> node.
 ***********************************************************************************************************/

Component::Component( Construction::Settings const &a_construction, std::string const &a_moniker, std::string const &a_keyName, 
                HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, 
                parseSuite a_parseSuite, Styles::Suite const *a_styles ) :
        Suite( a_construction, a_moniker, a_keyName, a_node, a_setupInfo, a_pops, a_internalPoPs, a_parseSuite, a_styles, true ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** moniker for the Suite instance.
 * @param a_keyName             [in]    The key name for elements of *this*.
 ***********************************************************************************************************/

Component::Component( std::string const &a_moniker, std::string const &a_keyName ) :
        Suite( a_moniker, a_keyName ) {

}

/*! \class LazyParsingHelperForm
 * This class stores information about a GNDS node so that it can be parsed at a later time if needed.
*/

/* *********************************************************************************************************//**
 * Constructor that stores information so the *a_node* can be parsed at a later time.
 *
 * @param a_construction            [in]    Used to pass user options for parsing.
 * @param a_parent                  [in]    The parent GIDI::Suite that the returned Form will be added to.
 * @param a_node                    [in]    The **HAPI::Node** to be parsed.
 * @param a_setupInfo               [in]    Information create my the Protare constructor to help in parsing.
 * @param a_pops                    [in]    A PoPs Database instance used to get particle indices and possibly other particle information.
 * @param a_internalPoPs            [in]    The *internal* PoPI::Database instance used to get particle indices and possibly other particle information.
 *                                          This is the <**PoPs**> node under the <**reactionSuite**> node.
 * @param a_name                    [in]    The moniker for the node to be parsed.
 * @param a_styles                  [in]    A pointer to the <**styles**> node.
 * @param a_parser                  [in]    The parser function for the suite the actual form will be inserted into.
 ***********************************************************************************************************/

LazyParsingHelperForm::LazyParsingHelperForm( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node,
                SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name,
                Styles::Suite const *a_styles, parseSuite a_parser ) :
        Form( a_node, a_setupInfo, FormType::lazyParsingHelperForm, a_parent ),
        m_construction( a_construction ),
        m_node( a_node ),
        m_setupInfo( a_setupInfo ),
        m_pops( &a_pops ),
        m_internalPoPs( &a_internalPoPs ),
        m_name( a_name ),
        m_styles( a_styles ),
        m_parser( a_parser ) {

    m_construction.setLazyParsing( false );
    m_setupInfo.m_protare->incrementNumberOfLazyParsingHelperForms( );
}

/* *********************************************************************************************************//**
 * Constructor that stores information so the *a_node* can be parsed at a later time.
 ***********************************************************************************************************/

Form *LazyParsingHelperForm::parse( ) {

    Form *form = m_parser( m_construction, parent( ), m_node, m_setupInfo, *m_pops, *m_internalPoPs, m_name, m_styles );
    m_setupInfo.m_protare->incrementNumberOfLazyParsingHelperFormsReplaced( );

    return( form );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LazyParsingHelperForm::~LazyParsingHelperForm( ) {

}

}
