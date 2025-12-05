/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include "GIDI.hpp"
#include <HAPI.hpp>

namespace GIDI {

/*! \class Form 
 * Base class inherited by most other GIDI classes. Mainly contains **label** and **type** members.
 */

/* *********************************************************************************************************//**
 *
 * @param a_type            [in]    The *FormType* the class represents.
 ***********************************************************************************************************/

Form::Form( FormType a_type ) :
        GUPI::Ancestry( "" ),
        m_parent( nullptr ),
        m_type( a_type ),
        m_keyName( GIDI_labelChars ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker         [in]    The moniker for *this*.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_label           [in]    The label for *this*.
 ***********************************************************************************************************/

Form::Form( std::string const &a_moniker, FormType a_type, std::string const &a_label ) :
        GUPI::Ancestry( a_moniker ),
        m_parent( nullptr ),
        m_type( a_type ),
        m_keyName( GIDI_labelChars ),
        m_label( a_label ) {

}

/* *********************************************************************************************************//**
 *
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Form::Form( HAPI::Node const &a_node, LUPI_maybeUnused SetupInfo &a_setupInfo, FormType a_type, Suite *a_parent ) :
        GUPI::Ancestry( a_node.name( ) ),
        m_parent( a_parent ),
        m_type( a_type ),
        m_keyName( GIDI_labelChars ),
        m_label( a_node.attribute_as_string( GIDI_labelChars ) ) {

}

/* *********************************************************************************************************//**
 * @param a_form            [in]    The Form to copy.
 ***********************************************************************************************************/

Form::Form( Form const &a_form ) :
        GUPI::Ancestry( a_form.moniker( ), a_form.attribute( ) ),
        m_parent( nullptr ),
        m_type( a_form.type( ) ),
        m_keyName( a_form.keyName( ) ),
        m_keyValue( a_form.keyValue( ) ),
        m_label( a_form.label( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Form::~Form( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for
 * the member *m_parent* which is set to **nullptr** and those not set by base classes.
 *
 * @param a_rhs                [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

Form &Form::operator=( Form const &a_rhs ) {

    if( this != &a_rhs ) {
        GUPI::Ancestry::operator=( a_rhs );

        m_parent = nullptr;
        m_type = a_rhs.type( );
        m_keyName = a_rhs.keyName( );
        m_keyValue = a_rhs.keyValue( );
        m_label = a_rhs.label( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * Set the *m_label* member per *a_label*. Also, if *m_keyName* is "label", calls **setKeyValue**.
 *
 * @param a_label               [in]    The value of the label.
 ***********************************************************************************************************/

void Form::setLabel( std::string const &a_label ) {

    m_label = a_label;
    if( m_keyName == GIDI_labelChars ) setKeyValue( GIDI_labelChars );
}

/* *********************************************************************************************************//**
 * Returns a *const* reference to the *m_keyName* member.
 *
 * @return                              The name of the key for *this*.
 ***********************************************************************************************************/

std::string const &Form::keyName( ) const {

    return( m_keyName );
}

/* *********************************************************************************************************//**
 * Set the *m_keyName* member to *a_keyName* and calls setKeyValue.
 *
 * @param a_keyName             [in]    The value of the key name.
 ***********************************************************************************************************/

void Form::setKeyName( std::string const &a_keyName ) {
    
    m_keyName = a_keyName;
    setKeyValue( m_keyName );
}

/* *********************************************************************************************************//**
 * Returns a *const* reference to the *m_keyValue* member.
 *
 * @return                              The value of the key for *this*.
 ***********************************************************************************************************/

std::string const &Form::keyValue( ) const {

    if( m_keyValue == "" ) setKeyValue( m_keyName );
    return( m_keyValue );
}

/* *********************************************************************************************************//**
 * Set the *m_keyValue* per the *a_keyName* name. This method assumes that *a_keyName* is "label". Otherwise, it executes a throw.
 *
 * @param a_keyName             [in]    The name of the key whose value is set.
 ***********************************************************************************************************/

void Form::setKeyValue( std::string const &a_keyName ) const {

    if( a_keyName != GIDI_labelChars ) throw Exception( "Form::setKeyValue: unsupported keyname \""  + a_keyName + "\"." );

    m_keyValue = m_label;
}

/* *********************************************************************************************************//**
 * Returns the sibling of *this* with label *a_label*.
 *
 * @param a_label           [in]    The label of the sibling to find.
 * @return                          The sibling with label *a_label*.
 ***********************************************************************************************************/

Form const *Form::sibling( std::string a_label ) const {

    Form *_form;

    try {
        _form = ((*parent( )).get<Form>( a_label ) ); }
    catch (...) {
        return( nullptr );
    }
    return( _form );
}

namespace Functions {

/*! \class FunctionForm
 * Base class inherited by other GIDI function classes.
 */

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_dimension           [in]    The dimension of the function.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

FunctionForm::FunctionForm( std::string const &a_moniker, FormType a_type, int a_dimension, ptwXY_interpolation a_interpolation, 
                int a_index, double a_outerDomainValue ) :
        Form( a_moniker, a_type, "" ),
        m_dimension( a_dimension ),
        m_interpolation( a_interpolation ),
        m_index( a_index ),
        m_outerDomainValue( a_outerDomainValue ) {

    m_interpolationString = ptwXY_interpolationToString( m_interpolation );
}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_dimension           [in]    The dimension of the function.
 * @param a_axes                [in]    The axes for the function.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

FunctionForm::FunctionForm( std::string const &a_moniker, FormType a_type, int a_dimension, Axes const &a_axes, ptwXY_interpolation a_interpolation, 
                int a_index, double a_outerDomainValue ) :
        Form( a_type ),
        m_dimension( a_dimension ), 
        m_axes( a_axes ),
        m_interpolation( a_interpolation ),
        m_index( a_index ),
        m_outerDomainValue( a_outerDomainValue ) {

    setMoniker( a_moniker );
    m_interpolationString = ptwXY_interpolationToString( m_interpolation );

    m_axes.setAncestor( this );
}

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the FunctionForm.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_dimension       [in]    The dimension of the function.
 * @param a_suite           [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

FunctionForm::FunctionForm( LUPI_maybeUnused Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		        FormType a_type, int a_dimension, Suite *a_suite ) :
        Form( a_node, a_setupInfo, a_type, a_suite ),
        m_dimension( a_dimension ),
        m_axes( a_node.child( GIDI_axesChars ), a_setupInfo, 0 ),
        m_interpolation( ptwXY_interpolationLinLin ),
        m_index( 0 ), 
        m_outerDomainValue( 0.0 ) {

    m_interpolationString = a_node.attribute_as_string( GIDI_interpolationChars );
    m_interpolation = ptwXY_stringToInterpolation( m_interpolationString.c_str( ) );
    if( m_interpolation != ptwXY_interpolationOther ) m_interpolationString = ptwXY_interpolationToString( m_interpolation );

    if( strcmp( a_node.attribute_as_string( GIDI_indexChars ).c_str( ), "" ) != 0 ) m_index = a_node.attribute_as_int( GIDI_indexChars );
    if( strcmp( a_node.attribute_as_string( GIDI_outerDomainValueChars ).c_str( ), "" ) != 0 ) m_outerDomainValue = a_node.attribute_as_double( GIDI_outerDomainValueChars );
}

/* *********************************************************************************************************//**
 * @param a_form            [in]    The FunctionForm to copy.
 ***********************************************************************************************************/

FunctionForm::FunctionForm( FunctionForm const &a_form ) :
        Form( a_form ),
        m_dimension( a_form.dimension( ) ),
        m_axes( a_form.axes( ) ),
        m_interpolation( a_form.interpolation( ) ),
        m_interpolationString( a_form.interpolationString( ) ),
        m_index( a_form.index( ) ),
        m_outerDomainValue( a_form.outerDomainValue( ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

FunctionForm::~FunctionForm( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

FunctionForm &FunctionForm::operator=( FunctionForm const &a_rhs ) {

    if( this != &a_rhs ) {
        Form::operator=( a_rhs );

        m_dimension = a_rhs.dimension( );
        m_axes = a_rhs.axes( );
        m_interpolation = a_rhs.interpolation( );
        m_interpolationString = a_rhs.interpolationString( );
        m_index = a_rhs.index( );
        m_outerDomainValue = a_rhs.outerDomainValue( );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * This method sets the integer (i.e., *ptwXY_interpolation*) interpolation value of **m_interpolation** to
 * **a_interpolation**. This method also sets the **m_interpolationString** member per **a_interpolation**.
 *
 * @param       a_interpolation     [in]        The *ptwXY_interpolation* integer value of the interpolaction.
 ***********************************************************************************************************/

void FunctionForm::setInterpolation( ptwXY_interpolation a_interpolation ) {

    m_interpolation = a_interpolation;
    m_interpolationString = ptwXY_interpolationToString( m_interpolation );
} 

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions container.
 ***********************************************************************************************************/

void FunctionForm::toXMLList_func( LUPI_maybeUnused GUPI::WriteInfo &a_writeInfo, LUPI_maybeUnused std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const {

    std::cout << "Node '" << moniker( ) << "' needs toXMLList methods." << std::endl;
}

/*! \class Function1dForm
 * Base class inherited by other GIDI 1d function classes.
 */

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Function1dForm::Function1dForm( std::string const &a_moniker, FormType a_type, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue ) :
        FunctionForm( a_moniker, a_type, 1, a_interpolation, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Function1dForm::Function1dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, 
                int a_index, double a_outerDomainValue ) :
        FunctionForm( a_moniker, a_type, 1, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the FunctionForm.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_suite           [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Function1dForm::Function1dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type,
		Suite *a_suite ) :
        FunctionForm( a_construction, a_node, a_setupInfo, a_type, 1, a_suite ) {

}

/* *********************************************************************************************************//**
 * @param a_form            [in]    The Function1dForm to copy.
 ***********************************************************************************************************/

Function1dForm::Function1dForm( Function1dForm const &a_form ) :
        FunctionForm( a_form ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Function1dForm::~Function1dForm( ) {

}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                 [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

Function1dForm &Function1dForm::operator=( Function1dForm const &a_rhs ) {

    if( this != &a_rhs ) {
        FunctionForm::operator=( a_rhs );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * This method executes a throw as the sub-class did not define it. Evaluates *this* at the X-values in *a_Xs*[*a_offset*:]
 * and adds the results to *a_results*[*a_offset*:].
 *
 * @param a_offset          [in]    The offset in *a_Xs* to start.
 * @param a_Xs              [in]    The list of domain values to evaluate *this* at.
 * @param a_results         [in]    The list whose values are added to by the Y-values of *this*.
 * @param a_scaleFactor     [in]    A factor applied to each evaluation before it is added to *a_results*. 
 ***********************************************************************************************************/

void Function1dForm::mapToXsAndAdd( LUPI_maybeUnused std::size_t a_offset, LUPI_maybeUnused std::vector<double> const &a_Xs, LUPI_maybeUnused std::vector<double> &a_results, LUPI_maybeUnused double a_scaleFactor ) const {

    throw Exception( "Function1dForm::mapToXsAndAdd: function " + moniker( ) + " not implemented." );
}
/* *********************************************************************************************************//**
 * This method writes a warning message stating that the write method has not been implemented for *this*.
 * 
 * @param       a_file              [in]    The C FILE instance to write the data to.
 * @param       a_format            [in]    The format string passed to each region's write method.
 ***********************************************************************************************************/

/* *********************************************************************************************************//**
 * This method returns a **nullptr** since this methods was not implemented in the the derived class.
 *
 * @param   a_asLinlin          [in]    If **true**, the inpolatation of the returned XYs1d instance will always be lin-lin. Otherwise,
 *                                      the interpolation depends on the child 1d functions. This argument is not needed or used for this class.
 * @param   a_accuracy          [in]    The accuracy use to convert the data to lin=lin interpolation if needed. This argument is not needed or used for this cl
 * @param   a_lowerEps          [in]    The relative domain ammount to put a point below a boundary between two regions. This argument is not needed or used for
 * @param   a_upperEps          [in]    The relative domain ammount to put a point above a boundary between two regions. This argument is not needed or used for
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Function1dForm::asXYs1d( LUPI_maybeUnused bool a_asLinlin, LUPI_maybeUnused double a_accuracy, LUPI_maybeUnused double a_lowerEps, LUPI_maybeUnused double a_upperEps ) const {

    return( nullptr );
}
 
void Function1dForm::write( FILE *a_file, LUPI_maybeUnused std::string const &a_format ) const {

    fprintf( a_file, "# Write method not implemented for %s.\n", moniker( ).c_str( ) );
}

/* *********************************************************************************************************//**
 * Calls the write method with stdout as the file stream.
 *
 * @param       a_format            [in]    The format string passed to the C printf function.
 ***********************************************************************************************************/

void Function1dForm::print( std::string const &a_format ) const {

    write( stdout, a_format.c_str( ) );
}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The **GNDS** node name for the 2d function.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    The GNDS **index** value for *this*.
 * @param a_outerDomainValue    [in]    The GNDS **outerDomainValue** value for *this*.
 ***********************************************************************************************************/

Function2dForm::Function2dForm( std::string const &a_moniker, FormType a_type, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue ) :
        FunctionForm( a_moniker, a_type, 2, a_interpolation, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Function2dForm::Function2dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, 
                int a_index, double a_outerDomainValue ) :
        FunctionForm( a_moniker, a_type, 2, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the FunctionForm.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_suite           [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Function2dForm::Function2dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		FormType a_type, Suite *a_suite ) :
        FunctionForm( a_construction, a_node, a_setupInfo, a_type, 2, a_suite ) {

}

/* *********************************************************************************************************//**
 * @param a_form            [in]    Function2dForm to copy.
 ***********************************************************************************************************/
Function2dForm::Function2dForm( Function2dForm const &a_form ) :
        FunctionForm( a_form ) {
        
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Function2dForm::~Function2dForm( ) {

}

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the FunctionForm.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_type            [in]    The *FormType* the class represents.
 * @param a_suite           [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Function3dForm::Function3dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
		FormType a_type, Suite *a_suite ) :
        FunctionForm( a_construction, a_node, a_setupInfo, a_type, 3, a_suite ) {

}

/* *********************************************************************************************************//**
 * @param a_moniker             [in]    The moniker for *this*.
 * @param a_type                [in]    The *FormType* the class represents.
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 ***********************************************************************************************************/

Function3dForm::Function3dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, 
                int a_index, double a_outerDomainValue ) :
        FunctionForm( a_moniker, a_type, 3, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Function3dForm::~Function3dForm( ) {

}

}               // End namespace Functions.

}               // End namespace GIDI.
