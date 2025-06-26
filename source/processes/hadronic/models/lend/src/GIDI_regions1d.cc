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

namespace Functions {

/* *********************************************************************************************************//**
 *
 *
 * @param   a_x         [in]    The 
 * @param   a_epsilon   [in]    The fractional amount to move the x-value.
 *
 * @return                      double.
 ***********************************************************************************************************/

static double getAdjustedX( double a_x, double a_epsilon ) {

    double xp = a_epsilon;

    if( a_x != 0.0 ) {
        if( a_x < 0 ) {
            xp = a_x * ( 1.0 - a_epsilon ); }
        else {
            xp = a_x * ( 1.0 + a_epsilon );
        }
    }

    return( xp );
}

/*! \class Regions1d
 * Class for the GNDS <**regions1d**> node.
 */

/* *********************************************************************************************************//**
 * @param a_construction    [in]    Used to pass user options to the constructor.
 * @param a_node            [in]    The **HAPI::Node** to be parsed and used to construct the XYs2d.
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent          [in]    The parent GIDI::Suite.
 ***********************************************************************************************************/

Regions1d::Regions1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::regions1d, a_parent ) {

    if( a_setupInfo.m_formatVersion.format( ) != GNDS_formatVersion_1_10Chars ) {
        data1dListParse( a_construction, a_node.child( GIDI_function1dsChars ), a_setupInfo, m_function1ds );
        checkSequentialDomainLimits1d( m_function1ds, m_Xs );
        return;                                     // Need to add uncertainty parsing.
    }

    for( HAPI::Node child = a_node.first_child( ); !child.empty( ); child.to_next_sibling( ) ) {
        std::string name( child.name( ) );

        if( name == GIDI_axesChars ) continue;
        if( name == GIDI_uncertaintyChars ) continue;

        Function1dForm *_form = data1dParse( a_construction, child, a_setupInfo, nullptr );
        if( _form == nullptr ) throw Exception( "Regions1d::Regions1d: data1dParse returned nullptr." );
        append( _form );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

Regions1d::~Regions1d( ) {

    for( std::vector<Function1dForm *>::iterator iter = m_function1ds.begin( ); iter < m_function1ds.end( ); ++iter ) delete *iter;
}

/* *********************************************************************************************************//**
 * Returns the domain minimum for the instance.
 *
 * @return          The domain minimum for the instance.
 ***********************************************************************************************************/

double Regions1d::domainMin( ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "Regions1d::domainMin: Regions1d has no regions" );
    return( m_Xs[0] );
}

/* *********************************************************************************************************//**
 * Returns the domain maximum for the instance.
 *
 * @return              The domain maximum for the instance.
 ***********************************************************************************************************/

double Regions1d::domainMax( ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "Regions1d::domainMax: Regions1d has not regions" );
    return( m_Xs[m_Xs.size( )-1] );
}

/* *********************************************************************************************************//**
 * Appends the 1d function *a_function* to the region.
 *
 * @param a_function            [in]    The 1d function (i.e., 1d region) to append to the Regions1d.
 ***********************************************************************************************************/
 
void Regions1d::append( Function1dForm *a_function ) {

    if( dimension( ) != a_function->dimension( ) ) throw Exception( "Regions1d::append: dimensions differ." );

    double _domainMin = a_function->domainMin( ), _domainMax = a_function->domainMax( );

    if( m_Xs.size( ) == 0 ) {
        m_Xs.push_back( _domainMin ); }
    else {
        if( m_Xs.back( ) != _domainMin ) throw Exception( "Regions1d::append: regions do not abut." );
    }

    m_Xs.push_back( _domainMax );
    m_function1ds.push_back( a_function );
}

/* *********************************************************************************************************//**
 * The value of *y(x1)* at the point *a_x1*.
 *
 * @param a_x1          [in]    Domain value to evaluate this at.
 * @return                      The value of this at the point *a_x1*.
 ***********************************************************************************************************/


double Regions1d::evaluate( double a_x1 ) const {

    if( m_Xs.size( ) == 0 ) throw Exception( "Regions1d::evaluate: Regions1d has not regions" );

    long iX1 = binarySearchVector( a_x1, m_Xs );

    if( iX1 < 0 ) {
        if( iX1 == -1 ) {       /* x1 > last value of Xs. */
            return( m_function1ds.back( )->evaluate( a_x1 ) );
        }
        iX1 = 0;                /* x1 < last value of Xs. */
    }
    return( m_function1ds[iX1]->evaluate( a_x1 ) );
}

/* *********************************************************************************************************//**
 * Evaluates *this* at the X-values in *a_Xs*[*a_offset*:] and adds the results to *a_results*[*a_offset*:].
 * *a_Xs* and *a_results* must be the same size otherwise a throw is executed.
 *
 * @param a_offset          [in]    The offset in *a_Xs* to start.
 * @param a_Xs              [in]    The list of domain values to evaluate *this* at.
 * @param a_results         [in]    The list whose values are added to by the Y-values of *this*.
 * @param a_scaleFactor     [in]    A factor applied to each evaluation before it is added to *a_results*.
 ***********************************************************************************************************/

void Regions1d::mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const {

    for( auto iter = m_function1ds.begin( ); iter < m_function1ds.end( ); ++iter ) {
        (*iter)->mapToXsAndAdd( a_offset, a_Xs, a_results, a_scaleFactor );
    }
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    If **true**, the inpolatation of the returned XYs1d instance will always be lin-lin. Otherwise, 
 *                                      the interpolation depends on the child 1d functions.
 * @param   a_accuracy          [in]    The accuracy use to convert the data to lin=lin interpolation if needed.
 * @param   a_lowerEps          [in]    The relative domain ammount to put a point below a boundary between two regions.
 * @param   a_upperEps          [in]    The relative domain ammount to put a point above a boundary between two regions.
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *Regions1d::asXYs1d( bool a_asLinlin, double a_accuracy, double a_lowerEps, double a_upperEps ) const {

    XYs1d *xys1d1 = nullptr, *xys1d2 = nullptr;

    if( a_lowerEps < 1e-14 ) a_lowerEps = 0.0;          // 1e-14 should be a global parameter.
    if( a_upperEps < 1e-14 ) a_upperEps = 0.0;

    if( !a_asLinlin ) {
        bool firstRegion = true;
        ptwXY_interpolation interpolation2 = ptwXY_interpolationLinLin;
        for( auto regionIter = m_function1ds.begin( ); regionIter != m_function1ds.end( ); ++regionIter ) {
            if( (*regionIter)->moniker( ) == GIDI_XYs1dChars ) {
                xys1d2 = static_cast<XYs1d *>( *regionIter );
                if( firstRegion ) interpolation2 = ptwXY_getInterpolation( xys1d2->ptwXY( ) );
                if( interpolation2 != ptwXY_getInterpolation( xys1d2->ptwXY( ) ) ) {
                    a_asLinlin = true;
                    break;
                } }
            else {
                a_asLinlin = true;
                break;
            }
        }
    }

    std::vector<double> xs;
    std::vector<double> ys;

    for( auto regionIter = m_function1ds.begin( ); regionIter != m_function1ds.end( ); ++regionIter ) {
        xys1d2 = (*regionIter)->asXYs1d( a_asLinlin, a_accuracy, a_lowerEps, a_upperEps );
        if( xys1d2 == nullptr ) {
            delete xys1d1;
            return( nullptr );
        }

        std::vector<double> xs2( xys1d2->xs( ) );
        if( xs2.size( ) < 2 ) continue;

        std::size_t xySize = xs.size( );
        std::vector<double> ys2( xys1d2->ys( ) );
        if( xySize > 0 ) {
            double y1u = ys.back( );
            double y2l = ys2[0];

            if( y1u == y2l ) {                              // Simple, just remove last point in xs and ys.
                xs.resize( xySize - 1 );
                ys.resize( xySize - 1 ); }
            else {
                bool addLower = false, addUpper = false;
                double x1l = xs[xySize-2];
                double x1u = xs.back( );
                double x2l = xs2[0];
                double x2u = xs2[1];
                double x1up = 0.0, x2lp= 0.0;

                if( a_lowerEps != 0.0 ) {
                    double x1lp = getAdjustedX( x1l, 0.8 *  a_lowerEps );
                    x1up        = getAdjustedX( x1u,       -a_lowerEps );       // Point to add in xs below xs.back() if room.

                    addLower = x1lp < x1up;                                     // Only add if relative spacing between existing point is greater than 1.8 * a_lowerEps.
                }

                if( a_upperEps > 0 ) {
                    x2lp        = getAdjustedX( x2l,        a_upperEps );       // Point to add in xs2 above xs2[0] if room.
                    double x2up = getAdjustedX( x2u, 0.8 * -a_upperEps );

                    addUpper = x2lp < x2up;                                     // Only add if relative spacing between existing point is greater than 1.8 * a_upperEps.
                }

                if( addLower ) {
                    if( addUpper ) {
                        xs.back( ) = x1up;
                        ys.back( ) = xys1d1->evaluate( x1up );
                        xs.push_back( x1u );
                        ys.push_back( 0.5 * ( y1u + y2l ) );
                        xs2[0] = x2lp;
                        ys2[0] = xys1d2->evaluate( x2lp ); }
                    else {
                        xs.back( ) = x1up;
                        ys.back( ) = xys1d1->evaluate( x1up );
                    } }
                else if( addUpper ) {
                    xs2[0] = x2lp;
                    ys2[0] = xys1d2->evaluate( x2lp ); }
                else {
                    ys2[0] = 0.5 * ( y1u + y2l );
                    xs.resize( xySize - 1 );
                    ys.resize( xySize - 1 );
                }
            }
        }
        xs.insert( xs.end( ),  xs2.begin( ), xs2.end( ) );
        ys.insert( ys.end( ),  ys2.begin( ), ys2.end( ) );
        xys1d1 = xys1d2;
    }
    delete xys1d1;

    return( new XYs1d( axes( ), ptwXY_interpolationLinLin, xs, ys ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void Regions1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
    std::string indent3 = a_writeInfo.incrementalIndent( indent2 );
    std::string attributes;

    if( a_embedded ) {
        attributes += a_writeInfo.addAttribute( GIDI_outerDomainValueChars, LUPI::Misc::doubleToShortestString( outerDomainValue( ) ) ); }
    else {
        if( a_inRegions ) {
            attributes = a_writeInfo.addAttribute( GIDI_indexChars, intToString( index( ) ) ); }
        else {
            if( label( ) != "" ) attributes = a_writeInfo.addAttribute( GIDI_labelChars, label( ) );
        }
    }

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    axes( ).toXMLList( a_writeInfo, indent2 );
    a_writeInfo.push_back( indent2 + "<function1ds>" );
    for( std::vector<Function1dForm *>::const_iterator iter = m_function1ds.begin( ); iter != m_function1ds.end( ); ++iter ) (*iter)->toXMLList_func( a_writeInfo, indent3, false, true );
    a_writeInfo.push_back( "</function1ds>" );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * This method calls the write method for each region of *this*.
 *
 * @param       a_file              [in]    The C FILE instance to write the data to.
 * @param       a_format            [in]    The format string passed to each region's write method.
 ***********************************************************************************************************/

void Regions1d::write( FILE *a_file, std::string const &a_format ) const {

    for( auto regionIter = m_function1ds.begin( ); regionIter != m_function1ds.end( ); ++regionIter ) {
        (*regionIter)->write( a_file, a_format );
    }
}

}               // End namespace Functions.

}               // End namespace GIDI.
