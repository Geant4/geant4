/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <LUPI.hpp>
#include <GIDI.hpp>
#include <HAPI.hpp>
#include <sstream>

static void mutualifyDomains( ptwXYPoints const *a_lhs, ptwXYPoints const *a_rhs, ptwXYPoints **ptwXY1, ptwXYPoints **ptwXY2 );

namespace GIDI {

namespace Functions {

/*! \class XYs1d
 * Class to store GNDS <**XYs1d**> node.
 */

/* *********************************************************************************************************//**
 * Constructor an empty XYs1d instance.
 *
 ***********************************************************************************************************/

XYs1d::XYs1d( ) :
        Function1dForm( GIDI_XYs1dChars, FormType::XYs1d, Axes(), ptwXY_interpolationLinLin, 0, 0.0 ) {

    double dummy[2];

    m_ptwXY = ptwXY_create2( nullptr, interpolation( ), 0, 0, 0, dummy, 0 );
}

/* *********************************************************************************************************//**
 * Constructor that creates an empty instance.
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 * @return
 ***********************************************************************************************************/

XYs1d::XYs1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_XYs1dChars, FormType::XYs1d, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

    double dummy[2];

    m_ptwXY = ptwXY_create2( nullptr, a_interpolation, 0, 0, 0, dummy, 0 );
}

/* *********************************************************************************************************//**
 * Constructor that create an instance from a list of doubles via **a_values** which must have an even number of values.
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_values              [in]    The data values as ( x_1, y_1, x_2, y_2, ..., x_n, y_n ). Must be an even number since they are n-pairs of xy values.
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 * @return
 ***********************************************************************************************************/

XYs1d::XYs1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, std::vector<double> const &a_values, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_XYs1dChars, FormType::XYs1d, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

    int64_t length = static_cast<int64_t>( a_values.size( ) ) / 2;

    m_ptwXY = ptwXY_create2( nullptr, a_interpolation, length, 0, length, a_values.data( ), 0 );
}

/* *********************************************************************************************************//**
 * Constructor that creates an instance from a list of x values and a list of y values. The size of **a_xs**
 * and **a_ys** must be the same.
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_interpolation       [in]    The interpolation along the outer most independent axis and the dependent axis.
 * @param a_xs                  [in]    The data values as ( x_1, x_2, ..., x_n ).
 * @param a_ys                  [in]    The data values as ( y_1, y_2, ..., y_n ).
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 * @return
 ***********************************************************************************************************/

XYs1d::XYs1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, std::vector<double> const &a_xs, std::vector<double> const &a_ys,
                int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_XYs1dChars, FormType::XYs1d, a_axes, a_interpolation, a_index, a_outerDomainValue ) {

    if( a_xs.size( ) != a_ys.size( ) ) throw Exception( "XYs1d::XYs1d: xs and ys not the same size" );
    int64_t length = static_cast<int64_t>( a_xs.size( ) );

    m_ptwXY = ptwXY_createFrom_Xs_Ys2( nullptr, a_interpolation, length, 0, length, a_xs.data( ), a_ys.data( ), 0 );
}

/* *********************************************************************************************************//**
 * Constructor that uses an existing **ptwXYPoints** instance. The **m_ptwXY** member is set to **ptwXYPoints** (i.e., this
 * XYs1d instance now owns the inputted **ptwXYPoints** instance).
 *
 * @param a_axes                [in]    The axes to copy for *this*.
 * @param a_ptwXY               [in]    The **ptwXYPoints** instance that *this* takes ownership of.*
 * @param a_index               [in]    Currently not used.
 * @param a_outerDomainValue    [in]    If embedded in a higher dimensional function, the value of the domain of the next higher dimension.
 * @return
 ***********************************************************************************************************/

XYs1d::XYs1d( Axes const &a_axes, ptwXYPoints *a_ptwXY, int a_index, double a_outerDomainValue ) :
        Function1dForm( GIDI_XYs1dChars, FormType::XYs1d, a_axes, ptwXY_getInterpolation( a_ptwXY ), a_index, a_outerDomainValue ),
        m_ptwXY( a_ptwXY ) {

}

/* *********************************************************************************************************//**
 * Constructs the instance from a HAPI::Node instance.
 *
 * @param a_construction        [in]    Used to pass user options for parsing.
 * @param a_node                [in]    The XYs1d HAPI::Node to be parsed and to construct the instance.
 * @param a_setupInfo           [in]    Information create my the Protare constructor to help in parsing.
 * @param a_parent              [in]    If imbedded in a two dimensional function, its pointers.
 ***********************************************************************************************************/

XYs1d::XYs1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent ) :
        Function1dForm( a_construction, a_node, a_setupInfo, FormType::XYs1d, a_parent ) {

    HAPI::Node values = a_node.child( GIDI_valuesChars );
    nf_Buffer<double> vals;
    parseValuesOfDoubles( a_construction, values, a_setupInfo, vals );

    int primarySize = vals.size() / 2, secondarySize = 0;
    double *dvals = new double[vals.size()];                  // Not sure we really need a copy here.
    for( size_t idx = 0; idx < vals.size(); idx++ ) dvals[idx] = vals[idx];
    m_ptwXY = ptwXY_create( NULL, interpolation( ), interpolationString( ).c_str( ), 12, 1e-3, primarySize, secondarySize, primarySize, dvals, 0 );
    delete[] dvals;
    if( m_ptwXY == nullptr ) throw Exception( "XYs1d::XYs1d: ptwXY_fromString failed" );
}

/* *********************************************************************************************************//**
 * The XYs1d copy constructor.
 *
 * @param a_XYs1d               [in]    The XYs1d instance to copy.
 ***********************************************************************************************************/

XYs1d::XYs1d( XYs1d const &a_XYs1d ) :
        Function1dForm( a_XYs1d ),
        m_ptwXY( nullptr ) {

    m_ptwXY = ptwXY_clone2( nullptr, a_XYs1d.ptwXY( ) );
    if( m_ptwXY == nullptr ) throw Exception( "XYs1d::XYs1d:2: ptwXY_clone2 failed" );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

XYs1d::~XYs1d( ) {

    ptwXY_free( m_ptwXY );
}

/* *********************************************************************************************************//**
 * The assignment operator. This method sets the members of *this* to those of *a_rhs* except for those
 * not set by base classes.
 *
 * @param a_rhs                     [in]    Instance whose member are used to set the members of *this*.
 ***********************************************************************************************************/

XYs1d &XYs1d::operator=( XYs1d const &a_rhs ) {

    if( this != &a_rhs ) {
        Function1dForm::operator=( a_rhs );

        LUPI::StatusMessageReporting smr;
        m_ptwXY = ptwXY_clone2( smr.smr( ), a_rhs.ptwXY( ) );
        if( m_ptwXY == nullptr ) throw Exception( smr.constructMessage( "XYs1d::operator=", -1, true ) );
    }

    return( *this );
}

/* *********************************************************************************************************//**
 * The element access methods that returns a point (i.e., an x1, y pair).
 *
 * @param a_index       [in]    The index of the element to access.
 * @return                      The x1, y values at **a_index**.
 ***********************************************************************************************************/

std::pair<double, double> XYs1d::operator[]( std::size_t a_index ) const {

    if( (int64_t) a_index >= ptwXY_length( nullptr, m_ptwXY ) ) throw Exception( "XYs1d::operator[]: index out of bounds." );

    ptwXYPoint *point = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, (int64_t) a_index );
    std::pair<double, double> CPPPoint( point->x, point->y );

    return( CPPPoint );
}

/* *********************************************************************************************************//**
 * Adds two **XYs1d** instances and returns the result.
 * 
 * @param a_rhs         [in]    The **XYs1d** instance to add to this instance.
 * @return                      An **XYs1d** instance that is the sum of this and *a_rhs*.
 ***********************************************************************************************************/

XYs1d XYs1d::operator+( XYs1d const &a_rhs ) const {

    XYs1d xys1d( *this );

    xys1d += a_rhs;
    return( xys1d );
}

/* *********************************************************************************************************//**
 * Adds an **XYs1d** instance to this.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to add to this instance.
 * @return                      This instance.
 ***********************************************************************************************************/

XYs1d &XYs1d::operator+=( XYs1d const &a_rhs ) {

    ptwXYPoints *sum, *ptwXY1, *ptwXY2;
    LUPI::StatusMessageReporting smr;

    mutualifyDomains( m_ptwXY, a_rhs.ptwXY( ), &ptwXY1, &ptwXY2 );
    sum = ptwXY_add_ptwXY( smr.smr( ), ptwXY1, ptwXY2 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    if( sum == nullptr ) 
        throw Exception( smr.constructMessage( "XYs1d::operator+=: ptwXY_add_ptwXY failed. ", -1, true ) );

    ptwXY_free( m_ptwXY );
    m_ptwXY = sum;

    return( *this );
}

/* *********************************************************************************************************//**
 * Subtracts two **XYs1d** instances and returns the result.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to substract from this instance.
 * @return                      An **XYs1d** instance that is the difference of this and *a_rhs*.
 ***********************************************************************************************************/

XYs1d XYs1d::operator-( XYs1d const &a_rhs ) const {

    XYs1d xys1d( *this );

    xys1d -= a_rhs;
    return( xys1d );
}

/* *********************************************************************************************************//**
 * Subtracts an **XYs1d** instance from this.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to subtract from this instance.
 * @return                      This instance.
 ***********************************************************************************************************/

XYs1d &XYs1d::operator-=( XYs1d const &a_rhs ) {

    ptwXYPoints *sum, *ptwXY1, *ptwXY2;
    LUPI::StatusMessageReporting smr;

    mutualifyDomains( m_ptwXY, a_rhs.ptwXY( ), &ptwXY1, &ptwXY2 );
    sum = ptwXY_sub_ptwXY( smr.smr( ), ptwXY1, ptwXY2 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    if( sum == nullptr ) 
        throw Exception( smr.constructMessage( "XYs1d::operator-=: ptwXY_sub_ptwXY failed. ", -1, true ) );

    ptwXY_free( m_ptwXY );
    m_ptwXY = sum;

    return( *this );
}

/* *********************************************************************************************************//**
 * Multiplies *this* by a double and returns the result.
 *
 * @param a_value       [in]    Number with this instance.
 *
 * @return                      An **XYs1d** instance that is the product of this and *a_rhs*.
 ***********************************************************************************************************/

XYs1d XYs1d::operator*( double a_value ) const {

    XYs1d xys1d( *this );

    xys1d *= a_value;
    return( xys1d );
}

/* *********************************************************************************************************//**
 * Multiplies *this* with another **XYs1d** instances and returns the result.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to multiply with this instance.
 *
 * @return                      An **XYs1d** instance that is the product of this and *a_rhs*.
 ***********************************************************************************************************/

XYs1d XYs1d::operator*( XYs1d const &a_rhs ) const {

    XYs1d xys1d( *this );

    xys1d *= a_rhs;
    return( xys1d );
}

/* *********************************************************************************************************//**
 * Multiplies a double with **this**.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to subtract from this instance.
 *
 * @return                      This instance.
 ***********************************************************************************************************/

XYs1d &XYs1d::operator*=( double a_value ) {

    LUPI::StatusMessageReporting smr;

    nfu_status status = ptwXY_mul_double( smr.smr( ), m_ptwXY, a_value );
    if( status != nfu_Okay ) 
        throw Exception( smr.constructMessage( "XYs1d::operator*=: ptwXY_mul_double failed. ", -1, true ) );

    return( *this );
}

/* *********************************************************************************************************//**
 * Multiplies **this** with an **XYs1d** instance.
 *
 * @param a_rhs         [in]    The **XYs1d** instance to subtract from this instance.
 *
 * @return                      This instance.
 ***********************************************************************************************************/

XYs1d &XYs1d::operator*=( XYs1d const &a_rhs ) {

    ptwXYPoints *mul, *ptwXY1, *ptwXY2;
    LUPI::StatusMessageReporting smr;

    mutualifyDomains( m_ptwXY, a_rhs.ptwXY( ), &ptwXY1, &ptwXY2 );
    mul = ptwXY_mul2_ptwXY( smr.smr( ), ptwXY1, ptwXY2 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    if( mul == nullptr )
        throw Exception( smr.constructMessage( "XYs1d::operator*=: ptwXY_mul2_ptwXY failed. ", -1, true ) );

    ptwXY_free( m_ptwXY );
    m_ptwXY = mul;

    return( *this );
}

/* *********************************************************************************************************//**
 * Returns the list of **x1** values of this.
 *
 * @return              The **x1** values of this.
 ***********************************************************************************************************/

std::vector<double> XYs1d::xs( ) const {

    int64_t n1 = size( );
    std::vector<double> _xs( n1, 0. );

    for( int64_t i1 = 0; i1 < n1; ++i1 ) {
        ptwXYPoint const *point = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, i1 );

        _xs[i1] = point->x;
    }
    return( _xs );
}

/* *********************************************************************************************************//**
 * Returns the list of **y** values of this.
 *
 * @return              The **y** values of this.
 ***********************************************************************************************************/

std::vector<double> XYs1d::ys( ) const {

    int64_t n1 = size( );
    std::vector<double> _ys( n1, 0. );

    for( int64_t i1 = 0; i1 < n1; ++i1 ) {
        ptwXYPoint const *point = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, i1 );

        _ys[i1] = point->y;
    }
    return( _ys );
}

/* *********************************************************************************************************//**
 * Returns a list of values that are this **y** mapped to the **x1** values in **a_xs**.
 *
 * @param a_xs              [in]    The list of **x1** values to map this' **y** values to.
 * @param a_offset          [out]   The index of the first value in **a_xs** where this starts.
 * @return                          The liist of **y** values.
 ***********************************************************************************************************/

std::vector<double> XYs1d::ysMappedToXs( std::vector<double> const &a_xs, std::size_t *a_offset ) const {

    int64_t n1 = size( ), i2, n2 = a_xs.size( );
    std::vector<double> _ys;

    *a_offset = 0;
    if( n1 == 0 ) return( _ys );

    ptwXYPoint const *point1 = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, 0 );
    for( i2 = 0; i2 < n2; ++i2 ) if( point1->x <= a_xs[i2] ) break;
    *a_offset = i2;
    if( i2 == n2 ) return( _ys );

    for( int64_t i1 = 1; i1 < n1; ++i1 ) {
        ptwXYPoint const *point2 = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, i1 );

        while( i2 < n2 ) {
            double x = a_xs[i2], y;
            if( x > point2->x ) break;           // Happens because of round off errors. Need to fix.

            ptwXY_interpolatePoint( nullptr, ptwXY_interpolationLinLin, x, &y, point1->x, point1->y, point2->x, point2->y );
            _ys.push_back( y );
            ++i2;
            if( x >= point2->x ) break;         // This check can fail hence check above.
        }
        point1 = point2;
    }

    return( _ys );
}

/* *********************************************************************************************************//**
 * Returns an **XYs1d** instance that is a domain slice of *this* from **domainMin** to **domainMax**.
 * If the minimum (maximum) of *this* is greater (less) than **domainMin** (**domainMax**) that domain is unchanged.
 * If *this*'s minimum is less than **domainMin**, **a_fill** is true and there is no point at **domainMin**
 * then an interpolated point is added at **domainMin**. Similarly for **domainMax**.
 *
 * @param a_domainMin           [in]    The minimum domain for the returned instance if *this*' minimum is less than **domainMin*.
 * @param a_domainMax           [in]    The maximum domain for the returned instance if *this*' maximum is less than **domainMax*.
 * @param a_fill                [in]    If true, values at **domainMin** and **domainMax** are added if needed.
 *
 * @return                              An **XYs1d** instance.
 ***********************************************************************************************************/

XYs1d XYs1d::domainSlice( double a_domainMin, double a_domainMax, bool a_fill ) const {

    LUPI::StatusMessageReporting smr;
    ptwXYPoints *ptwXY1 = ptwXY_clone2( smr.smr( ), m_ptwXY );
    if( ptwXY1 == nullptr ) throw Exception( smr.constructMessage( "XYs1d::domainSlice", -1, true ) );

    ptwXYPoints *ptwXYSliced = ptwXY_domainSlice( nullptr, ptwXY1, a_domainMin, a_domainMax, 10, a_fill ? 1 : 0 );
    ptwXY_free( ptwXY1 );
    if( ptwXYSliced == nullptr ) throw Exception( smr.constructMessage( "XYs1d::domainSlice", -1, true ) );

    return( XYs1d( axes( ), ptwXYSliced, index( ), outerDomainValue( ) ) );
}

/* *********************************************************************************************************//**
 * Returns an **XYs1d** instance that is this from its domain minimum to **domainMax**.
 *
 * @param a_domainMax           [in]    The maximum domain
 * @return                              An **XYs1d** instance.
 ***********************************************************************************************************/

XYs1d XYs1d::domainSliceMax( double a_domainMax ) const {

    ptwXYPoints *_ptwXY = ptwXY_clone2( nullptr, m_ptwXY );
    if( _ptwXY == nullptr ) throw Exception( "domainSliceMax: ptwXY_clone2 failed" );

    ptwXYPoints *ptwXYSliced = ptwXY_domainMaxSlice( nullptr, _ptwXY, a_domainMax, 10, 1 );
    ptwXY_free( _ptwXY );
    if( ptwXYSliced == nullptr ) throw Exception( "domainSliceMax: ptwXY_domainMaxSlice failed" );

    return( XYs1d( axes( ), ptwXYSliced ) );
}

/* *********************************************************************************************************//**
 * The **y** value of this at the domain value **a_x1**.
 *
 * @param a_x1          [in]    Domain value to evaluate this at.
 * @return                      The value of this at the domain value **a_x1**.
 ***********************************************************************************************************/

double XYs1d::evaluate( double a_x1 ) const {

    std::size_t length = ptwXY_length( nullptr, m_ptwXY );
    if( length == 0 ) throw Exception( "XYs1d::evaluate: XYs1d has no datum." );

    ptwXYPoint *point = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, 0 );
    if( point->x >= a_x1 ) return( point->y );

    point = ptwXY_getPointAtIndex_Unsafely( m_ptwXY, length - 1 );
    if( point->x <= a_x1 ) return( point->y );

    double y;
    nfu_status status = ptwXY_getValueAtX( nullptr, m_ptwXY, a_x1, &y );
    if( status != nfu_Okay ) throw Exception( "XYs1d::evaluate: status != nfu_Okay" );
    return( y );
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

void XYs1d::mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const {
    
    if( a_Xs.size( ) != a_results.size( ) ) throw Exception( "XYs1d::mapToXsAndAdd: a_Xs.size( ) != a_results.size( )" );
    if( a_offset < 0 ) throw Exception( "XYs1d::mapToXsAndAdd: a_offset < 0." );

    LUPI::StatusMessageReporting smr;
    int64_t length = static_cast<int64_t>( a_Xs.size( ) );

    nfu_status status = ptwXY_mapToXsAndAdd( smr.smr( ), m_ptwXY, a_offset, length, a_Xs.data( ), a_results.data( ), a_scaleFactor );
    if( ( status != nfu_Okay ) && ( status != nfu_tooFewPoints ) )
        throw Exception( smr.constructMessage( "XYs1d::mapToXsAndAdd", -1, true ) );
}

/* *********************************************************************************************************//**
 * This methods returns an XYs1d representation of *this*. The calling function owns the created instance and is responible
 * for freeing it.
 *
 * @param   a_asLinlin          [in]    If **true**, the inpolatation of the returned XYs1d instance will always be lin-lin. Otherwise,
 *                                      the interpolation depends on the child 1d functions. This argument is not needed or used for this class.
 * @param   a_accuracy          [in]    The accuracy use to convert the data to lin=lin interpolation if needed. This argument is not needed or used for this cl
 * @param   a_lowerEps          [in]    The dulling of the point at the domain minimum.
 * @param   a_upperEps          [in]    The dulling of the point at the domain maximum.
 *
 * @return                                  A pointer to an  XYs1d instance that must be freed by the calling function.
 ***********************************************************************************************************/

XYs1d *XYs1d::asXYs1d( LUPI_maybeUnused bool a_asLinlin, double a_accuracy, double a_lowerEps, double a_upperEps ) const {

    ptwXYPoints *ptwXY2 = nullptr;

    if( m_ptwXY->interpolation == ptwXY_interpolationFlat ) {
        ptwXY2 = ptwXY_flatInterpolationToLinear( nullptr, m_ptwXY, a_lowerEps, a_upperEps ); }
    else {
        ptwXY2 = ptwXY_toOtherInterpolation( nullptr, const_cast<ptwXYPoints *>( m_ptwXY ), ptwXY_interpolationLinLin, a_accuracy );
    }
    

    if( ptwXY2 == nullptr ) return( nullptr );

    return( new XYs1d( axes( ), ptwXY2 ) );
}

/* *********************************************************************************************************//**
 * This method returns **this** that is norimalized (i.e., its integral if 1.). The returned value is the integral
 * of **this** before being normalized.
 *
 * @return                          A double.
 ***********************************************************************************************************/

double XYs1d::integrate( double a_dommainMin, double a_dommainMax ) {

    double value = 0.0;
    LUPI::StatusMessageReporting smr;

    nfu_status status = ptwXY_integrate( smr.smr( ), m_ptwXY, a_dommainMin, a_dommainMax, &value );
    if( status != nfu_Okay ) throw Exception( smr.constructMessage( "XYs1d::normalize", -1, true ) );

    return( value );
}

/* *********************************************************************************************************//**
 * This method returns **this** that is norimalized (i.e., its integral if 1.). The returned value is the integral
 * of **this** before being normalized.
 *
 * @return                          A double.
 ***********************************************************************************************************/

double XYs1d::normalize( ) {

    double value = 0.0;
    LUPI::StatusMessageReporting smr;

    nfu_status status = ptwXY_integrateDomain( smr.smr( ), m_ptwXY, &value );
    if( status != nfu_Okay ) throw Exception( smr.constructMessage( "XYs1d::normalize", -1, true ) );

    status = ptwXY_normalize( smr.smr( ), m_ptwXY );
    if( status != nfu_Okay ) throw Exception( smr.constructMessage( "XYs1d::normalize", -1, true ) );

    return( value );
}

/* *********************************************************************************************************//**
 * This method returns a *Xs_pdf_cdf1d* representation of *this*.
 *
 * @return                          An instance of *Xs_pdf_cdf1d*.
 ***********************************************************************************************************/

Xs_pdf_cdf1d XYs1d::toXs_pdf_cdf1d( ) {
// FIXME, need to check that two consecutive cdf values are not the same.

    std::vector<double> xs1 = xs( );
    std::vector<double> pdf1 = ys( );

    LUPI::StatusMessageReporting smr;
    ptwXPoints *ptwX_cdf = ptwXY_runningIntegral( smr.smr( ), m_ptwXY );
    if( ptwX_cdf == nullptr ) throw Exception( smr.constructMessage( "XYs1d::toXs_pdf_cdf1d", -1, true ) );

    std::vector<double> cdf1( ptwX_cdf->length );
    for( int64_t index = 0; index < ptwX_cdf->length; ++index ) cdf1[index] = ptwX_cdf->points[index];
    ptwX_free( ptwX_cdf );

    return( Xs_pdf_cdf1d( GIDI::Axes( ), ptwXY_interpolationLinLin, xs1, pdf1, cdf1 ) );
}

/* *********************************************************************************************************//**
 * Fills the argument *a_writeInfo* with the XML lines that represent *this*. Recursively enters each sub-node.
 *
 * @param       a_writeInfo         [in/out]    Instance containing incremental indentation and other information and stores the appended lines.
 * @param       a_indent            [in]        The amount to indent *this* node.
 * @param       a_embedded          [in]        If *true*, *this* function is embedded in a higher dimensional function.
 * @param       a_inRegions         [in]        If *true*, *this* is in a Regions1d container.
 ***********************************************************************************************************/

void XYs1d::toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const {

    std::string indent2 = a_writeInfo.incrementalIndent( a_indent );
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

    if( interpolation( ) != ptwXY_interpolationLinLin ) attributes += a_writeInfo.addAttribute( GIDI_interpolationChars, interpolationString( ) );

    a_writeInfo.addNodeStarter( a_indent, moniker( ), attributes );
    axes( ).toXMLList( a_writeInfo, indent2 );

    std::vector<double> doubles( 2 * size( ) );
    for( std::size_t i1 = 0; i1 < size( ); ++i1 ) {
        std::pair<double, double> point = (*this)[i1];
        doubles[2*i1] = point.first;
        doubles[2*i1+1] = point.second;
    }

    doublesToXMLList( a_writeInfo, indent2, doubles );
    a_writeInfo.addNodeEnder( moniker( ) );
}

/* *********************************************************************************************************//**
 * Writes the (x,y) values to *a_file*. The format string must have two double conversion specifiers (e.g., "    %12.3e %.6f").
 *
 * @param       a_file              [in]    The C FILE instance to write the data to.
 * @param       a_format            [in]    The format string passed to the C printf function.
 ***********************************************************************************************************/

void XYs1d::write( FILE *a_file, std::string const &a_format ) const {

    ptwXY_simpleWrite( m_ptwXY, a_file, a_format.c_str( ) );
}

/* *********************************************************************************************************//**
 * This is a factory function for the XYs1d class that creates a XYs1d instance with the constant value
 * of *a_value* for the domain from *a_domainMin* to *a_domainMax*.
 *
 * @param       a_axes              [in]    The axes for the function.
 * @param       a_domainMin         [in]    The minimum of the domain for which the function is defined.
 * @param       a_domainMax         [in]    The maximum of the domain for which the function is defined.
 * @param       a_value             [in]    The y-value of the function over its domain.
 ***********************************************************************************************************/

XYs1d *XYs1d::makeConstantXYs1d( Axes const &a_axes, double a_domainMin, double a_domainMax, double a_value ) {

    std::vector<double> points( 4 );

    points[0] = a_domainMin;
    points[1] = a_value;
    points[2] = a_domainMax;
    points[3] = a_value;

     return( new XYs1d( a_axes, ptwXY_interpolationLinLin, points ) );
}

}               // End namespace Functions.

}               // End namespace GIDI.

/* *********************************************************************************************************//**
 * Returns to instances of **ptwXYPoints** that are the mutualified domains of **a_lhs** and **a_rhs**.
 *
 * @param a_lhs             [in]    One of the instances used to mutualify domains.
 * @param a_rhs             [in]    One of the instances used to mutualify domains.
 * @param a_ptwXY1          [out]   The mutualified domain for **a_lhs**.
 * @param a_ptwXY2          [out]   The mutualified domain for **a_rhs**.
 ***********************************************************************************************************/

static void mutualifyDomains( ptwXYPoints const *a_lhs, ptwXYPoints const *a_rhs, ptwXYPoints **a_ptwXY1, ptwXYPoints **a_ptwXY2 ) {

    double lowerEps = 1e-12, upperEps = 1e-12;

    *a_ptwXY1 = ptwXY_clone2( nullptr, a_lhs );
    if( *a_ptwXY1 == nullptr ) throw GIDI::Exception( "mutualifyDomains: ptwXY_clone2 failed for a_ptwXY1" );

    *a_ptwXY2 = ptwXY_clone2( nullptr, a_rhs );
    if( *a_ptwXY2 == nullptr ) {
        ptwXY_free( *a_ptwXY1 );
        throw GIDI::Exception( "mutualifyDomains: ptwXY_clone2 failed form a_ptwXY2" );
    }

    nfu_status status = ptwXY_mutualifyDomains( nullptr, *a_ptwXY1, lowerEps, upperEps, 1, *a_ptwXY2, lowerEps, upperEps, 1 );
    if( status != nfu_Okay ) {
        ptwXY_free( *a_ptwXY1 );
        ptwXY_free( *a_ptwXY2 );
        throw GIDI::Exception( "XYs1d::operator(+|-)=: mutualifyDomains in ptwXY_mutualifyDomains" );
    }
}
