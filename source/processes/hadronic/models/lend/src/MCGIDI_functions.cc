/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <typeinfo>

#include "MCGIDI.hpp"

namespace MCGIDI {

namespace Functions {

/*
============================================================
======================= FunctionBase =====================
============================================================
*/
LUPI_HOST_DEVICE FunctionBase::FunctionBase( ) :
        m_dimension( 0 ),
        m_domainMin( 0.0 ),
        m_domainMax( 0.0 ),
        m_interpolation( Interpolation::LINLIN ),
        m_outerDomainValue( 0.0 ) {

}
/*
============================================================
*/
LUPI_HOST FunctionBase::FunctionBase( GIDI::Functions::FunctionForm const &a_function ) :
        m_dimension( a_function.dimension( ) ),
        m_domainMin( a_function.domainMin( ) ),
        m_domainMax( a_function.domainMax( ) ),
        m_interpolation( GIDI2MCGIDI_interpolation( a_function.interpolation( ) ) ),
        m_outerDomainValue( a_function.outerDomainValue( ) ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE FunctionBase::FunctionBase( int a_dimension, double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue ) :
        m_dimension( a_dimension ),
        m_domainMin( a_domainMin ),
        m_domainMax( a_domainMax ),
        m_interpolation( a_interpolation ),
        m_outerDomainValue( a_outerDomainValue ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE FunctionBase::~FunctionBase( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void FunctionBase::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    DATA_MEMBER_INT( m_dimension, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_domainMin, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_domainMax, a_buffer, a_mode );

    int interpolation = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_interpolation ) {
        case Interpolation::FLAT :
            break;
        case Interpolation::LINLIN :
            interpolation = 1;
            break;
        case Interpolation::LINLOG :
            interpolation = 2;
            break;
        case Interpolation::LOGLIN :
            interpolation = 3;
            break;
        case Interpolation::LOGLOG :
            interpolation = 4;
            break;
        case Interpolation::OTHER :
            interpolation = 5;
            break;
        }
    }
    DATA_MEMBER_INT( interpolation, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( interpolation ) {
        case 0 :
            m_interpolation = Interpolation::FLAT;
            break;
        case 1 :
            m_interpolation = Interpolation::LINLIN;
            break;
        case 2 :
            m_interpolation = Interpolation::LINLOG;
            break;
        case 3 :
            m_interpolation = Interpolation::LOGLIN;
            break;
        case 4 :
            m_interpolation = Interpolation::LOGLOG;
            break;
        case 5 :
            m_interpolation = Interpolation::OTHER;
            break;
        }
    }

    DATA_MEMBER_DOUBLE( m_outerDomainValue, a_buffer, a_mode );
}

/*
============================================================
========================= Function1d =======================
============================================================
*/
LUPI_HOST_DEVICE Function1d::Function1d( ) :
        m_type( Function1dType::none ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE Function1d::Function1d( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue ) :
        FunctionBase( 1, a_domainMin, a_domainMax, a_interpolation, a_outerDomainValue ),
        m_type( Function1dType::none ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Function1d::~Function1d( ) {

}

/* *********************************************************************************************************//**
 * Returns a String representation of the **Function1d** type of *this*.
 *
 * @return                          A String instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE String Function1d::typeString( ) const {

    String typeStr( "Function1d::" );

    switch( m_type ) {
    case Function1dType::none :
        typeStr += "none";
        break;
    case Function1dType::constant :
        typeStr += "constant";
        break;
    case Function1dType::XYs :
        typeStr += "XYs";
        break;
    case Function1dType::polyomial :
        typeStr += "polyomial";
        break;
    case Function1dType::gridded :
        typeStr += "gridded";
        break;
    case Function1dType::regions :
        typeStr += "regions";
        break;
    case Function1dType::branching :
        typeStr += "branching";
        break;
    case Function1dType::TerrellFissionNeutronMultiplicityModel :
        typeStr += "TerrellFissionNeutronMultiplicityModel";
        break;
    }

    return( typeStr );
}

/* *********************************************************************************************************//**
 * Returns the value of the function at *a_x1*.
 *
 * @param a_x1                  [in]    The x-value to evaluate the function at.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Function1d::evaluate( double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case Function1dType::none:
        break;
    case Function1dType::constant:
        value = static_cast<Constant1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::XYs:
        value = static_cast<XYs1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::polyomial:
        value = static_cast<Polynomial1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::gridded:
        value = static_cast<Gridded1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::regions:
        value = static_cast<Regions1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::branching:
        value = static_cast<Branching1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::TerrellFissionNeutronMultiplicityModel:
        value = static_cast<TerrellFissionNeutronMultiplicityModel const *>( this )->evaluate( a_x1 );
        break;
    }

    return( value );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Function1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    FunctionBase::serialize( a_buffer, a_mode ); 

    int type = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_type ) {
        case Function1dType::none :
            break;
        case Function1dType::constant :
            type = 1;
            break;
        case Function1dType::XYs :
            type = 2;
            break;
        case Function1dType::polyomial :
            type = 3;
            break;
        case Function1dType::gridded :
            type = 4;
            break;
        case Function1dType::regions :
            type = 5;
            break;
        case Function1dType::branching :
            type = 6;
            break;
        case Function1dType::TerrellFissionNeutronMultiplicityModel :
            type = 7;
            break;
        }
    }
    DATA_MEMBER_INT( type, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            m_type = Function1dType::none;
            break;
        case 1 :
            m_type = Function1dType::constant;
            break;
        case 2 :
            m_type = Function1dType::XYs;
            break;
        case 3 :
            m_type = Function1dType::polyomial;
            break;
        case 4 :
            m_type = Function1dType::gridded;
            break;
        case 5 :
            m_type = Function1dType::regions;
            break;
        case 6 :
            m_type = Function1dType::branching;
            break;
        case 7 :
            m_type = Function1dType::TerrellFissionNeutronMultiplicityModel;
            break;
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the value of the function at *a_x1*.
 *
 * @param a_x1                  [in]    The x-value to evaluate the function at.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Function1d_d1::evaluate( double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case Function1dType::none:
        break;
    case Function1dType::constant:
        value = static_cast<Constant1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::XYs:
        value = static_cast<XYs1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::polyomial:
        value = static_cast<Polynomial1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::gridded:
        value = static_cast<Gridded1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::regions:
        value = static_cast<Regions1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::branching:
        value = static_cast<Branching1d const *>( this )->evaluate( a_x1 );
        break;
    default:
        String message( "Function1d_d1::evaluate: Unsupported Function1d_d1 " + typeString( ) );
        LUPI_THROW( message.c_str( ) );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * Returns the value of the function at *a_x1*.
 *
 * @param a_x1                  [in]    The x-value to evaluate the function at.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Function1d_d2::evaluate( double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case Function1dType::none:
        break;
    case Function1dType::constant:
        value = static_cast<Constant1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::XYs:
        value = static_cast<XYs1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::polyomial:
        value = static_cast<Polynomial1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::gridded:
        value = static_cast<Gridded1d const *>( this )->evaluate( a_x1 );
        break;
    case Function1dType::branching:
        value = static_cast<Branching1d const *>( this )->evaluate( a_x1 );
        break;
    default:                            // This should never happend.
        String message( "Function1d_d2::evaluate: Unsupported Function1d_d2 " + typeString( ) );
        LUPI_THROW( message.c_str( ) );
    }

    return( value );
}

/*
============================================================
======================== Constant1d ========================
============================================================
*/
LUPI_HOST_DEVICE Constant1d::Constant1d( ) :
        m_value( 0.0 ) {

    m_type = Function1dType::constant;
}
/*
============================================================
*/
LUPI_HOST_DEVICE Constant1d::Constant1d( double a_domainMin, double a_domainMax, double a_value, double a_outerDomainValue ) :
        Function1d_d2( a_domainMin, a_domainMax, Interpolation::FLAT, a_outerDomainValue ),
        m_value( a_value ) {

    m_type = Function1dType::constant;
}
/*
============================================================
*/
LUPI_HOST Constant1d::Constant1d( GIDI::Functions::Constant1d const &a_form1d ) :
        Function1d_d2( a_form1d.domainMin( ), a_form1d.domainMax( ), Interpolation::FLAT, a_form1d.outerDomainValue( ) ),
        m_value( a_form1d.value( ) ) {

    m_type = Function1dType::constant;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Constant1d::~Constant1d( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Constant1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_value, a_buffer, a_mode );
}

/*
============================================================
=========================== XYs1d ==========================
============================================================
*/
LUPI_HOST_DEVICE XYs1d::XYs1d( ) :
        m_Xs( ),
        m_Ys( ) {

    m_type = Function1dType::XYs;
}
/*
============================================================
*/
LUPI_HOST XYs1d::XYs1d( Interpolation a_interpolation, Vector<double> a_Xs, Vector<double> a_Ys, double a_outerDomainValue ) :
        Function1d_d2( a_Xs[0], a_Xs.back( ), a_interpolation, a_outerDomainValue ),
        m_Xs( a_Xs ),
        m_Ys( a_Ys ) {

    m_type = Function1dType::XYs;
}
/*
============================================================
*/
LUPI_HOST XYs1d::XYs1d( GIDI::Functions::XYs1d const &a_XYs1d ) :
        Function1d_d2( a_XYs1d.domainMin( ), a_XYs1d.domainMax( ), GIDI2MCGIDI_interpolation( a_XYs1d.interpolation( ) ), a_XYs1d.outerDomainValue( ) ) {

    m_type = Function1dType::XYs;
    std::size_t size = a_XYs1d.size( );

    m_Xs.resize( size );
    m_Ys.resize( size );
    for( std::size_t i1 = 0; i1 < size; ++i1 ) {
        std::pair<double, double> xy = a_XYs1d[i1];
        m_Xs[i1] = xy.first;
        m_Ys[i1] = xy.second;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE XYs1d::~XYs1d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double XYs1d::evaluate( double a_x1 ) const {

    int lower = binarySearchVector( a_x1, m_Xs );

    if( lower < 0 ) {
        if( lower == -2 ) return( m_Ys[0] );
        return( m_Ys.back( ) );
    }

    double evaluatedValue = 0.0;
    double y1 = m_Ys[lower];

    if( interpolation( ) == Interpolation::FLAT ) {
        evaluatedValue = y1; }
    else {
        double x1 = m_Xs[lower];
        double x2 = m_Xs[lower+1];
        double y2 = m_Ys[lower+1];

        if( interpolation( ) == Interpolation::LINLIN ) {
            double fraction = ( x2 - a_x1 ) / ( x2 - x1 );
            evaluatedValue = fraction * y1 + ( 1 - fraction ) * y2; }
        else if( interpolation( ) == Interpolation::LINLOG ) {
            double fraction = log( x2 / a_x1 ) / log( x2 / x1 );
            evaluatedValue = fraction * y1 + ( 1 - fraction ) * y2; }
        else if( interpolation( ) == Interpolation::LOGLIN ) {
            double fraction = ( x2 - a_x1 ) / ( x2 - x1 );
            evaluatedValue = exp( fraction * log( y1 ) + ( 1 - fraction ) * log( y2 ) ); }
        else if( interpolation( ) == Interpolation::LOGLOG ) {
            double fraction = log( x2 / a_x1 ) / log( x2 / x1 );
            evaluatedValue = exp( fraction * log( y1 ) + ( 1 - fraction ) * log( y2 ) ); }
        else {
            LUPI_THROW( "XYs1d::evaluate: unsupport interpolation." );
        }
    }

    return( evaluatedValue );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void XYs1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_Xs, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_Ys, a_buffer, a_mode );
}

/*
============================================================
======================= Polynomial1d =======================
============================================================
*/
LUPI_HOST_DEVICE Polynomial1d::Polynomial1d( ) :
        m_coefficients( ),
        m_coefficientsReversed( ) {

    m_type = Function1dType::polyomial;
}
/*
============================================================
*/
LUPI_HOST Polynomial1d::Polynomial1d( double a_domainMin, double a_domainMax, Vector<double> const &a_coefficients, double a_outerDomainValue ) :
        Function1d_d2( a_domainMin, a_domainMax, Interpolation::LINLIN, a_outerDomainValue ),
        m_coefficients( a_coefficients ) {

    m_type = Function1dType::polyomial;

    m_coefficientsReversed.reserve( m_coefficients.size( ) );
    for( Vector<double>::iterator iter = m_coefficients.begin( ); iter != m_coefficients.end( ); ++iter ) 
        m_coefficientsReversed.push_back( *iter );
}
/*
============================================================
*/
LUPI_HOST Polynomial1d::Polynomial1d( GIDI::Functions::Polynomial1d const &a_polynomial1d ) :
        Function1d_d2( a_polynomial1d.domainMin( ), a_polynomial1d.domainMax( ), Interpolation::LINLIN, a_polynomial1d.outerDomainValue( ) ) {

    m_type = Function1dType::polyomial;

    m_coefficients = a_polynomial1d.coefficients( );
    m_coefficientsReversed.reserve( m_coefficients.size( ) );
    for( Vector<double>::iterator iter = m_coefficients.begin( ); iter != m_coefficients.end( ); ++iter ) 
        m_coefficientsReversed.push_back( *iter );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Polynomial1d::~Polynomial1d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double Polynomial1d::evaluate( double a_x1 ) const {

    double d_value = 0;

    for( Vector<double>::const_iterator iter = m_coefficientsReversed.begin( ); iter != m_coefficientsReversed.end( ); ++iter ) {
        d_value = *iter + d_value * a_x1;
    }

    return( d_value );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Polynomial1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_coefficients, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_coefficientsReversed, a_buffer, a_mode );
}

/*
============================================================
========================= Gridded1d ========================
============================================================
*/
LUPI_HOST_DEVICE Gridded1d::Gridded1d( ) :
        m_grid( ),
        m_data( ) {

    m_type = Function1dType::gridded;
}
/*
============================================================
*/
LUPI_HOST Gridded1d::Gridded1d( GIDI::Functions::Gridded1d const &a_gridded1d ) :
        Function1d_d2( a_gridded1d.domainMin( ), a_gridded1d.domainMax( ), Interpolation::FLAT, a_gridded1d.outerDomainValue( ) ) {

    m_type = Function1dType::gridded;

    GIDI::Vector const &grid = a_gridded1d.grid( );
    m_grid.resize( grid.size( ) );
    for( std::size_t i1 = 0; i1 < grid.size( ); ++i1 ) m_grid[i1] = grid[i1];

    GIDI::Vector const &data = a_gridded1d.data( );
    m_data.resize( data.size( ) );
    for( std::size_t i1 = 0; i1 < data.size( ); ++i1 ) m_data[i1] = data[i1];
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Gridded1d::~Gridded1d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double Gridded1d::evaluate( double a_x1 ) const {

    return( m_data[binarySearchVector( a_x1, m_grid, true )] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Gridded1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_grid, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_data, a_buffer, a_mode );
}

/*
============================================================
========================= Regions1d ========================
============================================================
*/
LUPI_HOST_DEVICE Regions1d::Regions1d( ) :
        m_Xs( ),
        m_functions1d( ) {

    m_type = Function1dType::regions;
}
/*
============================================================
*/
LUPI_HOST Regions1d::Regions1d( GIDI::Functions::Regions1d const &a_regions1d ) :
        Function1d_d1( a_regions1d.domainMin( ), a_regions1d.domainMax( ), Interpolation::LINLIN, a_regions1d.outerDomainValue( ) ) {

    m_type = Function1dType::regions;

    m_Xs.reserve( a_regions1d.size( ) + 1 );
    m_functions1d.reserve( a_regions1d.size( ) );
    for( std::size_t i1 = 0; i1 < a_regions1d.size( ); ++i1 ) append( parseFunction1d_d2( a_regions1d[i1] ) );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Regions1d::~Regions1d( ) {

    for( std::size_t i1 = 0; i1 < m_functions1d.size( ); ++i1 ) delete m_functions1d[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE void Regions1d::append( Function1d_d2 *a_function1d ) {

    if( m_functions1d.size( ) == 0 ) m_Xs.push_back( a_function1d->domainMin( ) );
    m_Xs.push_back( a_function1d->domainMax( ) );

    m_functions1d.push_back( a_function1d );
}
/*
============================================================
*/
LUPI_HOST_DEVICE double Regions1d::evaluate( double a_x1 ) const {

    int lower = binarySearchVector( a_x1, m_Xs );

    if( lower < 0 ) {
        if( lower == -1 ) {                     // a_x1 > last value of m_Xs.
            return( m_functions1d.back( )->evaluate( a_x1 ) );
        }
        lower = 0;                              // a_x1 < last value of m_Xs.
    }

    return( m_functions1d[lower]->evaluate( a_x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Regions1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_Xs, a_buffer, a_mode );

    std::size_t vectorSize = m_functions1d.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_functions1d.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_functions1d.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_functions1d[vectorIndex] = serializeFunction1d_d2( a_buffer, a_mode, m_functions1d[vectorIndex] );
    }
}

/*
============================================================
======================== Branching1d =======================
============================================================
*/
LUPI_HOST_DEVICE Branching1d::Branching1d( ) {

}

/* *********************************************************************************************************//**
 * Function that parses a node one-d function node. Called from a Suite::parse instance.
 *
 * @param a_setupInfo       [in]    Information create my the Protare constructor to help in parsing.
 * @param a_form1d          [in]    The GIDI::Functions::Branching1d instance whose data is to be used to construct *this*.
 *
 * @return                          The parsed and constructed resonanceBackground region instance.
 ***********************************************************************************************************/

/*
============================================================
*/
LUPI_HOST Branching1d::Branching1d( SetupInfo &a_setupInfo, GIDI::Functions::Branching1d const &a_form1d ) :
        Function1d_d2( a_form1d.domainMin( ), a_form1d.domainMax( ), Interpolation::FLAT, 0.0 ),
        m_initialStateIndex( -1 ) {

    m_type = Function1dType::branching;

    std::map<std::string,int>::iterator iter = a_setupInfo.m_stateNamesToIndices.find( a_form1d.initialState( ) );
    if( iter == a_setupInfo.m_stateNamesToIndices.end( ) ) {
        std::string message( "Branching1d: initial state not found: pid = '" + a_form1d.initialState( ) + "'." );
        throw std::runtime_error( message.c_str( ) );
    }
    m_initialStateIndex = iter->second;

    a_setupInfo.m_initialStateIndex = m_initialStateIndex;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Branching1d::~Branching1d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double Branching1d::evaluate( LUPI_maybeUnused double a_x1 ) const {

    return( 0.0 );              // Returns 0 as needed by Product::sampleProducts.
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Branching1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_INT( m_initialStateIndex, a_buffer, a_mode );
}

/*
============================================================
========== TerrellFissionNeutronMultiplicityModel ==========
============================================================
*/

LUPI_HOST_DEVICE TerrellFissionNeutronMultiplicityModel::TerrellFissionNeutronMultiplicityModel( ) :
    m_multiplicity( nullptr ) {

    m_type = Function1dType::TerrellFissionNeutronMultiplicityModel;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST TerrellFissionNeutronMultiplicityModel::TerrellFissionNeutronMultiplicityModel( double a_width, Function1d_d1 *a_multiplicity ) :
        Function1d( a_multiplicity->domainMin( ), a_multiplicity->domainMax( ), a_multiplicity->interpolation( ), a_multiplicity->outerDomainValue( ) ),
        m_width( a_width ),
        m_multiplicity( a_multiplicity ) {

    m_type = Function1dType::TerrellFissionNeutronMultiplicityModel;
    if( a_width < 0.0 ) m_width = 1.079;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE TerrellFissionNeutronMultiplicityModel::~TerrellFissionNeutronMultiplicityModel( ) {

    delete m_multiplicity;
}

/* *********************************************************************************************************//**
 * Evaluated the *m_multiplicity* function at energy *a_energy*.
 *
 * @param a_energy              [in]        The energy of the projectile.
 *
 * @return                                  The number of emitted, prompt neutrons.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double TerrellFissionNeutronMultiplicityModel::evaluate( double a_energy ) const {

    return( m_multiplicity->evaluate( a_energy ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void TerrellFissionNeutronMultiplicityModel::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function1d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_width, a_buffer, a_mode );

    m_multiplicity = serializeFunction1d_d1( a_buffer, a_mode, m_multiplicity );
}

/*
============================================================
======================== Function2d ========================
============================================================
*/
LUPI_HOST_DEVICE Function2d::Function2d( ) :
        m_type( Function2dType::none ) {

}
/*
============================================================
*/
LUPI_HOST Function2d::Function2d( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue ) :
        FunctionBase( 2, a_domainMin, a_domainMax, a_interpolation, a_outerDomainValue ),
        m_type( Function2dType::none ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Function2d::~Function2d( ) {

}

/* *********************************************************************************************************//**
 * Returns a String representation of the **Function2d** type of *this*.
 *
 * @return                          A String instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE String Function2d::typeString( ) const {

    String typeStr( "Function2d::" );

    switch( m_type ) {
    case Function2dType::none :
        typeStr += "none";
        break;
    case Function2dType::XYs :
        typeStr += "XYs";
        break;
    }

    return( typeStr );
}

/* *********************************************************************************************************//**
 * Returns the value of the function f(x2, x1) at x2 = *a_x2* and x1 = *a_x1*.
 *
 * @param a_x2                  [in]    The x2 value.
 * @param a_x1                  [in]    The x1 value.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double Function2d::evaluate( double a_x2, double a_x1 ) const {

    return( static_cast<XYs2d const *>( this )->evaluate( a_x2, a_x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Function2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    FunctionBase::serialize( a_buffer, a_mode ); 

    int type = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_type ) {
        case Function2dType::none :
            break;
        case Function2dType::XYs :
            type = 1;
            break;
        }
    }
    DATA_MEMBER_INT( type, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            m_type = Function2dType::none;
            break;
        case 1 :
            m_type = Function2dType::XYs;
            break;
        }
    }
}

/*
============================================================
========================== XYs2d ===========================
============================================================
*/
LUPI_HOST_DEVICE XYs2d::XYs2d( ) :
        m_Xs( ),
        m_functions1d( ) {

    m_type = Function2dType::XYs;
}
/*
============================================================
*/
LUPI_HOST XYs2d::XYs2d( GIDI::Functions::XYs2d const &a_XYs2d ) :
        Function2d( a_XYs2d.domainMin( ), a_XYs2d.domainMax( ), GIDI2MCGIDI_interpolation( a_XYs2d.interpolation( ) ), a_XYs2d.outerDomainValue( ) ),
        m_Xs( a_XYs2d.Xs( ) ) {

    m_type = Function2dType::XYs;

    Vector<GIDI::Functions::Function1dForm *> const &function1ds = a_XYs2d.function1ds( );
    m_functions1d.resize( function1ds.size( ) );
    for( std::size_t i1 = 0; i1 < function1ds.size( ); ++i1 ) m_functions1d[i1] = parseFunction1d_d1( function1ds[i1] );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE XYs2d::~XYs2d( ) {

    for( std::size_t i1 = 0; i1 < m_functions1d.size( ); ++i1 ) delete m_functions1d[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE double XYs2d::evaluate( double a_x2, double a_x1 ) const {

    int lower = binarySearchVector( a_x2, m_Xs );
    double evaluatedValue = 0.0;

    if( lower < 0 ) {
        if( lower == -1 ) {               /* X2 > last value of Xs. */
            evaluatedValue = m_functions1d.back( )->evaluate( a_x1 ); }
        else {                          /* X2 < first value of Xs. */
            evaluatedValue = m_functions1d[0]->evaluate( a_x1 );
        } }
    else {
        double y1 = m_functions1d[lower]->evaluate( a_x1 );

        if( interpolation( ) == Interpolation::FLAT ) {
            evaluatedValue = y1; }
        else {
            double x1 = m_Xs[lower];
            double x2 = m_Xs[lower+1];
            double y2 = m_functions1d[lower+1]->evaluate( a_x1 );

            if( interpolation( ) == Interpolation::LINLIN ) {
                double fraction = ( x2 - a_x2 ) / ( x2 - x1 );
                evaluatedValue = fraction * y1 + ( 1 - fraction ) * y2; }
            else if( interpolation( ) == Interpolation::LINLOG ) {
                double fraction = log( x2 / a_x2 ) / log( x2 / x1 );
                evaluatedValue = fraction * y1 + ( 1 - fraction ) * y2; }
            else if( interpolation( ) == Interpolation::LOGLIN ) {
                double fraction = ( x2 - a_x2 ) / ( x2 - x1 );
                evaluatedValue = exp( fraction * log( y1 ) + ( 1 - fraction ) * log( y2 ) ); }
            else if( interpolation( ) == Interpolation::LOGLOG ) {
                double fraction = log( x2 / a_x2 ) / log( x2 / x1 );
                evaluatedValue = exp( fraction * log( y1 ) + ( 1 - fraction ) * log( y2 ) ); }
            else {
                LUPI_THROW( "XYs2d::evaluate: unsupport interpolation." );
            }
        }
    }

    return( evaluatedValue );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void XYs2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Function2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_Xs, a_buffer, a_mode );

    std::size_t vectorSize = m_functions1d.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_functions1d.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_functions1d.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_functions1d[vectorIndex] = serializeFunction1d_d1( a_buffer, a_mode, m_functions1d[vectorIndex] );
    }
}

/*
============================================================
========================== others ==========================
============================================================
*/
LUPI_HOST Function1d *parseMultiplicityFunction1d( SetupInfo &a_setupInfo, LUPI_maybeUnused Transporting::MC const &a_settings, GIDI::Suite const &a_suite ) {

    GIDI::Functions::Function1dForm const *form1d( a_suite.get<GIDI::Functions::Function1dForm>( 0 ) );

    if( form1d->type( ) == GIDI::FormType::branching1d ) return( new Branching1d( a_setupInfo, *static_cast<GIDI::Functions::Branching1d const *>( form1d ) ) );
    if( form1d->type( ) == GIDI::FormType::unspecified1d ) return( nullptr );

    return( parseFunction1d_d1( form1d ) );
}

/*
============================================================
*/
LUPI_HOST Function1d_d1 *parseFunction1d_d1( GIDI::Functions::Function1dForm const *a_form1d ) {

    GIDI::FormType type = a_form1d->type( );

    switch( type ) {
    case GIDI::FormType::constant1d :
        return( new Constant1d( *static_cast<GIDI::Functions::Constant1d const *>( a_form1d ) ) );
    case GIDI::FormType::XYs1d :
        return( new XYs1d( *static_cast<GIDI::Functions::XYs1d const *>( a_form1d ) ) );
    case GIDI::FormType::polynomial1d :
        return( new Polynomial1d( *static_cast<GIDI::Functions::Polynomial1d const *>( a_form1d ) ) );
    case GIDI::FormType::gridded1d :
        return( new Gridded1d( *static_cast<GIDI::Functions::Gridded1d const *>( a_form1d ) ) );
    case GIDI::FormType::regions1d :
        return( new Regions1d( *static_cast<GIDI::Functions::Regions1d const *>( a_form1d ) ) );
    case GIDI::FormType::reference1d :
        std::cout << "parseFunction1d_d1: Unsupported Function1d reference1d.";
        break;
    default :       // GIDI::FormTypes Legendre1d, reference1d, xs_pdFormType::cdf1d and resonancesWithBackground1d
        throw std::runtime_error( "Functions::parseFunction1d_d1: Unsupported Function1d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST Function1d_d2 *parseFunction1d_d2( GIDI::Functions::Function1dForm const *a_form1d ) {

    GIDI::FormType type = a_form1d->type( );

    switch( type ) {
    case GIDI::FormType::constant1d :
        return( new Constant1d( *static_cast<GIDI::Functions::Constant1d const *>( a_form1d ) ) );
    case GIDI::FormType::XYs1d :
        return( new XYs1d( *static_cast<GIDI::Functions::XYs1d const *>( a_form1d ) ) );
    case GIDI::FormType::polynomial1d :
        return( new Polynomial1d( *static_cast<GIDI::Functions::Polynomial1d const *>( a_form1d ) ) );
    case GIDI::FormType::gridded1d :
        return( new Gridded1d( *static_cast<GIDI::Functions::Gridded1d const *>( a_form1d ) ) );
    case GIDI::FormType::reference1d :
        std::cout << "Functions::parseFunction1d_d2: Unsupported Function1d reference1d.";
        break;
    default :
        throw std::runtime_error( "Functions::parseFunction1d_d2: Unsupported Function1d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST Function2d *parseFunction2d( GIDI::Functions::Function2dForm const *form2d ) {

    GIDI::FormType type = form2d->type( );

    switch( type ) {
    case GIDI::FormType::XYs2d :
        return( new XYs2d( *static_cast<GIDI::Functions::XYs2d const *>( form2d ) ) );
    default :
        throw std::runtime_error( "Functions::parseFunction2d: Unsupported Function2d" );
    }

    return( nullptr );
}

}           // End of namespace Functions.

/*
============================================================
============================================================
====================== probabilities =======================
============================================================
============================================================
*/
namespace Probabilities {

LUPI_HOST static nfu_status MCGIDI_NBodyPhaseSpacePDF_callback( statusMessageReporting *smr, double X, double *Y, void *argList );
LUPI_HOST static ProbabilityBase1d *ptwXY_To_Xs_pdf_cdf1d( ptwXYPoints *pdfXY );

/*
============================================================
===================== ProbabilityBase ======================
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase::ProbabilityBase( ) :
        m_Xs( ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase::ProbabilityBase( GIDI::Functions::FunctionForm const &a_probability ) :
        Functions::FunctionBase( a_probability ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase::ProbabilityBase( GIDI::Functions::FunctionForm const &a_probability, Vector<double> const &a_Xs ) :
        Functions::FunctionBase( a_probability ),
        m_Xs( a_Xs ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProbabilityBase::~ProbabilityBase( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProbabilityBase::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Functions::FunctionBase::serialize( a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_Xs, a_buffer, a_mode );
}

/*
============================================================
===================== ProbabilityBase1d ====================
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase1d::ProbabilityBase1d( ) :
        m_type( ProbabilityBase1dType::none ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase1d::ProbabilityBase1d( GIDI::Functions::FunctionForm const &a_probability, Vector<double> const &a_Xs ) :
        ProbabilityBase( a_probability, a_Xs ),
        m_type( ProbabilityBase1dType::none ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProbabilityBase1d::~ProbabilityBase1d( ) {

}

/* *********************************************************************************************************//**
 * Returns a String representation of the **ProbabilityBase1d** type of *this*.
 *
 * @return                          A String instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE String ProbabilityBase1d::typeString( ) const {
    
    String typeStr( "ProbabilityBase1d::" );
    
    switch( m_type ) {
    case ProbabilityBase1dType::none :
        typeStr += "none";
        break;
    case ProbabilityBase1dType::xs_pdf_cdf :
        typeStr += "xs_pdf_cdf";
        break;
    }
    
    return( typeStr );
}

/* *********************************************************************************************************//**
 * Returns the value of the function at *a_x1*.
 *
 * @param a_x1                  [in]    The x-value to evaluate the function at.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProbabilityBase1d::evaluate( double a_x1 ) const {

    return( static_cast<Xs_pdf_cdf1d const *>( this )->evaluate( a_x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProbabilityBase1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase::serialize( a_buffer, a_mode ); 

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_type ) {
        case ProbabilityBase1dType::none :
            break;
        case ProbabilityBase1dType::xs_pdf_cdf :
            type = 1;
            break;
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            m_type = ProbabilityBase1dType::none;
            break;
        case 1 :
            m_type = ProbabilityBase1dType::xs_pdf_cdf;
            break;
        }
    }
}

/*
============================================================
======================= Xs_pdf_cdf1d =======================
============================================================
*/
LUPI_HOST_DEVICE Xs_pdf_cdf1d::Xs_pdf_cdf1d( ) :
        m_pdf( ),
        m_cdf( ) {

    m_type = ProbabilityBase1dType::xs_pdf_cdf;
}
/*
============================================================
*/
LUPI_HOST Xs_pdf_cdf1d::Xs_pdf_cdf1d( GIDI::Functions::Xs_pdf_cdf1d const &a_xs_pdf_cdf1d ) :
        ProbabilityBase1d( a_xs_pdf_cdf1d, a_xs_pdf_cdf1d.Xs( ) ),
        m_pdf( a_xs_pdf_cdf1d.pdf( ) ),
        m_cdf( a_xs_pdf_cdf1d.cdf( ) ) {

    m_type = ProbabilityBase1dType::xs_pdf_cdf;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Xs_pdf_cdf1d::~Xs_pdf_cdf1d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double Xs_pdf_cdf1d::evaluate( double a_x1 ) const {

    int lower = binarySearchVector( a_x1, m_Xs );

    if( lower < 0 ) {
        if( lower == -2 ) return( m_pdf[0] );
        return( m_pdf.back( ) );
    }

    double fraction = ( a_x1 - m_Xs[lower] ) / ( m_Xs[lower+1] - m_Xs[lower] );
    return( ( 1. - fraction ) * m_pdf[lower] + fraction * m_pdf[lower+1] );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Xs_pdf_cdf1d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase1d::serialize( a_buffer, a_mode );

    DATA_MEMBER_VECTOR_DOUBLE( m_pdf, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_cdf, a_buffer, a_mode );
}

/*
============================================================
===================== ProbabilityBase2d ====================
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase2d::ProbabilityBase2d( ) :
        m_type( ProbabilityBase2dType::none ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase2d::ProbabilityBase2d( GIDI::Functions::FunctionForm const &a_probability ) :
        ProbabilityBase( a_probability ),
        m_type( ProbabilityBase2dType::none ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase2d::ProbabilityBase2d( GIDI::Functions::FunctionForm const &a_probability, Vector<double> const &a_Xs ) :
        ProbabilityBase( a_probability, a_Xs ),
        m_type( ProbabilityBase2dType::none ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProbabilityBase2d::~ProbabilityBase2d( ) {

}

/* *********************************************************************************************************//**
 * Returns a String representation of the **ProbabilityBase2d** type of *this*.
 *
 * @return                          A String instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE String ProbabilityBase2d::typeString( ) const {

    String typeStr( "ProbabilityBase2d::" );

    switch( m_type ) {
    case ProbabilityBase2dType::none :
        typeStr += "none";
        break;
    case ProbabilityBase2dType::XYs :
        typeStr += "XYs";
        break;
    case ProbabilityBase2dType::regions :
        typeStr += "regions";
        break;
    case ProbabilityBase2dType::isotropic :
        typeStr += "isotropic";
        break;
    case ProbabilityBase2dType::discreteGamma:
        typeStr += "discreteGamma";
        break;
    case ProbabilityBase2dType::primaryGamma:
        typeStr += "primaryGamma";
        break;
    case ProbabilityBase2dType::recoil:
        typeStr += "recoil";
        break;
    case ProbabilityBase2dType::NBodyPhaseSpace:
        typeStr += "NBodyPhaseSpace";
        break;
    case ProbabilityBase2dType::evaporation:
        typeStr += "evaporation";
        break;
    case ProbabilityBase2dType::generalEvaporation:
        typeStr += "generalEvaporation";
        break;
    case ProbabilityBase2dType::simpleMaxwellianFission:
        typeStr += "simpleMaxwellianFission";
        break;
    case ProbabilityBase2dType::Watt:
        typeStr += "Watt";
        break;
    case ProbabilityBase2dType::weightedFunctionals:
        typeStr += "weightedFunctionals";
        break;
    }

    return( typeStr );
}

/* *********************************************************************************************************//**
 * This method returns the value of pdf(x1|x2) at x1 of *a_x1* and x2 of *a_x2*.
 *
 * @param a_x2                  [in]        The value of x2.
 * @param a_x1                  [in]        The value of x1.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProbabilityBase2d::evaluate( double a_x2, double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::none:
        break;
    case ProbabilityBase2dType::weightedFunctionals:
        value = static_cast<WeightedFunctionals2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    default:
        value = static_cast<ProbabilityBase2d_d1 const *>( this )->evaluate( a_x2, a_x1 );
        break;
    }

    return( value );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProbabilityBase2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase::serialize( a_buffer, a_mode );
}

/* *********************************************************************************************************//**
 * Returns the value of the function evaluated at the point (*a_x2*, *a_x1*).
 *
 * @param a_x2                  [in]    Value of the outer most independent variable (i.e., *x2*).
 * @param a_x1                  [in]    Value of the inner most independent variable (i.e., *x1*).
 *
 * @return                              The value of the function at *a_x2* and *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProbabilityBase2d_d1::evaluate( double a_x2, double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::XYs:
    case ProbabilityBase2dType::isotropic:
    case ProbabilityBase2dType::discreteGamma:
    case ProbabilityBase2dType::primaryGamma:
    case ProbabilityBase2dType::recoil:
    case ProbabilityBase2dType::NBodyPhaseSpace:
    case ProbabilityBase2dType::evaporation:
    case ProbabilityBase2dType::generalEvaporation:
    case ProbabilityBase2dType::simpleMaxwellianFission:
    case ProbabilityBase2dType::Watt:
        value = static_cast<ProbabilityBase2d_d2 const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::regions:
        value = static_cast<Regions2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::none:
    case ProbabilityBase2dType::weightedFunctionals:
        LUPI_THROW( "ProbabilityBase2d_d1::evaluate: This should never happen." );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * Returns the value of the function evaluated at the point (*a_x2*, *a_x1*).
 *
 * @param a_x2                  [in]    Value of the outer most independent variable (i.e., *x2*).
 * @param a_x1                  [in]    Value of the inner most independent variable (i.e., *x1*).
 *
 * @return                              The value of the function at *a_x2* and *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProbabilityBase2d_d2::evaluate( double a_x2, double a_x1 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::XYs:
        value = static_cast<XYs2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::isotropic:
        value = static_cast<Isotropic2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::discreteGamma:
        value = static_cast<DiscreteGamma2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::primaryGamma:
        value = static_cast<PrimaryGamma2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::recoil:
        value = static_cast<Recoil2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::NBodyPhaseSpace:
        value = static_cast<NBodyPhaseSpace2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::evaporation:
        value = static_cast<Evaporation2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::generalEvaporation:
        value = static_cast<GeneralEvaporation2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::simpleMaxwellianFission:
        value = static_cast<SimpleMaxwellianFission2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::Watt:
        value = static_cast<Watt2d const *>( this )->evaluate( a_x2, a_x1 );
        break;
    case ProbabilityBase2dType::none:
    case ProbabilityBase2dType::weightedFunctionals:
    case ProbabilityBase2dType::regions:
        LUPI_THROW( "ProbabilityBase2d_d2::evaluate: This should never happen." );
    }

    return( value );
}

/*
============================================================
========================== XYs2d ===========================
============================================================
*/
LUPI_HOST_DEVICE XYs2d::XYs2d( ) :
        m_probabilities( ) {

    m_type = ProbabilityBase2dType::XYs;
}
/*
============================================================
*/
LUPI_HOST XYs2d::XYs2d( GIDI::Functions::XYs2d const &a_XYs2d ) :
        ProbabilityBase2d_d2( a_XYs2d, a_XYs2d.Xs( ) ) {

    m_type = ProbabilityBase2dType::XYs;

    Vector<GIDI::Functions::Function1dForm *> const &function1ds = a_XYs2d.function1ds( );
    m_probabilities.resize( function1ds.size( ) );
    for( std::size_t i1 = 0; i1 < function1ds.size( ); ++i1 ) m_probabilities[i1] = parseProbability1d( function1ds[i1] );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE XYs2d::~XYs2d( ) {

    for( std::size_t i1 = 0; i1 < m_probabilities.size( ); ++i1 ) delete m_probabilities[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE double XYs2d::evaluate( double a_x2, double a_x1 ) const {

    int lower = binarySearchVector( a_x2, m_Xs );

    if( lower < 0 ) {
        if( lower == -2 ) return( m_probabilities[0]->evaluate( a_x1 ) );
        return( m_probabilities.back( )->evaluate( a_x1 ) );
    }

    double fraction = ( a_x2 - m_Xs[lower] ) / ( m_Xs[lower+1] - m_Xs[lower] );
    double d_value = ( 1.0 - fraction ) * m_probabilities[lower]->evaluate( a_x1 ) + fraction * m_probabilities[lower+1]->evaluate( a_x1 );
    return( d_value );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void XYs2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );

    std::size_t vectorSize = m_probabilities.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_probabilities.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_probabilities.internalSize();

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_probabilities[vectorIndex] = serializeProbability1d( a_buffer, a_mode, m_probabilities[vectorIndex] );
    }
}

/*
============================================================
======================== Regions2d =========================
============================================================
*/
LUPI_HOST_DEVICE Regions2d::Regions2d( ) :
        m_probabilities( ) {

    m_type = ProbabilityBase2dType::regions;
}
/*
============================================================
*/
LUPI_HOST Regions2d::Regions2d( GIDI::Functions::Regions2d const &a_regions2d ) :
        ProbabilityBase2d_d1( a_regions2d, a_regions2d.Xs( ) ) {

    m_type = ProbabilityBase2dType::regions;

    Vector<GIDI::Functions::Function2dForm *> const &function2ds = a_regions2d.function2ds( );
    m_probabilities.resize( function2ds.size( ) );
    for( std::size_t i1 = 0; i1 < function2ds.size( ); ++i1 ) m_probabilities[i1] = parseProbability2d_d2( function2ds[i1], nullptr );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Regions2d::~Regions2d( ) {

    for( std::size_t i1 = 0; i1 < m_probabilities.size( ); ++i1 ) delete m_probabilities[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE double Regions2d::evaluate( double a_x2, double a_x1 ) const {

    int lower = binarySearchVector( a_x2, m_Xs );

    if( lower < 0 ) {
        if( lower == -1 ) {                         // a_x2 > last value of m_Xs.
            return( m_probabilities.back( )->evaluate( a_x2, a_x1 ) );
        }
        lower = 0;                                  // a_x2 < first value of m_Xs.
    }

    return( m_probabilities[lower]->evaluate( a_x2, a_x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Regions2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );

    std::size_t vectorSize = m_probabilities.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_probabilities.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_probabilities.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_probabilities[vectorIndex] = serializeProbability2d_d2( a_buffer, a_mode, m_probabilities[vectorIndex] );
    }
}


/*
============================================================
======================== Isotropic2d =======================
============================================================
*/
LUPI_HOST_DEVICE Isotropic2d::Isotropic2d( ) {

    m_type = ProbabilityBase2dType::isotropic;
}
/*
============================================================
*/
LUPI_HOST Isotropic2d::Isotropic2d( GIDI::Functions::Isotropic2d const &a_isotropic2d ) :
        ProbabilityBase2d_d2( a_isotropic2d ) {

    m_type = ProbabilityBase2dType::isotropic;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Isotropic2d::~Isotropic2d( ) {

}

/*
============================================================
====================== DiscreteGamma2d =====================
============================================================
*/
LUPI_HOST_DEVICE DiscreteGamma2d::DiscreteGamma2d( ) {

    m_type = ProbabilityBase2dType::discreteGamma;
}
/*
============================================================
*/
LUPI_HOST DiscreteGamma2d::DiscreteGamma2d( GIDI::Functions::DiscreteGamma2d const &a_discreteGamma2d ) :
        ProbabilityBase2d_d2( a_discreteGamma2d ),
        m_value( a_discreteGamma2d.value( ) ) {

    m_type = ProbabilityBase2dType::discreteGamma;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE DiscreteGamma2d::~DiscreteGamma2d( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void DiscreteGamma2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_value, a_buffer, a_mode );
}

/*
============================================================
====================== PrimaryGamma2d =====================
============================================================
*/
LUPI_HOST_DEVICE PrimaryGamma2d::PrimaryGamma2d( ) :
        m_primaryEnergy( 0.0 ),
        m_massFactor( 0.0 ),
        m_finalState( "" ),
        m_initialStateIndex( -1 ) {

    m_type = ProbabilityBase2dType::primaryGamma;
}
/*
============================================================
*/
LUPI_HOST PrimaryGamma2d::PrimaryGamma2d( GIDI::Functions::PrimaryGamma2d const &a_primaryGamma2d, SetupInfo *a_setupInfo ) :
        ProbabilityBase2d_d2( a_primaryGamma2d ),
        m_primaryEnergy( a_primaryGamma2d.value( ) ),
        m_massFactor( a_setupInfo->m_protare.targetMass( ) / ( a_setupInfo->m_protare.projectileMass( ) + a_setupInfo->m_protare.targetMass( ) ) ),
        m_finalState( a_primaryGamma2d.finalState( ).c_str( ) ),
        m_initialStateIndex( -1 ) {

    if( m_finalState.size( ) > 0 ) a_setupInfo->m_hasFinalStatePhotons = true;
    m_type = ProbabilityBase2dType::primaryGamma;

    if( a_primaryGamma2d.finalState( ) != "" ) {
        std::map<std::string,int>::iterator iter = a_setupInfo->m_stateNamesToIndices.find( a_primaryGamma2d.finalState( ) );
        if( iter == a_setupInfo->m_stateNamesToIndices.end( ) ) {
            std::string message( "Branching1d: final state not found: pid = '" + a_primaryGamma2d.finalState( ) + "'." );
            throw std::runtime_error( message.c_str( ) );
        }
        m_initialStateIndex = iter->second;

        a_setupInfo->m_initialStateIndex = m_initialStateIndex;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE PrimaryGamma2d::~PrimaryGamma2d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double PrimaryGamma2d::evaluate( double a_x2, double a_x1 ) const {

    double energy_out = m_primaryEnergy + a_x2 * m_massFactor;

// FIXME. I think this is correct but is it what we want.
    if( energy_out == a_x1 ) return( 1.0 );
    return( 0.0 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void PrimaryGamma2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_primaryEnergy, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_massFactor, a_buffer, a_mode );
    DATA_MEMBER_STRING( m_finalState, a_buffer, a_mode );
    DATA_MEMBER_INT( m_initialStateIndex, a_buffer, a_mode );
}

/*
============================================================
========================= Recoil2d =========================
============================================================
*/
LUPI_HOST_DEVICE Recoil2d::Recoil2d( ) :
        m_xlink( ) {

    m_type = ProbabilityBase2dType::recoil;
}
/*
============================================================
*/
LUPI_HOST Recoil2d::Recoil2d( GIDI::Functions::Recoil2d const &a_recoil2d ) :
        ProbabilityBase2d_d2( a_recoil2d ),
        m_xlink( a_recoil2d.xlink( ).c_str( ) ) {

    m_type = ProbabilityBase2dType::recoil;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Recoil2d::~Recoil2d( ) {

}
/*
============================================================
*/
LUPI_HOST_DEVICE double Recoil2d::evaluate( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const {

#if !defined(__NVCC__) && !defined(__HIP__)
    LUPI_THROW( "Recoil2d::evaluate: not implemented." );
#endif

    return( 0.0 );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Recoil2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_STRING( m_xlink, a_buffer, a_mode );
}

/*
============================================================
==================== NBodyPhaseSpace2d =====================
============================================================
*/
LUPI_HOST_DEVICE NBodyPhaseSpace2d::NBodyPhaseSpace2d( ) :
        m_numberOfProducts( 0 ),
        m_mass( 0.0 ),
        m_energy_in_COMFactor( 0.0 ),
        m_massFactor( 0.0 ),
        m_Q( 0.0 ),
        m_dist( nullptr ) {

    m_type = ProbabilityBase2dType::NBodyPhaseSpace;
}
/*
============================================================
*/
LUPI_HOST NBodyPhaseSpace2d::NBodyPhaseSpace2d( GIDI::Functions::NBodyPhaseSpace2d const &a_NBodyPhaseSpace2d, SetupInfo *a_setupInfo ) :
        ProbabilityBase2d_d2( a_NBodyPhaseSpace2d ),
        m_numberOfProducts( a_NBodyPhaseSpace2d.numberOfProducts( ) ),
        m_mass( PoPI_AMU2MeV_c2 * a_NBodyPhaseSpace2d.mass( ).value( ) ),
        m_energy_in_COMFactor( a_setupInfo->m_protare.targetMass( ) / ( a_setupInfo->m_protare.projectileMass( ) + a_setupInfo->m_protare.targetMass( ) ) ),
        m_massFactor( 1 - a_setupInfo->m_product1Mass / m_mass ),
        m_Q( a_setupInfo->m_Q ),
        m_dist( nullptr ) {

    m_type = ProbabilityBase2dType::NBodyPhaseSpace;
    double xs[2] = { 0.0, 1.0 };
    ptwXYPoints *pdf = nullptr;

    pdf = ptwXY_createFromFunction( nullptr, 2, xs, MCGIDI_NBodyPhaseSpacePDF_callback, (void *) &m_numberOfProducts, 1e-3, 0, 16 );
    if( pdf == nullptr ) throw std::runtime_error( "NBodyPhaseSpace2d::NBodyPhaseSpace2d: ptwXY_createFromFunction returned nullptr" );

    m_dist = ptwXY_To_Xs_pdf_cdf1d( pdf );
    ptwXY_free( pdf );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE NBodyPhaseSpace2d::~NBodyPhaseSpace2d( ) {

    delete m_dist;
}
/*
============================================================
*/
LUPI_HOST static nfu_status MCGIDI_NBodyPhaseSpacePDF_callback( LUPI_maybeUnused statusMessageReporting *smr, double X, double *Y, void *argList ) {

    int numberOfProducts = *((int *) argList);
    double exponent = 0.5 * ( 3 * numberOfProducts - 8 );

    *Y = sqrt( X ) * pow( 1.0 - X, exponent );
    return( nfu_Okay );
}
/*
============================================================
*/
LUPI_HOST_DEVICE double NBodyPhaseSpace2d::evaluate( double a_x2, double a_x1 ) const {

    double EMax = ( m_energy_in_COMFactor * a_x2 + m_Q ) * m_massFactor;
    double x1 = a_x1 / EMax;

    return( m_dist->evaluate( x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void NBodyPhaseSpace2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_INT( m_numberOfProducts, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_mass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_energy_in_COMFactor, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_massFactor, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_Q, a_buffer, a_mode );

    m_dist = serializeProbability1d( a_buffer, a_mode, m_dist );
}

/*
============================================================
====================== Evaporation2d =======================
============================================================
*/
LUPI_HOST_DEVICE Evaporation2d::Evaporation2d( ) :
        m_U( 0.0 ),
        m_theta( nullptr ) {

    m_type = ProbabilityBase2dType::evaporation;
}
/*
============================================================
*/
LUPI_HOST Evaporation2d::Evaporation2d( GIDI::Functions::Evaporation2d const &a_evaporation2d ) :
        ProbabilityBase2d_d2( a_evaporation2d ),
        m_U( a_evaporation2d.U( ) ),
        m_theta( Functions::parseFunction1d_d1( a_evaporation2d.theta( ) ) ) {

    m_type = ProbabilityBase2dType::evaporation;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Evaporation2d::~Evaporation2d( ) {

    delete m_theta;
}
/*
============================================================
*/
LUPI_HOST_DEVICE double Evaporation2d::evaluate( double a_x2, double a_x1 ) const {

    double theta = m_theta->evaluate( a_x2 );
    double E_U_theta = ( a_x2 - m_U ) / theta;
    double Ep_theta = a_x1 / theta;

    if( E_U_theta < 0 ) return( 0.0 );
    return( Ep_theta * exp( -Ep_theta ) / ( theta * ( 1.0 - exp( -E_U_theta ) * ( 1.0 + E_U_theta ) ) ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Evaporation2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_U, a_buffer, a_mode );

    m_theta = serializeFunction1d_d1( a_buffer, a_mode, m_theta );
}

/*
============================================================
=================== GeneralEvaporation2d ===================
============================================================
*/
LUPI_HOST_DEVICE GeneralEvaporation2d::GeneralEvaporation2d( ) :
        m_theta( nullptr ),
        m_g( nullptr ) {

    m_type = ProbabilityBase2dType::generalEvaporation;
}
/*
============================================================
*/
LUPI_HOST GeneralEvaporation2d::GeneralEvaporation2d( GIDI::Functions::GeneralEvaporation2d const &a_generalEvaporation2d ) :
        ProbabilityBase2d_d2( a_generalEvaporation2d ),
        m_theta( Functions::parseFunction1d_d1( a_generalEvaporation2d.theta( ) ) ),
        m_g( parseProbability1d( a_generalEvaporation2d.g( ) ) ) {

    m_type = ProbabilityBase2dType::generalEvaporation;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE GeneralEvaporation2d::~GeneralEvaporation2d( ) {

    delete m_theta;
    delete m_g;
}
/*
============================================================
*/
LUPI_HOST_DEVICE double GeneralEvaporation2d::evaluate( double a_x2, double a_x1 ) const {

    return( m_g->evaluate( a_x1 / m_theta->evaluate( a_x2 ) ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void GeneralEvaporation2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    m_theta = serializeFunction1d_d1( a_buffer, a_mode, m_theta );
    m_g = serializeProbability1d( a_buffer, a_mode, m_g );
}

/*
============================================================
================ SimpleMaxwellianFission2d =================
============================================================
*/
LUPI_HOST_DEVICE SimpleMaxwellianFission2d::SimpleMaxwellianFission2d( ) :
        m_U( 0.0 ),
        m_theta( nullptr ) {

    m_type = ProbabilityBase2dType::simpleMaxwellianFission;
}
/*
============================================================
*/
LUPI_HOST SimpleMaxwellianFission2d::SimpleMaxwellianFission2d( GIDI::Functions::SimpleMaxwellianFission2d const &a_simpleMaxwellianFission2d ) :
        ProbabilityBase2d_d2( a_simpleMaxwellianFission2d ),
        m_U( a_simpleMaxwellianFission2d.U( ) ),
        m_theta( Functions::parseFunction1d_d1( a_simpleMaxwellianFission2d.theta( ) ) ) {

    m_type = ProbabilityBase2dType::simpleMaxwellianFission;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE SimpleMaxwellianFission2d::~SimpleMaxwellianFission2d( ) {

    delete m_theta;
}
/*
============================================================
*/
LUPI_HOST_DEVICE double SimpleMaxwellianFission2d::evaluate( double a_x2, double a_x1 ) const {

    double theta = m_theta->evaluate( a_x2 );
    double E_U_theta = ( a_x2 - m_U ) / theta;

    if( E_U_theta < 0 ) return( 0.0 );

    double Ep_theta = a_x1 / theta;
    double sqrt_E_U_theta = sqrt( E_U_theta );

    return( sqrt( Ep_theta ) * exp( -Ep_theta ) / ( theta * ( erf( sqrt_E_U_theta ) / M_2_SQRTPI - sqrt_E_U_theta * exp( -E_U_theta ) ) ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void SimpleMaxwellianFission2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_U, a_buffer, a_mode );
    m_theta = serializeFunction1d_d1( a_buffer, a_mode, m_theta );
}

/*
============================================================
========================= Watt2d ===========================
============================================================
*/
LUPI_HOST_DEVICE Watt2d::Watt2d( ) :
        m_a( nullptr ),
        m_b( nullptr ) {

    m_type = ProbabilityBase2dType::Watt;
}
/*
============================================================
*/
LUPI_HOST Watt2d::Watt2d( GIDI::Functions::Watt2d const &a_Watt2d ) :
        ProbabilityBase2d_d2( a_Watt2d ),
        m_U( a_Watt2d.U( ) ),
        m_a( Functions::parseFunction1d_d1( a_Watt2d.a( ) ) ),
        m_b( Functions::parseFunction1d_d1( a_Watt2d.b( ) ) ) {

    m_type = ProbabilityBase2dType::Watt;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Watt2d::~Watt2d( ) {

    delete m_a;
    delete m_b;
}
/*
============================================================
*/
LUPI_HOST_DEVICE double Watt2d::evaluate( double a_x2, double a_x1 ) const {

    double Watt_a = m_a->evaluate( a_x2 );
    double E_U_a = ( a_x2 - m_U ) / Watt_a;

    if( E_U_a < 0 ) return( 0.0 );

    double Watt_b = m_b->evaluate( a_x2 );
    double sqrt_ab_4 = 0.5 * sqrt( Watt_a * Watt_b );
    double sqrt_E_U_a = sqrt( E_U_a );

    double I = 0.5 * sqrt_ab_4 * Watt_a * sqrt( M_PI ) * ( erf( sqrt_E_U_a - sqrt_ab_4 ) + erf( sqrt_E_U_a + sqrt_ab_4 ) )
                - Watt_a * exp( -E_U_a ) * sinh( 2.0 * sqrt_E_U_a * sqrt_ab_4 );
    return( exp( -a_x1 / Watt_a ) * sinh( sqrt( Watt_b * a_x1 ) ) / I );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Watt2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_U, a_buffer, a_mode );
    m_a = serializeFunction1d_d1( a_buffer, a_mode, m_a );
    m_b = serializeFunction1d_d1( a_buffer, a_mode, m_b );
}

/*
============================================================
=================== WeightedFunctionals2d ==================
============================================================
*/
LUPI_HOST_DEVICE WeightedFunctionals2d::WeightedFunctionals2d( ) :
        m_weight( ),
        m_energy( ) {

    m_type = ProbabilityBase2dType::weightedFunctionals;
}
/*
============================================================
*/
LUPI_HOST WeightedFunctionals2d::WeightedFunctionals2d( GIDI::Functions::WeightedFunctionals2d const &a_weightedFunctionals2d ) :
        ProbabilityBase2d( a_weightedFunctionals2d ) {

    m_type = ProbabilityBase2dType::weightedFunctionals;

    Vector<GIDI::Functions::Weighted_function2d *> const &weighted_function2d = a_weightedFunctionals2d.weighted_function2d( );
    m_weight.resize( weighted_function2d.size( ) );
    m_energy.resize( weighted_function2d.size( ) );
    for( std::size_t i1 = 0; i1 < weighted_function2d.size( ); ++i1 ) {
        m_weight[i1] = Functions::parseFunction1d_d1( weighted_function2d[i1]->weight( ) );
        m_energy[i1] = parseProbability2d_d1( weighted_function2d[i1]->energy( ), nullptr );
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE WeightedFunctionals2d::~WeightedFunctionals2d( ) {

    for( std::size_t i1 = 0; i1 < m_weight.size( ); ++i1 ) delete m_weight[i1];
    for( std::size_t i1 = 0; i1 < m_energy.size( ); ++i1 ) delete m_energy[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE double WeightedFunctionals2d::evaluate( double a_x2, double a_x1 ) const {

    std::size_t n1 = m_weight.size( );
    double evaluatedValue = 0;

    for( std::size_t i1 = 0; i1 < n1; ++i1 ) {
        evaluatedValue += m_weight[i1]->evaluate( a_x2 ) * m_energy[i1]->evaluate( a_x2, a_x1 );
    }
    return( evaluatedValue  );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void WeightedFunctionals2d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase2d::serialize( a_buffer, a_mode );

    std::size_t vectorSize = m_weight.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_weight.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_weight.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_weight[vectorIndex] = serializeFunction1d_d1( a_buffer, a_mode, m_weight[vectorIndex] );
    }

    vectorSize = m_energy.size( );
    vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_energy.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_energy.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_energy[vectorIndex] = serializeProbability2d_d1( a_buffer, a_mode, m_energy[vectorIndex] );
    }
}


/*
============================================================
===================== ProbabilityBase3d ====================
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase3d::ProbabilityBase3d( ) :
        m_type( ProbabilityBase3dType::none ) {

}
/*
============================================================
*/
LUPI_HOST ProbabilityBase3d::ProbabilityBase3d( GIDI::Functions::FunctionForm const &a_probability, Vector<double> const &a_Xs ) :
        ProbabilityBase( a_probability, a_Xs ),
        m_type( ProbabilityBase3dType::none ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE ProbabilityBase3d::~ProbabilityBase3d( ) {

}

/* *********************************************************************************************************//**
 * Returns a String representation of the **ProbabilityBase3d** type of *this*.
 *
 * @return                          A String instance.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE String ProbabilityBase3d::typeString( ) const {

    String typeStr( "ProbabilityBase3d::" );

    switch( m_type ) {
    case ProbabilityBase3dType::none :
        typeStr += "none";
        break;
    case ProbabilityBase3dType::XYs :
        typeStr += "XYs";
        break;
    }
 
    return( typeStr );
}

/* *********************************************************************************************************//**
 * Returns the value of the function f(x3, x2, x1) at x3 = *a_x3*, x2 = *a_x2* and x1 = *a_x1*.
 *
 * @param a_x3                  [in]    The x3 value.
 * @param a_x2                  [in]    The x2 value.
 * @param a_x1                  [in]    The x1 value.
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double ProbabilityBase3d::evaluate( double a_x3, double a_x2, double a_x1 ) const {

    return( static_cast<XYs3d const *>( this )->evaluate( a_x3, a_x2, a_x1 ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void ProbabilityBase3d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase::serialize( a_buffer, a_mode ); 

    int type = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_type ) {
        case ProbabilityBase3dType::none :
            break;
        case ProbabilityBase3dType::XYs :
            type = 1;
            break;
        }
    }
    DATA_MEMBER_INT( type, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            m_type = ProbabilityBase3dType::none;
            break;
        case 1 :
            m_type = ProbabilityBase3dType::XYs;
            break;
        }
    }
}

/*
============================================================
========================== XYs3d ===========================
============================================================
*/
LUPI_HOST_DEVICE XYs3d::XYs3d( ) :
        m_probabilities( ) {

    m_type = ProbabilityBase3dType::XYs;
}
/*
============================================================
*/
LUPI_HOST XYs3d::XYs3d( GIDI::Functions::XYs3d const &a_XYs3d ) :
        ProbabilityBase3d( a_XYs3d, a_XYs3d.Xs( ) ) {

    m_type = ProbabilityBase3dType::XYs;

    Vector<GIDI::Functions::Function2dForm *> const &functions2d = a_XYs3d.function2ds( );
    m_probabilities.resize( functions2d.size( ) );
    for( std::size_t i1 = 0; i1 < functions2d.size( ); ++i1 ) m_probabilities[i1] = parseProbability2d_d1( functions2d[i1], nullptr );
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE XYs3d::~XYs3d( ) {

    for( std::size_t i1 = 0; i1 < m_probabilities.size( ); ++i1 ) delete m_probabilities[i1];
}
/*
============================================================
*/
LUPI_HOST_DEVICE double XYs3d::evaluate( double a_x3, double a_x2, double a_x1 ) const {

    int lower = binarySearchVector( a_x3, m_Xs );
    double evaluatedValue;

    if( lower == -2 ) {                         // a_x3 < first value of Xs.
        evaluatedValue = m_probabilities[0]->evaluate( a_x2, a_x1 ); }
    else if( lower == -1 ) {                    // a_x3 > last value of Xs.
        evaluatedValue = m_probabilities.back( )->evaluate( a_x2, a_x1 ); }
    else {
        double value1 = m_probabilities[lower]->evaluate( a_x2, a_x1 );

        if( interpolation( ) == Interpolation::FLAT ) {
            evaluatedValue = value1; }
        else {
            double value2 = m_probabilities[lower+1]->evaluate( a_x2, a_x1 );

            if( interpolation( ) == Interpolation::LINLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x3 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                evaluatedValue = fraction * value1 + ( 1 - fraction ) * value2 ; }
            else if( interpolation( ) == Interpolation::LOGLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x3 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                evaluatedValue = value2 * pow( value2 / value1, fraction ); }
            else if( interpolation( ) == Interpolation::LINLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x3 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                evaluatedValue = fraction * value1 + ( 1 - fraction ) * value2; }
            else if( interpolation( ) == Interpolation::LOGLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x3 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                evaluatedValue = value2 * pow( value2 / value1, fraction ); }
            else {                                                              // This should never happen.
                LUPI_THROW( "XYs3d::evaluate: unsupported interpolation." );
            }
        }
    }

    return( evaluatedValue );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void XYs3d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    ProbabilityBase3d::serialize( a_buffer, a_mode );

    std::size_t vectorSize = m_probabilities.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_probabilities.resize( vectorSize, &a_buffer.m_placement );
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) a_buffer.m_placement += m_probabilities.internalSize();
    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        m_probabilities[vectorIndex] = serializeProbability2d_d1( a_buffer, a_mode, m_probabilities[vectorIndex] );
    }
}

/*
============================================================
========================== others ==========================
============================================================
*/
/*
============================================================
*/
LUPI_HOST ProbabilityBase1d *parseProbability1d( GIDI::Functions::Function1dForm const *a_form1d ) {

    GIDI::FormType type = a_form1d->type( );

    switch( type ) {
    case GIDI::FormType::xs_pdf_cdf1d :
        return( new Xs_pdf_cdf1d( *static_cast<GIDI::Functions::Xs_pdf_cdf1d const *>( a_form1d ) ) );
    default :
        throw std::runtime_error( "Probabilities::parseProbability1d: Unsupported Function1d with moniker " + a_form1d->moniker( ) 
                + " at " + a_form1d->toXLink( ) );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST ProbabilityBase2d *parseProbability2d( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo ) {

    GIDI::FormType type = form2d->type( );

    switch( type ) {
    case GIDI::FormType::XYs2d :
        return( new XYs2d( *static_cast<GIDI::Functions::XYs2d const *>( form2d ) ) );
    case GIDI::FormType::regions2d :
        return( new Regions2d( *static_cast<GIDI::Functions::Regions2d const *>( form2d ) ) );
    case GIDI::FormType::isotropic2d :
        return( new Isotropic2d( *static_cast<GIDI::Functions::Isotropic2d const *>( form2d ) ) );
    case GIDI::FormType::discreteGamma2d :
        return( new DiscreteGamma2d( *static_cast<GIDI::Functions::DiscreteGamma2d const *>( form2d ) ) );
    case GIDI::FormType::primaryGamma2d :
        return( new PrimaryGamma2d( *static_cast<GIDI::Functions::PrimaryGamma2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::recoil2d :
        return( new Recoil2d( *static_cast<GIDI::Functions::Recoil2d const *>( form2d ) ) );
    case GIDI::FormType::NBodyPhaseSpace2d :
        return( new NBodyPhaseSpace2d( *static_cast<GIDI::Functions::NBodyPhaseSpace2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::evaporation2d :
        return( new Evaporation2d( *static_cast<GIDI::Functions::Evaporation2d const *>( form2d ) ) );
    case GIDI::FormType::generalEvaporation2d :
        return( new GeneralEvaporation2d( *static_cast<GIDI::Functions::GeneralEvaporation2d const *>( form2d ) ) );
    case GIDI::FormType::simpleMaxwellianFission2d :
        return( new SimpleMaxwellianFission2d( *static_cast<GIDI::Functions::SimpleMaxwellianFission2d const *>( form2d ) ) );
    case GIDI::FormType::Watt2d :
        return( new Watt2d( *static_cast<GIDI::Functions::Watt2d const *>( form2d ) ) );
    case GIDI::FormType::weightedFunctionals2d :
        return( new WeightedFunctionals2d( *static_cast<GIDI::Functions::WeightedFunctionals2d const *>( form2d ) ) );
    default :
        throw std::runtime_error( "Probabilities::parseProbability2d: Unsupported Function2d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST ProbabilityBase2d_d1 *parseProbability2d_d1( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo ) {

    GIDI::FormType type = form2d->type( );

    switch( type ) {
    case GIDI::FormType::XYs2d :
        return( new XYs2d( *static_cast<GIDI::Functions::XYs2d const *>( form2d ) ) );
    case GIDI::FormType::regions2d :
        return( new Regions2d( *static_cast<GIDI::Functions::Regions2d const *>( form2d ) ) );
    case GIDI::FormType::isotropic2d :
        return( new Isotropic2d( *static_cast<GIDI::Functions::Isotropic2d const *>( form2d ) ) );
    case GIDI::FormType::discreteGamma2d :
        return( new DiscreteGamma2d( *static_cast<GIDI::Functions::DiscreteGamma2d const *>( form2d ) ) );
    case GIDI::FormType::primaryGamma2d :
        return( new PrimaryGamma2d( *static_cast<GIDI::Functions::PrimaryGamma2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::recoil2d :
        return( new Recoil2d( *static_cast<GIDI::Functions::Recoil2d const *>( form2d ) ) );
    case GIDI::FormType::NBodyPhaseSpace2d :
        return( new NBodyPhaseSpace2d( *static_cast<GIDI::Functions::NBodyPhaseSpace2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::evaporation2d :
        return( new Evaporation2d( *static_cast<GIDI::Functions::Evaporation2d const *>( form2d ) ) );
    case GIDI::FormType::generalEvaporation2d :
        return( new GeneralEvaporation2d( *static_cast<GIDI::Functions::GeneralEvaporation2d const *>( form2d ) ) );
    case GIDI::FormType::simpleMaxwellianFission2d :
        return( new SimpleMaxwellianFission2d( *static_cast<GIDI::Functions::SimpleMaxwellianFission2d const *>( form2d ) ) );
    case GIDI::FormType::Watt2d :
        return( new Watt2d( *static_cast<GIDI::Functions::Watt2d const *>( form2d ) ) );
    default :
        throw std::runtime_error( "Probabilities::parseProbability2d: Unsupported Function2d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST ProbabilityBase2d_d2 *parseProbability2d_d2( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo ) {

    GIDI::FormType type = form2d->type( );

    switch( type ) {
    case GIDI::FormType::XYs2d :
        return( new XYs2d( *static_cast<GIDI::Functions::XYs2d const *>( form2d ) ) );
    case GIDI::FormType::isotropic2d :
        return( new Isotropic2d( *static_cast<GIDI::Functions::Isotropic2d const *>( form2d ) ) );
    case GIDI::FormType::discreteGamma2d :
        return( new DiscreteGamma2d( *static_cast<GIDI::Functions::DiscreteGamma2d const *>( form2d ) ) );
    case GIDI::FormType::primaryGamma2d :
        return( new PrimaryGamma2d( *static_cast<GIDI::Functions::PrimaryGamma2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::NBodyPhaseSpace2d :
        return( new NBodyPhaseSpace2d( *static_cast<GIDI::Functions::NBodyPhaseSpace2d const *>( form2d ), a_setupInfo ) );
    case GIDI::FormType::evaporation2d :
        return( new Evaporation2d( *static_cast<GIDI::Functions::Evaporation2d const *>( form2d ) ) );
    case GIDI::FormType::generalEvaporation2d :
        return( new GeneralEvaporation2d( *static_cast<GIDI::Functions::GeneralEvaporation2d const *>( form2d ) ) );
    case GIDI::FormType::simpleMaxwellianFission2d :
        return( new SimpleMaxwellianFission2d( *static_cast<GIDI::Functions::SimpleMaxwellianFission2d const *>( form2d ) ) );
    case GIDI::FormType::Watt2d :
        return( new Watt2d( *static_cast<GIDI::Functions::Watt2d const *>( form2d ) ) );
    default :
        throw std::runtime_error( "Probabilities::parseProbability2d: Unsupported Function2d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST ProbabilityBase3d *parseProbability3d( GIDI::Functions::Function3dForm const *form3d ) {

    GIDI::FormType type = form3d->type( );

    switch( type ) {
    case GIDI::FormType::XYs3d :
        return( new XYs3d( *static_cast<GIDI::Functions::XYs3d const *>( form3d ) ) );
    default :
        throw std::runtime_error( "Probabilities::parseProbability3d: Unsupported Function3d" );
    }

    return( nullptr );
}

/*
============================================================
*/
LUPI_HOST static ProbabilityBase1d *ptwXY_To_Xs_pdf_cdf1d( ptwXYPoints *pdfXY ) {

    ptwXPoints *cdfX = nullptr;
    ptwXYPoint *point;
    std::size_t n1 = (std::size_t) ptwXY_length( nullptr, pdfXY );
    std::vector<double> Xs( n1 ), pdf( n1 ), cdf( n1 );

    if( ( cdfX = ptwXY_runningIntegral( nullptr, pdfXY ) ) == nullptr ) throw std::runtime_error( "ptwXY_To_Xs_pdf_cdf1d: ptwXY_runningIntegral returned error." );
    double norm = ptwX_getPointAtIndex_Unsafely( cdfX, n1 - 1 );
    if( norm <= 0 ) throw std::runtime_error( "ptwXY_To_Xs_pdf_cdf1d: norm <= 0." );

    norm = 1. / norm;
    for( std::size_t i1 = 0; i1 < n1; ++i1 ) {
        point = ptwXY_getPointAtIndex_Unsafely( pdfXY, i1 );
        Xs[i1] = point->x;
        pdf[i1] = norm * point->y;
        cdf[i1] = norm * ptwX_getPointAtIndex_Unsafely( cdfX, i1 );
    }
    cdf[n1-1] = 1.;

    ptwX_free( cdfX );

    GIDI::Axes axes;
    GIDI::Functions::Xs_pdf_cdf1d gidi_xs_pdf_cdf1d( axes, ptwXY_interpolationLinLin, Xs, pdf, cdf );

    Xs_pdf_cdf1d *xs_pdf_cdf1d = new Xs_pdf_cdf1d( gidi_xs_pdf_cdf1d );
    return( xs_pdf_cdf1d );
}

}           // End of namespace Probabilities.

/*
============================================================
========================== others ==========================
============================================================
*/
LUPI_HOST_DEVICE Interpolation GIDI2MCGIDI_interpolation( ptwXY_interpolation a_interpolation ) {

    if( a_interpolation == ptwXY_interpolationLinLin ) return( Interpolation::LINLIN );
    if( a_interpolation == ptwXY_interpolationLogLin ) return( Interpolation::LOGLIN );
    if( a_interpolation == ptwXY_interpolationLinLog ) return( Interpolation::LINLOG );
    if( a_interpolation == ptwXY_interpolationLogLog ) return( Interpolation::LOGLOG );
    if( a_interpolation == ptwXY_interpolationFlat ) return( Interpolation::FLAT );
    return( Interpolation::OTHER );
}
/*
============================================================
*/
LUPI_HOST_DEVICE Function1dType Function1dClass( Functions::Function1d *a_function ) {

    if( a_function == nullptr ) return( Function1dType::none );
    return( a_function->type( ) );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Functions::Function1d *serializeFunction1d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Functions::Function1d *a_function1d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        Function1dType fType = Function1dClass( a_function1d );

        switch( fType ) {
        case Function1dType::none :
            break;
        case Function1dType::constant :
            type = 1;
            break;
        case Function1dType::XYs :
            type = 2;
            break;
        case Function1dType::polyomial :
            type = 3;
            break;
        case Function1dType::gridded :
            type = 4;
            break;
        case Function1dType::regions :
            type = 5;
            break;
        case Function1dType::branching :
            type = 6;
            break;
        case Function1dType::TerrellFissionNeutronMultiplicityModel :
            type = 7;
            break;
        default:
            String message( "serializeFunction1d: Unsupported Function1d: " + a_function1d->typeString( ) );
            LUPI_THROW( message.c_str( ) );
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_function1d = nullptr;
        switch( type ) {
        case 0 :
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Constant1d;
                a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) ); }
            else {
                a_function1d = new Functions::Constant1d;
            }
            break;
        case 2 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::XYs1d;
                a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) ); }
            else {
                a_function1d = new Functions::XYs1d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Polynomial1d;
                a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) ); }
            else {
                a_function1d = new Functions::Polynomial1d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Gridded1d;
                a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) ); }
            else {
                a_function1d = new Functions::Gridded1d;
            }
            break;
        case 5 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Regions1d;
                a_buffer.incrementPlacement( sizeof( Functions::Regions1d ) ); }
            else {
                a_function1d = new Functions::Regions1d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Branching1d;
                a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) ); }
            else {
                a_function1d = new Functions::Branching1d;
            }
            break;
        case 7 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::TerrellFissionNeutronMultiplicityModel;
                a_buffer.incrementPlacement( sizeof( Functions::TerrellFissionNeutronMultiplicityModel ) ); }
            else {
                a_function1d = new Functions::TerrellFissionNeutronMultiplicityModel;
            }
            break;
        default:                                // This should never happen as Unpack should be called after Pack which checks type.
            LUPI_THROW( "serializeFunction1d: Unsupported Function1d:" );
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) );
            break;
        case 2 :
            a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) );
            break;
        case 5 :
            a_buffer.incrementPlacement( sizeof( Functions::Regions1d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) );
            break;
        case 7 :
            a_buffer.incrementPlacement( sizeof( Functions::TerrellFissionNeutronMultiplicityModel ) );
            break;
        default:
            LUPI_THROW( "serializeFunction1d: Unsupported Function1d:" );
        }
    }

    if( a_function1d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Functions::Constant1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 2 :
            static_cast<Functions::XYs1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Functions::Polynomial1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Functions::Gridded1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 5 :
            static_cast<Functions::Regions1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Functions::Branching1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 7 :
            static_cast<Functions::TerrellFissionNeutronMultiplicityModel *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        default:
            LUPI_THROW( "serializeFunction1d: Unsupported Function1d:" );
        }
    }

    return( a_function1d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Functions::Function1d_d1 *serializeFunction1d_d1( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Functions::Function1d_d1 *a_function1d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        Function1dType fType = Function1dClass( a_function1d );

        switch( fType ) {
        case Function1dType::none :
            break;
        case Function1dType::constant :
            type = 1;
            break;
        case Function1dType::XYs :
            type = 2;
            break;
        case Function1dType::polyomial :
            type = 3;
            break;
        case Function1dType::gridded :
            type = 4;
            break;
        case Function1dType::regions :
            type = 5;
            break;
        case Function1dType::branching :
            type = 6;
            break;
        default:
            String message( "serializeFunction1d_d1: Unsupported Function1d: " + a_function1d->typeString( ) );
            LUPI_THROW( message.c_str( ) );
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_function1d = nullptr;
        switch( type ) {
        case 0 :
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Constant1d;
                a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) ); }
            else {
                a_function1d = new Functions::Constant1d;
            }
            break;
        case 2 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::XYs1d;
                a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) ); }
            else {
                a_function1d = new Functions::XYs1d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Polynomial1d;
                a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) ); }
            else {
                a_function1d = new Functions::Polynomial1d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Gridded1d;
                a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) ); }
            else {
                a_function1d = new Functions::Gridded1d;
            }
            break;
        case 5 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Regions1d;
                a_buffer.incrementPlacement( sizeof( Functions::Regions1d ) ); }
            else {
                a_function1d = new Functions::Regions1d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Branching1d;
                a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) ); }
            else {
                a_function1d = new Functions::Branching1d;
            }
            break;
        default:                                // This should never happen as Unpack should be called after Pack which checks type.
            LUPI_THROW( "serializeFunction1d_d1: Unsupported Function1d:" );
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) );
            break;
        case 2 :
            a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) );
            break;
        case 5 :
            a_buffer.incrementPlacement( sizeof( Functions::Regions1d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) );
            break;
        default:
            LUPI_THROW( "serializeFunction1d_d1: Unsupported Function1d:" );
        }
    }

    if( a_function1d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Functions::Constant1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 2 :
            static_cast<Functions::XYs1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Functions::Polynomial1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Functions::Gridded1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 5 :
            static_cast<Functions::Regions1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Functions::Branching1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        default:
            LUPI_THROW( "serializeFunction1d_d1: Unsupported Function1d:" );
        }
    }

    return( a_function1d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Functions::Function1d_d2 *serializeFunction1d_d2( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Functions::Function1d_d2 *a_function1d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        Function1dType fType = Function1dClass( a_function1d );

        switch( fType ) {
        case Function1dType::none :
            break;
        case Function1dType::constant :
            type = 1;
            break;
        case Function1dType::XYs :
            type = 2;
            break;
        case Function1dType::polyomial :
            type = 3;
            break;
        case Function1dType::gridded :
            type = 4;
            break;
            break;
        case Function1dType::branching :
            type = 6;
            break;
        default:
            String message( "serializeFunction1d_d2: Unsupported Function1d: " + a_function1d->typeString( ) );
            LUPI_THROW( message.c_str( ) );
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_function1d = nullptr;
        switch( type ) {
        case 0 :
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Constant1d;
                a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) ); }
            else {
                a_function1d = new Functions::Constant1d;
            }
            break;
        case 2 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::XYs1d;
                a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) ); }
            else {
                a_function1d = new Functions::XYs1d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Polynomial1d;
                a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) ); }
            else {
                a_function1d = new Functions::Polynomial1d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Gridded1d;
                a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) ); }
            else {
                a_function1d = new Functions::Gridded1d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_function1d = new(a_buffer.m_placement) Functions::Branching1d;
                a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) ); }
            else {
                a_function1d = new Functions::Branching1d;
            }
            break;
        default:                                // This should never happen as Unpack should be called after Pack which checks type.
            LUPI_THROW( "serializeFunction1d_d2: Unsupported Function1d:" );
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Functions::Constant1d ) );
            break;
        case 2 :
            a_buffer.incrementPlacement( sizeof( Functions::XYs1d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Functions::Polynomial1d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Functions::Gridded1d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Functions::Branching1d ) );
            break;
        default:
            LUPI_THROW( "serializeFunction1d_d2: Unsupported Function1d:" );
        }
    }

    if( a_function1d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Functions::Constant1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 2 :
            static_cast<Functions::XYs1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Functions::Polynomial1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Functions::Gridded1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Functions::Branching1d *>( a_function1d )->serialize( a_buffer, a_mode );
            break;
        default:
            LUPI_THROW( "serializeFunction1d_d2: Unsupported Function1d:" );
        }
    }

    return( a_function1d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Function2dType Function2dClass( Functions::Function2d *a_function ) {

    if( a_function == nullptr ) return( Function2dType::none );
    return( a_function->type( ) );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Functions::Function2d *serializeFunction2d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Functions::Function2d *a_function2d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        Function2dType fType = Function2dClass( a_function2d );

        switch( fType ) {
        case Function2dType::none :
            break;
        case Function2dType::XYs :
            type = 1;
            break;
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_function2d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_function2d = new(a_buffer.m_placement) Functions::XYs2d;
                a_buffer.incrementPlacement( sizeof( Functions::XYs2d ) ); }
            else {
                a_function2d = new Functions::XYs2d;
            }
            break;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Functions::XYs2d ) );
            break;
        }
    }

    if( a_function2d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Functions::XYs2d *>( a_function2d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_function2d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase1dType ProbabilityBase1dClass( Probabilities::ProbabilityBase1d *a_function ) {

    if( a_function == nullptr ) return( ProbabilityBase1dType::none );
    return( a_function->type( ) );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Probabilities::ProbabilityBase1d *serializeProbability1d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Probabilities::ProbabilityBase1d *a_probability1d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        ProbabilityBase1dType pType = ProbabilityBase1dClass( a_probability1d );

        switch( pType ) {
        case ProbabilityBase1dType::none :
            break;
        case ProbabilityBase1dType::xs_pdf_cdf :
            type = 1;
            break;
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_probability1d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability1d = new(a_buffer.m_placement) Probabilities::Xs_pdf_cdf1d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Xs_pdf_cdf1d ) ); }
            else {
                a_probability1d = new Probabilities::Xs_pdf_cdf1d;
            }
            break;
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Xs_pdf_cdf1d ) );
            break;
        }
    }

    if( a_probability1d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Probabilities::Xs_pdf_cdf1d *>( a_probability1d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_probability1d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase2dType ProbabilityBase2dClass( Probabilities::ProbabilityBase2d *a_function ) {

    if( a_function == nullptr ) return( ProbabilityBase2dType::none );
    return( a_function->type( ) );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d *serializeProbability2d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Probabilities::ProbabilityBase2d *a_probability2d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        ProbabilityBase2dType pType = ProbabilityBase2dClass( a_probability2d );

        switch( pType ) {
        case ProbabilityBase2dType::none :
            break;
        case ProbabilityBase2dType::XYs :
            type = 1;
            break;
        case ProbabilityBase2dType::regions :
            type = 2;
            break;
        case ProbabilityBase2dType::isotropic :
            type = 3;
            break;
        case ProbabilityBase2dType::discreteGamma :
            type = 4;
            break;
        case ProbabilityBase2dType::primaryGamma :
            type = 5;
            break;
        case ProbabilityBase2dType::recoil :
            type = 6;
            break;
        case ProbabilityBase2dType::NBodyPhaseSpace :
            type = 7;
            break;
        case ProbabilityBase2dType::evaporation :
            type = 8;
            break;
        case ProbabilityBase2dType::generalEvaporation :
            type = 9;
            break;
        case ProbabilityBase2dType::simpleMaxwellianFission :
            type = 10;
            break;
        case ProbabilityBase2dType::Watt :
            type = 11;
            break;
        case ProbabilityBase2dType::weightedFunctionals :
            type = 12;
            break;
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_probability2d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::XYs2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) ); }
            else {
                a_probability2d = new Probabilities::XYs2d;
            }
            break;
        case 2 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Regions2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Regions2d ) ); }
            else {
                a_probability2d = new Probabilities::Regions2d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Isotropic2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) ); }
            else {
                a_probability2d = new Probabilities::Isotropic2d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::DiscreteGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::DiscreteGamma2d;
            }
            break;
        case 5 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::PrimaryGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::PrimaryGamma2d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Recoil2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) ); }
            else {
                a_probability2d = new Probabilities::Recoil2d;
            }
            break;
        case 7 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::NBodyPhaseSpace2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) ); }
            else {
                a_probability2d = new Probabilities::NBodyPhaseSpace2d;
            }
            break;
        case 8 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Evaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::Evaporation2d;
            }
            break;
        case 9 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::GeneralEvaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::GeneralEvaporation2d;
            }
            break;
        case 10 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::SimpleMaxwellianFission2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) ); }
            else {
                a_probability2d = new Probabilities::SimpleMaxwellianFission2d;
            }
            break;
        case 11 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Watt2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) ); }
            else {
                a_probability2d = new Probabilities::Watt2d;
            }
            break;
        case 12 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::WeightedFunctionals2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::WeightedFunctionals2d ) ); }
            else {
                a_probability2d = new Probabilities::WeightedFunctionals2d;
            }
            break;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) );
            break;
        case 2 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Regions2d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) );
            break;
        case 5 :
            a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) );
            break;
        case 7 :
            a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) );
            break;
        case 8 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) );
            break;
        case 9 :
            a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) );
            break;
        case 10 :
            a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) );
            break;
        case 11 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) );
            break;
        case 12 :
            a_buffer.incrementPlacement( sizeof( Probabilities::WeightedFunctionals2d ) );
            break;
        }
    }

    if( a_probability2d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Probabilities::XYs2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 2 :
            static_cast<Probabilities::Regions2d * >( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Probabilities::Isotropic2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Probabilities::DiscreteGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 5 :
            static_cast<Probabilities::PrimaryGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Probabilities::Recoil2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 7 :
            static_cast<Probabilities::NBodyPhaseSpace2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 8 :
            static_cast<Probabilities::Evaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 9 :
            static_cast<Probabilities::GeneralEvaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 10 :
            static_cast<Probabilities::SimpleMaxwellianFission2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 11 :
            static_cast<Probabilities::Watt2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 12 :
            static_cast<Probabilities::WeightedFunctionals2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_probability2d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *serializeProbability2d_d1( LUPI::DataBuffer &a_buffer, 
                LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase2d_d1 *a_probability2d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        ProbabilityBase2dType pType = ProbabilityBase2dClass( a_probability2d );

        switch( pType ) {
        case ProbabilityBase2dType::none :
            break;
        case ProbabilityBase2dType::XYs :
            type = 1;
            break;
        case ProbabilityBase2dType::regions :
            type = 2;
            break;
        case ProbabilityBase2dType::isotropic :
            type = 3;
            break;
        case ProbabilityBase2dType::discreteGamma :
            type = 4;
            break;
        case ProbabilityBase2dType::primaryGamma :
            type = 5;
            break;
        case ProbabilityBase2dType::recoil :
            type = 6;
            break;
        case ProbabilityBase2dType::NBodyPhaseSpace :
            type = 7;
            break;
        case ProbabilityBase2dType::evaporation :
            type = 8;
            break;
        case ProbabilityBase2dType::generalEvaporation :
            type = 9;
            break;
        case ProbabilityBase2dType::simpleMaxwellianFission :
            type = 10;
            break;
        case ProbabilityBase2dType::Watt :
            type = 11;
            break;
        default:
            LUPI_THROW( "Probabilities::ProbabilityBase2d_d1: Unsupported ProbabilityBase2d_d1." );
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_probability2d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::XYs2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) ); }
            else {
                a_probability2d = new Probabilities::XYs2d;
            }
            break;
        case 2 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Regions2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Regions2d ) ); }
            else {
                a_probability2d = new Probabilities::Regions2d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Isotropic2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) ); }
            else {
                a_probability2d = new Probabilities::Isotropic2d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::DiscreteGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::DiscreteGamma2d;
            }
            break;
        case 5 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::PrimaryGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::PrimaryGamma2d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Recoil2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) ); }
            else {
                a_probability2d = new Probabilities::Recoil2d;
            }
            break;
        case 7 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::NBodyPhaseSpace2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) ); }
            else {
                a_probability2d = new Probabilities::NBodyPhaseSpace2d;
            }
            break;
        case 8 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Evaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::Evaporation2d;
            }
            break;
        case 9 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::GeneralEvaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::GeneralEvaporation2d;
            }
            break;
        case 10 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::SimpleMaxwellianFission2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) ); }
            else {
                a_probability2d = new Probabilities::SimpleMaxwellianFission2d;
            }
            break;
        case 11 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Watt2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) ); }
            else {
                a_probability2d = new Probabilities::Watt2d;
            }
            break;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) );
            break;
        case 2 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Regions2d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) );
            break;
        case 5 :
            a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) );
            break;
        case 7 :
            a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) );
            break;
        case 8 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) );
            break;
        case 9 :
            a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) );
            break;
        case 10 :
            a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) );
            break;
        case 11 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) );
            break;
        }
    }

    if( a_probability2d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Probabilities::XYs2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 2 :
            static_cast<Probabilities::Regions2d * >( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Probabilities::Isotropic2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Probabilities::DiscreteGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 5 :
            static_cast<Probabilities::PrimaryGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Probabilities::Recoil2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 7 :
            static_cast<Probabilities::NBodyPhaseSpace2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 8 :
            static_cast<Probabilities::Evaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 9 :
            static_cast<Probabilities::GeneralEvaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 10 :
            static_cast<Probabilities::SimpleMaxwellianFission2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 11 :
            static_cast<Probabilities::Watt2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_probability2d );
}

/*  
============================================================
*/  
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d2 *serializeProbability2d_d2( LUPI::DataBuffer &a_buffer,
                LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase2d_d2 *a_probability2d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        ProbabilityBase2dType pType = ProbabilityBase2dClass( a_probability2d );

        switch( pType ) {
        case ProbabilityBase2dType::none :
            break;
        case ProbabilityBase2dType::XYs :
            type = 1;
            break;
        case ProbabilityBase2dType::isotropic :
            type = 3;
            break;
        case ProbabilityBase2dType::discreteGamma :
            type = 4;
            break;
        case ProbabilityBase2dType::primaryGamma :
            type = 5;
            break;
        case ProbabilityBase2dType::recoil :
            type = 6;
            break;
        case ProbabilityBase2dType::NBodyPhaseSpace :
            type = 7;
            break;
        case ProbabilityBase2dType::evaporation :
            type = 8;
            break;
        case ProbabilityBase2dType::generalEvaporation :
            type = 9;
            break;
        case ProbabilityBase2dType::simpleMaxwellianFission :
            type = 10;
            break;
        case ProbabilityBase2dType::Watt :
            type = 11;
            break;
        default:
            LUPI_THROW( "Probabilities::ProbabilityBase2d_d1: Unsupported ProbabilityBase2d_d2" );
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_probability2d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::XYs2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) ); }
            else {
                a_probability2d = new Probabilities::XYs2d;
            }
            break;
        case 3 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Isotropic2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) ); }
            else {
                a_probability2d = new Probabilities::Isotropic2d;
            }
            break;
        case 4 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::DiscreteGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::DiscreteGamma2d;
            }
            break;
        case 5 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::PrimaryGamma2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) ); }
            else {
                a_probability2d = new Probabilities::PrimaryGamma2d;
            }
            break;
        case 6 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Recoil2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) ); }
            else {
                a_probability2d = new Probabilities::Recoil2d;
            }
            break;
        case 7 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::NBodyPhaseSpace2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) ); }
            else {
                a_probability2d = new Probabilities::NBodyPhaseSpace2d;
            }
            break;
        case 8 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Evaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::Evaporation2d;
            }
            break;
        case 9 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::GeneralEvaporation2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) ); }
            else {
                a_probability2d = new Probabilities::GeneralEvaporation2d;
            }
            break;
        case 10 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::SimpleMaxwellianFission2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) ); }
            else {
                a_probability2d = new Probabilities::SimpleMaxwellianFission2d;
            }
            break;
        case 11 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability2d = new(a_buffer.m_placement) Probabilities::Watt2d;
                a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) ); }
            else {
                a_probability2d = new Probabilities::Watt2d;
            }
            break;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Probabilities::XYs2d ) );
            break;
        case 3 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Isotropic2d ) );
            break;
        case 4 :
            a_buffer.incrementPlacement( sizeof( Probabilities::DiscreteGamma2d ) );
            break;
        case 5 :
            a_buffer.incrementPlacement( sizeof( Probabilities::PrimaryGamma2d ) );
            break;
        case 6 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Recoil2d ) );
            break;
        case 7 :
            a_buffer.incrementPlacement( sizeof( Probabilities::NBodyPhaseSpace2d ) );
            break;
        case 8 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Evaporation2d ) );
            break;
        case 9 :
            a_buffer.incrementPlacement( sizeof( Probabilities::GeneralEvaporation2d ) );
            break;
        case 10 :
            a_buffer.incrementPlacement( sizeof( Probabilities::SimpleMaxwellianFission2d ) );
            break;
        case 11 :
            a_buffer.incrementPlacement( sizeof( Probabilities::Watt2d ) );
            break;
        }
    }

    if( a_probability2d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Probabilities::XYs2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 3 :
            static_cast<Probabilities::Isotropic2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 4 :
            static_cast<Probabilities::DiscreteGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 5 :
            static_cast<Probabilities::PrimaryGamma2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 6 :
            static_cast<Probabilities::Recoil2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 7 :
            static_cast<Probabilities::NBodyPhaseSpace2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 8 :
            static_cast<Probabilities::Evaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 9 :
            static_cast<Probabilities::GeneralEvaporation2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 10 :
            static_cast<Probabilities::SimpleMaxwellianFission2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        case 11 :
            static_cast<Probabilities::Watt2d *>( a_probability2d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_probability2d );
}

/*
============================================================
*/
LUPI_HOST_DEVICE ProbabilityBase3dType ProbabilityBase3dClass( Probabilities::ProbabilityBase3d *a_function ) {

    if( a_function == nullptr ) return( ProbabilityBase3dType::none );
    return( a_function->type( ) );
}

/*
============================================================
*/
LUPI_HOST_DEVICE Probabilities::ProbabilityBase3d *serializeProbability3d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase3d *a_probability3d ) {

    int type = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        ProbabilityBase3dType pType = ProbabilityBase3dClass( a_probability3d );

        switch( pType ) {
        case ProbabilityBase3dType::none :
            break;
        case ProbabilityBase3dType::XYs :
            type = 1;
            break;
        }
    }

    DATA_MEMBER_INT( type, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case 0 :
            a_probability3d = nullptr;
            break;
        case 1 :
            if( a_buffer.m_placement != nullptr ) {
                a_probability3d = new(a_buffer.m_placement) Probabilities::XYs3d;
                a_buffer.incrementPlacement( sizeof( Probabilities::XYs3d ) ); }
            else {
                a_probability3d = new Probabilities::XYs3d;
            }
            break;
        }
    }
    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            a_buffer.incrementPlacement( sizeof( Probabilities::XYs3d ) );
            break;
        }
    }

    if( a_probability3d != nullptr ) {
        switch( type ) {
        case 0 :
            break;
        case 1 :
            static_cast<Probabilities::XYs3d *>( a_probability3d )->serialize( a_buffer, a_mode );
            break;
        }
    }

    return( a_probability3d );
}

}           // End of namespace MCGIDI.
