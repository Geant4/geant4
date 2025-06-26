/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef GIDI_data_hpp_included
#define GIDI_data_hpp_included 1

#include <stdio.h>
#include <string>
#include <vector>

namespace GIDI {

/*
============================================================
=========================== Data1d =========================
============================================================
*/
class Data1d {     // BRB: currently not used.

    private:
        std::vector<double> m_xs;
        std::vector<double> m_ys;

    public:
        Data1d( std::size_t a_number, double const *const a_xs );
        Data1d( std::size_t a_number, double const *const a_xs, double const *const a_ys );
        Data1d( std::vector<double> const &a_xs );
        Data1d( std::vector<double> const &a_xs, std::vector<double> const &a_ys );
        Data1d( Data1d const &a_1dData );
        ~Data1d( );

        std::size_t size( ) const { return( m_xs.size( ) ); }

        Data1d operator+( double a_value ) const ;
        Data1d &operator+=( double a_value );
        Data1d operator+( Data1d const &a_rhs ) const ;
        Data1d &operator+=( Data1d const &a_rhs );

        Data1d operator-( double a_value ) const ;
        Data1d &operator-=( double a_value );
        Data1d operator-( Data1d const &a_rhs ) const ;
        Data1d &operator-=( Data1d const &a_rhs );

        Data1d operator*( double a_value ) const ;
        Data1d &operator*=( double a_value );
        Data1d operator/( double a_value ) const ;
        Data1d &operator/=( double a_value );

        void print( std::string const &a_prefix ) const ;
};

/*
============================================================
========================== Vector ==========================
============================================================
*/
class Vector {

    private:
        std::vector<double> m_vector;                       /**< The list of elements, each is a double instance. */

        void writeWithBoundaries2( FILE *a_file, char const *a_format, std::vector<double> const &a_boundaries, double a_epsilon ) const ;

    public:
        Vector( std::size_t a_number = 0 );
        Vector( std::vector<double> const &a_values );
        Vector( std::size_t a_number, double const *a_values );
        Vector( Vector const &a_vector );
        ~Vector( );

        Vector &operator=( Vector const &a_rhs );

        std::size_t size( ) const { return( m_vector.size( ) ); }                                   /**< Returns a number of elements of *this*. */
        void resize( std::size_t a_number, double a_value = 0.0 ) { m_vector.resize( a_number, a_value ); }                     /**< Resizes *this* to *a_number* elements. For details, see std::vector.resize. */
        std::vector<double> &data( ) { return( m_vector ); }

        double &operator[]( std::size_t a_index ) { return( m_vector[a_index] ); }           /**< Returns a reference to the (*a_index*-1)th element. */
        double operator[]( std::size_t a_index ) const { return( m_vector[a_index] ); }      /**< Returns a reference to the (*a_index*-1)th element. */

        Vector operator+( double a_value ) const ;
        Vector &operator+=( double a_value );
        Vector operator+( Vector const &a_rhs ) const ;
        Vector &operator+=( Vector const &a_rhs );

        Vector operator-( double a_value ) const ;
        Vector &operator-=( double a_value );
        Vector operator-( Vector const &a_rhs ) const ;
        Vector &operator-=( Vector const &a_rhs );

        Vector operator*( double a_value ) const ;
        Vector &operator*=( double a_value );

        Vector operator/( double a_value ) const ;
        Vector &operator/=( double a_value );

        void reverse( );

        void setToValueInFlatRange( std::size_t a_start, std::size_t a_end, double a_value );
        double sum( );
        void print( std::string const &a_prefix ) const ;
        void write( FILE *a_file, std::string const &a_prefix ) const ;
        void writeWithBoundaries( FILE *a_file, char const *a_format, std::vector<double> const &a_boundaries, double a_epsilon ) const ;
};

/*
============================================================
========================== Matrix ==========================
============================================================
*/
class Matrix {

    private:
        std::vector<Vector> m_matrix;                       /**< The list of rows, each is a Vector instance. */

    public:
        Matrix( std::size_t a_rows, std::size_t a_columns );
        Matrix( Matrix const &a_gidi_matrix );
        ~Matrix( );
        Matrix &operator=( Matrix const &a_rhs );

        std::size_t size( ) const { return( m_matrix.size( ) ); }                                       /**< Returns the number of rows or *this*. */

        Vector       &operator[]( std::size_t a_index )       { return( m_matrix[a_index] ); }   /**< Returns a reference to the (*a_index*-1)th row. */
        Vector const &operator[]( std::size_t a_index ) const { return( m_matrix[a_index] ); }   /**< Returns a reference to the (*a_index*-1)th row. */

                                                            /** Sets the cell at row **a_row** and column **a_column** to **a_value**. */
        void operator()( std::size_t a_row           /**< The cell's row. */,
                                std::size_t a_column        /**< The cell's row. */,
                                double a_value              /**< The value to put in the cell. */ )
                                        { m_matrix[a_row][a_column] = a_value; }

        std::vector<Vector> const &matrix( ) const { return( m_matrix ); }

        Matrix operator+( double a_value ) const ;
        Matrix &operator+=( double a_value );
        Matrix operator+( Matrix const &a_rhs ) const ;
        Matrix &operator+=( Matrix const &a_rhs );

        Matrix operator-( double a_value ) const ;
        Matrix &operator-=( double a_value );
        Matrix operator-( Matrix const &a_rhs ) const ;
        Matrix &operator-=( Matrix const &a_rhs );

        Matrix operator*( double a_value ) const ;
        Matrix &operator*=( double a_value );

        Matrix operator/( double a_value ) const ;
        Matrix &operator/=( double a_value );

        std::size_t numberOfColumns( ) const ;
                                                        /** Sets the cell at row **a_row** and column **a_column** to **a_value**. */
        void set( std::size_t a_row                     /**< The cell's row. */,
                  std::size_t a_column                  /**< The cell's row. */,
                  double a_value                        /**< The value to put in the cell. */ ) 
                                        {  m_matrix[a_row][a_column] = a_value; }
                                                        /** Sets the row at **a_row** to **a_vector**. */
        void set( std::size_t a_row                     /**< The row to set. */,
                  Vector const &a_vector                /**< The Vector to set at row **a_row**. */ )
                                        {  m_matrix[a_row] = a_vector; }
        void push_back( Vector const &a_vector );
        Matrix transpose( );
        void reverse( );
        void print( std::string const &a_prefixForRow ) const ;
};

}

#endif      // End of GIDI_data_hpp_included
