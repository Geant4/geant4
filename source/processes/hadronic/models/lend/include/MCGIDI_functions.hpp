/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_functions_hpp_included
#define MCGIDI_functions_hpp_included 1

#include <nf_utilities.h>
#include <ptwXY.h>
#include <LUPI_dataBuffer.hpp>

namespace MCGIDI {

enum class Interpolation { LINLIN, LINLOG, LOGLIN, LOGLOG, FLAT, OTHER };
enum class Function1dType { none, constant, XYs, polyomial, gridded, regions, branching, TerrellFissionNeutronMultiplicityModel };
enum class Function2dType { none, XYs };
enum class ProbabilityBase1dType { none, xs_pdf_cdf };
enum class ProbabilityBase2dType { none, XYs, regions, isotropic, discreteGamma, primaryGamma, recoil, NBodyPhaseSpace, evaporation, 
        generalEvaporation, simpleMaxwellianFission, Watt, weightedFunctionals };

enum class ProbabilityBase3dType { none, XYs };

namespace Functions {

/*
============================================================
====================== FunctionBase ========================
============================================================
*/
class FunctionBase {

    private:
        int m_dimension;
        double m_domainMin;
        double m_domainMax;
        Interpolation m_interpolation;
        double m_outerDomainValue;

    public:
        LUPI_HOST_DEVICE FunctionBase( );
        LUPI_HOST   FunctionBase( GIDI::Functions::FunctionForm const &a_function );
        LUPI_HOST_DEVICE FunctionBase( int a_dimension, double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue = 0 );
        LUPI_HOST_DEVICE virtual ~FunctionBase( ) = 0;

        LUPI_HOST_DEVICE Interpolation interpolation( ) const { return( m_interpolation ); }
        LUPI_HOST_DEVICE double domainMin( ) const { return( m_domainMin ); }
        LUPI_HOST_DEVICE double domainMax( ) const { return( m_domainMax ); }
        LUPI_HOST_DEVICE double outerDomainValue( ) const { return( m_outerDomainValue ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Function1d ========================
============================================================
*/
class Function1d : public FunctionBase {

    protected:
        Function1dType m_type;

    public:
        LUPI_HOST_DEVICE Function1d( );
        LUPI_HOST_DEVICE Function1d( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue = 0 );
        LUPI_HOST_DEVICE ~Function1d( );

        LUPI_HOST_DEVICE Function1dType type( ) const { return( m_type ); }
        LUPI_HOST_DEVICE String typeString( ) const ;

        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION int sampleBoundingInteger( double a_x1, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double evaluate( double a_x1 ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== Function1d_d1 =======================
============================================================
*/
class Function1d_d1 : public Function1d {

    public:
        LUPI_HOST_DEVICE Function1d_d1( ) :
            Function1d( ) { }
        LUPI_HOST_DEVICE Function1d_d1( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue = 0 ) :
                Function1d( a_domainMin, a_domainMax, a_interpolation, a_outerDomainValue ) { }

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
};

/*
============================================================
====================== Function1d_d2 =======================
============================================================
*/
class Function1d_d2 : public Function1d_d1 {

    public:
        LUPI_HOST_DEVICE Function1d_d2( ) :
            Function1d_d1( ) { }
        LUPI_HOST_DEVICE Function1d_d2( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue = 0 ) :
                Function1d_d1( a_domainMin, a_domainMax, a_interpolation, a_outerDomainValue ) { }

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
};

/*
============================================================
======================== Constant1d ========================
============================================================
*/
class Constant1d : public Function1d_d2 {

    private:
        double m_value;

    public:
        LUPI_HOST_DEVICE Constant1d( );
        LUPI_HOST_DEVICE Constant1d( double a_domainMin, double a_domainMax, double a_value, double a_outerDomainValue = 0 );
        LUPI_HOST Constant1d( GIDI::Functions::Constant1d const &a_constant1d );
        LUPI_HOST_DEVICE ~Constant1d( );

        LUPI_HOST_DEVICE double evaluate( LUPI_maybeUnused double a_x1 ) const { return( m_value ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=========================== XYs1d ==========================
============================================================
*/
class XYs1d : public Function1d_d2 {

    private:
        Vector<double> m_Xs;
        Vector<double> m_Ys;

    public:
        LUPI_HOST_DEVICE XYs1d( );
        LUPI_HOST XYs1d( Interpolation a_interpolation, Vector<double> a_Xs, Vector<double> a_Ys, double a_outerDomainValue = 0 );
        LUPI_HOST XYs1d( GIDI::Functions::XYs1d const &a_XYs1d );
        LUPI_HOST_DEVICE ~XYs1d( );

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================= Polynomial1d =======================
============================================================
*/          
class Polynomial1d : public Function1d_d2 {

    private:
        Vector<double> m_coefficients;
        Vector<double> m_coefficientsReversed;

    public:
        LUPI_HOST_DEVICE Polynomial1d( );
        LUPI_HOST Polynomial1d( double a_domainMin, double a_domainMax, Vector<double> const &a_coefficients, double a_outerDomainValue = 0 );
        LUPI_HOST Polynomial1d( GIDI::Functions::Polynomial1d const &a_polynomial1d );
        LUPI_HOST_DEVICE ~Polynomial1d( );

        LUPI_HOST_DEVICE Vector<double> const &coefficients( ) const { return( m_coefficients ); }
        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================= Gridded1d ========================
============================================================
*/
class Gridded1d : public Function1d_d2 {

    private:
        Vector<double> m_grid;
        Vector<double> m_data;

    public:
        LUPI_HOST_DEVICE Gridded1d( );
        LUPI_HOST Gridded1d( GIDI::Functions::Gridded1d const &a_gridded1d );
        LUPI_HOST_DEVICE ~Gridded1d( );

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================= Regions1d ========================
============================================================
*/
class Regions1d : public Function1d_d1 {

    private:
        Vector<double> m_Xs;
        Vector<Function1d_d2 *> m_functions1d;

    public:
        LUPI_HOST_DEVICE Regions1d( );
        LUPI_HOST Regions1d( GIDI::Functions::Regions1d const &a_regions1d );
        LUPI_HOST_DEVICE ~Regions1d( );

        LUPI_HOST_DEVICE void append( Function1d_d2 *a_function1d );
        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Branching1d =======================
============================================================
*/
class Branching1d : public Function1d_d2 {

    private:
        int m_initialStateIndex;

    public:
        LUPI_HOST_DEVICE Branching1d( );
        LUPI_HOST Branching1d( SetupInfo &a_setupInfo, GIDI::Functions::Branching1d const &a_branching1d );
        LUPI_HOST_DEVICE ~Branching1d( );

        LUPI_HOST_DEVICE int initialStateIndex( ) const { return( m_initialStateIndex ); }

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========== TerrellFissionNeutronMultiplicityModel ==========
============================================================
*/
class TerrellFissionNeutronMultiplicityModel : public Function1d {

    private:
        double m_width;
        Function1d_d1 *m_multiplicity;

    public:
        LUPI_HOST_DEVICE TerrellFissionNeutronMultiplicityModel( );
        LUPI_HOST TerrellFissionNeutronMultiplicityModel( double a_width, Function1d_d1 *a_multiplicity );
        LUPI_HOST_DEVICE ~TerrellFissionNeutronMultiplicityModel( );

        template <typename RNG>
        LUPI_HOST_DEVICE int sampleBoundingInteger( double a_energy, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE double evaluate( double a_energy ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Function2d ========================
============================================================
*/
class Function2d : public FunctionBase {

    protected:
        Function2dType m_type;

    public:
        LUPI_HOST_DEVICE Function2d( );
        LUPI_HOST Function2d( double a_domainMin, double a_domainMax, Interpolation a_interpolation, double a_outerDomainValue = 0 );
        LUPI_HOST_DEVICE ~Function2d( );

        LUPI_HOST_DEVICE Function2dType type( ) const { return m_type; }
        LUPI_HOST_DEVICE String typeString( ) const ;

        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double evaluate( double a_x2, double a_x1 ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=========================== XYs2d ==========================
============================================================
*/
class XYs2d : public Function2d {

    private:
        Vector<double> m_Xs;
        Vector<Function1d_d1 *> m_functions1d;

    public:
        LUPI_HOST_DEVICE XYs2d( );
        LUPI_HOST XYs2d( GIDI::Functions::XYs2d const &a_XYs2d );
        LUPI_HOST_DEVICE ~XYs2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== others ==========================
============================================================
*/
LUPI_HOST Function1d *parseMultiplicityFunction1d( SetupInfo &a_setupInfo, Transporting::MC const &a_settings, GIDI::Suite const &a_suite );
LUPI_HOST Function1d_d1 *parseFunction1d_d1( Transporting::MC const &a_settings, GIDI::Suite const &a_suite );
LUPI_HOST Function1d_d1 *parseFunction1d_d1( GIDI::Functions::Function1dForm const *form1d );
LUPI_HOST Function1d_d2 *parseFunction1d_d2( GIDI::Functions::Function1dForm const *form1d );
LUPI_HOST Function2d *parseFunction2d( Transporting::MC const &a_settings, GIDI::Suite const &a_suite );
LUPI_HOST Function2d *parseFunction2d( GIDI::Functions::Function2dForm const *form2d );

}           // End of namespace Functions.

/*
============================================================
============================================================
================== namespace Probabilities ==================
============================================================
============================================================
*/
namespace Probabilities {

/*
============================================================
===================== ProbabilityBase ======================
============================================================
*/
class ProbabilityBase : public Functions::FunctionBase {

    protected:
        Vector<double> m_Xs;

    public:

        LUPI_HOST_DEVICE ProbabilityBase( );
        LUPI_HOST ProbabilityBase( GIDI::Functions::FunctionForm const &a_probabilty );
        LUPI_HOST ProbabilityBase( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs );
        LUPI_HOST_DEVICE ~ProbabilityBase( );
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
===================== ProbabilityBase1d ====================
============================================================
*/
class ProbabilityBase1d : public ProbabilityBase {

    protected:
        ProbabilityBase1dType m_type;

    public:
        LUPI_HOST_DEVICE ProbabilityBase1d( );
        LUPI_HOST ProbabilityBase1d( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs );
        LUPI_HOST_DEVICE ~ProbabilityBase1d( );

        LUPI_HOST_DEVICE ProbabilityBase1dType type( ) const { return m_type; }
        LUPI_HOST_DEVICE String typeString( ) const ;

        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double evaluate( double a_x1 ) const MCGIDI_TRUE_VIRTUAL;
        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double sample( double a_rngValue, RNG && a_rng ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================= Xs_pdf_cdf1d =======================
============================================================
*/
class Xs_pdf_cdf1d : public ProbabilityBase1d {

    private:
        Vector<double> m_pdf;
        Vector<double> m_cdf;

    public:
        LUPI_HOST_DEVICE Xs_pdf_cdf1d( );
        LUPI_HOST Xs_pdf_cdf1d( GIDI::Functions::Xs_pdf_cdf1d const &a_xs_pdf_cdf1d );
        LUPI_HOST_DEVICE ~Xs_pdf_cdf1d( );

        LUPI_HOST_DEVICE double evaluate( double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
===================== ProbabilityBase2d ====================
============================================================
*/
class ProbabilityBase2d : public ProbabilityBase {

    protected:
        ProbabilityBase2dType m_type;

    public:
        LUPI_HOST_DEVICE ProbabilityBase2d( );
        LUPI_HOST ProbabilityBase2d( GIDI::Functions::FunctionForm const &a_probabilty );
        LUPI_HOST ProbabilityBase2d( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs );
        LUPI_HOST_DEVICE ~ProbabilityBase2d( );

        LUPI_HOST_DEVICE ProbabilityBase2dType type( ) const { return m_type; }
        LUPI_HOST_DEVICE String typeString( ) const ;

        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double evaluate( double a_x2, double a_x1 ) const MCGIDI_TRUE_VIRTUAL;
        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double sample( double a_x2, double a_rngValue, RNG && a_rng ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=================== ProbabilityBase2d_d1 ===================
============================================================
*/
class ProbabilityBase2d_d1 : public ProbabilityBase2d {

    public:
        LUPI_HOST_DEVICE ProbabilityBase2d_d1( ) :
                ProbabilityBase2d( ) { }
        LUPI_HOST ProbabilityBase2d_d1( GIDI::Functions::FunctionForm const &a_probabilty ) :
                ProbabilityBase2d( a_probabilty ) { }
        LUPI_HOST ProbabilityBase2d_d1( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs ) :
                ProbabilityBase2d( a_probabilty, a_Xs ) { }

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng, double *a_x1_1, double *a_x1_2 ) const ;
};

/*
============================================================
=================== ProbabilityBase2d_d2 ===================
============================================================
*/
class ProbabilityBase2d_d2 : public ProbabilityBase2d_d1 {

    public:
        LUPI_HOST_DEVICE ProbabilityBase2d_d2( ) :
                ProbabilityBase2d_d1( ) { }
        LUPI_HOST ProbabilityBase2d_d2( GIDI::Functions::FunctionForm const &a_probabilty ) :
                ProbabilityBase2d_d1( a_probabilty ) { }
        LUPI_HOST ProbabilityBase2d_d2( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs ) :
                ProbabilityBase2d_d1( a_probabilty, a_Xs ) { }

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng, double *a_x1_1, double *a_x1_2 ) const ;
};

/*
============================================================
========================== XYs2d ===========================
============================================================
*/
class XYs2d : public ProbabilityBase2d_d2 {

    private:
        Vector<ProbabilityBase1d *> m_probabilities;

    public:
        LUPI_HOST_DEVICE XYs2d( );
        LUPI_HOST XYs2d( GIDI::Functions::XYs2d const &a_XYs2d );
        LUPI_HOST_DEVICE ~XYs2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng, double *a_x1_1, double *a_x1_2 ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Regions2d =========================
============================================================
*/
class Regions2d : public ProbabilityBase2d_d1 {

    private:
        Vector<ProbabilityBase2d_d2 *> m_probabilities;

    public:
        LUPI_HOST_DEVICE Regions2d( );
        LUPI_HOST Regions2d( GIDI::Functions::Regions2d const &a_regions2d );
        LUPI_HOST_DEVICE ~Regions2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Isotropic2d =======================
============================================================
*/
class Isotropic2d : public ProbabilityBase2d_d2 {

    public:
        LUPI_HOST_DEVICE Isotropic2d( );
        LUPI_HOST Isotropic2d( GIDI::Functions::Isotropic2d const &a_isotropic2d );
        LUPI_HOST_DEVICE ~Isotropic2d( );

        LUPI_HOST_DEVICE double evaluate( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const { return( 0.5 ); }
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( LUPI_maybeUnused double a_x2, double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const { return( 1. - 2. * a_rngValue ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) { 
            ProbabilityBase2d::serialize( a_buffer, a_mode ); }
};

/*
============================================================
====================== DiscreteGamma2d =====================
============================================================
*/
class DiscreteGamma2d : public ProbabilityBase2d_d2 {

    private:
        double m_value;

    public:
        LUPI_HOST_DEVICE DiscreteGamma2d( );
        LUPI_HOST DiscreteGamma2d( GIDI::Functions::DiscreteGamma2d const &a_discreteGamma2d );
        LUPI_HOST_DEVICE ~DiscreteGamma2d( );

        LUPI_HOST_DEVICE double evaluate( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const { return( m_value ); }        // FIXME This is wrong, should be something like 1 when domainMin <= a_x1 <= domainMax ), I think. I.e., should be a probability.
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const { return( m_value ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== PrimaryGamma2d =====================
============================================================
*/
class PrimaryGamma2d : public ProbabilityBase2d_d2 {

    private:
        double m_primaryEnergy;
        double m_massFactor;
        String m_finalState;
        int m_initialStateIndex;

    public:
        LUPI_HOST_DEVICE PrimaryGamma2d( );
        LUPI_HOST PrimaryGamma2d( GIDI::Functions::PrimaryGamma2d const &a_primaryGamma2d, SetupInfo *a_setupInfo );
        LUPI_HOST_DEVICE ~PrimaryGamma2d( );

        double primaryEnergy( ) const { return( m_primaryEnergy ); }                /**< Returns the value of the *m_primaryEnergy* member. */
        double massFactor( ) const { return( m_massFactor ); }                      /**< Returns the value of the *m_massFactor* member. */
        String const &finalState( ) const { return( m_finalState ); }               /**< Returns a const reference to the *m_finalState* member. */
        int initialStateIndex( ) const { return( m_initialStateIndex ); }           /**< Returns the value of the *m_initialStateIndex* member. */

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, LUPI_maybeUnused double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const { return( m_primaryEnergy + a_x2 * m_massFactor ); }
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================= Recoil2d =========================
============================================================
*/
class Recoil2d: public ProbabilityBase2d_d2 {

    private:
        String m_xlink;

    public:
        LUPI_HOST_DEVICE Recoil2d( );
        LUPI_HOST Recoil2d( GIDI::Functions::Recoil2d const &a_recoil2d );
        LUPI_HOST_DEVICE ~Recoil2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
==================== NBodyPhaseSpace2d =====================
============================================================
*/
class NBodyPhaseSpace2d : public ProbabilityBase2d_d2 {

    private:
        int m_numberOfProducts;
        double m_mass;
        double m_energy_in_COMFactor;
        double m_massFactor;
        double m_Q;
        ProbabilityBase1d *m_dist;

    public:
        LUPI_HOST_DEVICE NBodyPhaseSpace2d( );
        LUPI_HOST NBodyPhaseSpace2d( GIDI::Functions::NBodyPhaseSpace2d const &a_NBodyPhaseSpace2d, SetupInfo *a_setupInfo );
        LUPI_HOST_DEVICE ~NBodyPhaseSpace2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== Evaporation2d =======================
============================================================
*/
class Evaporation2d: public ProbabilityBase2d_d2 {

    private:
        double m_U;
        Functions::Function1d_d1 *m_theta;

    public:
        LUPI_HOST_DEVICE Evaporation2d( );
        LUPI_HOST Evaporation2d( GIDI::Functions::Evaporation2d const &a_generalEvaporation2d );
        LUPI_HOST_DEVICE ~Evaporation2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=================== GeneralEvaporation2d ===================
============================================================
*/
class GeneralEvaporation2d: public ProbabilityBase2d_d2 {

    private:
        Functions::Function1d_d1 *m_theta;
        ProbabilityBase1d *m_g;

    public:
        LUPI_HOST_DEVICE GeneralEvaporation2d( );
        LUPI_HOST GeneralEvaporation2d( GIDI::Functions::GeneralEvaporation2d const &a_generalEvaporation2d );
        LUPI_HOST_DEVICE ~GeneralEvaporation2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
================= SimpleMaxwellianFission2d ================
============================================================
*/
class SimpleMaxwellianFission2d: public ProbabilityBase2d_d2 {

    private:
        double m_U;
        Functions::Function1d_d1 *m_theta;

    public:
        LUPI_HOST_DEVICE SimpleMaxwellianFission2d( );
        LUPI_HOST SimpleMaxwellianFission2d( GIDI::Functions::SimpleMaxwellianFission2d const &a_simpleMaxwellianFission2d );
        LUPI_HOST_DEVICE ~SimpleMaxwellianFission2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== Watt2d ==========================
============================================================
*/
class Watt2d : public ProbabilityBase2d_d2 {

    private:
        double m_U;
        Functions::Function1d_d1 *m_a;
        Functions::Function1d_d1 *m_b;

    public:
        LUPI_HOST_DEVICE Watt2d( );
        LUPI_HOST Watt2d( GIDI::Functions::Watt2d const &a_Watt2d );
        LUPI_HOST_DEVICE ~Watt2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
=================== WeightedFunctionals2d ==================
============================================================
*/
class WeightedFunctionals2d: public ProbabilityBase2d {

    private:
      Vector<Functions::Function1d_d1 *> m_weight;
      Vector<ProbabilityBase2d_d1 *> m_energy;

    public:
        LUPI_HOST_DEVICE WeightedFunctionals2d( );
        LUPI_HOST WeightedFunctionals2d( GIDI::Functions::WeightedFunctionals2d const &a_weightedFunctionals2d );
        LUPI_HOST_DEVICE ~WeightedFunctionals2d( );

        LUPI_HOST_DEVICE double evaluate( double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
===================== ProbabilityBase3d ====================
============================================================
*/
class ProbabilityBase3d : public ProbabilityBase {

    protected:
        ProbabilityBase3dType m_type;

    public:
        LUPI_HOST_DEVICE ProbabilityBase3d( );
        LUPI_HOST ProbabilityBase3d( GIDI::Functions::FunctionForm const &a_probabilty, Vector<double> const &a_Xs );
        LUPI_HOST_DEVICE ~ProbabilityBase3d( );

        LUPI_HOST_DEVICE ProbabilityBase3dType type( ) const { return m_type; }
        LUPI_HOST_DEVICE String typeString( ) const ;

        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double evaluate( double a_x3, double a_x2, double a_x1 ) const MCGIDI_TRUE_VIRTUAL;
        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double sample( double a_x3, double a_x2_1, double a_x2_2, double a_rngValue, RNG && a_rng ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== XYs3d ===========================
============================================================
*/
class XYs3d : public ProbabilityBase3d {

    private:
        Vector<ProbabilityBase2d_d1 *> m_probabilities;

    public:
        LUPI_HOST_DEVICE XYs3d( );
        LUPI_HOST XYs3d( GIDI::Functions::XYs3d const &a_XYs3d );
        LUPI_HOST_DEVICE ~XYs3d( );

        LUPI_HOST_DEVICE double evaluate( double a_x3, double a_x2, double a_x1 ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double sample( double a_x3, double a_x2_1, double a_x2_2, double a_rngValue, RNG && a_rng ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== others ==========================
============================================================
*/
LUPI_HOST ProbabilityBase1d *parseProbability1d( Transporting::MC const &a_settings, GIDI::Suite const &a_suite );
LUPI_HOST ProbabilityBase1d *parseProbability1d( GIDI::Functions::Function1dForm const *form1d );
LUPI_HOST ProbabilityBase2d *parseProbability2d( Transporting::MC const &a_settings, GIDI::Suite const &a_suite, SetupInfo *a_setupInfo );
LUPI_HOST ProbabilityBase2d *parseProbability2d( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo );
LUPI_HOST ProbabilityBase2d_d1 *parseProbability2d_d1( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo );
LUPI_HOST ProbabilityBase2d_d2 *parseProbability2d_d2( GIDI::Functions::Function2dForm const *form2d, SetupInfo *a_setupInfo );
LUPI_HOST ProbabilityBase3d *parseProbability3d( Transporting::MC const &a_settings, GIDI::Suite const &a_suite );
LUPI_HOST ProbabilityBase3d *parseProbability3d( GIDI::Functions::Function3dForm const *form3d );


}           // End of namespace Probabilities.

/*
============================================================
========================== others ==========================
============================================================
*/
LUPI_HOST_DEVICE Interpolation GIDI2MCGIDI_interpolation( ptwXY_interpolation a_interpolation );

LUPI_HOST_DEVICE Function1dType Function1dClass( Functions::Function1d *funct );
LUPI_HOST_DEVICE Functions::Function1d *serializeFunction1d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Functions::Function1d *a_function1d );
LUPI_HOST_DEVICE Functions::Function1d_d1 *serializeFunction1d_d1( LUPI::DataBuffer &a_buffer, 
                LUPI::DataBuffer::Mode a_mode, Functions::Function1d_d1 *a_function1d );
LUPI_HOST_DEVICE Functions::Function1d_d2 *serializeFunction1d_d2( LUPI::DataBuffer &a_buffer, 
                LUPI::DataBuffer::Mode a_mode, Functions::Function1d_d2 *a_function1d );

LUPI_HOST_DEVICE Function2dType Function2dClass( Functions::Function2d *funct );
LUPI_HOST_DEVICE Functions::Function2d *serializeFunction2d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Functions::Function2d *a_function2d );

LUPI_HOST_DEVICE ProbabilityBase1dType ProbabilityBase1dClass( Probabilities::ProbabilityBase1d *funct );
LUPI_HOST_DEVICE Probabilities::ProbabilityBase1d *serializeProbability1d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase1d *a_probability1d );

LUPI_HOST_DEVICE ProbabilityBase2dType ProbabilityBase2dClass( Probabilities::ProbabilityBase2d *funct );
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d *serializeProbability2d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Probabilities::ProbabilityBase2d *a_probability2d );
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *serializeProbability2d_d1( LUPI::DataBuffer &a_buffer, 
                LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase2d_d1 *a_probability2d );
LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d2 *serializeProbability2d_d2( LUPI::DataBuffer &a_buffer, 
                LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase2d_d2 *a_probability2d );

LUPI_HOST_DEVICE ProbabilityBase3dType ProbabilityBase3dClass( Probabilities::ProbabilityBase3d *funct );
LUPI_HOST_DEVICE Probabilities::ProbabilityBase3d *serializeProbability3d( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Probabilities::ProbabilityBase3d *a_probability3d );

}           // End of namespace MCGIDI.

#endif      // End of MCGIDI_functions_hpp_included
