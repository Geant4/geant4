/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include "MCGIDI.hpp"

namespace MCGIDI {

namespace Distributions {

LUPI_HOST_DEVICE static double coherentPhotoAtomicScatteringIntegrateSub( int a_n, double a_a, double a_logX, double a_energy1, double a_y1, double a_energy2, double a_y2 );
static LUPI_HOST Distribution *parseGIDI2( GIDI::Distributions::Distribution const &a_GIDI_distribution, SetupInfo &a_setupInfo, 
                Transporting::MC const &a_settings );

/*! \class Distribution
 * This class is the base class for all distribution forms.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Distribution::Distribution( ) :
        m_type( Type::none ),
        m_productFrame( GIDI::Frame::lab ),
        m_projectileMass( 0.0 ),
        m_targetMass( 0.0 ),
        m_productMass( 0.0 ) {

}

/* *********************************************************************************************************//**
 * @param a_type                [in]    The Type of the distribution.
 * @param a_distribution        [in]    The GIDI::Distributions::Distribution instance whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST Distribution::Distribution( Type a_type, GIDI::Distributions::Distribution const &a_distribution, SetupInfo &a_setupInfo ) :
        m_type( a_type ),
        m_productFrame( a_distribution.productFrame( ) ),
        m_projectileMass( a_setupInfo.m_protare.projectileMass( ) ),
        m_targetMass( a_setupInfo.m_protare.targetMass( ) ),
        m_productMass( a_setupInfo.m_product1Mass ) {                           // Includes nuclear excitation energy.

}

/* *********************************************************************************************************//**
 * @param a_type                [in]    The Type of the distribution.
 * @param a_productFrame        [in]    The frame of the product's data for distribution.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST Distribution::Distribution( Type a_type, GIDI::Frame a_productFrame, SetupInfo &a_setupInfo ) :
        m_type( a_type ),
        m_productFrame( a_productFrame ),
        m_projectileMass( a_setupInfo.m_protare.projectileMass( ) ),
        m_targetMass( a_setupInfo.m_protare.targetMass( ) ),
        m_productMass( a_setupInfo.m_product1Mass ) {                           // Includes nuclear excitation energy.

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Distribution::~Distribution( ) {

}

/* *********************************************************************************************************//**
 * This method calls the **setModelDBRC_data2* method if the distribution is AngularTwoBody, otherwise it * executes a throw.
 *
 * @param a_modelDBRC_data      [in]    The instance storing data needed to treat the DRRC upscatter mode.
 ***********************************************************************************************************/

LUPI_HOST void Distribution::setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data ) {

    if( type( ) != Type::angularTwoBody ) throw std::runtime_error( "Setting ModelDBRC_data for non-two body distribution is not allowed." );

    static_cast<AngularTwoBody *>( this )->setModelDBRC_data2( a_modelDBRC_data );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Distribution::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int distributionType = distributionTypeToInt( m_type );
    DATA_MEMBER_INT( distributionType, a_buffer, a_mode );
    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) m_type = intToDistributionType( distributionType );

    int frame = 0;
    if( m_productFrame == GIDI::Frame::centerOfMass ) frame = 1;
    DATA_MEMBER_INT( frame, a_buffer, a_mode );
    m_productFrame = GIDI::Frame::lab;
    if( frame == 1 ) m_productFrame = GIDI::Frame::centerOfMass;

    DATA_MEMBER_DOUBLE( m_projectileMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_targetMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_productMass, a_buffer, a_mode );
}

/*! \class AngularTwoBody
 * This class represents the distribution for an outgoing product for a two-body interaction.
 */

/* *********************************************************************************************************//**
 * Base contructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE AngularTwoBody::AngularTwoBody( ) :
        m_residualMass( 0.0 ),
        m_Q( 0.0 ),
        m_twoBodyThreshold( 0.0 ),
        m_Upscatter( false ),
        m_angular( nullptr ),
        m_modelDBRC_data( nullptr )  {

}

/* *********************************************************************************************************//**
 * @param a_angularTwoBody          [in]    The GIDI::Distributions::AngularTwoBody instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST AngularTwoBody::AngularTwoBody( GIDI::Distributions::AngularTwoBody const &a_angularTwoBody, SetupInfo &a_setupInfo ) :
        Distribution( Type::angularTwoBody, a_angularTwoBody, a_setupInfo ),
        m_residualMass( a_setupInfo.m_product2Mass ),                           // Includes nuclear excitation energy.
        m_Q( a_setupInfo.m_Q ),
        m_twoBodyThreshold( a_setupInfo.m_reaction->twoBodyThreshold( ) ),
        m_Upscatter( false ),
        m_angular( Probabilities::parseProbability2d_d1( a_angularTwoBody.angular( ), &a_setupInfo ) ),
        m_modelDBRC_data( nullptr )  {

    if( a_setupInfo.m_protare.projectileIntid( ) == PoPI::Intids::neutron ) {
        m_Upscatter = a_setupInfo.m_reaction->ENDF_MT( ) == 2;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE AngularTwoBody::~AngularTwoBody( ) {

    delete m_angular;
    delete m_modelDBRC_data;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void AngularTwoBody::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_residualMass, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_Q, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_twoBodyThreshold, a_buffer, a_mode );
    DATA_MEMBER_INT( m_Upscatter, a_buffer, a_mode );

    m_angular = serializeProbability2d_d1( a_buffer, a_mode, m_angular );
    m_modelDBRC_data = serializeModelDBRC_data( a_buffer, a_mode, m_modelDBRC_data );
}

/* *********************************************************************************************************//**
 * This method sets *this* *m_modelDBRC_data* to *a_modelDBRC_data*. It also deletes the current *m_modelDBRC_data* member.
 *
 * @param a_modelDBRC_data      [in]    The instance storing data needed to treat the DRRC upscatter mode.
 ***********************************************************************************************************/

LUPI_HOST void AngularTwoBody::setModelDBRC_data2( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data ) {

    delete m_modelDBRC_data;
    m_modelDBRC_data = a_modelDBRC_data;
}

/*! \class Uncorrelated
 * This class represents the distribution for an outgoing product for which the distribution is the product of uncorrelated
 * angular (i.e., P(mu|E)) and energy (i.e., P(E'|E)) distributions.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Uncorrelated::Uncorrelated( ) :
        m_angular( nullptr ),
        m_energy( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_uncorrelated            [in]    The GIDI::Distributions::Uncorrelated instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST Uncorrelated::Uncorrelated( GIDI::Distributions::Uncorrelated const &a_uncorrelated, SetupInfo &a_setupInfo ) :
        Distribution( Type::uncorrelated, a_uncorrelated, a_setupInfo ),
        m_angular( Probabilities::parseProbability2d_d1( a_uncorrelated.angular( ), nullptr ) ),
        m_energy( Probabilities::parseProbability2d( a_uncorrelated.energy( ), &a_setupInfo ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Uncorrelated::~Uncorrelated( ) {

    delete m_angular;
    delete m_energy;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Uncorrelated::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    m_angular = serializeProbability2d_d1( a_buffer, a_mode, m_angular );
    m_energy = serializeProbability2d( a_buffer, a_mode, m_energy );
}

/*! \class Branching3d
 * This class represents the distribution for an outgoing product for which the distribution is the product of uncorrelated
 * angular (i.e., P(mu|E)) and energy (i.e., P(E'|E)) distributions.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Branching3d::Branching3d( ) :
        m_initialStateIndex( -1 ) {

}

/* *********************************************************************************************************//**
 * @param a_branching3d             [in]    The GIDI::Distributions::Branching3dinstance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST Branching3d::Branching3d( GIDI::Distributions::Branching3d const &a_branching3d, SetupInfo &a_setupInfo ) :
        Distribution( Type::branching3d, a_branching3d, a_setupInfo ),
        m_initialStateIndex( -1 ) {

    auto iter = a_setupInfo.m_stateNamesToIndices.find( a_branching3d.initialState( ) );
    if( iter == a_setupInfo.m_stateNamesToIndices.end( ) ) {
        std::string message( "Branching3d: initial state not found: pid = '" + a_branching3d.initialState( ) + "'." );
        throw std::runtime_error( message.c_str( ) );
    }
    m_initialStateIndex = iter->second;

    a_setupInfo.m_initialStateIndex = m_initialStateIndex;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Branching3d::~Branching3d( ) {

}


/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Branching3d::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    DATA_MEMBER_INT( m_initialStateIndex, a_buffer, a_mode );
}

/*! \class EnergyAngularMC
 * This class represents the distribution for an outgoing particle where the distribution is give as 
 * P(E'|E) * P(mu|E,E') where E is the projectile's energy, E' is the product's outgoing energy, mu is the 
 * cosine of the product's outgoing angle relative to the projectile's velocity, P(E'|E) is the probability for E' given E
 * and (P(mu|E,E') is the probability for mu given E and E'.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE EnergyAngularMC::EnergyAngularMC( ) :
        m_energy( nullptr ),
        m_angularGivenEnergy( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_energyAngularMC         [in]    The GIDI::Distributions::EnergyAngularMC instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST EnergyAngularMC::EnergyAngularMC( GIDI::Distributions::EnergyAngularMC const &a_energyAngularMC, SetupInfo &a_setupInfo ) :
        Distribution( Type::energyAngularMC, a_energyAngularMC, a_setupInfo ),
        m_energy( Probabilities::parseProbability2d_d1( a_energyAngularMC.energy( ), nullptr ) ),
        m_angularGivenEnergy( Probabilities::parseProbability3d( a_energyAngularMC.energyAngular( ) ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE EnergyAngularMC::~EnergyAngularMC( ) {

    delete m_energy;
    delete m_angularGivenEnergy;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void EnergyAngularMC::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    m_energy = serializeProbability2d_d1( a_buffer, a_mode, m_energy );
    m_angularGivenEnergy = serializeProbability3d( a_buffer, a_mode, m_angularGivenEnergy );
}

/*! \class AngularEnergyMC
 * This class represents the distribution for an outgoing particle where the distribution is give as 
 * P(mu|E) * P(E'|E,mu) where E is the projectile's energy, E' is the product's outgoing energy, mu is the 
 * cosine of the product's outgoing angle relative to the projectile's velocity, P(mu|E) is the probability for mu given E
 * and (P(E'|E,mu) is the probability for E' given E and mu.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE AngularEnergyMC::AngularEnergyMC( ) :
        m_angular( nullptr ),
        m_energyGivenAngular( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_angularEnergyMC         [in]    The GIDI::Distributions::AngularEnergyMC instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST AngularEnergyMC::AngularEnergyMC( GIDI::Distributions::AngularEnergyMC const &a_angularEnergyMC, SetupInfo &a_setupInfo ) :
        Distribution( Type::angularEnergyMC, a_angularEnergyMC, a_setupInfo ),
        m_angular( Probabilities::parseProbability2d_d1( a_angularEnergyMC.angular( ), nullptr ) ),
        m_energyGivenAngular( Probabilities::parseProbability3d( a_angularEnergyMC.angularEnergy( ) ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE AngularEnergyMC::~AngularEnergyMC( ) {

    delete m_angular;
    delete m_energyGivenAngular;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void AngularEnergyMC::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    m_angular = serializeProbability2d_d1( a_buffer, a_mode, m_angular );
    m_energyGivenAngular = serializeProbability3d( a_buffer, a_mode, m_energyGivenAngular );
}

/*! \class KalbachMann
 * This class represents the distribution for an outgoing product whose distribution is represented by Kalbach-Mann systematics.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE KalbachMann::KalbachMann( ) :
        m_energyToMeVFactor( 0.0 ),
        m_eb_massFactor( 0.0 ),
        m_f( nullptr ),
        m_r( nullptr ),
        m_a( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_KalbachMann             [in]    The GIDI::Distributions::KalbachMann instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST KalbachMann::KalbachMann( GIDI::Distributions::KalbachMann const &a_KalbachMann, SetupInfo &a_setupInfo ) :
        Distribution( Type::KalbachMann, a_KalbachMann, a_setupInfo ),
        m_energyToMeVFactor( 1 ),                                           // FIXME.
        m_eb_massFactor( 1 ),                                               // FIXME.
        m_f( Probabilities::parseProbability2d_d1( a_KalbachMann.f( ), nullptr ) ),
        m_r( Functions::parseFunction2d( a_KalbachMann.r( ) ) ),
        m_a( Functions::parseFunction2d( a_KalbachMann.a( ) ) ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE KalbachMann::~KalbachMann( ) {

    delete m_f;
    delete m_r;
    delete m_a;
}

/* *********************************************************************************************************//**
 * This method evaluates the Kalbach-Mann formalism at the projectile energy a_energy, and outgoing product energy a_energyOut and a_mu.
 *
 * @param a_energy                  [in]    The energy of the projectile in the lab frame.
 * @param a_energyOut               [in]    The energy of the product in the center-of-mass frame.
 * @param a_mu                      [in]    The mu of the product in the center-of-mass frame.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double KalbachMann::evaluate( double a_energy, double a_energyOut, double a_mu ) {

//    double f_0 = m_f->evaluate( a_energy, a_energyOut );
    double rValue = m_r->evaluate( a_energy, a_energyOut );
    double aValue = m_a->evaluate( a_energy, a_energyOut );
//    double pdf_val = aValue * f_0 / 2.0 / sinh( aValue ) * (cosh(aValue * a_mu) + rValue * sinh( aValue * a_mu ) ); // double-differential PDF for a_energyOut and a_mu (Eq. 6.4 in ENDF-102, 2012)
    double pdf_val = aValue / ( 2.0 * sinh( aValue ) ) * ( cosh( aValue * a_mu ) + rValue * sinh( aValue * a_mu ) ); // double-differential PDF for a_energyOut and mu (Eq. 6.4 in ENDF-102, 2012)
    return pdf_val;
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void KalbachMann::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_energyToMeVFactor, a_buffer, a_mode );
    DATA_MEMBER_DOUBLE( m_eb_massFactor, a_buffer, a_mode );

    m_f = serializeProbability2d_d1( a_buffer, a_mode, m_f );
    m_r = serializeFunction2d( a_buffer, a_mode, m_r );
    m_a = serializeFunction2d( a_buffer, a_mode, m_a );
}

/*! \class CoherentPhotoAtomicScattering
 * This class represents the distribution for an outgoing photon via coherent photo-atomic elastic scattering.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering( ) :
        m_realAnomalousFactor( nullptr ),
        m_imaginaryAnomalousFactor( nullptr ) {

}

/* *********************************************************************************************************//**
 * @param a_coherentPhotoAtomicScattering   [in]    GIDI::Distributions::CoherentPhotoAtomicScattering instance whose data is to be used to construct *this*.
 * @param a_setupInfo                       [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering( GIDI::Distributions::CoherentPhotoAtomicScattering const &a_coherentPhotoAtomicScattering, SetupInfo &a_setupInfo ) :
        Distribution( Type::coherentPhotoAtomicScattering, a_coherentPhotoAtomicScattering, a_setupInfo ),
        m_anomalousDataPresent( false ),
        m_realAnomalousFactor( nullptr ),
        m_imaginaryAnomalousFactor( nullptr ) {

    GUPI::Ancestry const *link = a_coherentPhotoAtomicScattering.findInAncestry( a_coherentPhotoAtomicScattering.href( ) );
    GIDI::DoubleDifferentialCrossSection::CoherentPhotoAtomicScattering const &coherentPhotoAtomicScattering = 
            *static_cast<GIDI::DoubleDifferentialCrossSection::CoherentPhotoAtomicScattering const *>( link );

    std::string domainUnit;
    GIDI::Functions::XYs1d const *xys1d0, *xys1d1;
    std::size_t dataSize = 0, offset = 0;

    GIDI::Functions::Function1dForm const *formFactor = coherentPhotoAtomicScattering.formFactor( );
    if( formFactor->type( ) == GIDI::FormType::XYs1d ) {
        xys1d0 = static_cast<GIDI::Functions::XYs1d const *>( formFactor );
        xys1d1 = xys1d0;

        domainUnit = xys1d0->axes( )[0]->unit( );

        dataSize = xys1d1->size( );
        offset = 1; }
    else if( formFactor->type( ) == GIDI::FormType::regions1d ) {
        GIDI::Functions::Regions1d const *regions1d = static_cast<GIDI::Functions::Regions1d const *>( formFactor );
        if( regions1d->size( ) != 2 ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor size." );

        domainUnit = regions1d->axes( )[0]->unit( );

        GIDI::Functions::Function1dForm const *region0 = (*regions1d)[0];
        if( region0->type( ) != GIDI::FormType::XYs1d ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor for region 0." );
        xys1d0 = static_cast<GIDI::Functions::XYs1d const *>( region0 );
        if( xys1d0->size( ) != 2 ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported size of region 1 of form factor." );

        GIDI::Functions::Function1dForm const *region1 = (*regions1d)[1];
        if( region1->type( ) != GIDI::FormType::XYs1d ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor for region 1." );
        xys1d1 = static_cast<GIDI::Functions::XYs1d const *>( region1 );

        dataSize = xys1d1->size( ) + 1; }
    else {
        throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor. Must be XYs1d or regions1d." );
    }

    double domainFactor = 1.0;
    if( domainUnit == "1/Ang" ) {
        domainFactor = 0.012398419739640716; }                      // Converts 'h * c /Ang' to MeV.
    else if( domainUnit == "1/cm" ) {
        domainFactor = 0.012398419739640716 * 1e-8; }               // Converts 'h * c /cm' to MeV.
    else {
        throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported domain unit" );
    }

    m_energies.resize( dataSize );
    m_formFactor.resize( dataSize );
    m_a.resize( dataSize );
    m_integratedFormFactor.resize( dataSize );
    m_integratedFormFactorSquared.resize( dataSize );
    m_probabilityNorm1_1.resize( dataSize );
    m_probabilityNorm1_3.resize( dataSize );
    m_probabilityNorm1_5.resize( dataSize );
    m_probabilityNorm2_1.resize( dataSize );
    m_probabilityNorm2_3.resize( dataSize );
    m_probabilityNorm2_5.resize( dataSize );

    std::pair<double, double> xy = (*xys1d0)[0];
    m_energies[0] = 0.0;
    m_formFactor[0] = xy.second;
    m_a[0] = 0.0;
    m_integratedFormFactor[0] = 0.0;
    m_integratedFormFactorSquared[0] = 0.0;

    xy = (*xys1d1)[offset];
    double energy1 = domainFactor * xy.first;
    double y1 = xy.second;
    m_energies[1] = energy1;
    m_formFactor[1] = y1;
    m_integratedFormFactor[1] = 0.5 * energy1 * energy1 * y1;
    m_integratedFormFactorSquared[1] = 0.5 * energy1 * energy1 * y1 * y1;

    double sum1 = m_integratedFormFactor[1];
    double sum2 = m_integratedFormFactorSquared[1];
    for( std::size_t i1 = 1 + offset; i1 < xys1d1->size( ); ++i1 ) {
        xy = (*xys1d1)[i1];
        double energy2 = domainFactor * xy.first;
        double y2 = xy.second;

        double logEs = log( energy2 / energy1 );
        double _a = log( y2 / y1 ) / logEs;

        m_energies[i1+1-offset] = energy2;
        m_formFactor[i1+1-offset] = y2;
        m_a[i1-offset] = _a;

        sum1 += coherentPhotoAtomicScatteringIntegrateSub( 1,       _a, logEs, energy1,      y1, energy2,      y2 );
        m_integratedFormFactor[i1+1-offset] = sum1;

        sum2 += coherentPhotoAtomicScatteringIntegrateSub( 1, 2.0 * _a, logEs, energy1, y1 * y1, energy2, y2 * y2 );
        m_integratedFormFactorSquared[i1+1-offset] = sum2;

        energy1 = energy2;
        y1 = y2;
    }

    m_a[m_a.size()-1] = 0.0;

    if( coherentPhotoAtomicScattering.realAnomalousFactor( ) != nullptr ) {
        m_anomalousDataPresent = true;
        m_realAnomalousFactor = Functions::parseFunction1d_d1( coherentPhotoAtomicScattering.realAnomalousFactor( ) );
        m_imaginaryAnomalousFactor = Functions::parseFunction1d_d1( coherentPhotoAtomicScattering.imaginaryAnomalousFactor( ) );
    }

    m_probabilityNorm1_1[0] = 0.0;
    m_probabilityNorm1_3[0] = 0.0;
    m_probabilityNorm1_5[0] = 0.0;
    m_probabilityNorm2_1[0] = 0.0;
    m_probabilityNorm2_3[0] = 0.0;
    m_probabilityNorm2_5[0] = 0.0;
    energy1 = m_energies[1];
    y1 = m_formFactor[0];
    for( std::size_t i1 = 1; i1 < m_probabilityNorm1_1.size( ); ++i1 ) {
        double energy2 = m_energies[i1];
        double y2 = m_formFactor[i1];
        double logEs = log( energy2 / energy1 );

        m_probabilityNorm1_1[i1] = m_probabilityNorm1_1[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 1,       m_a[i1-1], logEs, energy1,      y1, energy2,      y2 );
        m_probabilityNorm1_3[i1] = m_probabilityNorm1_3[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 3,       m_a[i1-1], logEs, energy1,      y1, energy2,      y2 );
        m_probabilityNorm1_5[i1] = m_probabilityNorm1_5[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 5,       m_a[i1-1], logEs, energy1,      y1, energy2,      y2 );

        m_probabilityNorm2_1[i1] = m_probabilityNorm2_1[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 1, 2.0 * m_a[i1-1], logEs, energy1, y1 * y1, energy2, y2 * y2 );
        m_probabilityNorm2_3[i1] = m_probabilityNorm2_3[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 3, 2.0 * m_a[i1-1], logEs, energy1, y1 * y1, energy2, y2 * y2 );
        m_probabilityNorm2_5[i1] = m_probabilityNorm2_5[i1-1] + coherentPhotoAtomicScatteringIntegrateSub( 5, 2.0 * m_a[i1-1], logEs, energy1, y1 * y1, energy2, y2 * y2 );

        energy1 = energy2;
        y1 = y2;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE CoherentPhotoAtomicScattering::~CoherentPhotoAtomicScattering( ) {

    delete m_realAnomalousFactor;
    delete m_imaginaryAnomalousFactor;
}

/* *********************************************************************************************************//**
 * This method evaluates the coherent photo-atomic scattering double differentil at the projectile energy a_energy and product cosine of angle a_mu.
 *
 * @param a_energyIn                [in]    The energy of the projectile in the lab frame.
 * @param a_mu                      [in]    The mu of the product in the center-of-mass frame.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double CoherentPhotoAtomicScattering::evaluate( double a_energyIn, double a_mu ) const {

    double probability;
    int intLowerIndexEnergy = binarySearchVector( a_energyIn, m_energies, true );      // FIXME - need to handle case where lowerIndexEnergy = 0 like in evaluateScatteringFactor.
    std::size_t lowerIndexEnergy = static_cast<std::size_t>( intLowerIndexEnergy );
    double _a = m_a[lowerIndexEnergy];
    double _a_2 = _a * _a;
    double X1 = m_energies[lowerIndexEnergy];
    double logEs = log( a_energyIn / X1 );
    double formFactor_1 = m_formFactor[lowerIndexEnergy];
    double formFactor_2 = formFactor_1 * formFactor_1;
    double formFactorEnergyIn_1 = formFactor_1 * pow( a_energyIn / X1, _a );
    double formFactorEnergyIn_2 = formFactorEnergyIn_1 * formFactorEnergyIn_1;
    double inverseEnergyIn_1 = 1.0 / a_energyIn;
    double inverseEnergyIn_2 = inverseEnergyIn_1 * inverseEnergyIn_1;
    double inverseEnergyIn_3 = inverseEnergyIn_1 * inverseEnergyIn_2;
    double inverseEnergyIn_4 = inverseEnergyIn_2 * inverseEnergyIn_2;
    double inverseEnergyIn_5 = inverseEnergyIn_1 * inverseEnergyIn_4;
    double inverseEnergyIn_6 = inverseEnergyIn_2 * inverseEnergyIn_4;

    double norm = 0.5 * inverseEnergyIn_2 * ( m_probabilityNorm2_1[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 1, _a_2, logEs, X1, formFactor_2, a_energyIn, formFactorEnergyIn_2 ) )
                      - inverseEnergyIn_4 * ( m_probabilityNorm2_3[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 3, _a_2, logEs, X1, formFactor_2, a_energyIn, formFactorEnergyIn_2 ) )
                      + inverseEnergyIn_6 * ( m_probabilityNorm2_5[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 5, _a_2, logEs, X1, formFactor_2, a_energyIn, formFactorEnergyIn_2 ) );

    double realAnomalousFactor = 0.0;
    double imaginaryAnomalousFactor = 0.0;
    if( m_anomalousDataPresent ) {
        realAnomalousFactor = m_realAnomalousFactor->evaluate( a_energyIn );
        imaginaryAnomalousFactor = m_imaginaryAnomalousFactor->evaluate( a_energyIn );
        norm += realAnomalousFactor * (         inverseEnergyIn_1 * ( m_probabilityNorm1_1[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 1, _a, logEs, X1, formFactor_1, a_energyIn, formFactorEnergyIn_1 ) )
                                        - 2.0 * inverseEnergyIn_3 * ( m_probabilityNorm1_3[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 3, _a, logEs, X1, formFactor_1, a_energyIn, formFactorEnergyIn_1 ) )
                                        + 2.0 * inverseEnergyIn_5 * ( m_probabilityNorm1_5[lowerIndexEnergy] + coherentPhotoAtomicScatteringIntegrateSub( 5, _a, logEs, X1, formFactor_1, a_energyIn, formFactorEnergyIn_1 ) ) );
    }
    norm *= 16.0;
    norm += 8.0 / 3.0 * ( realAnomalousFactor * realAnomalousFactor + imaginaryAnomalousFactor * imaginaryAnomalousFactor );

    double _formFactor = evaluateFormFactor( a_energyIn, a_mu );
    probability = ( 1.0 + a_mu * a_mu ) * ( ( _formFactor + realAnomalousFactor ) * ( _formFactor + realAnomalousFactor ) + imaginaryAnomalousFactor * imaginaryAnomalousFactor ) / norm;

    return( probability );
}

/* *********************************************************************************************************//**
 * This method evaluates the coherent photo-atomic form factor at the projectile energy a_energy and product cosine of angle a_mu.
 *
 * @param a_energyIn                [in]    The energy of the projectile in the lab frame.
 * @param a_mu                      [in]    The mu of the product in the center-of-mass frame.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double CoherentPhotoAtomicScattering::evaluateFormFactor( double a_energyIn, double a_mu ) const {

    double X = a_energyIn * sqrt( 0.5 * ( 1 - a_mu ) );
    int intLowerIndex = binarySearchVector( X, m_energies );

    if( intLowerIndex < 1 ) {
        if( intLowerIndex == 0 ) return( m_formFactor[0] );
        if( intLowerIndex == -2 ) return( m_formFactor[0] );               // This should never happend for proper a_energyIn and a_mu.
        return( m_formFactor.back( ) );
    }

    std::size_t lowerIndex = static_cast<std::size_t>( intLowerIndex );

    return( m_formFactor[lowerIndex] * pow( X / m_energies[lowerIndex] , m_a[lowerIndex] ) );
}

/* *********************************************************************************************************//**
 * FIX ME.
 *
 * @param a_Z                       [in]    
 * @param a_a                       [in]    
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double CoherentPhotoAtomicScattering::Z_a( double a_Z, double a_a ) const {

    if( fabs( a_a ) < 1e-3 ) {
        double logZ = log( a_Z );
        double a_logZ = a_a * logZ;
        return( logZ * ( 1.0 + 0.5 * a_logZ * ( 1.0 + a_logZ / 3.0 * ( 1.0 + 0.25 * a_logZ ) ) ) );
    }
    return( ( pow( a_Z, a_a ) - 1.0 ) / a_a );
}

/* *********************************************************************************************************//**
 * FIX ME.
 *
 * @param a_n                       [in]    
 * @param a_a                       [in]    
 * @param a_logX                    [in]    
 * @param a_energy1                 [in]    
 * @param a_y1                      [in]    
 * @param a_energy2                 [in]    
 * @param a_y2                      [in]    
 ***********************************************************************************************************/

LUPI_HOST_DEVICE static double coherentPhotoAtomicScatteringIntegrateSub( int a_n, double a_a, double a_logX, double a_energy1, double a_y1, double a_energy2, double a_y2 ) {

    double epsilon = a_a + a_n + 1.0;
    double integral = 0.0;

    if( fabs( epsilon ) < 1e-3 ) {
        double epsilon_logX = epsilon * a_logX;
        integral = a_y1 * pow( a_energy1, a_n + 1.0 ) * a_logX * ( 1.0 + 0.5 * epsilon_logX * ( 1.0 + epsilon_logX / 3.0 * ( 1.0 + 0.25 * epsilon_logX ) ) ); }
    else {
        integral = ( a_y2 * pow( a_energy2, a_n + 1.0 ) - a_y1 * pow( a_energy1, a_n + 1.0 ) ) / epsilon;
    }

    return( integral );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void CoherentPhotoAtomicScattering::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_INT( m_anomalousDataPresent, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_formFactor, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_a, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_integratedFormFactor, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_integratedFormFactorSquared, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm1_1, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm1_3, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm1_5, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm2_1, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm2_3, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_probabilityNorm2_5, a_buffer, a_mode );

    if( m_anomalousDataPresent ) {
        m_realAnomalousFactor = serializeFunction1d_d1( a_buffer, a_mode, m_realAnomalousFactor );
        m_imaginaryAnomalousFactor = serializeFunction1d_d1( a_buffer, a_mode, m_imaginaryAnomalousFactor );
    }
}

/*! \class IncoherentPhotoAtomicScattering
 * This class represents the distribution for an outgoing photon via incoherent photo-atomic elastic scattering.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentPhotoAtomicScattering::IncoherentPhotoAtomicScattering( ) {

}

/* *********************************************************************************************************//**
 * @param a_incoherentPhotoAtomicScattering     [in]    The GIDI::Distributions::IncoherentPhotoAtomicScattering instance whose data is to be used to construct *this*.
 * @param a_setupInfo                           [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST IncoherentPhotoAtomicScattering::IncoherentPhotoAtomicScattering( GIDI::Distributions::IncoherentPhotoAtomicScattering const &a_incoherentPhotoAtomicScattering, 
                SetupInfo &a_setupInfo ) :
        Distribution( Type::incoherentPhotoAtomicScattering, a_incoherentPhotoAtomicScattering, a_setupInfo ) {

    GUPI::Ancestry const *link = a_incoherentPhotoAtomicScattering.findInAncestry( a_incoherentPhotoAtomicScattering.href( ) );
    GIDI::DoubleDifferentialCrossSection::IncoherentPhotoAtomicScattering const &incoherentPhotoAtomicScattering = 
            *static_cast<GIDI::DoubleDifferentialCrossSection::IncoherentPhotoAtomicScattering const *>( link );

    std::string domainUnit;
    GIDI::Functions::XYs1d const *xys1d0, *xys1d1;
    std::size_t dataSize = 0, offset = 0;

    GIDI::Functions::Function1dForm const *scatteringFactor = incoherentPhotoAtomicScattering.scatteringFactor( );
    if( scatteringFactor->type( ) == GIDI::FormType::XYs1d ) {
        xys1d0 = static_cast<GIDI::Functions::XYs1d const *>( scatteringFactor );
        xys1d1 = xys1d0;

        domainUnit = xys1d0->axes( )[0]->unit( );

        dataSize = xys1d1->size( );
        offset = 1; }
    else if( scatteringFactor->type( ) == GIDI::FormType::regions1d ) {
        GIDI::Functions::Regions1d const *regions1d = static_cast<GIDI::Functions::Regions1d const *>( scatteringFactor );
        if( regions1d->size( ) != 2 ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor size." );

        domainUnit = regions1d->axes( )[0]->unit( );

        GIDI::Functions::Function1dForm const *region0 = (*regions1d)[0];
        if( region0->type( ) != GIDI::FormType::XYs1d ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor for region 0." );
        xys1d0 = static_cast<GIDI::Functions::XYs1d const *>( region0 );
        if( xys1d0->size( ) != 2 ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported size of region 1 of form factor." );

        GIDI::Functions::Function1dForm const *region1 = (*regions1d)[1];
        if( region1->type( ) != GIDI::FormType::XYs1d ) throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor for region 1." );
        xys1d1 = static_cast<GIDI::Functions::XYs1d const *>( region1 );

        dataSize = xys1d1->size( ) + 1; }
    else {
        throw std::runtime_error( "MCGIDI::CoherentPhotoAtomicScattering::CoherentPhotoAtomicScattering: unsupported form factor. Must be XYs1d or regions1d." );
    }

    double domainFactor = 1.0;
    if( domainUnit == "1/Ang" ) {
        domainFactor = 0.012398419739640716; }                      // Converts 'h * c /Ang' to MeV.
    else if( domainUnit == "1/cm" ) {
        domainFactor = 0.012398419739640716 * 1e-8; }               // Converts 'h * c /cm' to MeV.
    else {
        throw std::runtime_error( "MCGIDI::IncoherentPhotoAtomicScattering::IncoherentPhotoAtomicScattering: unsupported domain unit" );
    }

    m_energies.resize( dataSize );
    m_scatteringFactor.resize( dataSize );
    m_a.resize( dataSize );

    std::pair<double, double> xy = (*xys1d0)[0];
    m_energies[0] = domainFactor * xy.first;
    m_scatteringFactor[0] = xy.second;
    m_a[0] = 1.0;

    xy = (*xys1d1)[offset];
    double energy1 = domainFactor * xy.first;
    double y1 = xy.second;

    m_energies[1] = energy1;
    m_scatteringFactor[1] = y1;

    for( std::size_t i1 = 1 + offset; i1 < xys1d1->size( ); ++i1 ) {
        xy = (*xys1d1)[i1];
        double energy2 = domainFactor * xy.first;
        double y2 = xy.second;

        m_energies[i1+1-offset] = energy2;
        m_scatteringFactor[i1+1-offset] = y2;

        double _a = log( y2 / y1 ) / log( energy2 / energy1 );
        m_a[i1-offset] = _a;

        energy1 = energy2;
        y1 = y2;
    }
    m_a[m_a.size()-1] = 0.0;
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentPhotoAtomicScattering::~IncoherentPhotoAtomicScattering( ) {

}

/* *********************************************************************************************************//**
 * FIX ME.
 *
 * @param a_energyIn                [in]    
 * @param a_mu                      [in]    
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentPhotoAtomicScattering::energyRatio( double a_energyIn, double a_mu ) const {

    double relativeEnergy = a_energyIn / PoPI_electronMass_MeV_c2;

    return( 1.0 / ( 1.0 + relativeEnergy * ( 1.0 - a_mu ) ) );
}

/* *********************************************************************************************************//**
 * This method evaluates the Klein-Nishina. FIX ME. This should be a function as it does not use member data.
 *
 * @param a_energyIn                [in]    
 * @param a_mu                      [in]    
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentPhotoAtomicScattering::evaluateKleinNishina( double a_energyIn, double a_mu ) const {

    double relativeEnergy = a_energyIn / PoPI_electronMass_MeV_c2;
    double _energyRatio = energyRatio( a_energyIn, a_mu );
    double one_minus_mu = 1.0 - a_mu;

    double norm = ( 1.0 + 2.0 * relativeEnergy );
    norm = 2.0 * relativeEnergy * ( 2.0 + relativeEnergy * ( 1.0 + relativeEnergy ) * ( 8.0 + relativeEnergy ) ) / ( norm * norm );
    norm += ( ( relativeEnergy - 2.0 ) * relativeEnergy - 2.0 ) * log( 1.0 + 2.0 * relativeEnergy );
    norm /= relativeEnergy * relativeEnergy * relativeEnergy;

    return( _energyRatio * _energyRatio * ( _energyRatio + a_mu * a_mu + relativeEnergy * one_minus_mu * one_minus_mu ) / norm );
}

/* *********************************************************************************************************//**
 * This method evaluates the Klein-Nishina.
 *
 * @param a_energyIn                [in]    The energy of the projectile.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentPhotoAtomicScattering::evaluateScatteringFactor( double a_energyIn ) const {

    int intLowerIndex = binarySearchVector( a_energyIn, m_energies );

    if( intLowerIndex < 1 ) {
        if( intLowerIndex == -1 ) return( m_scatteringFactor.back( ) );
        return( m_scatteringFactor[1] * a_energyIn / m_energies[1] );
    }

    std::size_t lowerIndex = static_cast<std::size_t>( intLowerIndex );
    return( m_scatteringFactor[lowerIndex] * pow( a_energyIn / m_energies[lowerIndex], m_a[lowerIndex] ) );
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void IncoherentPhotoAtomicScattering::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_scatteringFactor, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_a, a_buffer, a_mode );
}

/*! \class IncoherentBoundToFreePhotoAtomicScattering
 * This class represents the distribution for an outgoing photon via incoherent photo-atomic elastic scattering.
 */

/* *********************************************************************************************************//**
 * Default constructor used when broadcasting a Protare as needed by MPI or GPUs.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentBoundToFreePhotoAtomicScattering::IncoherentBoundToFreePhotoAtomicScattering( ) {

}

/* *********************************************************************************************************//**
 * @param a_incoherentBoundToFreePhotoAtomicScattering  [in]    The GIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering instance whose data is to be used to construct *this*.
 * @param a_setupInfo                                   [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST IncoherentBoundToFreePhotoAtomicScattering::IncoherentBoundToFreePhotoAtomicScattering( 
                GIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering const &a_incoherentBoundToFreePhotoAtomicScattering,
                SetupInfo &a_setupInfo ) :
        Distribution( Type::incoherentBoundToFreePhotoAtomicScattering, a_incoherentBoundToFreePhotoAtomicScattering, a_setupInfo ),
        m_bindingEnergy( 0.0 ) {

    GIDI::ProtareSingle const &GIDI_protare = a_setupInfo.m_GIDI_protare;
    auto monikers = GIDI_protare.styles( ).findAllOfMoniker( GIDI_MonteCarlo_cdfStyleChars );
    std::string MonteCarlo_cdf = "";
    if( monikers.size( ) == 1 ) {
        MonteCarlo_cdf = monikers[0][0]->label( ); 
    }

    std::string Compton_href = a_incoherentBoundToFreePhotoAtomicScattering.href( );

    if( Compton_href.find( MonteCarlo_cdf ) != std::string::npos ) {
        const GUPI::Ancestry *link = a_incoherentBoundToFreePhotoAtomicScattering.findInAncestry( Compton_href );
        GIDI::DoubleDifferentialCrossSection::IncoherentBoundToFreePhotoAtomicScattering const &dd = *static_cast<GIDI::DoubleDifferentialCrossSection::IncoherentBoundToFreePhotoAtomicScattering const *>( link );
        GIDI::Functions::Xs_pdf_cdf1d const *xpcCompton;
        GIDI::Functions::Function1dForm  const *ComptonProfile = dd.ComptonProfile( );
        xpcCompton = static_cast<GIDI::Functions::Xs_pdf_cdf1d const *>( ComptonProfile );

        std::vector<double> occupationNumbers = xpcCompton->cdf( );
        std::vector<double> pz_grid = xpcCompton->Xs( );
        std::size_t dataSize = pz_grid.size( );
        m_occupationNumber.resize( dataSize );
        m_pz.resize( dataSize );
        for( std::size_t index = 0; index < occupationNumbers.size( ); ++index ) {
            m_occupationNumber[index] = occupationNumbers[index];
            m_pz[index] = pz_grid[index];
        }

        m_bindingEnergy = a_setupInfo.m_Q;
    }
}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentBoundToFreePhotoAtomicScattering::~IncoherentBoundToFreePhotoAtomicScattering( ) {

}

/* *********************************************************************************************************//**
 * FIX ME.
 *
 * @param a_energyIn                [in]
 * @param a_mu                      [in]
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentBoundToFreePhotoAtomicScattering::energyRatio( double a_energyIn, double a_mu ) const {

    double relativeEnergy = a_energyIn / PoPI_electronMass_MeV_c2;

    return( 1.0 / ( 1.0 + relativeEnergy * ( 1.0 - a_mu ) ) );
}

/* *********************************************************************************************************//**
 * This method evaluates the Klein-Nishina. FIX ME. This should be a function as it does not use member data.
 *
 * @param a_energyIn                [in]
 * @param a_mu                      [in]
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentBoundToFreePhotoAtomicScattering::evaluateKleinNishina( double a_energyIn, double a_mu ) const {

    double relativeEnergy = a_energyIn / PoPI_electronMass_MeV_c2;
    double _energyRatio = energyRatio( a_energyIn, a_mu );
    double one_minus_mu = 1.0 - a_mu;

    double norm = ( 1.0 + 2.0 * relativeEnergy );
    norm = 2.0 * relativeEnergy * ( 2.0 + relativeEnergy * ( 1.0 + relativeEnergy ) * ( 8.0 + relativeEnergy ) ) / ( norm * norm );
    norm += ( ( relativeEnergy - 2.0 ) * relativeEnergy - 2.0 ) * log( 1.0 + 2.0 * relativeEnergy );
    norm /= relativeEnergy * relativeEnergy * relativeEnergy;

    return( _energyRatio * _energyRatio * ( _energyRatio + a_mu * a_mu + relativeEnergy * one_minus_mu * one_minus_mu ) / norm );
}


/* *********************************************************************************************************//**
 * This method evaluates the occupation number cdf.
 *
 * @param a_energyIn                [in]
 * @param a_mu                      [in]
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double IncoherentBoundToFreePhotoAtomicScattering::evaluateOccupationNumber( double a_energyIn, double a_mu ) const {

    const double alpha_in = a_energyIn / PoPI_electronMass_MeV_c2;
    //const double mec = 2.7309245307378233e-22; // m_e * c in SI units
    const double alpha_binding = -m_bindingEnergy/PoPI_electronMass_MeV_c2;  // BE [MeV] / 0.511 [MeV]
    const double pzmax = ( -alpha_binding + alpha_in*(alpha_in - alpha_binding)*(1-a_mu) )/( sqrt( 2*alpha_in*(alpha_in-alpha_binding)*(1-a_mu) + alpha_binding*alpha_binding ) ); // *mec

    int intLowerIndex = binarySearchVector( pzmax, m_pz );
    std::size_t lowerIndex = static_cast<std::size_t>( intLowerIndex );
    int size1 = static_cast<int>( m_occupationNumber.size( ) );

    if( intLowerIndex == -1 || intLowerIndex == ( size1 - 1 ) ) {
         return( m_occupationNumber.back( ) );
    }
    if( intLowerIndex == -2 ){
        return( m_occupationNumber[0] );
    }

    return(m_occupationNumber[lowerIndex] + (pzmax-m_pz[lowerIndex])*(m_occupationNumber[lowerIndex+1]-m_occupationNumber[lowerIndex])/(m_pz[lowerIndex+1]-m_pz[lowerIndex]) );

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void IncoherentBoundToFreePhotoAtomicScattering::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_VECTOR_DOUBLE( m_pz, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_occupationNumber, a_buffer, a_mode );

}

/*
======================================================================================================
========== IncoherentPhotoAtomicScatteringElectron                                          ==========
======================================================================================================
*/

/*! \class IncoherentPhotoAtomicScatteringElectron
 * This class represents the distribution for the outgoing electron for incoherent photo-atomic scattering.
 */

/* *********************************************************************************************************//**
 * Plain constructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentPhotoAtomicScatteringElectron::IncoherentPhotoAtomicScatteringElectron( ) {

}

/* *********************************************************************************************************//**
 * Constructor.
 *
 * @param a_setupInfo                           [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST IncoherentPhotoAtomicScatteringElectron::IncoherentPhotoAtomicScatteringElectron( SetupInfo &a_setupInfo ) :
        Distribution( Type::incoherentPhotoAtomicScatteringElectron, GIDI::Frame::lab, a_setupInfo ) {

}

/* *********************************************************************************************************//**
 * Destructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentPhotoAtomicScatteringElectron::~IncoherentPhotoAtomicScatteringElectron( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void IncoherentPhotoAtomicScatteringElectron::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
}

/*
======================================================================================================
========== PairProductionGamma                                                              ==========
======================================================================================================
*/

/*! \class PairProductionGamma
 * This class represents the distribution for an outgoing photon the is the result of an electron annihilating with a positron.
 */

/* *********************************************************************************************************//**
 * Basic constructor.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE PairProductionGamma::PairProductionGamma( ) :
        m_firstSampled( false ) {

}

/* *********************************************************************************************************//**
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_firstSampled            [in]    FIX ME
 ***********************************************************************************************************/

LUPI_HOST PairProductionGamma::PairProductionGamma( SetupInfo &a_setupInfo, bool a_firstSampled ) :
        Distribution( Type::pairProductionGamma, GIDI::Frame::lab, a_setupInfo ),
        m_firstSampled( a_firstSampled ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE PairProductionGamma::~PairProductionGamma( ) {

}


/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void PairProductionGamma::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_INT( m_firstSampled, a_buffer, a_mode );
}

/*! \class CoherentElasticTNSL
 * This class represents the distribution for an outgoing product whose distribution is TNSL coherent elastic scattering.
 * This class samples directly from the Debye/Waller function.
 */

/* *********************************************************************************************************//**
 * Constructor for the CoherentElasticTNSL class.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE CoherentElasticTNSL::CoherentElasticTNSL( ) :
        m_temperatureInterpolation( Interpolation::LINLIN ) {

}

/* *********************************************************************************************************//**
 * Constructor for the CoherentElasticTNSL class.
 *
 * @param a_coherentElasticTNSL     [in]    GIDI::CoherentElastic instance containing the Debye/Waller data.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST CoherentElasticTNSL::CoherentElasticTNSL( GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::CoherentElastic const *a_coherentElasticTNSL,
                SetupInfo &a_setupInfo ) :
        Distribution( Type::coherentElasticTNSL, GIDI::Frame::lab, a_setupInfo ),
        m_temperatureInterpolation( Interpolation::LINLIN ) {

    GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::S_table const &s_table = a_coherentElasticTNSL->s_table( );
    GIDI::Functions::Gridded2d const *gridded2d = dynamic_cast<GIDI::Functions::Gridded2d const *>( s_table.function2d( ) );
    GIDI::Axes const &axes = gridded2d->axes( );

    GIDI::Grid const *axis = dynamic_cast<GIDI::Grid const *>( axes[0] );
    m_temperatureInterpolation = GIDI2MCGIDI_interpolation( ptwXY_stringToInterpolation( axis->interpolation( ).c_str( ) ) );

    double temperatureToMeV_K = 1.0;
    if( axis->unit( ) == "K" ) temperatureToMeV_K = 8.617330337217212e-11;           // This is a kludge until units are properly supported.
    nf_Buffer<double> const *grid = &axis->values( );
    m_temperatures.resize( grid->size( ) );
    for( std::size_t index = 0; index < grid->size( ); ++index ) m_temperatures[index] = temperatureToMeV_K * (*grid)[index];

    axis = dynamic_cast<GIDI::Grid const *>( axes[1] );
    grid = &axis->values( );
    m_energies.resize( grid->size( ) );
    for( std::size_t index = 0; index < grid->size( ); ++index ) m_energies[index] = (*grid)[index];

    GIDI::Array::FullArray fullArray = gridded2d->array( ).constructArray( );
    m_S_table.resize( fullArray.size( ) );
    for( std::size_t index = 0; index < fullArray.m_flattenedValues.size( ); ++index ) m_S_table[index] = fullArray.m_flattenedValues[index];
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void CoherentElasticTNSL::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_VECTOR_DOUBLE( m_temperatures, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_energies, a_buffer, a_mode );
    DATA_MEMBER_VECTOR_DOUBLE( m_S_table, a_buffer, a_mode );

    int interpolation = 0;
    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( m_temperatureInterpolation ) {
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
            m_temperatureInterpolation = Interpolation::FLAT;
            break;
        case 1 :
            m_temperatureInterpolation = Interpolation::LINLIN;
            break;
        case 2 :
            m_temperatureInterpolation = Interpolation::LINLOG;
            break;
        case 3 :
            m_temperatureInterpolation = Interpolation::LOGLIN;
            break;
        case 4 :
            m_temperatureInterpolation = Interpolation::LOGLOG;
            break;
        case 5 :
            m_temperatureInterpolation = Interpolation::OTHER;
            break;
        }
    }
}

/*! \class IncoherentElasticTNSL
 * This class represents the distribution for an outgoing product whose distribution is TNSL incoherent elastic scattering.
 * This class samples directly from the Debye/Waller function.
 */

/* *********************************************************************************************************//**
 * Constructor for the IncoherentElasticTNSL class.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE IncoherentElasticTNSL::IncoherentElasticTNSL( ) :
        m_temperatureToMeV_K( 1.0 ),
        m_DebyeWallerIntegral( nullptr ) {

}

/* *********************************************************************************************************//**
 * Constructor for the IncoherentElasticTNSL class.
 *
 * @param a_incoherentElasticTNSL   [in]    GIDI::IncoherentElastic instance containing the Debye/Waller data.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST IncoherentElasticTNSL::IncoherentElasticTNSL( GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentElastic const *a_incoherentElasticTNSL,
                SetupInfo &a_setupInfo ) :
        Distribution( Type::incoherentElasticTNSL, GIDI::Frame::lab, a_setupInfo ),
        m_temperatureToMeV_K( 1.0 ),
        m_DebyeWallerIntegral( nullptr ) {

    GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::DebyeWallerIntegral const &debyeWallerIntegral = a_incoherentElasticTNSL->debyeWallerIntegral( );
    m_DebyeWallerIntegral = Functions::parseFunction1d_d1( debyeWallerIntegral.function1d( ) );
    GIDI::Axes const &axes = debyeWallerIntegral.function1d( )->axes( );
    GIDI::Axis const *axis = dynamic_cast<GIDI::Axis const *>( axes[0] );
    if( axis->unit( ) == "K" ) m_temperatureToMeV_K = 8.617330337217212e-11;        // This is a kludge until units are properly supported.
}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void IncoherentElasticTNSL::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );

    DATA_MEMBER_DOUBLE( m_temperatureToMeV_K, a_buffer, a_mode );
    m_DebyeWallerIntegral = serializeFunction1d_d1( a_buffer, a_mode, m_DebyeWallerIntegral );
}

/*! \class Unspecified
 * This class represents the distribution for an outgoing product whose distribution is not specified.
 */

LUPI_HOST_DEVICE Unspecified::Unspecified( ) {

}

/* *********************************************************************************************************//**
 * @param a_distribution            [in]    The GIDI::Distributions::Unspecified instance whose data is to be used to construct *this*.
 * @param a_setupInfo               [in]    Used internally when constructing a Protare to pass information to other constructors.
 ***********************************************************************************************************/

LUPI_HOST Unspecified::Unspecified( GIDI::Distributions::Distribution const &a_distribution, SetupInfo &a_setupInfo ) :
        Distribution( Type::unspecified, a_distribution, a_setupInfo ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Unspecified::~Unspecified( ) {

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void Unspecified::serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    Distribution::serialize( a_buffer, a_mode );
}


/* *********************************************************************************************************//**
 * This function is used to call the proper distribution constructor for *a_distribution*.
 *
 * @param a_distribution        [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 ***********************************************************************************************************/

LUPI_HOST Distribution *parseGIDI( GIDI::Suite const &a_distribution, SetupInfo &a_setupInfo, Transporting::MC const &a_settings ) {

    if( a_setupInfo.m_protare.projectileIntid( ) == PoPI::Intids::neutron ) {
        if( a_settings.wantRawTNSL_distributionSampling( ) ) {
            if( a_setupInfo.m_reaction->doubleDifferentialCrossSection( ).size( ) > 0 ) {
                GIDI::Form const *form = a_setupInfo.m_reaction->doubleDifferentialCrossSection( ).get<GIDI::Form>( 0 );

                if( form->type( ) == GIDI::FormType::coherentElastic ) {
                    GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::CoherentElastic const *coherentElasticTNSL 
                            = static_cast<GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::CoherentElastic const *>( form );
                    return( new CoherentElasticTNSL( coherentElasticTNSL, a_setupInfo ) ); }
                else if( form->type( ) == GIDI::FormType::incoherentElastic ) {
                    GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentElastic const *incoherentElasticTNSL 
                            = static_cast<GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentElastic const *>( form );
                    return( new IncoherentElasticTNSL( incoherentElasticTNSL, a_setupInfo ) );
                }
            }
        }
    }

    std::string const *label = a_settings.styles( )->findLabelInLineage( a_distribution, a_setupInfo.m_distributionLabel );
    GIDI::Distributions::Distribution const &GIDI_distribution = *a_distribution.get<GIDI::Distributions::Distribution>( *label );

    return parseGIDI2( GIDI_distribution, a_setupInfo, a_settings );
}

/* *********************************************************************************************************//**
 * This function is used to convert the GIDI distribution *a_GIDI_distribution* into an MCGIDI distribution. It was split off of
 * **parseGIDI** to be able to handle referenced distribution by calling itself.
 *
 * @param a_distribution        [in]    The GIDI::Protare whose data is to be used to construct *this*.
 * @param a_setupInfo           [in]    Used internally when constructing a Protare to pass information to other constructors.
 * @param a_settings            [in]    Used to pass user options to the *this* to instruct it which data are desired.
 ***********************************************************************************************************/

static LUPI_HOST Distribution *parseGIDI2( GIDI::Distributions::Distribution const &a_GIDI_distribution, SetupInfo &a_setupInfo, 
                Transporting::MC const &a_settings ) {

    Distribution *distribution = nullptr;

    GIDI::FormType type = a_GIDI_distribution.type( );

    switch( type ) {
    case GIDI::FormType::angularTwoBody :
        distribution = new AngularTwoBody( static_cast<GIDI::Distributions::AngularTwoBody const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::uncorrelated :
        distribution = new Uncorrelated( static_cast<GIDI::Distributions::Uncorrelated const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::KalbachMann :
        distribution = new KalbachMann( static_cast<GIDI::Distributions::KalbachMann const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::energyAngularMC :
        distribution = new EnergyAngularMC( static_cast<GIDI::Distributions::EnergyAngularMC const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::angularEnergyMC :
        distribution = new AngularEnergyMC( static_cast<GIDI::Distributions::AngularEnergyMC const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::coherentPhotonScattering :
        distribution = new CoherentPhotoAtomicScattering( static_cast<GIDI::Distributions::CoherentPhotoAtomicScattering const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::incoherentPhotonScattering :
        distribution = new IncoherentPhotoAtomicScattering( static_cast<GIDI::Distributions::IncoherentPhotoAtomicScattering const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::incoherentBoundToFreePhotonScattering :
        distribution = new IncoherentBoundToFreePhotoAtomicScattering( static_cast<GIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::branching3d :
        distribution = new Branching3d( static_cast<GIDI::Distributions::Branching3d const &>( a_GIDI_distribution ), a_setupInfo );
        break;
    case GIDI::FormType::unspecified :
        distribution = new Unspecified( a_GIDI_distribution, a_setupInfo );
        break;
    case GIDI::FormType::reference3d : {
        GIDI::Distributions::Reference3d const *reference3d = static_cast<GIDI::Distributions::Reference3d const *>( &a_GIDI_distribution );
        GIDI::Distributions::Distribution const *linkedForm = static_cast<GIDI::Distributions::Distribution const *>( reference3d->findInAncestry( reference3d->href( ) ) );
        if( linkedForm == nullptr ) 
            throw std::runtime_error( "MCGIDI::Distributions::parseGIDI: could not find link '" + a_GIDI_distribution.toXLink( ) + "." );
        distribution = parseGIDI2( *linkedForm, a_setupInfo, a_settings ); }
        break;
    default :
        throw std::runtime_error( "MCGIDI::Distributions::parseGIDI: unsupported distribution: " + a_GIDI_distribution.toXLink( ) + "." );
    }

    return( distribution );
}

/* *********************************************************************************************************//**
 * @param a_distribution        [in]    The GIDI::Protare whose data is to be used to construct *this*.
 *
 * @return                              The type of the distribution or Distributions::Type::none if *a_distribution* is a *nullptr* pointer.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Type DistributionType( Distribution const *a_distribution ) {

    if( a_distribution == nullptr ) return( Type::none );
    return( a_distribution->type( ) );
}

}

/* *********************************************************************************************************//**
 * This method serializes *this* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *this* or unpack *this* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.A
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Distributions::Distribution *serializeDistribution( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, 
                Distributions::Distribution *a_distribution ) {

    Distributions::Type type = Distributions::Type::none;
    if( a_distribution != nullptr ) type = a_distribution->type( );
    int distributionType = distributionTypeToInt( type );
    DATA_MEMBER_INT( distributionType, a_buffer, a_mode );
    type = intToDistributionType( distributionType );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( type ) {
        case Distributions::Type::none :
            a_distribution = nullptr;
            break;
        case Distributions::Type::unspecified :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::Unspecified;
                a_buffer.incrementPlacement( sizeof( Distributions::Unspecified ) ); }
            else {
                a_distribution = new Distributions::Unspecified;
            }
            break;
        case Distributions::Type::angularTwoBody :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::AngularTwoBody;
                a_buffer.incrementPlacement( sizeof( Distributions::AngularTwoBody ) ); }
            else {
                a_distribution = new Distributions::AngularTwoBody;
            }
            break;
        case Distributions::Type::KalbachMann :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::KalbachMann;
                a_buffer.incrementPlacement( sizeof( Distributions::KalbachMann ) ); }
            else {
                a_distribution = new Distributions::KalbachMann;
            }
            break;
        case Distributions::Type::uncorrelated :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::Uncorrelated;
                a_buffer.incrementPlacement( sizeof( Distributions::Uncorrelated ) ); }
            else {
                a_distribution = new Distributions::Uncorrelated;
            }
            break;
        case Distributions::Type::branching3d:
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::Branching3d;
                a_buffer.incrementPlacement( sizeof( Distributions::Branching3d ) ); }
            else {
                a_distribution = new Distributions::Branching3d;
            }
            break;
        case Distributions::Type::energyAngularMC :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::EnergyAngularMC;
                a_buffer.incrementPlacement( sizeof( Distributions::EnergyAngularMC ) ); }
            else {
                a_distribution = new Distributions::EnergyAngularMC;
            }
            break;
        case Distributions::Type::angularEnergyMC :
            if (a_buffer.m_placement != nullptr) {
                a_distribution = new(a_buffer.m_placement) Distributions::AngularEnergyMC;
                a_buffer.incrementPlacement( sizeof( Distributions::AngularEnergyMC ) ); }
            else {
                a_distribution = new Distributions::AngularEnergyMC;
            }
            break;
        case Distributions::Type::coherentPhotoAtomicScattering :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::CoherentPhotoAtomicScattering;
                 a_buffer.incrementPlacement( sizeof( Distributions::CoherentPhotoAtomicScattering ) ); }
             else {
                 a_distribution = new Distributions::CoherentPhotoAtomicScattering;
             }
             break;
        case Distributions::Type::incoherentPhotoAtomicScattering :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::IncoherentPhotoAtomicScattering;
                 a_buffer.incrementPlacement( sizeof( Distributions::IncoherentPhotoAtomicScattering ) ); }
             else {
                 a_distribution = new Distributions::IncoherentPhotoAtomicScattering;
             }
             break;
        case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::IncoherentBoundToFreePhotoAtomicScattering;
                 a_buffer.incrementPlacement( sizeof( Distributions::IncoherentBoundToFreePhotoAtomicScattering ) ); }
             else {
                 a_distribution = new Distributions::IncoherentBoundToFreePhotoAtomicScattering;
             }
             break;
        case Distributions::Type::pairProductionGamma :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::PairProductionGamma;
                 a_buffer.incrementPlacement( sizeof( Distributions::PairProductionGamma ) ); }
             else {
                 a_distribution = new Distributions::PairProductionGamma;
             }
             break;
        case Distributions::Type::coherentElasticTNSL :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::CoherentElasticTNSL;
                 a_buffer.incrementPlacement( sizeof( Distributions::CoherentElasticTNSL ) ); }
             else {
                 a_distribution = new Distributions::CoherentElasticTNSL;
             }
             break;
        case Distributions::Type::incoherentElasticTNSL :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::IncoherentElasticTNSL;
                 a_buffer.incrementPlacement( sizeof( Distributions::IncoherentElasticTNSL ) ); }
             else {
                 a_distribution = new Distributions::IncoherentElasticTNSL;
             }
             break;
        case Distributions::Type::incoherentPhotoAtomicScatteringElectron :
            if( a_buffer.m_placement != nullptr ) {
                 a_distribution = new(a_buffer.m_placement) Distributions::IncoherentPhotoAtomicScatteringElectron;
                 a_buffer.incrementPlacement( sizeof( Distributions::IncoherentPhotoAtomicScatteringElectron ) ); }
             else {
                 a_distribution = new Distributions::IncoherentPhotoAtomicScatteringElectron;
             }
             break;
        }
    }

    if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        switch( type ) {
        case Distributions::Type::none :
            break;
        case Distributions::Type::unspecified :
            a_buffer.incrementPlacement( sizeof( Distributions::Unspecified ) );
            break;
        case Distributions::Type::angularTwoBody :
            a_buffer.incrementPlacement( sizeof( Distributions::AngularTwoBody ) );
            break;
        case Distributions::Type::KalbachMann :
            a_buffer.incrementPlacement( sizeof( Distributions::KalbachMann ) );
            break;
        case Distributions::Type::uncorrelated :
            a_buffer.incrementPlacement( sizeof( Distributions::Uncorrelated ) );
            break;
        case Distributions::Type::branching3d:
            a_buffer.incrementPlacement( sizeof( Distributions::Branching3d ) );
            break;
        case Distributions::Type::energyAngularMC :
            a_buffer.incrementPlacement( sizeof( Distributions::EnergyAngularMC ) );
            break;
        case Distributions::Type::angularEnergyMC :
            a_buffer.incrementPlacement( sizeof( Distributions::AngularEnergyMC ) );
            break;
        case Distributions::Type::coherentPhotoAtomicScattering :
             a_buffer.incrementPlacement( sizeof( Distributions::CoherentPhotoAtomicScattering ) );
             break;
        case Distributions::Type::incoherentPhotoAtomicScattering :
             a_buffer.incrementPlacement( sizeof( Distributions::IncoherentPhotoAtomicScattering ) );
             break;
        case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering :
             a_buffer.incrementPlacement( sizeof( Distributions::IncoherentBoundToFreePhotoAtomicScattering ) );
             break;
        case Distributions::Type::pairProductionGamma :
             a_buffer.incrementPlacement( sizeof( Distributions::PairProductionGamma ) );
             break;
        case Distributions::Type::coherentElasticTNSL :
             a_buffer.incrementPlacement( sizeof( Distributions::CoherentElasticTNSL ) );
             break;
        case Distributions::Type::incoherentElasticTNSL :
             a_buffer.incrementPlacement( sizeof( Distributions::IncoherentElasticTNSL ) );
             break;
        case Distributions::Type::incoherentPhotoAtomicScatteringElectron :
             a_buffer.incrementPlacement( sizeof( Distributions::IncoherentPhotoAtomicScatteringElectron ) );
             break;
        }
    }

    switch( type ) {
    case Distributions::Type::none :
        break;
    case Distributions::Type::unspecified :
        static_cast<Distributions::Unspecified *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::angularTwoBody :
        static_cast<Distributions::AngularTwoBody *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::KalbachMann :
        static_cast<Distributions::KalbachMann *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::uncorrelated :
        static_cast<Distributions::Uncorrelated *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::branching3d:
        static_cast<Distributions::Branching3d *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::energyAngularMC :
        static_cast<Distributions::EnergyAngularMC *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::angularEnergyMC :
        static_cast<Distributions::AngularEnergyMC *>( a_distribution )->serialize( a_buffer, a_mode );
        break;
    case Distributions::Type::coherentPhotoAtomicScattering :
        static_cast<Distributions::CoherentPhotoAtomicScattering *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::incoherentPhotoAtomicScattering :
        static_cast<Distributions::IncoherentPhotoAtomicScattering *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering :
        static_cast<Distributions::IncoherentBoundToFreePhotoAtomicScattering *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::pairProductionGamma :
        static_cast<Distributions::PairProductionGamma *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::coherentElasticTNSL :
        static_cast<Distributions::CoherentElasticTNSL *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::incoherentElasticTNSL :
        static_cast<Distributions::IncoherentElasticTNSL *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    case Distributions::Type::incoherentPhotoAtomicScatteringElectron :
        static_cast<Distributions::IncoherentPhotoAtomicScatteringElectron *>( a_distribution )->serialize( a_buffer, a_mode );
         break;
    }

    return( a_distribution );
}

}
