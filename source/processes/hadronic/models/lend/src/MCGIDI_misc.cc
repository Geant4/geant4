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

/*! \class SetupInfo
 * This class is used internally when constructing a Protare to pass internal information to other constructors.
 */

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST SetupInfo::SetupInfo( ProtareSingle &a_protare, GIDI::ProtareSingle const &a_GIDI_protare, PoPI::Database const &a_popsUser, 
                PoPI::Database const &a_pops ) :
        m_protare( a_protare ),
        m_GIDI_protare( a_GIDI_protare ),
        m_popsUser( a_popsUser ),
        m_pops( a_pops ),
        m_neutronIndex( MCGIDI_popsIndex( a_popsUser, PoPI::IDs::neutron ) ),
        m_photonIndex( MCGIDI_popsIndex( a_popsUser, PoPI::IDs::photon ) ),
        m_initialStateIndex( -1 ),
        m_GRIN_continuumGammas( nullptr ) {

}

/* *********************************************************************************************************//**
 ***********************************************************************************************************/

LUPI_HOST SetupInfo::~SetupInfo( ) {

    for( auto iter = m_ACE_URR_probabilityTablesFromGIDI.begin( ); iter != m_ACE_URR_probabilityTablesFromGIDI.end( ); ++iter ) delete (*iter).second;
}

/* *********************************************************************************************************//**
 * This function returns the intid for particle *a_id* or -1 if *a_id* not in *a_pops*.
 *
 * @param       a_pops      [in]    A PoPI::Database to retrived the particle's intid from.
 * @param       a_id        [in]    The GNDS PoPs id of the particle whose intid is requested.
 *
 * @return                          The *intid*.
 ***********************************************************************************************************/

LUPI_HOST int MCGIDI_popsIntid( PoPI::Database const &a_pops, std::string const &a_id ) {

    if( a_id == PoPI::IDs::FissionProductENDL99120 ) return( PoPI::Intids::FissionProductENDL99120 );
    if( a_id == PoPI::IDs::FissionProductENDL99125 ) return( PoPI::Intids::FissionProductENDL99125 );
    int intid =  a_pops.intid( a_id );

    if( intid < 0 ) {
        if(      a_id == PoPI::IDs::neutron ) {
            intid = PoPI::Intids::neutron; }
        else if( a_id == PoPI::IDs::photon ) {
            intid = PoPI::Intids::photon ; }
        else {
            PoPI::ParseIdInfo parseIdInfo( a_id );
            if( parseIdInfo.isSupported( ) ) {
                if( parseIdInfo.isNuclear( ) ) {
                    intid = 1000 * parseIdInfo.Z( ) + parseIdInfo.A( );
                }
            }
        }
    }
    return( intid );
}

/* *********************************************************************************************************//**
 * This function returns the index in *a_pops* for particle *a_id* or -1 if *a_id* not in *a_pops*.
 *
 * @param       a_pops      [in]    A PoPI::Database to retrived the particle's index from.
 * @param       a_id        [in]    The GNDS PoPs id of the particle whose index is requested.
 *
 * @return                          The *index*.
 ***********************************************************************************************************/

LUPI_HOST int MCGIDI_popsIndex( PoPI::Database const &a_pops, std::string const &a_id ) {

    if( !a_pops.exists( a_id ) ) return( -1 );
    return( a_pops[a_id] );
}

/* *********************************************************************************************************//**
 * @param           a_vector    [in]    The GIDI::Vector whose contents are coped to a MCGIGI::Vector.
 *
 * @return                              The MCGIGI::Vector.
 ***********************************************************************************************************/


LUPI_HOST Vector<double> GIDI_VectorDoublesToMCGIDI_VectorDoubles( GIDI::Vector a_vector ) {

    Vector<double> vector( a_vector.size( ) );

    for( std::size_t i1 = 0; i1 < a_vector.size( ); ++i1 ) vector[i1] = a_vector[i1];

    return( vector );
}

/* *********************************************************************************************************//**
 * Adds the items in Vector *a_from* to the set *a_to*.
 *
 * @param a_to              [in]    The list of ints to add to the set.
 * @param a_from            [in]    The set to add the ints to.
 ***********************************************************************************************************/

LUPI_HOST void addVectorItemsToSet( Vector<int> const &a_from, std::set<int> &a_to ) {

    for( Vector<int>::const_iterator iter = a_from.begin( ); iter != a_from.end( ); ++iter ) a_to.insert( *iter );
}

/* *********************************************************************************************************//**
 * This function returns a particle kinetic energy from its mass and beta (i.e., v/c) using a relativistic formula.
 *
 * @param a_mass_unitOfEnergy   [in]    The particle's mass in units of energy.
 * @param a_particleBeta        [in]    The particle's velocity divided by the speed of light (i.e., beta = v/c).
 *
 * @return                              The relativistic kinetic energy of the particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double particleKineticEnergy( double a_mass_unitOfEnergy, double a_particleBeta ) {

    if( a_particleBeta < 1e-4 ) return( 0.5 * a_mass_unitOfEnergy * a_particleBeta * a_particleBeta );

    return( a_mass_unitOfEnergy * ( 1.0 / sqrt( 1.0 - a_particleBeta * a_particleBeta ) - 1.0 ) );
}


/* *********************************************************************************************************//**
 * This function is like particleKineticEnergy except that *a_particleBeta2* is beta squared (i.e., (v/c)^2).
 *
 * @param a_mass_unitOfEnergy   [in]    The particle's mass in units of energy.
 * @param a_particleBeta2       [in]    The square of beta (i.e., beta^2 where beta = v/c).
 *
 * @return                              The relativistic kinetic energy of the particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double particleKineticEnergyFromBeta2( double a_mass_unitOfEnergy, double a_particleBeta2 ) {

    if( a_particleBeta2 < 1e-8 ) return( 0.5 * a_mass_unitOfEnergy * a_particleBeta2 );

    return( a_mass_unitOfEnergy * ( 1.0 / sqrt( 1.0 - a_particleBeta2 ) - 1.0 ) );
}

/* *********************************************************************************************************//**
 * This function returns the boost speed required to boost to the center-of-mass for a projectile hitting a target.
 *
 * @param a_massProjectile              [in]    The mass of the projectile in energy units.
 * @param a_kineticEnergyProjectile     [in]    The kinetic energy of the projectile.
 * @param a_massTarget                  [in]    The mass of the target in energy units.
 *
 * @return                              The relativistic kinetic energy of the particle.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE double boostSpeed( double a_massProjectile, double a_kineticEnergyProjectile, double a_massTarget ) {

    double betaProjectile = MCGIDI_particleBeta( a_massProjectile, a_kineticEnergyProjectile );

    return( betaProjectile / ( 1.0 + a_massTarget / ( a_massProjectile + a_kineticEnergyProjectile ) ) );
}

/* *********************************************************************************************************//**
 * This function determines the mu value(s) in the center-of-mass frame for a specified mu value in the lab frame for a
 * product with speed *a_productBeta* and a boost speed of *a_productBeta*. The returned value is the number of mu values
 * in the center-of-mass frame. The return value can be 0, 1 or 2. *a_muMinus* and *a_JacobianMinus* are undefined when the
 * returned value is less than 2. *a_muPlus* and *a_JacobianPlus* are undefined when the returned value is 0.
 *
 * @param a_muLab                       [in]    The mu specified mu value in the lab frame.
 * @param a_boostBeta                   [in]    The boost speed from the lab from to the center-of-mass frame in units of the speed-of-light.
 * @param a_productBeta                 [in]    The speed of the product in the center-of-mass frame in units of the speed-of-light.
 * @param a_muPlus                      [in]    The first mu value if the returned is greater than 0.
 * @param a_JacobianPlus                [in]    The partial derivative of mu_com with respect to mu_lab at a_muPlus.
 * @param a_muMinus                     [in]    The second mu value if the returned value is 2.
 * @param a_JacobianMinus               [in]    The partial derivative of mu_com with respect to mu_lab at a_muMinus.
 *
 * @return                                      The number of returned center-of-mass frame mu values. Can be 0, 1 or 2.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE int muCOM_From_muLab( double a_muLab, double a_boostBeta, double a_productBeta, double &a_muPlus, double &a_JacobianPlus, 
                double &a_muMinus, double &a_JacobianMinus ) {

    int numberOfSolutions = 0;
    double boostBeta2 = a_boostBeta * a_boostBeta;
    double productBeta2 = a_productBeta * a_productBeta;
    double muLab2 = a_muLab * a_muLab;
    double oneMinusMuLab2 = 1.0 - muLab2;
    double oneMinusBoostBeta2 = 1.0 - boostBeta2;
    double oneMinusBoostBeta2MuLab2 = 1.0 - boostBeta2 * muLab2;

    a_muPlus = a_muLab;                                 // Handles case when a_productBeta is 0.0 or a_muLab is +/-1.0. Actually, when a_productBeta is 0.0 it is not defined.
    a_muMinus = 0.0;

    if( ( a_productBeta == 0.0 ) || ( a_muLab == 1.0 ) ) return( 1 );

    if( a_productBeta >= a_boostBeta ) {                // Intentionally treating case where a_productBeta == a_boostBeta as one solution even though is it
        numberOfSolutions = 1; }
    else {
        if( a_muLab > 0.0 ) {                           // Only have solutions for positive mu. The next expression only test mu^2 and therefore treats negative mu like positive mu.
            if( productBeta2 * oneMinusBoostBeta2MuLab2 > boostBeta2 * oneMinusMuLab2 ) numberOfSolutions = 2;      // This ignores the case for numberOfSolutions = 1 as it probabily will never happen.
        }
    }

    if( numberOfSolutions == 0 ) return( 0 );

    double sqrt_b2minus4ac = sqrt( oneMinusBoostBeta2 * ( productBeta2 * oneMinusBoostBeta2MuLab2 - boostBeta2 * oneMinusMuLab2 ) );
    double minusbTerm = a_boostBeta * oneMinusMuLab2;
    double inv2a = 1.0 / ( a_productBeta * oneMinusBoostBeta2MuLab2 );

    a_muPlus =  ( a_muLab * sqrt_b2minus4ac - minusbTerm ) * inv2a;
    a_muMinus = ( -a_muLab * sqrt_b2minus4ac - minusbTerm ) * inv2a;      // This is meaningless when numberOfSolutions is not 1, but why add an if test.

    double JacobianTerm1 = 2.0 * boostBeta2 * a_muLab / oneMinusBoostBeta2MuLab2;
    double JacobianTerm2 = 2.0 * a_muLab * a_boostBeta / ( a_productBeta * oneMinusBoostBeta2MuLab2 );
    double JacobianTerm3 = productBeta2 * ( 1.0 - 2.0 * boostBeta2 * muLab2 ) - boostBeta2 * ( 1.0 - 2.0 * muLab2 );
    JacobianTerm3 *= oneMinusBoostBeta2 / ( a_productBeta * oneMinusBoostBeta2MuLab2 * sqrt_b2minus4ac );

    a_JacobianPlus  = fabs( a_muPlus  * JacobianTerm1 + JacobianTerm2 + JacobianTerm3 );
    a_JacobianMinus = fabs( a_muMinus * JacobianTerm1 + JacobianTerm2 - JacobianTerm3 );

    return( numberOfSolutions );
}

/* *********************************************************************************************************//**
 * This function returns a unique integer for the **Distributions::Type**. For internal use when broadcasting a
 * distribution for MPI and GPUs needs.
 *              
 * @param a_type                [in]    The distribution's type.
 *
 * @return                              Returns a unique integer for the distribution type.
 ***********************************************************************************************************/
            
LUPI_HOST_DEVICE int distributionTypeToInt( Distributions::Type a_type ) {

    int distributionType = 0;

    switch( a_type ) {
    case Distributions::Type::none :
        distributionType = 0;
        break;
    case Distributions::Type::unspecified :
        distributionType = 1;
        break;
    case Distributions::Type::angularTwoBody :
        distributionType = 2;
        break;
    case Distributions::Type::KalbachMann :
        distributionType = 3;
        break;
    case Distributions::Type::uncorrelated :
        distributionType = 4;
        break;
    case Distributions::Type::energyAngularMC :
        distributionType = 5;
        break;
    case Distributions::Type::angularEnergyMC :
        distributionType = 6;
        break;
    case Distributions::Type::coherentPhotoAtomicScattering :
        distributionType = 7;
        break;
    case Distributions::Type::incoherentPhotoAtomicScattering :
        distributionType = 8;
        break;
    case Distributions::Type::pairProductionGamma :
        distributionType = 9;
        break;
    case Distributions::Type::coherentElasticTNSL :
        distributionType = 10;
        break;
    case Distributions::Type::incoherentElasticTNSL :
        distributionType = 11;
        break;
    case Distributions::Type::incoherentPhotoAtomicScatteringElectron :
        distributionType = 12;
        break;
    case Distributions::Type::branching3d :
        distributionType = 13;
        break;
    case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering :
        distributionType = 14;
        break;
    }

    return( distributionType );
}

/* *********************************************************************************************************//**
 * This function returns the **Distributions::Type** corresponding to the integer returned by **distributionTypeToInt**.
 *
 * @param a_type                [in]    The value returned by **distributionTypeToInt**.
 *
 * @return                              The **Distributions::Type** corresponding to *a_type*.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE Distributions::Type intToDistributionType( int a_type ) {

    Distributions::Type type = Distributions::Type::none;

    switch( a_type ) {
    case 0 :
        type = Distributions::Type::none;
        break;
    case 1 :
        type = Distributions::Type::unspecified;
        break;
    case 2 :
        type = Distributions::Type::angularTwoBody;
        break;
    case 3 :
        type = Distributions::Type::KalbachMann;
        break;
    case 4 :
        type = Distributions::Type::uncorrelated;
        break;
    case 5 :
        type = Distributions::Type::energyAngularMC;
        break;
    case 6 :
        type = Distributions::Type::angularEnergyMC;
        break;
    case 7 :
        type = Distributions::Type::coherentPhotoAtomicScattering;
        break;
    case 8 :
        type = Distributions::Type::incoherentPhotoAtomicScattering;
        break;
    case 9 :
        type = Distributions::Type::pairProductionGamma;
        break;
    case 10 :
        type = Distributions::Type::coherentElasticTNSL;
        break;
    case 11 :
        type = Distributions::Type::incoherentElasticTNSL;
        break;
    case 12 :
        type = Distributions::Type::incoherentPhotoAtomicScatteringElectron;
        break;
    case 13 :
        type = Distributions::Type::branching3d;
        break;
    case 14 :
        type = Distributions::Type::incoherentBoundToFreePhotoAtomicScattering;
        break;
    default:
        LUPI_THROW( "intToDistributionType: unsupported distribution type." );
    }

    return( type );
}

/* *********************************************************************************************************//**
 * This method serializes *a_products* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *a_products* or unpack *a_products* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 * @param a_products            [in]    The products to serialize.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void serializeProducts( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<Product *> &a_products ) {

    std::size_t vectorSize = a_products.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_products.resize( vectorSize, &a_buffer.m_placement );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( a_buffer.m_placement != nullptr ) {
                a_products[vectorIndex] = new(a_buffer.m_placement) Product;
                a_buffer.incrementPlacement( sizeof( Product ) );
            }
            else {
               a_products[vectorIndex] = new Product;
            }
        } }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += a_products.internalSize( );
        a_buffer.incrementPlacement( sizeof( Product ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        a_products[vectorIndex]->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * This method serializes *a_delayedNeutrons* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *a_delayedNeutrons* or unpack *a_delayedNeutrons* depending on *a_mode*.
 *
 * @param a_delayedNeutrons     [in]    The delayed neutrons to serialize.
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/
        
LUPI_HOST_DEVICE void serializeDelayedNeutrons( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<DelayedNeutron *> &a_delayedNeutrons ) {

    std::size_t vectorSize = a_delayedNeutrons.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_delayedNeutrons.resize( vectorSize, &a_buffer.m_placement );
        for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
            if( a_buffer.m_placement != nullptr ) {
                a_delayedNeutrons[vectorIndex] = new(a_buffer.m_placement) DelayedNeutron;
                a_buffer.incrementPlacement( sizeof( DelayedNeutron ) );
            }
            else {
                a_delayedNeutrons[vectorIndex] = new DelayedNeutron;
            }
        } }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += a_delayedNeutrons.internalSize( );
        a_buffer.incrementPlacement( sizeof( DelayedNeutron ) * vectorSize );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        a_delayedNeutrons[vectorIndex]->serialize( a_buffer, a_mode );
    }
}

/* *********************************************************************************************************//**
 * This method serializes *a_Qs* for broadcasting as needed for MPI and GPUs. The method can count the number of required
 * bytes, pack *a_Qs* or unpack *a_Qs* depending on *a_mode*.
 *
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 * @param a_Qs                  [in]    The Q functions to serialize.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void serializeQs( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, Vector<Functions::Function1d_d1 *> &a_Qs ) {

    std::size_t vectorSize = a_Qs.size( );
    int vectorSizeInt = (int) vectorSize;
    DATA_MEMBER_INT( vectorSizeInt, a_buffer, a_mode );
    vectorSize = (std::size_t) vectorSizeInt;

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        a_Qs.resize( vectorSize, &a_buffer.m_placement ); }
    else if( a_mode == LUPI::DataBuffer::Mode::Memory ) {
        a_buffer.m_placement += a_Qs.internalSize( );
    }

    for( std::size_t vectorIndex = 0; vectorIndex < vectorSize; ++vectorIndex ) {
        a_Qs[vectorIndex] = serializeFunction1d_d1( a_buffer, a_mode, a_Qs[vectorIndex] );
    }
}


/* *********************************************************************************************************//**
 * 
 * @param a_fissionResiduals    [in]    A reference to the GIDI::Construction::FissionResiduals reference serialize.
 * @param a_buffer              [in]    The buffer to read or write data to depending on *a_mode*.
 * @param a_mode                [in]    Specifies the action of this method.
 ***********************************************************************************************************/

LUPI_HOST_DEVICE void serializeFissionResiduals( GIDI::Construction::FissionResiduals &a_fissionResiduals, 
                LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode ) {

    int fissionResidualsInt = 0;

    if( a_mode != LUPI::DataBuffer::Mode::Unpack ) {
        switch( a_fissionResiduals ) {
        case GIDI::Construction::FissionResiduals::none :
            break;
        case GIDI::Construction::FissionResiduals::ENDL99120 :
            fissionResidualsInt = 1;
            break;
        case GIDI::Construction::FissionResiduals::ENDL99125 :
            fissionResidualsInt = 2;
            break;
        }
    }

    DATA_MEMBER_INT( fissionResidualsInt, a_buffer, a_mode );

    if( a_mode == LUPI::DataBuffer::Mode::Unpack ) {
        switch( fissionResidualsInt ) {
        case 0 :
            a_fissionResiduals = GIDI::Construction::FissionResiduals::none;
            break;
        case 1 :
            a_fissionResiduals = GIDI::Construction::FissionResiduals::ENDL99120;
            break;
        case 2 :
            a_fissionResiduals = GIDI::Construction::FissionResiduals::ENDL99125;
            break;
        }
    }
}

/* *********************************************************************************************************//**
 * This function returns a std::vector<double> that represents **a_input**.
 *
 * @param a_input               [in]    The input for the returned std::vector<double> instance.
 *
 * @returns                             A std::vector<double> instance.
 ***********************************************************************************************************/

LUPI_HOST std::vector<double> vectorToSTD_vector( Vector<double> a_input ) {

    std::vector<double> vector( a_input.size( ) );

    std::size_t index = 0;
    for( auto iter = a_input.begin( ); iter != a_input.end( ); ++iter, ++index ) vector[index] = *iter;

    return( vector );
}

/* *********************************************************************************************************//**
 * This function returns a std::vector<double> that represents **a_input**.
 *
 * @param a_input               [in]    The input for the returned std::vector<double> instance.
 *
 * @returns                             A std::vector<double> instance.
 ***********************************************************************************************************/

LUPI_HOST std::vector<double> vectorToSTD_vector( Vector<float> a_input ) {

    std::vector<double> vector( a_input.size( ) );

    std::size_t index = 0;
    for( auto iter = a_input.begin( ); iter != a_input.end( ); ++iter, ++index ) vector[index] = *iter;

    return( vector );
}

}
