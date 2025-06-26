/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_headerSource_hpp_included
#define MCGIDI_headerSource_hpp_included 1

#include <climits>


#ifdef MCGIDI_USE_DOUBLES
    #define crossSectionSumError 1e-8
#else
    #define crossSectionSumError 1e-6
#endif

// From file: MCGIDI_URR.cpp

/* *********************************************************************************************************//**
 * Updates *this* if *a_protare* has a non-negative *URR_index*.
 *
 * @param a_protare             [in]    The protare whose *URR_index* is used to see if *this* needs updating.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::URR_protareInfos::updateProtare( MCGIDI::Protare const *a_protare, double a_energy, RNG && a_rng ) {

    for( std::size_t i1 = 0; i1 < a_protare->numberOfProtares( ); ++i1 ) {
        ProtareSingle *protareSingle = const_cast<ProtareSingle *>( a_protare->protare( i1 ) );

        if( protareSingle->URR_index( ) >= 0 ) {
            URR_protareInfo &URR_protare_info = m_URR_protareInfos[protareSingle->URR_index( )];

            URR_protare_info.m_inURR = protareSingle->inURR( a_energy );
            if( URR_protare_info.inURR( ) ) URR_protare_info.m_rng_Value = a_rng( );
        }
    }
}


/* *********************************************************************************************************//**
 * This function samples an energy and cosine of the angle for a photon for Klein Nishina scattering (i.e, incoherent photo-atomic scattering).
 *
 * @param a_energyIn            [in]    The energy of the incoming photon.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energyOut           [in]    The energy of the scattered photon.
 * @param a_mu                  [in]    The cosine of the angle of the scattered photon's z-axis and the incoming photon's z-axis.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI_sampleKleinNishina( double a_energyIn, RNG && a_rng, double *a_energyOut, double *a_mu ) {
/*
  Description
    Sample the Klein-Nishina distribution.
      The unit of energy is the rest mass of the electron.
      Reference: R. N. Blomquist and E. N. Gelbard, Nuclear Science
      and Engineering, 83, 380-384 (1983)

   This routine was taken from MCAPM which was from MCNP with only cosmetic changes.

   Input
     a_energyIn     - incident photon energy ( in electron rest mass units )
     *rng           - user supplied random number generator
   Output
     *a_energyOut   - exiting photon energy ( in electron rest mass units )
     *a_mu          - exiting photon cosine
*/

    double a1, b1, t1, s1, r1, mu, energyOut;

    a1 = 1.0 / a_energyIn;
    b1 = 1.0 / ( 1.0 + 2.0 * a_energyIn );

    if( a_energyIn < 3.0 ) {                      // Kahn''s method ( e < 1.5 MeV ) AECU-3259.
        bool reject = true;

        t1 = 1.0 / ( 1.0 + 8.0 * b1 );
        do {
            if( a_rng( ) <= t1 ) {
                r1 = 2.0 * a_rng( );
                s1 = 1.0 / ( 1.0 + a_energyIn * r1 );
                mu = 1.0 - r1;
                reject = a_rng( ) > 4.0 * s1 * ( 1.0 - s1 ); }
            else {
                s1 = ( 1.0 + 2.0 * a_energyIn * a_rng( ) ) * b1;
                mu = 1.0 + a1 * ( 1.0 - 1.0 / s1 );
                reject = a_rng( ) > 0.5 * ( mu * mu + s1 );
            }
        } while( reject );
        energyOut = a_energyIn / ( 1 + a_energyIn * ( 1 - mu ) ); }
    else {                                        // Koblinger''s method ( e > 1.5 MeV ) NSE 56, 218 ( 1975 ).
        t1 = a_rng( ) * ( 4.0 * a1 + 0.5 * ( 1.0 - b1 * b1 ) - ( 1.0 - 2.0 * ( 1.0 + a_energyIn ) * ( a1 * a1 ) ) * log( b1 ) );
        if( t1 > 2.0 * a1 ) {
            if( t1 > 4.0 * a1 ) {
                if( t1 > 4.0 * a1 + 0.5 * ( 1.0 - b1 * b1 ) ) {
                    energyOut = a_energyIn * pow( b1, a_rng( ) );
                    mu = 1.0 + a1 - 1.0 / energyOut; }
                else {
                    energyOut = a_energyIn * sqrt( 1.0 - a_rng( ) * ( 1.0 - b1 * b1 ) );
                    mu = 1.0 + a1 - 1.0 / energyOut;
                  } }
            else {
                energyOut = a_energyIn * ( 1.0 + a_rng( ) * ( b1 - 1.0 ) );
                mu =  1.0 + a1 - 1.0 / energyOut; } }
        else {
            r1 = 2.0 * a_rng( );
            mu = 1.0 - r1;
            energyOut = 1.0 / ( a1 + r1 );
          }
    }

    *a_mu = mu;
    *a_energyOut = energyOut;

    return;
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data for the two outgoing particles in a two-body outgoing channel.
 * First, is samples *mu*, the cosine of the product's outgoing angle, since this is for two-body interactions, *mu*
 * is in the center-of-mass frame. It then calls kinetics_COMKineticEnergy2LabEnergyAndMomentum.
 *
 * @param a_X                       [in]    The energy of the projectile in the lab frame.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/
template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::Distribution::sample( double a_X, MCGIDI::Sampling::Input &a_input, RNG && a_rng ) const {

    switch( type( ) ) {
    case Distributions::Type::none:
        break;
    case Distributions::Type::unspecified:
        static_cast<Distributions::Unspecified const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::angularTwoBody:
        static_cast<Distributions::AngularTwoBody const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::KalbachMann:
        static_cast<Distributions::KalbachMann const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::uncorrelated:
        static_cast<Distributions::Uncorrelated const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::branching3d:
        static_cast<Distributions::Branching3d const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::energyAngularMC:
        static_cast<Distributions::EnergyAngularMC const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::angularEnergyMC:
        static_cast<Distributions::AngularEnergyMC const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::coherentPhotoAtomicScattering:
        static_cast<Distributions::CoherentPhotoAtomicScattering const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::incoherentPhotoAtomicScattering:
        static_cast<Distributions::IncoherentPhotoAtomicScattering const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering:
        static_cast<Distributions::IncoherentBoundToFreePhotoAtomicScattering const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::incoherentPhotoAtomicScatteringElectron:
        static_cast<Distributions::IncoherentPhotoAtomicScatteringElectron const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::pairProductionGamma:
        static_cast<Distributions::PairProductionGamma const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::coherentElasticTNSL:
        static_cast<Distributions::CoherentElasticTNSL const *>( this )->sample( a_X, a_input, a_rng );
        break;
    case Distributions::Type::incoherentElasticTNSL:
        static_cast<Distributions::IncoherentElasticTNSL const *>( this )->sample( a_X, a_input, a_rng );
        break;
    }
}

// From file: MCGIDI_distributions.cpp

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [out]   The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::Distribution::angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab,
                RNG && a_rng, double &a_energy_out ) const {

    double probability = 0.0;
    a_energy_out = 0.0;

    switch( type( ) ) {
    case Distributions::Type::none:
        break;
    case Distributions::Type::unspecified:
        probability = static_cast<Distributions::Unspecified const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::angularTwoBody:
        probability = static_cast<Distributions::AngularTwoBody const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::KalbachMann:
        probability = static_cast<Distributions::KalbachMann const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::uncorrelated:
        probability = static_cast<Distributions::Uncorrelated const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::branching3d:
        probability = static_cast<Distributions::Branching3d const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::energyAngularMC:
        probability = static_cast<Distributions::EnergyAngularMC const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::angularEnergyMC:
        probability = static_cast<Distributions::AngularEnergyMC const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::coherentPhotoAtomicScattering:
        probability = static_cast<Distributions::CoherentPhotoAtomicScattering const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::incoherentPhotoAtomicScattering:
        probability = static_cast<Distributions::IncoherentPhotoAtomicScattering const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::incoherentBoundToFreePhotoAtomicScattering:
        probability = static_cast<Distributions::IncoherentBoundToFreePhotoAtomicScattering const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::incoherentPhotoAtomicScatteringElectron:
        probability = static_cast<Distributions::IncoherentPhotoAtomicScatteringElectron const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::pairProductionGamma:
        probability = static_cast<Distributions::PairProductionGamma const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::coherentElasticTNSL:
        probability = static_cast<Distributions::CoherentElasticTNSL const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    case Distributions::Type::incoherentElasticTNSL:
        probability = static_cast<Distributions::IncoherentElasticTNSL const *>( this )->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab,
                a_rng, a_energy_out );
        break;
    }

    return( probability );
}

/* *********************************************************************************************************//**
 * This function calculates the products outgoing data (i.e., energy, velocity/momentum) for the two products of
 * a two-body interaction give the cosine of the first product's outgoing angle.
 *
 * @param a_beta                    [in]    The velocity/speedOflight of the com frame relative to the lab frame.
 * @param a_kinetic_com             [in]    Total kinetic energy (K1 + K2) in the COM frame.
 * @param a_m3cc                    [in]    The mass of the first product.
 * @param a_m4cc                    [in]    The mass of the second product.
 * @param a_input                   [in]    Sample options requested by user and where the products' outgoing data are returned.
 ***********************************************************************************************************/

inline LUPI_HOST_DEVICE void kinetics_COMKineticEnergy2LabEnergyAndMomentum( double a_beta, double a_kinetic_com,
        double a_m3cc, double a_m4cc, MCGIDI::Sampling::Input &a_input ) {
/*
    Relativity:
        E = K + m, E^2 = K^2 + 2 K m + m^2, E^2 - m^2 = p^2 = K^2 + 2 K m

         pc          p     v
        ---- = v,   --- = --- = beta = b
         E           E     c

           K ( K + 2 m )
    b^2 = ---------------
            ( K + m )^2
*/
    double x, v_p, p, pp3, pp4, px3, py3, pz3, pz4, pz, p_perp2, E3, E4, gamma, m3cc2 = a_m3cc * a_m3cc, m4cc2 = a_m4cc * a_m4cc;

    p = sqrt( a_kinetic_com * ( a_kinetic_com + 2. * a_m3cc ) * ( a_kinetic_com + 2. * a_m4cc )  *
            ( a_kinetic_com + 2. * ( a_m3cc + a_m4cc ) ) ) / ( 2. * ( a_kinetic_com + a_m3cc + a_m4cc ) );

    py3 = p * sqrt( 1 - a_input.m_mu * a_input.m_mu );
    px3 = py3 * cos( a_input.m_phi );
    py3 *= sin( a_input.m_phi );
    pz = p * a_input.m_mu;
    if( 1 ) {                           // FIXME Assuming the answer is wanted in the lab frame for now.
        a_input.m_frame = GIDI::Frame::lab;
        E3 = sqrt( p * p + m3cc2 );
        E4 = sqrt( p * p + m4cc2 );
        gamma = sqrt( 1. / ( 1. - a_beta * a_beta ) );
        pz3 = gamma * (  pz + a_beta * E3 );
        pz4 = gamma * ( -pz + a_beta * E4 ); }
    else {                              // COM frame.
        a_input.m_frame = GIDI::Frame::centerOfMass;
        pz3 = pz;
        pz4 = -pz;
    }

    p_perp2 = px3 * px3 + py3 * py3;

    a_input.m_px_vx1 = px3;
    a_input.m_py_vy1 = py3;
    a_input.m_pz_vz1 = pz3;
    pp3 = p_perp2 + pz3 * pz3;
    x = ( a_m3cc > 0 ) ? pp3 / ( 2 * m3cc2 ) : 1.;
    if( x < 1e-5 ) {
        a_input.m_energyOut1 = a_m3cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        a_input.m_energyOut1 = sqrt( m3cc2 + pp3 ) - a_m3cc;
    }
    a_input.m_px_vx2 = -px3;
    a_input.m_py_vy2 = -py3;
    a_input.m_pz_vz2 = pz4;
    pp4 = p_perp2 + pz4 * pz4;
    x = ( a_m4cc > 0 ) ? pp4 / ( 2 * m4cc2 ) : 1.;
    if( x < 1e-5 ) {
        a_input.m_energyOut2 = a_m4cc * x  * ( 1 - 0.5 * x * ( 1 - x ) ); }
    else {
        a_input.m_energyOut2 = sqrt( m4cc2 + pp4 ) - a_m4cc;
    }

    if( a_input.wantVelocity( ) ) {
        v_p = MCGIDI_speedOfLight_cm_sec / sqrt( pp3 + m3cc2 );
        a_input.m_px_vx1 *= v_p;
        a_input.m_py_vy1 *= v_p;
        a_input.m_pz_vz1 *= v_p;

        v_p = MCGIDI_speedOfLight_cm_sec / sqrt( pp4 + m4cc2 );
        a_input.m_px_vx2 *= v_p;
        a_input.m_py_vy2 *= v_p;
        a_input.m_pz_vz2 *= v_p;
    }
}


/* *********************************************************************************************************//**
 * This method samples the outgoing product data for the two outgoing particles in a two-body outgoing channel.
 * First, is samples *mu*, the cosine of the product's outgoing angle, since this is for two-body interactions, *mu*
 * is in the center-of-mass frame. It then calls kinetics_COMKineticEnergy2LabEnergyAndMomentum.
 *
 * @param a_X                       [in]    The energy of the projectile in the lab frame.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/
template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::AngularTwoBody::sample( double a_X, MCGIDI::Sampling::Input &a_input, RNG && a_rng ) const {

    double initialMass = projectileMass( ) + targetMass( ), finalMass = productMass( ) + m_residualMass;
    double beta = sqrt( a_X * ( a_X + 2. * projectileMass( ) ) ) / ( a_X + initialMass );      // beta = v/c.
    double _x = targetMass( ) * ( a_X - m_twoBodyThreshold ) / ( finalMass * finalMass );
    double Kp;                          // Kp is the total kinetic energy for m3 and m4 in the COM frame.

    a_input.m_sampledType = Sampling::SampledType::firstTwoBody;

    if( m_Upscatter ) {
       if( ( a_input.m_upscatterModel == Sampling::Upscatter::Model::B ) || ( a_input.m_upscatterModel == Sampling::Upscatter::Model::BSnLimits )
                || ( a_input.m_upscatterModel == Sampling::Upscatter::Model::DBRC ) ) {
            if( upscatterModelB( a_X, a_input, a_rng ) ) return;
        }
    }

    if( _x < 2e-5 ) {
        Kp = finalMass * _x * ( 1 - 0.5 * _x * ( 1 - _x ) ); }
    else {          // This is the relativistic formula derived from E^2 - (pc)^2 is frame independent.
        Kp = sqrt( finalMass * finalMass + 2 * targetMass( ) * ( a_X - m_twoBodyThreshold ) ) - finalMass;
    }
    if( Kp < 0 ) Kp = 0.;           // FIXME There needs to be a better test here.

    a_input.m_mu = m_angular->sample( a_X, a_rng( ), a_rng );
    a_input.m_phi = 2. * M_PI * a_rng( );
    kinetics_COMKineticEnergy2LabEnergyAndMomentum( beta, Kp, productMass( ), m_residualMass, a_input );
}


/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [out]   The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::AngularTwoBody::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
                double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = 0.0;

    double initialMass = projectileMass( ) + targetMass( ), finalMass = productMass( ) + m_residualMass;
    double _x = targetMass( ) * ( a_energy_in - m_twoBodyThreshold ) / ( finalMass * finalMass );
    double Kp;                      // Total kinetic energy of products in the center-of-mass.

    if( _x < 2e-5 ) {
        Kp = finalMass * _x * ( 1 - 0.5 * _x * ( 1 - _x ) ); }
    else {          // This is the relativistic formula derived from E^2 - (pc)^2 which is frame independent (i.e., an invariant).
        Kp = sqrt( finalMass * finalMass + 2.0 * targetMass( ) * ( a_energy_in - m_twoBodyThreshold ) ) - finalMass;
    }
    if( Kp < 0 ) Kp = 0.;           // FIXME There needs to be a better test here.

    double energy_product_com = 0.5 * Kp * ( Kp + 2.0 * m_residualMass ) / ( Kp + productMass( ) + m_residualMass );

    if( productMass( ) == 0.0 ) {
        double boostBeta = sqrt( a_energy_in * ( a_energy_in + 2. * projectileMass( ) ) ) / ( a_energy_in + initialMass );  // Good, even for projectileMass = 0.
        double one_mu_beta = 1.0 - a_mu_lab * boostBeta;
        double mu_com = ( a_mu_lab - boostBeta ) / one_mu_beta;
        double Jacobian = ( 1.0 - boostBeta * boostBeta ) / ( one_mu_beta * one_mu_beta );

        a_energy_out = sqrt( 1.0 - boostBeta * boostBeta ) * energy_product_com * ( 1.0 + mu_com * boostBeta );

        return( Jacobian * m_angular->evaluate( a_energy_in, mu_com ) );
    }

    double productBeta = MCGIDI_particleBeta( productMass( ), energy_product_com );
    double boostBeta = sqrt( a_energy_in * ( a_energy_in + 2. * projectileMass( ) ) ) / ( a_energy_in + initialMass );      // beta = v/c.
    double muPlus = 0.0, JacobianPlus = 0.0, muMinus = 0.0, JacobianMinus = 0.0;

    int numberOfMus = muCOM_From_muLab( a_mu_lab, boostBeta, productBeta, muPlus, JacobianPlus, muMinus, JacobianMinus );

    if( numberOfMus == 0 ) return( 0.0 );

    double probability = JacobianPlus * m_angular->evaluate( a_energy_in, muPlus );

    if( numberOfMus == 2 ) {
        double probabilityMinus = JacobianMinus * m_angular->evaluate( a_energy_in, muMinus );
        probability += probabilityMinus;
        if( probabilityMinus > a_rng( ) * probability ) {
            muPlus = muMinus;
        }
    }

    double productBeta2 = productBeta * productBeta;
    double productBetaLab2 = productBeta2 + boostBeta * boostBeta * ( 1.0 - productBeta2 * ( 1.0 - muPlus * muPlus ) ) + 2.0 * muPlus * productBeta * boostBeta;
    productBetaLab2 /= 1.0 - muPlus * productBeta * boostBeta;
    a_energy_out = MCGIDI::particleKineticEnergyFromBeta2( productMass( ), productBetaLab2 );

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples a targets velocity for elastic upscattering for upscatter model B and then calculates the outgoing
 * product data for the projectile and target.
 *
 * @param a_kineticLab              [in]    The kinetic energy of the projectile in the lab frame.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE bool MCGIDI::Distributions::AngularTwoBody::upscatterModelB( double a_kineticLab, Sampling::Input &a_input, RNG && a_rng ) const {

    const double Two_sqrtPi = 1.1283791670955125739;
    const double C0 = 1.0410423479, C1 = 3.9626339162e-4, C2 =-1.8654539193e-3, C3 = 1.0264818153e-4;
    double neutronMass = projectileMass( );                             // Mass are in incident energy unit / c**2.
    double _targetMass = targetMass( );
    double temperature = 1e-3 * a_input.m_temperature;                  // Assumes m_temperature is in keV/K.
    double kineticLabMax = 1e4 * temperature;

    if( a_input.m_upscatterModel == Sampling::Upscatter::Model::BSnLimits ) {
        double kineticLabMax200 = 200.0 * temperature;

        kineticLabMax = 1e3 * temperature * neutronMass / _targetMass;
        if( kineticLabMax < kineticLabMax200 ) kineticLabMax = kineticLabMax200;
        if( a_kineticLab >= 0.1 ) kineticLabMax = 0.9 * a_kineticLab; }
    else {
        if( kineticLabMax > 1e-2 ) {
            kineticLabMax = 1e-2;
            if( kineticLabMax < 100.0 * temperature ) {
                kineticLabMax = 100.0 * temperature;
                if( kineticLabMax > 10.0 ) kineticLabMax = 10.0;        // Assumes energy is in MeV.
            }
        }
    }

    if( a_kineticLab > kineticLabMax ) return( false );                   // Only for low neutron energy.

    a_input.m_frame = GIDI::Frame::lab;

    double muProjectileTarget, relativeBeta, targetBeta;
    double targetThermalBeta = MCGIDI_particleBeta( _targetMass, temperature );
    double neutronBeta = MCGIDI_particleBeta( neutronMass, a_kineticLab );

    a_input.m_numberOfDBRC_rejections = 0;
    bool continueLoop = false;
    double crossSectionMax = 0.0;
    if( a_input.m_upscatterModel == Sampling::Upscatter::Model::DBRC ) {
        double targetThermalSpeed = m_modelDBRC_data->targetThermalSpeed( temperature );          // Non-relativistic calculations, unlike targetThermalBeta.
        crossSectionMax = m_modelDBRC_data->crossSectionMax( a_kineticLab, targetThermalSpeed );
    }
    do {
        continueLoop = false;
        ++a_input.m_numberOfDBRC_rejections;
        do {
            int MethodP1orP2 = 0;         /* Assume P2 */
            if( a_rng( ) * ( neutronBeta + Two_sqrtPi * targetThermalBeta ) < neutronBeta ) MethodP1orP2 = 1;
            muProjectileTarget = 1.0 - 2.0 * a_rng( );
            if( MethodP1orP2 == 0 ) {                                       // x Exp( -x ) term.
                targetBeta = targetThermalBeta * sqrt( -log( ( 1.0 - a_rng( ) ) * ( 1.0 - a_rng( ) ) ) ); }
            else {                                                          // x^2 Exp( -x^2 ) term.
                double x1;
                do {
                    x1 = a_rng( );
                    x1 = sqrt( -log( ( 1.0 - a_rng( ) ) * ( 1.0 - x1 * x1 ) ) );
                    x1 = x1 / ( ( ( C3 * x1 + C2 ) * x1 + C1 ) * x1 + C0 );
                } while( x1 > 4.0 );
                targetBeta = targetThermalBeta * x1;
            }
            relativeBeta = sqrt( targetBeta * targetBeta + neutronBeta * neutronBeta - 2 * muProjectileTarget * targetBeta * neutronBeta );
        } while( relativeBeta < ( targetBeta + neutronBeta ) * a_rng( ) );
       if( a_input.m_upscatterModel == Sampling::Upscatter::Model::DBRC ) {
            double relativeNeutronEnergy = 0.5 * productMass( ) * relativeBeta * relativeBeta;
            if( m_modelDBRC_data->evaluate( relativeNeutronEnergy ) < a_rng( ) * crossSectionMax ) continueLoop = true;
        }
    } while( continueLoop );

    double m1_12 = neutronMass / ( neutronMass + _targetMass );
    double m2_12 = _targetMass / ( neutronMass + _targetMass );

    double cosRelative = 0.0;                                           // Cosine of angle between projectile velocity and relative velocity.
    if( relativeBeta != 0.0 ) cosRelative = ( neutronBeta - muProjectileTarget * targetBeta ) / relativeBeta;
    if( cosRelative > 1.0 ) {
            cosRelative = 1.0; }
    else if( cosRelative < -1.0 ) {
            cosRelative = -1.0;
    }
    double sinRelative = sqrt( 1.0 - cosRelative * cosRelative );       // Sine of angle between projectile velocity and relative velocity.

    a_input.m_relativeMu = cosRelative;
    a_input.m_targetBeta = targetBeta;
    a_input.m_relativeBeta = relativeBeta;

    double betaNeutronOut = m2_12 * relativeBeta;
    double kineticEnergyRelative = particleKineticEnergy( neutronMass, betaNeutronOut );
    double muCOM = m_angular->sample( kineticEnergyRelative, a_rng( ), a_rng );
    double phiCOM = 2.0 * M_PI * a_rng( );
    double SCcom = sqrt( 1.0 - muCOM * muCOM );
    double SScom = SCcom * sin( phiCOM );
    SCcom *= cos( phiCOM );

    a_input.m_pz_vz1 = betaNeutronOut * ( muCOM * cosRelative - SCcom * sinRelative );
    a_input.m_px_vx1 = betaNeutronOut * ( muCOM * sinRelative + SCcom * cosRelative );
    a_input.m_py_vy1 = betaNeutronOut * SScom;

    double massRatio = -neutronMass / _targetMass;
    a_input.m_pz_vz2 = massRatio * a_input.m_pz_vz1;
    a_input.m_px_vx2 = massRatio * a_input.m_px_vx1;
    a_input.m_py_vy2 = massRatio * a_input.m_py_vy1;

    double vCOMz = m1_12 * neutronBeta + m2_12 * muProjectileTarget * targetBeta;                   // Boost from center-of-mass to lab frame.
    double vCOMx = m2_12 * sqrt( 1.0 - muProjectileTarget * muProjectileTarget ) * targetBeta;
    a_input.m_pz_vz1 += vCOMz;
    a_input.m_px_vx1 += vCOMx;
    a_input.m_pz_vz2 += vCOMz;
    a_input.m_px_vx2 += vCOMx;

    double vx2_vy2 = a_input.m_px_vx1 * a_input.m_px_vx1 + a_input.m_py_vy1 * a_input.m_py_vy1;
    double v2 = a_input.m_pz_vz1 * a_input.m_pz_vz1 + vx2_vy2;
    a_input.m_mu = 0.0;
    if( v2 != 0.0 ) a_input.m_mu = a_input.m_pz_vz1 / sqrt( v2 );
    a_input.m_phi = atan2( a_input.m_py_vy1, a_input.m_px_vx1 );

    a_input.m_energyOut1 = MCGIDI::particleKineticEnergyFromBeta2( neutronMass, v2 );
    a_input.m_energyOut2 = MCGIDI::particleKineticEnergyFromBeta2( _targetMass, a_input.m_px_vx2 * a_input.m_px_vx2 + a_input.m_py_vy2 * a_input.m_py_vy2 + a_input.m_pz_vz2 * a_input.m_pz_vz2 );

    a_input.m_px_vx1 *= MCGIDI_speedOfLight_cm_sec;
    a_input.m_py_vy1 *= MCGIDI_speedOfLight_cm_sec;
    a_input.m_pz_vz1 *= MCGIDI_speedOfLight_cm_sec;

    a_input.m_px_vx2 *= MCGIDI_speedOfLight_cm_sec;
    a_input.m_py_vy2 *= MCGIDI_speedOfLight_cm_sec;
    a_input.m_pz_vz2 *= MCGIDI_speedOfLight_cm_sec;

    if( !a_input.wantVelocity( ) ) {                // Return momenta.
        a_input.m_px_vx1 *= neutronMass;            // Non-relativistic.
        a_input.m_py_vy1 *= neutronMass;
        a_input.m_pz_vz1 *= neutronMass;

        a_input.m_px_vx2 *= _targetMass;
        a_input.m_py_vy2 *= _targetMass;
        a_input.m_pz_vz2 *= _targetMass;
    }

    double phi = 2.0 * M_PI * a_rng( );
    double sine = sin( phi );
    double cosine = cos( phi );

    double saved = a_input.m_px_vx1;
    a_input.m_px_vx1 = cosine * a_input.m_px_vx1 - sine   * a_input.m_py_vy1;
    a_input.m_py_vy1 = sine   * saved            + cosine * a_input.m_py_vy1;

    return( true );
}


/* *********************************************************************************************************//**
 * This method samples the outgoing product data by sampling the outgoing energy E' and mu from the uncorrelated
 * E and mu probabilities. It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::Uncorrelated::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_mu = m_angular->sample( a_X, a_rng( ), a_rng );
    a_input.m_energyOut1 = m_energy->sample( a_X, a_rng( ), a_rng );
    a_input.m_phi = 2. * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::Uncorrelated::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
                double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    if( productFrame( ) != GIDI::Frame::lab ) {
        a_energy_out = 0.0;

        double initialMass = projectileMass( ) + targetMass( );
        double boostBeta = sqrt( a_energy_in * ( a_energy_in + 2. * projectileMass( ) ) ) / ( a_energy_in + initialMass );  // Good, even for projectileMass = 0.
        double energy_out_com = m_energy->sample( a_energy_in, a_rng( ), a_rng );

        if( productMass( ) == 0.0 ) {
            double one_mu_beta = 1.0 - a_mu_lab * boostBeta;
            double mu_com = ( a_mu_lab - boostBeta ) / one_mu_beta;
            double Jacobian = ( 1.0 - boostBeta * boostBeta ) / ( one_mu_beta * one_mu_beta );

            a_energy_out = sqrt( 1.0 - boostBeta * boostBeta ) * energy_out_com * ( 1.0 + mu_com * boostBeta );

            return( Jacobian * m_angular->evaluate( a_energy_in, mu_com ) );
        }

        double productBeta = MCGIDI_particleBeta( productMass( ), energy_out_com );
        double muPlus = 0.0, JacobianPlus = 0.0, muMinus = 0.0, JacobianMinus = 0.0;
        int numberOfMus = muCOM_From_muLab( a_mu_lab, boostBeta, productBeta, muPlus, JacobianPlus, muMinus, JacobianMinus );

        if( numberOfMus == 0 ) return( 0.0 );

        double probability = JacobianPlus * m_angular->evaluate( a_energy_in, muPlus );

        if( numberOfMus == 2 ) {
            double probabilityMinus = JacobianMinus * m_angular->evaluate( a_energy_in, muMinus );

            probability += probabilityMinus;
            if( probabilityMinus > a_rng( ) * probability ) muPlus = muMinus;
        }

        double productBeta2 = productBeta * productBeta;
        double productBetaLab2 = productBeta2 + boostBeta * boostBeta * ( 1.0 - productBeta2 * ( 1.0 - muPlus * muPlus ) ) + 2.0 * muPlus * productBeta * boostBeta;
        productBetaLab2 /= 1.0 - muPlus * productBeta * boostBeta;
        a_energy_out = MCGIDI::particleKineticEnergyFromBeta2( productMass( ), productBetaLab2 );

        return( probability );
    }

    a_energy_out = m_energy->sample( a_energy_in, a_rng( ), a_rng );
    return( m_angular->evaluate( a_energy_in, a_mu_lab ) );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing branching photons.
 *
 * @param a_X                       [in]    The energy of the projectile in the lab frame.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::Branching3d::sample( LUPI_maybeUnused double a_X, LUPI_maybeUnused Sampling::Input &a_input, LUPI_maybeUnused RNG && a_rng ) const {

}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::Branching3d::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		LUPI_maybeUnused double a_energy_in, LUPI_maybeUnused double a_mu_lab, LUPI_maybeUnused RNG && a_rng, LUPI_maybeUnused double &a_energy_out ) const {

    double probability = 0.0;

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data by sampling the outgoing energy E' from the probability P(E'|E) and then samples mu from
 * the probability P(mu|E,E'). It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::EnergyAngularMC::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    double energyOut_1, energyOut_2;

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = m_energy->sample2dOf3d( a_X, a_rng( ), a_rng, &energyOut_1, &energyOut_2 );
    a_input.m_mu = m_angularGivenEnergy->sample( a_X, energyOut_1, energyOut_2, a_rng( ), a_rng );
    a_input.m_phi = 2. * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::EnergyAngularMC::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    double probability = 0.0;

    if( productFrame( ) == GIDI::Frame::centerOfMass ) {
        a_energy_out = m_energy->sample( a_energy_in, a_rng( ), a_rng );

        double initialMass = projectileMass( ) + targetMass( );
        double boostBeta = sqrt( a_energy_in * ( a_energy_in + 2. * projectileMass( ) ) ) / ( a_energy_in + initialMass );  // Good, even for projectileMass = 0.
        double energy_out_com = m_energy->sample( a_energy_in, a_rng( ), a_rng );

        if( productMass( ) == 0.0 ) {
            double one_mu_beta = 1.0 - a_mu_lab * boostBeta;
            double mu_com = ( a_mu_lab - boostBeta ) / one_mu_beta;
            double Jacobian = ( 1.0 - boostBeta * boostBeta ) / ( one_mu_beta * one_mu_beta );

            a_energy_out = sqrt( 1.0 - boostBeta * boostBeta ) * energy_out_com * ( 1.0 + mu_com * boostBeta );

            return( Jacobian * m_angularGivenEnergy->evaluate( a_energy_in, energy_out_com, mu_com ) );
        }

        double productBeta = MCGIDI_particleBeta( productMass( ), energy_out_com );
        double muPlus = 0.0, JacobianPlus = 0.0, muMinus = 0.0, JacobianMinus = 0.0;
        int numberOfMus = muCOM_From_muLab( a_mu_lab, boostBeta, productBeta, muPlus, JacobianPlus, muMinus, JacobianMinus );

        if( numberOfMus == 0 ) return( 0.0 );

        probability = JacobianPlus * m_angularGivenEnergy->evaluate( a_energy_in, energy_out_com, muPlus );

        if( numberOfMus == 2 ) {
            double probabilityMinus = JacobianMinus * m_angularGivenEnergy->evaluate( a_energy_in, energy_out_com, muMinus );

            probability += probabilityMinus;
            if( probabilityMinus > a_rng( ) * probability ) muPlus = muMinus;
        }

        double productBeta2 = productBeta * productBeta;
        double productBetaLab2 = productBeta2 + boostBeta * boostBeta * ( 1.0 - productBeta2 * ( 1.0 - muPlus * muPlus ) ) + 2.0 * muPlus * productBeta * boostBeta;
        productBetaLab2 /= 1.0 - muPlus * productBeta * boostBeta;
        a_energy_out = MCGIDI::particleKineticEnergyFromBeta2( productMass( ), productBetaLab2 ); }
    else {
        a_energy_out = m_energy->sample( a_energy_in, a_rng( ), a_rng );
        probability =  m_angularGivenEnergy->evaluate( a_energy_in, a_energy_out, a_mu_lab );
    }

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data by sampling the outgoing mu from the probability P(mu|E) and then samples E' from
 * the probability P(E'|E,mu). It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::AngularEnergyMC::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    double mu_1, mu_2;

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_mu = m_angular->sample2dOf3d( a_X, a_rng( ), a_rng, &mu_1, &mu_2 );
    a_input.m_energyOut1 = m_energyGivenAngular->sample( a_X, mu_1, mu_2, a_rng( ), a_rng );
    a_input.m_phi = 2. * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::AngularEnergyMC::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    if( productFrame( ) != GIDI::Frame::lab ) LUPI_THROW( "AngularEnergyMC::angleBiasing: center-of-mass not supported." );

    a_energy_out = m_energyGivenAngular->sample( a_energy_in, a_mu_lab, a_mu_lab, a_rng( ), a_rng );
    return( m_angular->evaluate( a_energy_in, a_mu_lab ) );
}


/* *********************************************************************************************************//**
 * This method samples the outgoing product data using the Kalbach-Mann formalism.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::KalbachMann::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = m_f->sample( a_X, a_rng( ), a_rng );
    double rValue = m_r->evaluate( a_X, a_input.m_energyOut1 );
    double aValue = m_a->evaluate( a_X, a_input.m_energyOut1 );

        // In the following: Cosh[ a mu ] + r Sinh[ a mu ] = ( 1 - r ) Cosh[ a mu ] + r ( Cosh[ a mu ] + Sinh[ a mu ] ).
    if( a_rng( ) >= rValue ) { // Sample the '( 1 - r ) Cosh[ a mu ]' term.
        double T = ( 2. * a_rng( ) - 1. ) * sinh( aValue );

        a_input.m_mu = log( T + sqrt( T * T + 1. ) ) / aValue; }
    else {                                                                  // Sample the 'r ( Cosh[ a mu ] + Sinh[ a mu ] )' term.
        double rng1 = a_rng( ), exp_a = exp( aValue );

        a_input.m_mu = log( rng1 * exp_a + ( 1. - rng1 ) / exp_a ) / aValue;
    }
    if( a_input.m_mu < -1 ) a_input.m_mu = -1;
    if( a_input.m_mu >  1 ) a_input.m_mu = 1;

    a_input.m_phi = 2. * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::KalbachMann::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = 0.0;

    double initialMass = projectileMass( ) + targetMass( );
    double energy_out_com = m_f->sample( a_energy_in, a_rng( ), a_rng );
    double productBeta = MCGIDI_particleBeta( productMass( ), energy_out_com );
    double boostBeta = sqrt( a_energy_in * ( a_energy_in + 2. * projectileMass( ) ) ) / ( a_energy_in + initialMass );      // beta = v/c.

    double muPlus = 0.0, JacobianPlus = 0.0, muMinus = 0.0, JacobianMinus = 0.0;

    int numberOfMus = muCOM_From_muLab( a_mu_lab, boostBeta, productBeta, muPlus, JacobianPlus, muMinus, JacobianMinus );

    if( numberOfMus == 0 ) return( 0.0 );

    double rAtEnergyEnergyPrime = m_r->evaluate( a_energy_in, energy_out_com );
    double aAtEnergyEnergyPrime = m_a->evaluate( a_energy_in, energy_out_com );
    double aMu = aAtEnergyEnergyPrime * muPlus;

    double probability = 0.5 * JacobianPlus;
    if( productMass( ) == 0.0 ) {
        probability *= 1.0 - rAtEnergyEnergyPrime + rAtEnergyEnergyPrime * aAtEnergyEnergyPrime * exp( aMu ) / sinh( aAtEnergyEnergyPrime ); }
    else {
        probability *= aAtEnergyEnergyPrime * ( cosh( aMu ) + rAtEnergyEnergyPrime * cosh( aMu ) ) / sinh( aAtEnergyEnergyPrime );
    }

    if( numberOfMus == 2 ) {
        aMu = aAtEnergyEnergyPrime * muMinus;

        double probabilityMinus = 0.5 * JacobianMinus;
        if( productMass( ) == 0.0 ) {
            probabilityMinus *= 1.0 - rAtEnergyEnergyPrime + rAtEnergyEnergyPrime * aAtEnergyEnergyPrime * exp( aMu ) / sinh( aAtEnergyEnergyPrime ); }
        else {
            probabilityMinus *= aAtEnergyEnergyPrime * ( cosh( aMu ) + rAtEnergyEnergyPrime * cosh( aMu ) ) / sinh( aAtEnergyEnergyPrime );
        }
        probability += probabilityMinus;

        if( probabilityMinus > a_rng( ) * probability ) muPlus = muMinus;
    }

    double productBeta2 = productBeta * productBeta;
    double productBetaLab2 = productBeta2 + boostBeta * boostBeta * ( 1.0 - productBeta2 * ( 1.0 - muPlus * muPlus ) ) + 2.0 * muPlus * productBeta * boostBeta;
    productBetaLab2 /= 1.0 - muPlus * productBeta * boostBeta;
    a_energy_out = MCGIDI::particleKineticEnergyFromBeta2( productMass( ), productBetaLab2 );

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data from the coherent photo-atomic scattering law.
 * It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile in the lab frame.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::CoherentPhotoAtomicScattering::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    a_input.m_energyOut1 = a_X;

    int lowerIndex = binarySearchVector( a_X, m_energies );

    if( lowerIndex < 1 ) {
        do {
            a_input.m_mu = 1.0 - 2.0 * a_rng( );
        } while( ( 1.0 + a_input.m_mu * a_input.m_mu ) < 2.0 * a_rng( ) ); }
    else {
        double _a = m_a[lowerIndex];
        double X_i = m_energies[lowerIndex];
        double formFactor_i = m_formFactor[lowerIndex];
        double formFactor_X_i = formFactor_i * X_i;
        double Z = a_X / X_i;
        double realAnomalousFactor = 0.0;
        double imaginaryAnomalousFactor = 0.0;

        if( m_anomalousDataPresent ) {
            realAnomalousFactor = m_realAnomalousFactor->evaluate( a_X );
            imaginaryAnomalousFactor = m_imaginaryAnomalousFactor->evaluate( a_X );
        }

        double anomalousFactorSquared = realAnomalousFactor * realAnomalousFactor + imaginaryAnomalousFactor * imaginaryAnomalousFactor;
        double normalization = m_integratedFormFactorSquared[lowerIndex] + formFactor_X_i * formFactor_X_i * Z_a( Z, 2.0 * _a + 2.0 );

        anomalousFactorSquared = 0.0;
        if( anomalousFactorSquared != 0.0 ) {
            double integratedFormFactor_i = m_integratedFormFactor[lowerIndex] + formFactor_X_i * Z_a( Z, _a + 2.0 );

            normalization += 2.0  * integratedFormFactor_i * realAnomalousFactor + 0.5 * anomalousFactorSquared * a_X * a_X;
        }

        do {
            double partialIntegral = a_rng( ) * normalization;
            double X;
            if( anomalousFactorSquared == 0.0 ) {
                lowerIndex = binarySearchVector( partialIntegral, m_integratedFormFactorSquared );

                if( lowerIndex == 0 ) {
                    X = sqrt( 2.0 * partialIntegral ) / m_formFactor[0]; }
                else {
                    double remainer = partialIntegral - m_integratedFormFactorSquared[lowerIndex];
                    double epsilon = 2.0 * m_a[lowerIndex] + 2.0;

                    X_i = m_energies[lowerIndex];
                    formFactor_i = m_formFactor[lowerIndex];
                    formFactor_X_i = formFactor_i * X_i;

                    remainer /= formFactor_X_i * formFactor_X_i;
                    if( fabs( epsilon ) < 1e-6 ) {
                        X = X_i * exp( remainer ); }
                    else {
                        X = X_i * pow( 1.0 + epsilon * remainer, 1.0 / epsilon );
                    }
                } }
            else {                                  // Currently not implemented.
                X = 0.5 * a_X;
            }
            double X_E = X / a_X;
            a_input.m_mu = 1.0 - 2.0 * X_E * X_E;
        } while( ( 1.0 + a_input.m_mu * a_input.m_mu ) < 2.0 * a_rng( ) );
    }

    a_input.m_phi = 2.0 * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::CoherentPhotoAtomicScattering::angleBiasing( Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = a_energy_in;

    URR_protareInfos URR_protareInfos1;
    double sigma = a_reaction->protareSingle( )->reactionCrossSection( a_reaction->reactionIndex( ), URR_protareInfos1, 0.0, a_energy_in );
    double formFactor = evaluateFormFactor( a_energy_in, a_mu_lab );
    double imaginaryAnomalousFactor = 0.0;

    if( m_anomalousDataPresent ) {
        formFactor += m_realAnomalousFactor->evaluate( a_energy_in );
        imaginaryAnomalousFactor = m_imaginaryAnomalousFactor->evaluate( a_energy_in );
    }

    double probability = M_PI * MCGIDI_classicalElectronRadius * MCGIDI_classicalElectronRadius * ( 1.0 + a_mu_lab * a_mu_lab )
                        * ( formFactor * formFactor + imaginaryAnomalousFactor * imaginaryAnomalousFactor ) / sigma;

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data by sampling the outgoing energy E' from the probability P(E'|E) and then samples mu from
 * the probability P(mu|E,E'). It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::IncoherentPhotoAtomicScattering::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    double k1 = a_X / PoPI_electronMass_MeV_c2;
    double energyOut, mu, scatteringFactor;

    if( a_X >= m_energies.back( ) ) {
        MCGIDI_sampleKleinNishina( k1, a_rng, &energyOut, &mu ); }
    else {
        double scatteringFactorMax = evaluateScatteringFactor( a_X );
        do {
            MCGIDI_sampleKleinNishina( k1, a_rng, &energyOut, &mu );
            scatteringFactor = evaluateScatteringFactor( a_X * sqrt( 0.5 * ( 1.0 - mu ) ) );
        } while( scatteringFactor < a_rng( ) * scatteringFactorMax );
    }

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = energyOut * PoPI_electronMass_MeV_c2;
    a_input.m_mu = mu;
    a_input.m_phi = 2.0 * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing product data by sampling the outgoing energy E' from the probability P(E'|E) and then samples mu from
 * the probability P(mu|E,E'). It also samples the outgoing phi uniformly between 0 and 2 pi.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_userrng                 [in]    A random number generator that takes the state *a_rngState* and returns a double in the range [0.0, 1.0).
 * @param a_rngState                [in]    The current state for the random number generator.
 ***********************************************************************************************************/
template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering::sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    double energyOut, mu, occupationNumber;
    // Convert incident photon energy [MeV] to units of rest mass energy of the electron
    const double alpha_in = a_X / PoPI_electronMass_MeV_c2;
    double alpha_ratio, occupation_pz, occupationNumberMax;
    double quad_a = 0, quad_b = 0, quad_c = 0, pz = 0;  // Initialize with dummy values to silence compiler warnings

    bool energetically_possible = false;
    int ep_it = 0;
    while( energetically_possible == false && ep_it < 1000 ){
        // Sample outgoing angle
        occupationNumberMax = evaluateOccupationNumber( a_X, -1.0 );
        if( a_X >= 10.0 ) {  // This condition is not yet correct
            MCGIDI_sampleKleinNishina( alpha_in, a_rng, &energyOut, &mu ); }
        else {
            do {
                MCGIDI_sampleKleinNishina( alpha_in, a_rng, &energyOut, &mu );
                occupationNumber = evaluateOccupationNumber( a_X, mu );
            } while( occupationNumber < occupationNumberMax * a_rng( ) );
        }

        // Sample electron momentum projection, pz
        occupation_pz = occupationNumberMax*a_rng();
        int lowerIndex = binarySearchVector( occupation_pz, m_occupationNumber );
        if( lowerIndex == -1 ){
            pz = m_pz.back();
        }
        else{
            pz = m_pz[lowerIndex] + (occupation_pz-m_occupationNumber[lowerIndex])*(m_pz[lowerIndex+1]-m_pz[lowerIndex])/(m_occupationNumber[lowerIndex+1]-m_occupationNumber[lowerIndex]);
        }

        // Convert pz to outgoing photon energy
        alpha_ratio = energyRatio(a_X, mu);
        quad_a = pz*pz - (1/alpha_ratio)*(1/alpha_ratio);
        quad_b = -2*alpha_in*( pz*pz * mu - (1/alpha_ratio));
        quad_c = alpha_in*alpha_in*( pz*pz - 1 );

        if(quad_b*quad_b - 4*quad_a*quad_c > 0){
            energetically_possible = true;
        }
        ep_it = ep_it + 1;
    }

    const double quad_1 = -quad_b/(2*quad_a) + sqrt( quad_b*quad_b - 4*quad_a*quad_c )/( 2*quad_a );
    const double quad_2 = -quad_b/(2*quad_a) - sqrt( quad_b*quad_b - 4*quad_a*quad_c )/( 2*quad_a );

    // Select the correct outgoing energy based on the pz value
    if(pz >= 0){
        if(quad_1 >= quad_2){
            energyOut = quad_1;
        }
        else{
            energyOut = quad_2;
        }
    }
    else{
        if(quad_1 >= quad_2){
            energyOut = quad_2;
        }
        else{
            energyOut = quad_1;
        }
    }

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = energyOut * PoPI_electronMass_MeV_c2;
    a_input.m_mu = mu;
    a_input.m_phi = 2.0 * M_PI * a_rng( );
    a_input.m_frame = productFrame( );

}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::IncoherentPhotoAtomicScattering::angleBiasing( Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    URR_protareInfos URR_protareInfos1;
    double sigma = a_reaction->protareSingle( )->reactionCrossSection( a_reaction->reactionIndex( ), URR_protareInfos1, 0.0, a_energy_in );

    double norm = M_PI * MCGIDI_classicalElectronRadius * MCGIDI_classicalElectronRadius / sigma;

    double one_minus_mu = 1.0 - a_mu_lab;
    double k_in = a_energy_in / PoPI_electronMass_MeV_c2;
    a_energy_out = a_energy_in / ( 1.0 + k_in * one_minus_mu );
    double k_out = a_energy_out / PoPI_electronMass_MeV_c2;

    double k_ratio = k_out / k_in;
    double probability = evaluateScatteringFactor( a_energy_in * sqrt( 0.5 * one_minus_mu ) );
    probability *= k_ratio * k_ratio * ( 1.0 + a_mu_lab * a_mu_lab + k_in * k_out * one_minus_mu * one_minus_mu ) * norm;

    return( probability );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_userrng                 [in]    A random number generator that takes the state *a_rngState* and returns a double in the range [0.0, 1.0).
 * @param a_rngState                [in]    The current state for the random number generator.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering::angleBiasing( Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, double a_mu_lab, RNG && a_rng, double &a_energy_out ) const {

    URR_protareInfos URR_protareInfos1;
    double sigma = a_reaction->protareSingle( )->reactionCrossSection( a_reaction->reactionIndex( ), URR_protareInfos1, 0.0, a_energy_in );

    double norm = M_PI * MCGIDI_classicalElectronRadius * MCGIDI_classicalElectronRadius / sigma;

    double one_minus_mu = 1.0 - a_mu_lab;
    double alpha_in = a_energy_in / PoPI_electronMass_MeV_c2;

    double quad_a, quad_b, quad_c, alpha_ratio, pz, occupation_pz, occupationNumberMax;

    bool energetically_possible = false;
    int ep_it = 0;
    while( energetically_possible == false && ep_it < 1000 ){

        // Sample electron momentum projection, pz
        occupationNumberMax = evaluateOccupationNumber( a_energy_in, -1.0 );
        occupation_pz = occupationNumberMax*a_rng();
        int lowerIndex = binarySearchVector( occupation_pz, m_occupationNumber );
        pz = 0;
        if( lowerIndex == -1 ){
            pz = m_pz.back();
        }
        else{
            pz = m_pz[lowerIndex] + (occupation_pz-m_occupationNumber[lowerIndex])*(m_pz[lowerIndex+1]-m_pz[lowerIndex])/(m_occupationNumber[lowerIndex+1]-m_occupationNumber[lowerIndex]);
        }

        // Convert pz to outgoing photon energy
        alpha_ratio = energyRatio(a_energy_in, a_mu_lab);
        quad_a = pz*pz - (1/alpha_ratio)*(1/alpha_ratio);
        quad_b = -2*alpha_in*( pz*pz * a_mu_lab - (1/alpha_ratio));
        quad_c = alpha_in*alpha_in*( pz*pz - 1 );
        if(quad_b*quad_b - 4*quad_a*quad_c > 0){
            energetically_possible = true;
        }
        ep_it = ep_it + 1;
    }

    const double quad_1 = -quad_b/(2*quad_a) + sqrt( quad_b*quad_b - 4*quad_a*quad_c )/( 2*quad_a );
    const double quad_2 = -quad_b/(2*quad_a) - sqrt( quad_b*quad_b - 4*quad_a*quad_c )/( 2*quad_a );

    // Select the correct outgoing energy based on the pz value
    double alpha_out = 0;
    if(pz >= 0){
        if(quad_1 >= quad_2){
            alpha_out = quad_1;
        }
        else{
            alpha_out = quad_2;
        }
    }
    else{
        if(quad_1 >= quad_2){
            alpha_out = quad_2;
        }
        else{
            alpha_out = quad_1;
        }
    }
    alpha_ratio = alpha_out / alpha_in;
    a_energy_out = alpha_out * PoPI_electronMass_MeV_c2;
    double probability = evaluateOccupationNumber( a_energy_in, a_mu_lab );
    probability *= alpha_ratio * alpha_ratio * ( 1.0 + a_mu_lab * a_mu_lab + alpha_in * alpha_out * one_minus_mu * one_minus_mu ) * norm;

    return( probability );
}

/* *********************************************************************************************************//**
 * This method returns the outgoing electron energy and angle given that the photon when out at an angle of *a_input.m_mu*.
 * Ergo, this method must be called directly after the photon has been sampled.
 *
 * @param a_energy                  [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::IncoherentPhotoAtomicScatteringElectron::sample( double a_energy, Sampling::Input &a_input, RNG && a_rng ) const {

    double halfTheta = 0.5 * acos( a_input.m_mu );
    double cot_psi = ( 1.0 + a_energy / PoPI_electronMass_MeV_c2 ) * tan( halfTheta );
    double psi = atan( 1.0 / cot_psi );

    double deltaE_photon = a_energy - a_input.m_energyOut1;
    double electronMomentum2 = deltaE_photon * ( deltaE_photon + 2.0 * PoPI_electronMass_MeV_c2 );  // Square of the electron outlgoing momentum.

    a_input.m_energyOut1 = electronMomentum2 / ( sqrt( electronMomentum2 + PoPI_electronMass_MeV_c2 * PoPI_electronMass_MeV_c2 ) + PoPI_electronMass_MeV_c2 );
    a_input.m_mu = cos( psi );
    a_input.m_phi = 2.0 * M_PI * a_rng( );
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* causing a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 * Currently, this method only returns 0.0 for the probability and outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::IncoherentPhotoAtomicScatteringElectron::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
                LUPI_maybeUnused double a_energy_in, LUPI_maybeUnused double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = 0;
    return( 0.0 );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing photon by assigning the electron rest mass energy as the photon's energy and,
 * if m_firstSampled is true, randomly picking mu and phi. If m_firstSampled is false, the previous sampled particle
 * that filled in a_input must be the other sampled photon, then, the mu and phi for the second-sampled photon is such that 
 * it is back-to-back with the other photon.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::PairProductionGamma::sample( LUPI_maybeUnused double a_X, Sampling::Input &a_input, RNG && a_rng ) const {

    if( m_firstSampled ) {
        a_input.m_mu = 1.0 - 2.0 * a_rng( );
        a_input.m_phi = M_PI * a_rng( ); }
    else {
        a_input.m_mu *= -1.0;
        a_input.m_phi += M_PI;
    }
    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = PoPI_electronMass_MeV_c2;
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::PairProductionGamma::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		LUPI_maybeUnused double a_energy_in, LUPI_maybeUnused double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = PoPI_electronMass_MeV_c2;
    return( 1.0 );                          // 1.0 as there are two photons, each with 1/2 probability.
}

/* *********************************************************************************************************//**
 * This method samples the outgoing neutron data for coherent elastic TSNL from the Debye/Waller function.
 *
 * @param a_energy                  [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::CoherentElasticTNSL::sample( double a_energy, Sampling::Input &a_input,
                RNG && a_rng ) const {

    if( a_energy <= m_energies[0] ) {
        a_input.m_mu = 1.0; }
    else {
        double temperature = 1e-3 * a_input.m_temperature;                  // Assumes m_temperature is in keV/K.
        if( temperature < m_temperatures[0] ) temperature = m_temperatures[0];
        if( temperature > m_temperatures.back( ) ) temperature = m_temperatures.back( );
        std::size_t temperatureIndex = (std::size_t) MCGIDI::binarySearchVector( temperature, m_temperatures, true );
        double const *pointer1 = &m_S_table[temperatureIndex * m_energies.size( )];

        double const *pointer2 = pointer1;
        double fractionFirstTemperature = 1.0;
        if( temperatureIndex != ( m_temperatures.size( ) - 1 ) ) {
            fractionFirstTemperature = ( m_temperatures[temperatureIndex+1] - temperature ) / ( m_temperatures[temperatureIndex+1] -  m_temperatures[temperatureIndex] );
            pointer2 += m_energies.size( );
        }
        double fractionSecondTemperature = 1.0 - fractionFirstTemperature;

        int energyIndexMax = MCGIDI::binarySearchVector( a_energy, m_energies, true );
        if( a_energy == m_energies[energyIndexMax] ) --energyIndexMax;

        double randomTotal = a_rng( ) * ( fractionFirstTemperature * pointer1[energyIndexMax] + fractionSecondTemperature * pointer2[energyIndexMax] );
        int energyIndex = 0;
        for( ; energyIndex < energyIndexMax; ++energyIndex ) {
            if( randomTotal <= fractionFirstTemperature * pointer1[energyIndex] + fractionSecondTemperature * pointer2[energyIndex] ) break;
        }
        a_input.m_mu = 1.0 - 2.0 * m_energies[energyIndex] / a_energy;
    }

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = a_energy;
    a_input.m_phi = 2.0 * M_PI * a_rng( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::CoherentElasticTNSL::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		double a_energy_in, LUPI_maybeUnused double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

//    double temperature = 1e-3 * a_temperature;                          // Assumes a_temperature is in keV/K.
    double probability = 0.0;

    a_energy_out = a_energy_in;

    return( probability );
}

/* *********************************************************************************************************//**
 * This method samples the outgoing neutron data for incoherent elastic TSNL from the Debye/Waller function.
 *
 * @param a_energy                  [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::IncoherentElasticTNSL::sample( double a_energy, Sampling::Input &a_input, 
                RNG && a_rng ) const {

    double temperature = 1e-3 * a_input.m_temperature / m_temperatureToMeV_K;                  // Assumes m_temperature is in keV/K.
    double W_prime = m_DebyeWallerIntegral->evaluate( temperature );
    double twoEW = 2 * a_energy * W_prime;
    double expOfTwice_twoEW = exp( -2 * twoEW );
    double sampled_cdf = a_rng( );

    if( sampled_cdf > ( 1 - 1e-5 ) ) {
        double Variable = ( 1.0 - sampled_cdf ) * ( 1.0 - expOfTwice_twoEW );
        a_input.m_mu = 1.0 - Variable * ( 1.0 + 0.5 * Variable ) / twoEW; }
    else if( sampled_cdf < expOfTwice_twoEW ) {
        a_input.m_mu = -1.0 + log( sampled_cdf / expOfTwice_twoEW * ( 1.0 - expOfTwice_twoEW ) + 1.0 ) / twoEW; }
    else {
        a_input.m_mu = 1.0 + log( expOfTwice_twoEW + sampled_cdf * ( 1.0 - expOfTwice_twoEW ) ) / twoEW;
    }

    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_energyOut1 = a_energy;
    a_input.m_phi = 2.0 * M_PI * a_rng( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    The temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::IncoherentElasticTNSL::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, double a_temperature,
		double a_energy_in, double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    double temperature = 1e-3 * a_temperature / m_temperatureToMeV_K;                          // Assumes a_temperature is in keV/K.
    double W_prime = m_DebyeWallerIntegral->evaluate( temperature );
    double twoEW = 2 * a_energy_in * W_prime;
    double probability = exp( -twoEW * ( 1.0 - a_mu_lab ) ) * twoEW / ( 1.0 - exp( -2 * twoEW ) );

    a_energy_out = a_energy_in;

    return( probability );
}

/* *********************************************************************************************************//**
 * The method sets all outgoing product data to 0.0 and set the sampledType to Sampling::unspecified.
 *
 * @param a_X                       [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Distributions::Unspecified::sample( LUPI_maybeUnused double a_X, Sampling::Input &a_input, LUPI_maybeUnused RNG && a_rng ) const {

    a_input.m_sampledType = Sampling::SampledType::unspecified;
    a_input.m_energyOut1 = 0.;
    a_input.m_mu = 0.;
    a_input.m_phi = 0.;
    a_input.m_frame = productFrame( );
}

/* *********************************************************************************************************//**
 * Returns the probability for a projectile with energy *a_energy_in* to cause a particle to be emitted 
 * at angle *a_mu_lab* as seen in the lab frame. *a_energy_out* is the sampled outgoing energy. This one should never
 * be called. If called, returns 0.0 for a probability.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 *
 * @return                                  The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Distributions::Unspecified::angleBiasing( LUPI_maybeUnused Reaction const *a_reaction, LUPI_maybeUnused double a_temperature,
		LUPI_maybeUnused double a_energy_in, LUPI_maybeUnused double a_mu_lab, LUPI_maybeUnused RNG && a_rng, double &a_energy_out ) const {

    a_energy_out = 0.0;
    return( 0.0 );
}


// From file: MCGIDI_functions.cpp

/*
============================================================
*/
template <typename RNG >
LUPI_HOST_DEVICE int MCGIDI::Functions::Function1d::sampleBoundingInteger( double a_x1, RNG && a_rng ) const {

    if( type( ) == Function1dType::TerrellFissionNeutronMultiplicityModel ) 
        return( static_cast<TerrellFissionNeutronMultiplicityModel const *>( this )->sampleBoundingInteger( a_x1, a_rng ) );

    double d_value = evaluate( a_x1 );
    int iValue = (int) d_value;
    if( iValue == d_value ) return( iValue );
    if( d_value - iValue > a_rng( ) ) ++iValue;

    return( iValue );
}

/* *********************************************************************************************************//**
 * Sample the number of fission prompt neutrons using Terrell's modified Gaussian distribution.
 * Method uses Red Cullen's algoritm (see UCRL-TR-222526).
 *
 * @param a_energy              [in]        The energy of the projectile.
 * @param a_rng                 [in]        The random number generator function the uses *a_rngState* to generator a double in the range [0, 1.0).
 * @param a_rngState            [in/out]    The random number generator state.
 *
 * @return                                  The sampled number of emitted, prompt neutrons for fission.
 ***********************************************************************************************************/
template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::Functions::TerrellFissionNeutronMultiplicityModel::sampleBoundingInteger( double a_energy, RNG && a_rng ) const {

    const double Terrell_BSHIFT = -0.43287;
    double width = M_SQRT2 * m_width;
    double temp1 = m_multiplicity->evaluate( a_energy ) + 0.5;
    double temp2 = temp1 / width;
    double expo = exp( -temp2 * temp2 );
    double cshift = temp1 + Terrell_BSHIFT * m_width * expo / ( 1.0 - expo );

    double multiplicity = 1.0;
    do {
      double rw = sqrt( -log( a_rng( ) ) );
      double theta = ( 2.0 * M_PI ) * a_rng();

      multiplicity = width * rw * cos( theta ) + cshift;
    } while ( multiplicity < 0.0 );

    return( static_cast<int>( floor( multiplicity ) ) );
}

/* *********************************************************************************************************//**
 * Returns the x-value corresponding cumulative probability *a_rngValue*.
 *
 * @param a_rngValue            [in]    The x-value to evaluate the function at.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                              The value of the function at *a_x1*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase1d::sample( double a_rngValue, RNG && a_rng ) const {

    return( static_cast<Xs_pdf_cdf1d const *>( this )->sample( a_rngValue, a_rng ) );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::Xs_pdf_cdf1d::sample( double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const {

    int lower = binarySearchVector( a_rngValue, m_cdf );
    double domainValue = 0;


    if( lower < 0 ) {                                   // This should never happen.
        LUPI_THROW( "Xs_pdf_cdf1d::sample: lower < 0." );
    }

    if( interpolation( ) == Interpolation::FLAT ) {
        double fraction = ( m_cdf[lower+1] - a_rngValue ) / ( m_cdf[lower+1] - m_cdf[lower] );
        domainValue = fraction * m_Xs[lower] + ( 1 - fraction ) * m_Xs[lower+1]; }
    else {                                              // Assumes lin-lin interpolation.
        double slope = m_pdf[lower+1] - m_pdf[lower];

        if( slope == 0.0 ) {
            if( m_pdf[lower] == 0.0 ) {
                domainValue = m_Xs[lower];
                if( lower == 0 ) domainValue = m_Xs[1]; }
            else {
                double fraction = ( m_cdf[lower+1] - a_rngValue ) / ( m_cdf[lower+1] - m_cdf[lower] );
                domainValue = fraction * m_Xs[lower] + ( 1 - fraction ) * m_Xs[lower+1];
            } }
        else {
            double d1, d2;

            slope = slope / ( m_Xs[lower+1] - m_Xs[lower] );
            d1 = a_rngValue - m_cdf[lower];
            d2 = m_cdf[lower+1] - a_rngValue;
            if( d2 > d1 ) {                         // Closer to lower.
                domainValue = m_Xs[lower] + ( sqrt( m_pdf[lower] * m_pdf[lower] + 2. * slope * d1 ) - m_pdf[lower] ) / slope; }
            else {                                  // Closer to lower + 1.
                domainValue = m_Xs[lower+1] - ( m_pdf[lower+1] - sqrt( m_pdf[lower+1] * m_pdf[lower+1] - 2. * slope * d2 ) ) / slope;
            }
        }
    }
    return( domainValue );
}

/* *********************************************************************************************************//**
 * This method samples an x1 from a pdf(x1|x2) given x2 and the cumulative value of the pdf as *a_rngValue*.
 *
 * @param a_x2                  [in]        The value of x2.
 * @param a_rngValue            [in]        The value of the cumulative used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::none:
        break;
    case ProbabilityBase2dType::weightedFunctionals:
        value = static_cast<WeightedFunctionals2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    default:
        value = static_cast<ProbabilityBase2d_d1 const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    }

    return( value );
}

/* *********************************************************************************************************//**
 * Returns the value of x1, given x2 and the cumulative probability *a_rngValue*.
 *
 * @param a_x2                  [in]        Value of the outer most independent variable (i.e., *x2*).
 * @param a_rngValue            [in]        The value of the cumulative probability used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                                  The value *x1* where the cumulative probability is *a_rngValue* for x2 = *a_x2*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase2d_d1::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

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
        value = static_cast<ProbabilityBase2d_d2 const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::regions:
        value = static_cast<Regions2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::none:
    case ProbabilityBase2dType::weightedFunctionals:
        LUPI_THROW( "ProbabilityBase2d_d1::sample: This should never happen." );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * This method returns two x1 values for use with ProbabilityBase3d functions.
 *
 * @param a_x2                  [in]        The value of x2.
 * @param a_rngValue            [in]        The value of the cumulative value used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 * @param a_x1_1                [in]        The lower value of the x1 value.
 * @param a_x1_2                [in]        The upper value of the x1 value.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase2d_d1::sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng,
                double *a_x1_1, double *a_x1_2 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::XYs:
        value = static_cast<XYs2d const *>( this )->sample2dOf3d( a_x2, a_rngValue, a_rng, a_x1_1, a_x1_2 );
        break;
    default:
        LUPI_THROW( "ProbabilityBase2d_d1::sample2dOf3d: not implemented." );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * Returns the value of x1, given x2 and the cumulative probability *a_rngValue*.
 *
 * @param a_x2                  [in]        Value of the outer most independent variable (i.e., *x2*).
 * @param a_rngValue            [in]        The value of the cumulative probability used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                                  The value *x1* where the cumulative probability is *a_rngValue* for x2 = *a_x2*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase2d_d2::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::XYs:
        value = static_cast<XYs2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::isotropic:
        value = static_cast<Isotropic2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::discreteGamma:
        value = static_cast<DiscreteGamma2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::primaryGamma:
        value = static_cast<PrimaryGamma2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::recoil:
        value = static_cast<Recoil2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::NBodyPhaseSpace:
        value = static_cast<NBodyPhaseSpace2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::evaporation:
        value = static_cast<Evaporation2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::generalEvaporation:
        value = static_cast<GeneralEvaporation2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::simpleMaxwellianFission:
        value = static_cast<SimpleMaxwellianFission2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::Watt:
        value = static_cast<Watt2d const *>( this )->sample( a_x2, a_rngValue, a_rng );
        break;
    case ProbabilityBase2dType::none:
    case ProbabilityBase2dType::weightedFunctionals:
    case ProbabilityBase2dType::regions:
        LUPI_THROW( "ProbabilityBase2d_d2::sample: This should never happen." );
    }

    return( value );
}

/* *********************************************************************************************************//**
 * This method returns two x1 values for use with ProbabilityBase3d functions.
 *
 * @param a_x2                  [in]        The value of x2.
 * @param a_rngValue            [in]        The value of the cumulative value used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 * @param a_x1_1                [in]        The lower value of the x1 value.
 * @param a_x1_2                [in]        The upper value of the x1 value.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase2d_d2::sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng,
                double *a_x1_1, double *a_x1_2 ) const {

    double value = 0.0;

    switch( type( ) ) {
    case ProbabilityBase2dType::XYs:
        value = static_cast<XYs2d const *>( this )->sample2dOf3d( a_x2, a_rngValue, a_rng, a_x1_1, a_x1_2 );
        break;
    default:
        LUPI_THROW( "ProbabilityBase2d_d2::sample2dOf3d: not implemented." );
    }

    return( value );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::XYs2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {
/*
C    Samples from a pdf(x1|x2). First determine which pdf(s) to sample from given x2.
C    Then use rngValue to sample from pdf1(x1) and maybe pdf2(x1) and interpolate to
C    determine x1.
*/
    double sampledValue = 0;
    int lower = binarySearchVector( a_x2, m_Xs );

    if( lower == -2 ) {
        sampledValue = m_probabilities[0]->sample( a_rngValue, a_rng ); }
    else if( lower == -1 ) {
        sampledValue = m_probabilities.back( )->sample( a_rngValue, a_rng ); }
    else {
        double sampled1 = m_probabilities[lower]->sample( a_rngValue, a_rng );

        if( interpolation( ) == Interpolation::FLAT ) {
            sampledValue = sampled1; }
        else {
            double sampled2 = m_probabilities[lower+1]->sample( a_rngValue, a_rng );

            if( interpolation( ) == Interpolation::LINLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x2 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = fraction * sampled1 + ( 1 - fraction ) * sampled2; }
            else if( interpolation( ) == Interpolation::LOGLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x2 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = sampled2 * pow( sampled2 / sampled1, fraction ); }
            else if( interpolation( ) == Interpolation::LINLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x2 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = fraction * sampled1 + ( 1 - fraction ) * sampled2; }
            else if( interpolation( ) == Interpolation::LOGLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x2 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = sampled2 * pow( sampled2 / sampled1, fraction ); }
            else {                                                              // This should never happen.
                LUPI_THROW( "XYs2d::sample: unsupported interpolation." );
            }
        }
    }
    return( sampledValue );
}

/* *********************************************************************************************************//**
 * This method returns two x1 values for use with ProbabilityBase3d functions.
 *
 * @param a_x2                  [in]        The value of x2.
 * @param a_rngValue            [in]        The value of the cumulative value used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 * @param a_x1_1                [in]        The lower value of the x1 value.
 * @param a_x1_2                [in]        The upper value of the x1 value.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::XYs2d::sample2dOf3d( double a_x2, double a_rngValue, RNG && a_rng, 
                double *a_x1_1, double *a_x1_2 ) const {
/*
C   Samples from a pdf(x1|x2). First determine which pdf(s) to sample from given x2. Then use rngValue to sample from pdf1(x1) 
C   and maybe pdf2(x1) and interpolate to determine x1.
*/
    double sampledValue = 0;
    int lower = binarySearchVector( a_x2, m_Xs );

    if( lower == -2 ) {
        sampledValue = m_probabilities[0]->sample( a_rngValue, a_rng );
        *a_x1_2 = *a_x1_1 = sampledValue; }
    else if( lower == -1 ) {
        sampledValue = m_probabilities.back( )->sample( a_rngValue, a_rng );
        *a_x1_2 = *a_x1_1 = sampledValue; }
    else {
        *a_x1_1 = m_probabilities[lower]->sample( a_rngValue, a_rng );

        if( interpolation( ) == Interpolation::FLAT ) {
            sampledValue = *a_x1_2 = *a_x1_1; }
        else {
            *a_x1_2 = m_probabilities[lower+1]->sample( a_rngValue, a_rng );

            if( interpolation( ) == Interpolation::LINLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x2 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = fraction * *a_x1_1 + ( 1 - fraction ) * *a_x1_2; }
            else if( interpolation( ) == Interpolation::LOGLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x2 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = *a_x1_2 * pow( *a_x1_2 / *a_x1_1, fraction ); }
            else if( interpolation( ) == Interpolation::LINLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x2 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = fraction * *a_x1_1 + ( 1 - fraction ) * *a_x1_2; }
            else if( interpolation( ) == Interpolation::LOGLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x2 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = *a_x1_2 * pow( *a_x1_2 / *a_x1_1 , fraction ); }
            else {                                                              // This should never happen.
                LUPI_THROW( "XYs2d::sample: unsupported interpolation." );
            }
        }
    }
    return( sampledValue );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::Regions2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

    int lower = binarySearchVector( a_x2, m_Xs );

    if( lower < 0 ) {
        if( lower == -1 ) {                         // a_x2 > last value of m_Xs.
            return( m_probabilities.back( )->sample( a_x2, a_rngValue, a_rng ) );
        }
        lower = 0;                                  // a_x2 < first value of m_Xs.
    }

    return( m_probabilities[lower]->sample( a_x2, a_rngValue, a_rng ) );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::Recoil2d::sample( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const {

#if !defined(__NVCC__) && !defined(__HIP__)
    LUPI_THROW( "Recoil2d::sample: not implemented." );
#endif

    return( 0.0 );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::NBodyPhaseSpace2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

    return( ( m_energy_in_COMFactor * a_x2 + m_Q ) * m_massFactor * m_dist->sample( a_rngValue, a_rng ) );
}

inline LUPI_HOST_DEVICE static double MCGIDI_sampleEvaporation( double a_xMax, double a_rngValue ) {

    double b1, c1, xMid, norm, xMin = 0.;

    norm = 1 - ( 1 + a_xMax ) * exp( -a_xMax );
    b1 = 1. - norm * a_rngValue;
    for( int i1 = 0; i1 < 16; i1++ ) {
        xMid = 0.5 * ( xMin + a_xMax );
        c1 = ( 1 + xMid ) * exp( -xMid );
        if( b1 > c1 ) {
            a_xMax = xMid; }
        else {
            xMin = xMid;
        }
    }
    return( xMid );
}


template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::Evaporation2d::sample( double a_x2, double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const {

    double theta = m_theta->evaluate( a_x2 );

    return( theta * MCGIDI_sampleEvaporation( ( a_x2 - m_U ) / theta, a_rngValue ) );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::GeneralEvaporation2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {

    return( m_theta->evaluate( a_x2 ) * m_g->sample( a_rngValue, a_rng ) );
}

inline LUPI_HOST_DEVICE static double MCGIDI_sampleSimpleMaxwellianFission( double a_xMax, double a_rngValue ) {

    double b1, c1, xMid, norm, xMin = 0., sqrt_xMid, sqrt_pi_2 = 0.5 * sqrt( M_PI );

    sqrt_xMid = sqrt( a_xMax );
    norm = sqrt_pi_2 * erf( sqrt_xMid ) - sqrt_xMid * exp( -a_xMax );
    b1 = norm * a_rngValue;
    for( int i1 = 0; i1 < 16; i1++ ) {
        xMid = 0.5 * ( xMin + a_xMax );
        sqrt_xMid = sqrt( xMid );
        c1 = sqrt_pi_2 * erf( sqrt_xMid ) - sqrt_xMid * exp( -xMid );
        if( b1 < c1 ) {
            a_xMax = xMid; }
        else {
            xMin = xMid;
        }
    }
    return( xMid );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::SimpleMaxwellianFission2d::sample( double a_x2, double a_rngValue, LUPI_maybeUnused RNG && a_rng ) const {

    double theta = m_theta->evaluate( a_x2 );

    return( theta * MCGIDI_sampleSimpleMaxwellianFission( ( a_x2 - m_U ) / theta, a_rngValue ) );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::Watt2d::sample( double a_x2, LUPI_maybeUnused double a_rngValue, RNG && a_rng ) const {
/*
*   From MCAPM via Sample Watt Spectrum as in TART ( Kalos algorithm ).
*/
    double WattMin = 0., WattMax = a_x2 - m_U, x, y, z, energyOut, rand1, rand2;
    double Watt_a = 1./m_a->evaluate( a_x2 );  // Kalos algorithm uses the inverse of the ‘a’ parameter stored in GNDS
    double Watt_b = m_b->evaluate( a_x2 );

    x = 1. + ( Watt_b / ( 8. * Watt_a ) );
    y = ( x + sqrt( x * x - 1. ) ) / Watt_a;
    z = Watt_a * y - 1.;
    do {
        rand1 = -log( a_rng( ) );
        rand2 = -log( a_rng( ) );
        energyOut = y * rand1;
    } while( ( ( rand2 - z * ( rand1 + 1. ) ) * ( rand2 - z * ( rand1 + 1. ) ) > Watt_b * y * rand1 ) || 
             ( energyOut < WattMin ) || ( energyOut > WattMax ) );
    return( energyOut );
}


template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::WeightedFunctionals2d::sample( double a_x2, double a_rngValue, RNG && a_rng ) const {
/*
c   This routine assumes that the weights sum to 1.
*/
    std::size_t i1;
    std::size_t n1 = m_weight.size( ) - 1;      // Take last point if others do not add to randomWeight.
    double randomWeight = a_rng( ), cumulativeWeight = 0.;

    for( i1 = 0; i1 < n1; ++i1 ) {
        cumulativeWeight += m_weight[i1]->evaluate( a_x2 );
        if( cumulativeWeight >= randomWeight ) break;
    }
    return( m_energy[i1]->sample( a_x2, a_rngValue, a_rng) );
}

/* *********************************************************************************************************//**
 * This method samples an x1 from a pdf(x1|x2) given x2 and the cumulative value of the pdf as *a_rngValue*.
 *
 * @param a_x3                  [in]        The value of x3.
 * @param a_x2_1                [in]        The value of ?.
 * @param a_x2_2                [in]        The value of ?.
 * @param a_rngValue            [in]        The value of the cumulative used to determine the x1 value.
 * @param a_rng                 [in]        The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::ProbabilityBase3d::sample( double a_x3, double a_x2_1, double a_x2_2, double a_rngValue, RNG && a_rng ) const {

    return( static_cast<XYs3d const *>( this )->sample( a_x3, a_x2_1, a_x2_2, a_rngValue, a_rng ) );
}

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Probabilities::XYs3d::sample( double a_x3, double a_x2_1, double a_x2_2, double a_rngValue, RNG && a_rng ) const {
/*
C    Samples from a pdf(x1|x3,x2). First determine which pdf(s) to sample from given x3
C    Then use rngValue to sample from pdf2_1(x2) and maybe pdf2_2(x2) and interpolate to
C    determine x1.
*/
    double sampledValue = 0;
    int lower = binarySearchVector( a_x3, m_Xs );

    if( lower == -2 ) {                         // x3 < first value of Xs.
        sampledValue = m_probabilities[0]->sample( a_x2_1, a_rngValue, a_rng ); }
    else if( lower == -1 ) {                    // x3 > last value of Xs.
        sampledValue = m_probabilities.back( )->sample( a_x2_1, a_rngValue, a_rng ); }
    else {
        double sampled1 = m_probabilities[lower]->sample( a_x2_1, a_rngValue, a_rng );

        if( interpolation( ) == Interpolation::FLAT ) {
            sampledValue = sampled1; }
        else {
            double sampled2 = m_probabilities[lower+1]->sample( a_x2_2, a_rngValue, a_rng );

            if( interpolation( ) == Interpolation::LINLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x3 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = fraction * sampled1 + ( 1 - fraction ) * sampled2; }
            else if( interpolation( ) == Interpolation::LOGLIN ) {
                double fraction = ( m_Xs[lower+1] - a_x3 ) / ( m_Xs[lower+1] - m_Xs[lower] );
                sampledValue = sampled2 * pow( sampled2 / sampled1, fraction ); }
            else if( interpolation( ) == Interpolation::LINLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x3 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = fraction * sampled1 + ( 1 - fraction ) * sampled2; }
            else if( interpolation( ) == Interpolation::LOGLOG ) {
                double fraction = log( m_Xs[lower+1] / a_x3 ) / log( m_Xs[lower+1] / m_Xs[lower] );
                sampledValue = sampled2 * pow( sampled2 / sampled1, fraction ); }
            else {                                                              // This should never happen.
                LUPI_THROW( "XYs3d::sample: unsupported interpolation." );
            }
        }
    }

    return( sampledValue );
}


// From file: MCGIDI_heatedCrossSections.cpp


template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::HeatedCrossSectionsContinuousEnergy::sampleReaction( URR_protareInfos const &a_URR_protareInfos, int a_URR_index, 
                int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng ) const {

    int i1, sampled_reaction_index, temperatureIndex1, temperatureIndex2, number_of_temperatures = static_cast<int>( m_temperatures.size( ) );
    double sampleCrossSection = a_crossSection * a_rng( );

    if( a_temperature <= m_temperatures[0] ) {
        temperatureIndex1 = 0;
        temperatureIndex2 = temperatureIndex1; }
    else if( a_temperature >= m_temperatures.back( ) ) {
        temperatureIndex1 = static_cast<int>( m_temperatures.size( ) ) - 1;
        temperatureIndex2 = temperatureIndex1; }
    else {
        for( i1 = 0; i1 < number_of_temperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        temperatureIndex1 = i1 - 1;
        temperatureIndex2 = i1;
    }

    int numberOfReactions = m_heatedCrossSections[0]->numberOfReactions( );
    double energyFraction1, energyFraction2, crossSectionSum = 0.0;

    HeatedCrossSectionContinuousEnergy &heatedCrossSection1 = *m_heatedCrossSections[temperatureIndex1];
    int energyIndex1 = heatedCrossSection1.evaluationInfo( a_hashIndex, a_energy, &energyFraction1 );

    if( temperatureIndex1 == temperatureIndex2 ) {
        for( sampled_reaction_index = 0; sampled_reaction_index < numberOfReactions; ++sampled_reaction_index ) {
            crossSectionSum += heatedCrossSection1.reactionCrossSection2( sampled_reaction_index, a_URR_protareInfos, a_URR_index, a_energy, 
                    energyIndex1, energyFraction1 );
            if( crossSectionSum >= sampleCrossSection ) break;
        } }
    else {
        double temperatureFraction2 = ( a_temperature - m_temperatures[temperatureIndex1] ) 
                / ( m_temperatures[temperatureIndex2] - m_temperatures[temperatureIndex1] );
        double temperatureFraction1 = 1.0 - temperatureFraction2;
        HeatedCrossSectionContinuousEnergy &heatedCrossSection2 = *m_heatedCrossSections[temperatureIndex2];
        int energyIndex2 = heatedCrossSection2.evaluationInfo( a_hashIndex, a_energy, &energyFraction2 );

        for( sampled_reaction_index = 0; sampled_reaction_index < numberOfReactions; ++sampled_reaction_index ) {
            if( m_thresholds[sampled_reaction_index] >= a_energy ) continue;
            crossSectionSum += temperatureFraction1 * heatedCrossSection1.reactionCrossSection2( sampled_reaction_index, a_URR_protareInfos, 
                    a_URR_index, a_energy, energyIndex1, energyFraction1 );
            crossSectionSum += temperatureFraction2 * heatedCrossSection2.reactionCrossSection2( sampled_reaction_index, a_URR_protareInfos, 
                    a_URR_index, a_energy, energyIndex2, energyFraction2 );
            if( crossSectionSum >= sampleCrossSection ) break;
        }
    }

    if( sampled_reaction_index == numberOfReactions ) {
        if( crossSectionSum < ( 1.0 - crossSectionSumError ) * a_crossSection ) {
#if LUPI_ON_GPU
            MCGIDI_PRINTF( "HeatedCrossSectionsContinuousEnergy::sampleReaction: crossSectionSum %.17e less than a_crossSection =  %.17e.", 
                    crossSectionSum, a_crossSection );
#else
            std::string errorString = "HeatedCrossSectionsContinuousEnergy::sampleReaction: crossSectionSum " 
                    + LUPI::Misc::doubleToString3( "%.17e", crossSectionSum ) + " less than a_crossSection = " 
                    + LUPI::Misc::doubleToString3( "%.17e", a_crossSection ) + ".";
            LUPI_THROW( errorString.c_str( ) );
#endif
        }
        for( sampled_reaction_index = 0; sampled_reaction_index < numberOfReactions; ++sampled_reaction_index ) {   // This should rarely happen so just pick the first reaction with non-zero cross section.
            if( heatedCrossSection1.reactionCrossSection2( sampled_reaction_index, a_URR_protareInfos, a_URR_index, a_energy, energyIndex1, 
                    energyFraction1, true ) > 0 ) break;
        }
    }

    return( sampled_reaction_index );
}

/* *********************************************************************************************************//**
 * Returns the requested reaction's multi-group cross section for target temperature *a_temperature* and projectile multi-group *a_hashIndex*.
 *
 * @param a_hashIndex           [in]    The multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_crossSection        [in]    The index of the reaction.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::HeatedCrossSectionsMultiGroup::sampleReaction( int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, 
                RNG && a_rng ) const {

    int i1, sampled_reaction_index, temperatureIndex1, temperatureIndex2, numberOfTemperatures = static_cast<int>( m_temperatures.size( ) );
    double sampleCrossSection = a_crossSection * a_rng( );

    if( a_temperature <= m_temperatures[0] ) {
        temperatureIndex1 = 0;
        temperatureIndex2 = temperatureIndex1; }
    else if( a_temperature >= m_temperatures.back( ) ) {
        temperatureIndex1 = static_cast<int>( m_temperatures.size( ) ) - 1;
        temperatureIndex2 = temperatureIndex1; }
    else {
        for( i1 = 0; i1 < numberOfTemperatures; ++i1 ) if( a_temperature < m_temperatures[i1] ) break;
        temperatureIndex1 = i1 - 1;
        temperatureIndex2 = i1;
    }

    int numberOfReactions = m_heatedCrossSections[0]->numberOfReactions( );
    double crossSectionSum = 0;
    HeatedCrossSectionMultiGroup &heatedCrossSection1 = *m_heatedCrossSections[temperatureIndex1];

    if( temperatureIndex1 == temperatureIndex2 ) {
        for( sampled_reaction_index = 0; sampled_reaction_index < numberOfReactions; ++sampled_reaction_index ) {
            crossSectionSum += heatedCrossSection1.reactionCrossSection( sampled_reaction_index, a_hashIndex, true );
            if( crossSectionSum >= sampleCrossSection ) break;
        } }
    else {
        double temperatureFraction2 = ( a_temperature - m_temperatures[temperatureIndex1] ) / ( m_temperatures[temperatureIndex2] - m_temperatures[temperatureIndex1] );
        double temperatureFraction1 = 1.0 - temperatureFraction2;
        HeatedCrossSectionMultiGroup &heatedCrossSection2 = *m_heatedCrossSections[temperatureIndex2];

        for( sampled_reaction_index = 0; sampled_reaction_index < numberOfReactions; ++sampled_reaction_index ) {
            if( m_thresholds[sampled_reaction_index] >= a_energy ) continue;
            crossSectionSum += temperatureFraction1 * heatedCrossSection1.reactionCrossSection( sampled_reaction_index, a_hashIndex, true );
            crossSectionSum += temperatureFraction2 * heatedCrossSection2.reactionCrossSection( sampled_reaction_index, a_hashIndex, true );
            if( crossSectionSum >= sampleCrossSection ) break;
        }
    }

    if( sampled_reaction_index == numberOfReactions ) return( MCGIDI_nullReaction );

    if( m_multiGroupThresholdIndex[sampled_reaction_index] == a_hashIndex ) {
        double energyAboveThreshold = a_energy - m_thresholds[sampled_reaction_index];

        if( energyAboveThreshold <= ( a_rng( ) * ( m_projectileMultiGroupBoundariesCollapsed[a_hashIndex+1] - m_thresholds[sampled_reaction_index] ) ) )
            return( MCGIDI_nullReaction );
    }

    return( sampled_reaction_index );
}

// From file: MCGIDI_misc.cpp

/* *********************************************************************************************************//**
 * The function returns a normalized Maxwellian speed (i.e., v = |velocity|) in 3d (i.e., v^2 Exp( -v^2 )).
 * Using formula in https://link.springer.com/content/pdf/10.1007%2Fs10955-011-0364-y.pdf.
 * Author Nader M.A. Mohamed, title "Efficient Algorithm for Generating Maxwell Random Variables".
 *
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                              The sampled normalized Maxwellian speed.
 ***********************************************************************************************************/

template <typename RNG>
inline LUPI_HOST_DEVICE double sampleBetaFromMaxwellian( RNG && a_rng ) {

    double _g = 2.0 / ( 1.37 * 0.5 * 1.772453850905516 );      // 1.772453850905516 = sqrt( pi ).
    double beta, r1;

    do {
        r1 = a_rng( );
        beta = sqrt( -2.0 * log( r1 ) );
    } while( _g * r1 * beta < a_rng( ) );

    return( beta );
}

/* *********************************************************************************************************//**
 * This function is used internally to sample a target's velocity (speed and cosine of angle relative to projectile)
 * for a heated target using zero temperature, multi-grouped cross sections.
 *
 * @param a_protare             [in]    The Protare instance for the projectile and target.
 * @param a_projectileEnergy    [in]    The energy of the projectile in the lab frame of the target.
 * @param a_input               [in]    Contains needed input like the targets temperature. Also will have the target sampled velocity on return if return value is *true*.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                              Returns *true* if target velocity is sampled and false otherwise.
 ***********************************************************************************************************/
namespace MCGIDI {
template <typename RNG>
inline LUPI_HOST_DEVICE bool sampleTargetBetaForUpscatterModelA( Protare const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng ) {

    double projectileBeta = MCGIDI_particleBeta( a_protare->projectileMass( ), a_projectileEnergy );

    double temperature = a_input.m_temperature * 1e-3;                   // FIXME Assumes m_temperature is in keV/k for now.
    double targetThermalBeta = MCGIDI_particleBeta( a_protare->targetMass( ), temperature );

    if( targetThermalBeta < 1e-4 * projectileBeta ) return( false );

    double relativeBetaMin = projectileBeta - 2.0 * targetThermalBeta;
    double relativeBetaMax = projectileBeta + 2.0 * targetThermalBeta;

    Vector<double> const &upscatterModelAGroupVelocities = a_protare->upscatterModelAGroupVelocities( );
    int maxIndex = (int) upscatterModelAGroupVelocities.size( ) - 2;
    int relativeBetaMinIndex = binarySearchVector( relativeBetaMin, upscatterModelAGroupVelocities, true );
    int relativeBetaMaxIndex = binarySearchVector( relativeBetaMax, upscatterModelAGroupVelocities, true );
    double targetBeta, relativeBeta, mu;

    if( relativeBetaMinIndex >= maxIndex ) relativeBetaMinIndex = maxIndex;
    if( relativeBetaMaxIndex >= maxIndex ) relativeBetaMaxIndex = maxIndex;

    if( relativeBetaMinIndex == relativeBetaMaxIndex ) {
        targetBeta = targetThermalBeta * sampleBetaFromMaxwellian( a_rng );
        mu = 1.0 - 2.0 * a_rng( );
        relativeBeta = sqrt( targetBeta * targetBeta + projectileBeta * projectileBeta - 2.0 * mu * targetBeta * projectileBeta ); }
    else {

        Vector<double> const &upscatterModelACrossSection = a_input.m_reaction->upscatterModelACrossSection( );
        double reactionRate;
        double reactionRateMax = 0;
        for( int i1 = relativeBetaMinIndex; i1 <= relativeBetaMaxIndex; ++i1 ) {
            reactionRate = upscatterModelACrossSection[i1] * upscatterModelAGroupVelocities[i1+1];
            if( reactionRate > reactionRateMax ) reactionRateMax = reactionRate;
        }

        do {
            targetBeta = targetThermalBeta * sampleBetaFromMaxwellian( a_rng );
            mu = 1.0 - 2.0 * a_rng( );
            relativeBeta = sqrt( targetBeta * targetBeta + projectileBeta * projectileBeta - 2.0 * mu * targetBeta * projectileBeta );

            int index = binarySearchVector( relativeBeta, upscatterModelAGroupVelocities, true );
            if( index > maxIndex ) index = maxIndex;
            reactionRate = upscatterModelACrossSection[index] * relativeBeta;
        } while( reactionRate <  a_rng( ) * reactionRateMax );
    }

    a_input.m_projectileBeta = projectileBeta;
    a_input.m_relativeMu = mu;
    a_input.m_targetBeta = targetBeta;
    a_input.m_relativeBeta = relativeBeta;
    a_input.m_projectileEnergy = particleKineticEnergy( a_protare->projectileMass( ), relativeBeta );

    return( true );
}

/* *********************************************************************************************************//**
 * This function boost a particle from one frame to another frame. The frames have a relative speed *a_boostSpeed*
 * and cosine of angle *a_boostMu* between their z-axes. BRB FIXME, currently it is the x-axis.
 *
 * @param a_input                   [in]    Instance containing a random number generator that returns a double in the range [0, 1).
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_product                 [in]    The particle to boost.
 ***********************************************************************************************************/

template <typename RNG>
inline LUPI_HOST_DEVICE void upScatterModelABoostParticle( Sampling::Input &a_input, RNG && a_rng, Sampling::Product &a_product ) {

    double C_rel = 1.0;
    if( a_input.m_relativeBeta != 0.0 ) {
        C_rel = ( a_input.m_projectileBeta - a_input.m_relativeMu * a_input.m_targetBeta ) / a_input.m_relativeBeta;
        if( C_rel > 1.0 ) C_rel = 1.0;                  // Handle round-off issue. Probably should check who big the issue is.
    }
    double S_rel = sqrt( 1.0 - C_rel * C_rel );

    double pz_vz = a_product.m_pz_vz;
    a_product.m_pz_vz =  C_rel * a_product.m_pz_vz + S_rel * a_product.m_px_vx;
    a_product.m_px_vx = -S_rel * pz_vz             + C_rel * a_product.m_px_vx;

    double targetSpeed = MCGIDI_speedOfLight_cm_sec * a_input.m_targetBeta;
    a_product.m_pz_vz += a_input.m_relativeMu * targetSpeed;
    a_product.m_px_vx += sqrt( 1.0 - a_input.m_relativeMu * a_input.m_relativeMu ) * targetSpeed;

    double phi = 2.0 * M_PI * a_rng( );
    double sine = sin( phi );
    double cosine = cos( phi );
    double px_vx = a_product.m_px_vx;
    a_product.m_px_vx = cosine * a_product.m_px_vx - sine   * a_product.m_py_vy;
    a_product.m_py_vy = sine   * px_vx             + cosine * a_product.m_py_vy;

    double speed2 = a_product.m_px_vx * a_product.m_px_vx + a_product.m_py_vy * a_product.m_py_vy + a_product.m_pz_vz * a_product.m_pz_vz;
    speed2 /= MCGIDI_speedOfLight_cm_sec * MCGIDI_speedOfLight_cm_sec;

    a_product.m_kineticEnergy = particleKineticEnergyFromBeta2( a_product.m_productMass, speed2 );
}
}


// From file: MCGIDI_outputChannel.cpp


/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*.
 *
 * @param a_protare                 [in]    The Protare this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::OutputChannel::sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

    if( m_hasFinalStatePhotons ) {
        double random = a_rng( );
        double cumulative = 0.0;
        bool sampled = false;
        for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
            cumulative += (*productIter)->multiplicity( )->evaluate( a_projectileEnergy );
            if( cumulative >= random ) {
                (*productIter)->sampleFinalState( a_protare, a_projectileEnergy, a_input, a_rng, a_push_back, a_products );
                sampled = true;
                break;
            }
        }
        if( !sampled ) {     // BRB: FIXME: still need to code for continuum photon.
        } }
    else {
        for( Vector<Product *>::const_iterator iter = m_products.begin( ); iter != m_products.end( ); ++iter )
            (*iter)->sampleProducts( a_protare, a_projectileEnergy, a_input, a_rng, a_push_back, a_products );
    }

    if( m_totalDelayedNeutronMultiplicity != nullptr ) {
        double totalDelayedNeutronMultiplicity = m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );

        if( a_rng( ) < totalDelayedNeutronMultiplicity ) {       // Assumes that totalDelayedNeutronMultiplicity < 1.0, which it is.
            double sum = 0.0;

            totalDelayedNeutronMultiplicity *= a_rng( );
            for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
                DelayedNeutron const *delayedNeutron1( delayedNeutron( i1 ) );
                Product const &product = delayedNeutron1->product( );

                sum += product.multiplicity( )->evaluate( a_projectileEnergy );
                if( sum >= totalDelayedNeutronMultiplicity ) {
                    product.distribution( )->sample( a_projectileEnergy, a_input, a_rng );
                    a_input.m_delayedNeutronIndex = delayedNeutron1->delayedNeutronIndex( );
                    a_input.m_delayedNeutronDecayRate = delayedNeutron1->rate( );
                    a_products.add( a_projectileEnergy, product.intid( ), product.index( ), product.userParticleIndex( ), product.mass( ), 
                            a_input, a_rng, a_push_back, false );
                    break;
                }
            }
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the probability for a project with energy *a_energy_in* to cause this channel to emitted a particle of index
 * *a_index* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_index                   [in]    The index of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_weight                  [in]    The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The sum of the multiplicity for other outgoing particles with index *a_index*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::OutputChannel::angleBiasing( Reaction const *a_reaction, int a_index, double a_temperature, double a_energy_in, double a_mu_lab, 
                double &a_weight, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const {

    for( Vector<Product *>::const_iterator iter = m_products.begin( ); iter != m_products.end( ); ++iter )
        (*iter)->angleBiasing( a_reaction, a_index, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight );

    if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_index == m_neutronIndex ) ) {
        for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
            DelayedNeutron const *delayedNeutron1( delayedNeutron( i1 ) );
            Product const &product = delayedNeutron1->product( );

            product.angleBiasing( a_reaction, a_index, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight );
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the probability for a project with energy *a_energy_in* to cause this channel to emitted a particle of intid
 * *a_intid* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_intid                   [in]    The intid of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_weight                  [in]    The probability of emitting outgoing particle into lab angle *a_mu_lab*.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The sum of the multiplicity for other outgoing particles with intid *a_intid*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::OutputChannel::angleBiasingViaIntid( Reaction const *a_reaction, int a_intid, double a_temperature, double a_energy_in, double a_mu_lab,
                double &a_weight, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const {

    for( Vector<Product *>::const_iterator iter = m_products.begin( ); iter != m_products.end( ); ++iter )
        (*iter)->angleBiasingViaIntid( a_reaction, a_intid, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight );

    if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_intid == PoPI::Intids::neutron ) ) {
        for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
            DelayedNeutron const *delayedNeutron1( delayedNeutron( i1 ) );
            Product const &product = delayedNeutron1->product( );

            product.angleBiasingViaIntid( a_reaction, a_intid, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight );
        }
    }
}


// From file: MCGIDI_product.cpp

/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*.
 *
 * @param a_protare                 [in]    The Protare this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::Product::sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) {
        m_outputChannel->sampleProducts( a_protare, a_projectileEnergy, a_input, a_rng, a_push_back, a_products ); }
    else {
#endif
        if( m_twoBodyOrder == TwoBodyOrder::secondParticle ) {
            a_products.add( a_projectileEnergy, intid( ), index( ), userParticleIndex( ), mass( ), a_input, a_rng, a_push_back, m_intid == PoPI::Intids::photon ); }
        else {
            int _multiplicity = m_multiplicity->sampleBoundingInteger( a_projectileEnergy, a_rng );
            int __multiplicity = _multiplicity;

            for( ; _multiplicity > 0; --_multiplicity ) {
                m_distribution->sample( a_projectileEnergy, a_input, a_rng );
                a_input.m_delayedNeutronIndex = -1;
                a_input.m_delayedNeutronDecayRate = 0.0;
                a_products.add( a_projectileEnergy, intid( ), index( ), userParticleIndex( ), mass( ), a_input, a_rng, a_push_back, m_intid == PoPI::Intids::photon );
            }
            if( m_initialStateIndex >= 0 ) {
                if( __multiplicity == 0 ) {
                    a_protare->sampleBranchingGammas( a_input, a_projectileEnergy, m_initialStateIndex, a_rng, a_push_back, a_products );
                }
            }
        }
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    }
#endif
}

/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*. In particular, the product is a capture reaction 
 * primary gamma what has a finalState attribute. This gamma is added as well as the gammas from the
 * gamma cascade.
 *
 * @param a_protare                 [in]    The Protare this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/
template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::Product::sampleFinalState( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

    m_distribution->sample( a_projectileEnergy, a_input, a_rng );
    a_input.m_delayedNeutronIndex = -1;
    a_input.m_delayedNeutronDecayRate = 0.0;
    a_products.add( a_projectileEnergy, m_intid, m_index, m_userParticleIndex, mass( ), a_input, a_rng, a_push_back, m_intid == PoPI::Intids::photon );

    if( m_initialStateIndex >= 0 ) {
        a_protare->sampleBranchingGammas( a_input, a_projectileEnergy, m_initialStateIndex, a_rng, a_push_back, a_products );
    }
}

/* *********************************************************************************************************//**
 * Returns the weight for a projectile with energy *a_energy_in* to cause this channel to emitted a particle of index
 * *a_pid* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_index                   [in]    The index of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_weight                  [in]    The weight of emitting outgoing particle into lab angle *a_mu_lab*.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The sum of the multiplicity for other outgoing particles with index *a_index*.
 ***********************************************************************************************************/
template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Product::angleBiasing( Reaction const *a_reaction, int a_index, double a_temperature, double a_energy_in, double a_mu_lab, 
                double &a_weight, double &a_energy_out, RNG && a_rng, double &a_cumulative_weight ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) {
        m_outputChannel->angleBiasing( a_reaction, a_index, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight ); }
    else {
#endif
        if( m_index != a_index ) return;

        double probability = 0.0;
        double energy_out = 0.0;

        if( a_cumulative_weight == 0.0 ) a_energy_out = 0.0;

        if( m_multiplicity->type( ) == Function1dType::branching ) { // Needs to handle F1_Branching.
            }
        else {
            probability = m_distribution->angleBiasing( a_reaction, a_temperature, a_energy_in, a_mu_lab, a_rng, energy_out );
        }

        double weight = m_multiplicity->evaluate( a_energy_in ) * probability;
        a_cumulative_weight += weight;
        if( weight > a_rng( ) * a_cumulative_weight ) {
            a_weight = weight;
            a_energy_out = energy_out;
        }
#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    }
#endif
}

/* *********************************************************************************************************//**
 * Returns the weight for a projectile with energy *a_energy_in* to cause this channel to emitted a particle of intid
 * *a_intid* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy.
 *
 * @param a_reaction                [in]    The reaction containing the particle which this distribution describes.
 * @param a_intid                   [in]    The intid of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_weight                  [in]    The weight of emitting outgoing particle into lab angle *a_mu_lab*.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The sum of the multiplicity for other outgoing particles with intid *a_intid*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE void MCGIDI::Product::angleBiasingViaIntid( Reaction const *a_reaction, int a_intid, double a_temperature, double a_energy_in, double a_mu_lab,
                double &a_weight, double &a_energy_out, RNG && a_rng , double &a_cumulative_weight ) const {

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    if( m_outputChannel != nullptr ) {
        m_outputChannel->angleBiasingViaIntid( a_reaction, a_intid, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight ); }
    else {
#endif
        if( m_intid != a_intid ) return;

        angleBiasing( a_reaction, m_index, a_temperature, a_energy_in, a_mu_lab, a_weight, a_energy_out, a_rng, a_cumulative_weight );

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    }
#endif
}



// From file: MCGIDI_protare.cpp


/* *********************************************************************************************************//**
 * Samples a reaction of *this* and returns its index.
 *
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_hashIndex         [in]    The cross section hash index.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 * @param   a_crossSection      [in]    The total cross section at *a_temperature* and *a_energy*.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                          The index of the sampled reaction.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::Protare::sampleReaction( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, 
                double a_energy, double a_crossSection, RNG && a_rng ) const {

    int reactionIndex = -1;

    switch( protareType( ) ) {
    case ProtareType::single: 
        reactionIndex = static_cast<ProtareSingle const *>( this )->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, 
                a_crossSection, a_rng );
        break;
    case ProtareType::composite:
        reactionIndex = static_cast<ProtareComposite const *>( this )->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, 
                a_crossSection, a_rng );
        break;
    case ProtareType::TNSL:
        reactionIndex = static_cast<ProtareTNSL const *>( this )->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, 
                a_crossSection, a_rng );
        break;
    }

    return( reactionIndex );
}

/* *********************************************************************************************************//**
 * Samples gammas from a nuclide electro-magnetic decay.
 *
 * @param a_input               [in]    Sample options requested by user.
 * @param a_projectileEnergy    [in]    The energy of the projectile.
 * @param a_initialStateIndex   [in]    The index in *m_nuclideGammaBranchStateInfos* whose nuclide data are used for sampling.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products            [in]    The object to add all sampled gammas to.
 ***********************************************************************************************************/
template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::ProtareSingle::sampleBranchingGammas( Sampling::Input &a_input, double a_projectileEnergy, int a_initialStateIndex, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

    int initialStateIndex = a_initialStateIndex;
    double energyLevelSampleWidthUpper = 0.0;           // Used for GRIN continuum levels to add variaction to outgoing photons.

    NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = nullptr;
    if( initialStateIndex >= 0 ) nuclideGammaBranchStateInfo = m_nuclideGammaBranchStateInfos[initialStateIndex];
// std::cout << initialStateIndex << " " << nuclideGammaBranchStateInfo->nuclearLevelEnergy( ) << std::endl;
    while( initialStateIndex >= 0 ) {
        Vector<int> const &branchIndices = nuclideGammaBranchStateInfo->branchIndices( );

        double random = a_rng( );
        double sum = 0.0;
        initialStateIndex = -1;             // Just in case the for loop never has "sum >= random".
        for( std::size_t i1 = 0; i1 < branchIndices.size( ); ++i1 ) {
            NuclideGammaBranchInfo *nuclideGammaBranchInfo = m_branches[branchIndices[i1]];

            sum += nuclideGammaBranchInfo->probability( );
            if( sum >= random ) {
                double energyLevelSampleWidthLower = 0.0;
                initialStateIndex = nuclideGammaBranchInfo->residualStateIndex( );
                if( initialStateIndex >= 0 ) {
                    nuclideGammaBranchStateInfo = m_nuclideGammaBranchStateInfos[initialStateIndex];
                    energyLevelSampleWidthLower = a_rng( ) * nuclideGammaBranchStateInfo->nuclearLevelEnergyWidth( );
                }
                if( nuclideGammaBranchInfo->photonEmissionProbability( ) > a_rng( ) ) {
                    a_input.m_sampledType = Sampling::SampledType::photon;
                    a_input.m_dataInTargetFrame = false;
                    a_input.m_frame = GIDI::Frame::lab;

                    a_input.m_energyOut1 = nuclideGammaBranchInfo->gammaEnergy( ) + energyLevelSampleWidthUpper - energyLevelSampleWidthLower;
// std::cout << a_input.m_energyOut1 << " " << nuclideGammaBranchInfo->gammaEnergy( ) << " " << energyLevelSampleWidthUpper << " " << energyLevelSampleWidthLower << std::endl;
                    a_input.m_mu = 1.0 - 2.0 * a_rng( );
                    a_input.m_phi = 2.0 * M_PI * a_rng( );

                    a_products.add( a_projectileEnergy, PoPI::Intids::photon, m_photonIndex, userPhotonIndex( ), 0.0, a_input, a_rng, a_push_back, true );
                }
                energyLevelSampleWidthUpper = energyLevelSampleWidthLower;
                break;
            }
        }
    }
}

/* *********************************************************************************************************//**
 * Returns the index of a sampled reaction for a target with termpature *a_temperature*, a projectile with energy *a_energy* and total cross section 
 * *a_crossSection*. Random numbers are obtained via *a_rng*.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    Specifies the continuous energy or multi-group index.
 * @param a_temperature         [in]    The temperature of the target.
 * @param a_energy              [in]    The energy of the projectile.
 * @param a_crossSection        [in]    The total cross section.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::ProtareSingle::sampleReaction( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, 
                double a_crossSection, RNG && a_rng ) const {

    if( m_continuousEnergy ) return( m_heatedCrossSections.sampleReaction( a_URR_protareInfos, m_URR_index, a_hashIndex, a_temperature, a_energy, 
            a_crossSection, a_rng ) );

    return( m_heatedMultigroupCrossSections.sampleReaction( a_hashIndex, a_temperature, a_energy, a_crossSection, a_rng ) );
}


// From file: MCGIDI_protareComposite.cpp

/* *********************************************************************************************************//**
 * Samples a reaction of *this* and returns its index.
 *
 * @param   a_URR_protareInfos  [in]    URR information.
 * @param   a_hashIndex         [in]    The cross section hash index.
 * @param   a_temperature       [in]    The target temperature.
 * @param   a_energy            [in]    The projectile energy.
 * @param   a_crossSection      [in]    The total cross section at *a_temperature* and *a_energy*.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                          The index of the sampled reaction.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::ProtareComposite::sampleReaction( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng ) const {

    std::size_t length = static_cast<std::size_t>( m_protares.size( ) );
    int reaction_index = 0;
    double cross_section_sum = 0.0;
    double cross_section_rng = a_rng( ) * a_crossSection;

    for( std::size_t i1 = 0; i1 < length; ++i1 ) {
        double cross_section = m_protares[i1]->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, true );

        cross_section_sum += cross_section;
        if( cross_section_sum > cross_section_rng ) {
            int reaction_index2 = m_protares[i1]->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, cross_section, a_rng );

            reaction_index += reaction_index2;
            if( reaction_index2  == MCGIDI_nullReaction ) reaction_index = MCGIDI_nullReaction;
            break;
        }
        reaction_index += m_protares[i1]->numberOfReactions( );
    }

    return( reaction_index );
}


// From file: MCGIDI_protareTNSL.cpp

/* *********************************************************************************************************//**
 * Returns the total cross section.
 *
 * @param a_URR_protareInfos    [in]    URR information.
 * @param a_hashIndex           [in]    The cross section hash index.
 * @param a_temperature         [in]    The target temperature.
 * @param a_energy              [in]    The projectile energy.
 * @param a_crossSection        [in]    The total cross section at *a_temperature* and *a_energy*.
 * @param a_rng                 [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                              The index of the sampled reaction.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::ProtareTNSL::sampleReaction( URR_protareInfos const &a_URR_protareInfos, int a_hashIndex, double a_temperature, double a_energy, double a_crossSection, RNG && a_rng ) const {

    int reactionIndex = 0;

    if( ( a_energy < m_TNSL_maximumEnergy ) && ( a_temperature <= m_TNSL_maximumTemperature ) ) {
        double TNSL_crossSection = m_TNSL->crossSection( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, true );

        if( TNSL_crossSection > a_rng( ) * a_crossSection ) {
            reactionIndex = m_TNSL->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, TNSL_crossSection, a_rng ); }
        else { 
            reactionIndex = m_protareWithoutElastic->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_crossSection - TNSL_crossSection, a_rng );
            if( reactionIndex != MCGIDI_nullReaction ) reactionIndex += m_numberOfTNSLReactions + 1;
        } }
    else {
        reactionIndex = m_protareWithElastic->sampleReaction( a_URR_protareInfos, a_hashIndex, a_temperature, a_energy, a_crossSection, a_rng );
        if( reactionIndex != MCGIDI_nullReaction ) reactionIndex += m_numberOfTNSLReactions;
    }

    return( reactionIndex );
}


// From file: MCGIDI_reaction.cpp

/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*.
 *
 * @param a_protare                 [in]    The Protare this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 * @param a_checkOrphanProducts     [in]    If true, associated orphan products are also sampled.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::Reaction::sampleProducts( Protare const *a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products, bool a_checkOrphanProducts ) const {

    double projectileEnergy = a_projectileEnergy;

    a_input.m_GRIN_intermediateResidual = -1;
    a_input.m_reaction = this;
    a_input.m_projectileMass = m_projectileMass;
    a_input.m_targetMass = m_targetMass;
    a_input.m_relativeMu = 0.0;
    a_input.m_targetBeta = 0.0;
    a_input.m_relativeBeta = MCGIDI_particleBeta( m_projectileMass, a_projectileEnergy );

    a_input.m_dataInTargetFrame = false;
    if( upscatterModelASupported( ) && ( a_input.m_upscatterModel == Sampling::Upscatter::Model::A ) ) {

        a_input.m_dataInTargetFrame = sampleTargetBetaForUpscatterModelA( m_protareSingle, a_projectileEnergy, a_input, a_rng );
        if( a_input.m_dataInTargetFrame ) projectileEnergy = a_input.m_projectileEnergy;
    }

    if( m_GRIN_specialSampleProducts ) {
        if( m_GRIN_capture != nullptr ) {
            if( projectileEnergy < m_GRIN_maximumCaptureIncidentEnergy ) {
                if( m_GRIN_capture->sampleProducts( (ProtareSingle const *) a_protare, projectileEnergy, a_input, a_rng, a_push_back, a_products ) ) {
                    return;
                }
            } }
        else if( m_GRIN_inelastic != nullptr ) {
            if( m_GRIN_inelastic->sampleProducts( (ProtareSingle const *) a_protare, projectileEnergy, a_input, a_rng, a_push_back, a_products ) ) {
                return;
            }
        }
    }

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->sampleProducts( m_protareSingle, projectileEnergy, a_input, a_rng, a_push_back, a_products );
#else
    if( m_hasFinalStatePhotons ) {
        double random = a_rng( );
        double cumulative = 0.0;
        bool sampled = false;
        for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
            cumulative += (*productIter)->multiplicity( )->evaluate( projectileEnergy );
            if( cumulative >= random ) {
                (*productIter)->sampleFinalState( m_protareSingle, projectileEnergy, a_input, a_rng, a_push_back, a_products );
                sampled = true;
                break;
            }
        }
        if( !sampled ) {     // BRB: FIXME: still need to code for continuum photon.
        } }
    else {
        for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
            (*productIter)->sampleProducts( m_protareSingle, projectileEnergy, a_input, a_rng, a_push_back, a_products );
        }
    }

    if( m_totalDelayedNeutronMultiplicity != nullptr ) {
        double totalDelayedNeutronMultiplicity = m_totalDelayedNeutronMultiplicity->evaluate( a_projectileEnergy );

        if( a_rng( ) < totalDelayedNeutronMultiplicity ) {       // Assumes that totalDelayedNeutronMultiplicity < 1.0, which it is.
            double sum = 0.0;

            totalDelayedNeutronMultiplicity *= a_rng( );
            for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
                DelayedNeutron const *delayedNeutron1 = m_delayedNeutrons[i1];
                Product const &product = delayedNeutron1->product( );

                sum += product.multiplicity( )->evaluate( a_projectileEnergy );
                if( sum >= totalDelayedNeutronMultiplicity ) {
                    product.distribution( )->sample( a_projectileEnergy, a_input, a_rng );
                    a_input.m_delayedNeutronIndex = delayedNeutron1->delayedNeutronIndex( );
                    a_input.m_delayedNeutronDecayRate = delayedNeutron1->rate( );
                    a_products.add( a_projectileEnergy, product.intid( ), product.index( ), product.userParticleIndex( ), product.mass( ), a_input, a_rng, a_push_back, false );
                    break;
                }
            }
        }
    }

    if( m_fissionResiduaIntid != -1 ) {             // Special treatment to add 2 ENDL 99120 or 99125 products.
        a_input.m_sampledType = MCGIDI::Sampling::SampledType::unspecified;
        a_input.m_frame = GIDI::Frame::lab;
        a_input.m_energyOut1 = 0.0;
        a_input.m_mu = 0.0;
        a_input.m_phi = 0.0;
        a_input.m_delayedNeutronIndex = -1;
        a_input.m_delayedNeutronDecayRate = 0.0;
        a_products.add( 0.0, m_fissionResiduaIntid, m_fissionResiduaIndex, m_fissionResiduaUserIndex, m_fissionResidualMass, a_input, a_rng, a_push_back, false );
        a_products.add( 0.0, m_fissionResiduaIntid, m_fissionResiduaIndex, m_fissionResiduaUserIndex, m_fissionResidualMass, a_input, a_rng, a_push_back, false );
    }

#endif

    if( a_checkOrphanProducts ) {
        for( auto productIter = m_associatedOrphanProducts.begin( ); productIter != m_associatedOrphanProducts.end( ); ++productIter ) {
            (*productIter)->sampleProducts( m_protareSingle, projectileEnergy, a_input, a_rng, a_push_back, a_products );
        }
    }
}

/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*.
 *
 * @param a_protare                 [in]    The ProtareSingle this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE bool MCGIDI::GRIN_capture::sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

    std::size_t index = 0;
    double random = a_rng( );
    for( ; index < m_summedProbabilities.size( ) - 1; ++index ) {
        if( random < m_summedProbabilities[index] ) break;
    }
    GRIN_captureLevelProbability *GRIN_captureLevelProbability1 = m_captureLevelProbabilities[index];

    double availableEnergy = m_captureNeutronSeparationEnergy + a_projectileEnergy;
    int primaryCaptureLevelIndex = GRIN_captureLevelProbability1->sampleCaptureLevel( a_protare, availableEnergy, a_rng );
    NuclideGammaBranchStateInfo const *nuclideGammaBranchStateInfo = a_protare->nuclideGammaBranchStateInfos( )[primaryCaptureLevelIndex];

    a_input.m_GRIN_intermediateResidual = nuclideGammaBranchStateInfo->intid( );

    a_input.m_sampledType = Sampling::SampledType::photon;
    a_input.m_dataInTargetFrame = false;
    a_input.m_frame = GIDI::Frame::lab;

    a_input.m_energyOut1 = availableEnergy - nuclideGammaBranchStateInfo->nuclearLevelEnergy( );
    a_input.m_mu = 2 * a_rng( ) - 1.0;
    a_input.m_phi = 2.0 * M_PI * a_rng( );

    a_products.add( a_projectileEnergy, PoPI::Intids::photon, a_protare->photonIndex( ), a_protare->userPhotonIndex( ), 
            0.0, a_input, a_rng, a_push_back, true );

    a_protare->sampleBranchingGammas( a_input, a_projectileEnergy, primaryCaptureLevelIndex, a_rng, a_push_back, a_products );

    if( m_residualIntid != -1 ) {
        a_input.m_energyOut1 = 0.0;
        a_input.m_mu = 0.0;
        a_input.m_phi = 0.0;
        a_products.add( a_projectileEnergy, m_residualIntid, m_residualIndex, m_residualUserIndex, m_residualMass, a_input, a_rng, a_push_back, false );
    }

    return( true );
}

/* *********************************************************************************************************//**
 * This method adds sampled products to *a_products*.
 *
 * @param a_protare                 [in]    The ProtareSingle this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE bool MCGIDI::GRIN_inelastic::sampleProducts( ProtareSingle const *a_protare, double a_projectileEnergy, Sampling::Input &a_input,
                RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) const {

    std::size_t index = 1;
    for( ; index < m_energies.size( ); ++index ) {
        if( m_energies[index] > a_projectileEnergy ) break;
    }
    --index;

    GRIN_inelasticForEnergy *inelasticForEnergy = m_inelasticForEnergy[index];
    int levelIndex = inelasticForEnergy->sampleLevelIndex( a_projectileEnergy, a_rng( ) );
    if( levelIndex < 0 ) return( false );

    NuclideGammaBranchStateInfo const *nuclideGammaBranchStateInfo = a_protare->nuclideGammaBranchStateInfos( )[levelIndex];

    a_input.m_GRIN_intermediateResidual = nuclideGammaBranchStateInfo->intid( );

    double residualMass = m_targetMass + nuclideGammaBranchStateInfo->nuclearLevelEnergy( );
    double initialMass = m_neutronMass + m_targetMass;
    double finalMass = m_neutronMass + residualMass;
    double twoBodyThreshold = 0.5 * ( finalMass * finalMass - initialMass * initialMass ) / m_targetMass;
    double betaBoast = sqrt( a_projectileEnergy * ( a_projectileEnergy + 2. * m_neutronMass ) ) 
            / ( a_projectileEnergy + m_neutronMass + m_targetMass );      // betaBoast = v/c.
    double _x = m_targetMass * ( a_projectileEnergy - twoBodyThreshold ) / ( finalMass * finalMass );
    if( _x < 0 ) _x = 0.;           // FIXME There needs to be a better test here.
    double Kp;
    if( _x < 2e-5 ) {
        Kp = finalMass * _x * ( 1 - 0.5 * _x * ( 1 - _x ) ); }
    else {          // This is the relativistic formula derived from E^2 - (pc)^2 is frame independent.
        Kp = sqrt( finalMass * finalMass + 2 * m_targetMass * ( a_projectileEnergy - twoBodyThreshold ) ) - finalMass;
    }
    if( Kp < 0 ) Kp = 0.;           // FIXME There needs to be a better test here.

    a_input.m_sampledType = Sampling::SampledType::firstTwoBody;
    a_input.m_mu = 1.0 - 2.0 * a_rng( );
    a_input.m_phi = 2. * M_PI * a_rng( );
    kinetics_COMKineticEnergy2LabEnergyAndMomentum( betaBoast, Kp, m_neutronMass, residualMass, a_input );

    a_input.m_delayedNeutronIndex = -1;
    a_input.m_delayedNeutronDecayRate = 0.0;
    a_products.add( a_projectileEnergy, PoPI::Intids::neutron, m_neutronIndex, m_neutronUserParticleIndex, m_neutronMass, a_input, a_rng, 
            a_push_back, false );
    a_products.add( a_projectileEnergy, m_targetIntid, m_targetIndex, m_targetUserParticleIndex, 
            m_targetMass, a_input, a_rng, a_push_back, false );

    a_protare->sampleBranchingGammas( a_input, a_projectileEnergy, levelIndex, a_rng, a_push_back, a_products );

    return( true );
}

/* *********************************************************************************************************//**
 * This method samples a capture state level and returns an index into the a_protare->m_nuclideGammaBranchStateInfos vector
 * of the sampled state level.
 *
 * @param a_energy                  [in]    The neutron separation energy plus the projectile energy.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                                  An integer of the sampled state in a_protare->m_nuclideGammaBranchStateInfos.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::GRIN_captureToCompound::sampleCaptureLevel( ProtareSingle const *a_protare, double a_energy, 
                RNG && a_rng, bool a_checkEnergy ) const {

    if( a_checkEnergy ) {
        NuclideGammaBranchStateInfo *nuclideGammaBranchStateInfo = a_protare->nuclideGammaBranchStateInfos( )[m_index];
        if( nuclideGammaBranchStateInfo->nuclearLevelEnergy( ) < a_energy ) return( -1 );
    }

    double random = a_rng( );
    std::size_t index = 0;
    for( ; index < m_continuumIndices.m_levels.size( ) - 1; ++index ) {
        if( m_continuumIndices.m_summedProbabilities[index] >= random ) break;
    }
    return( m_continuumIndices.m_levels[index] );
}

/* *********************************************************************************************************//**
 * This method samples a capture state level and returns an index into the a_protare->m_nuclideGammaBranchStateInfos vector
 * of the sampled state level.
 *
 * @param a_protare                 [in]    The ProtareSingle this Reaction belongs to.
 * @param a_energy                  [in]    The neutron separation energy plus the projectile energy.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 *
 * @return                                  An integer of the sampled state in a_protare->m_nuclideGammaBranchStateInfos.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE int MCGIDI::GRIN_captureLevelProbability::sampleCaptureLevel( ProtareSingle const *a_protare, double a_energy, RNG && a_rng ) {

    double random = a_rng( );
    for( std::size_t index = 0; index < m_knownLevelsAndProbabilities.m_levels.size( ); ++index ) {
        if( m_knownLevelsAndProbabilities.m_summedProbabilities[index] >= random ) {
            return( m_knownLevelsAndProbabilities.m_levels[index] );
        }
    }

    for( std::size_t i1 = 0; i1 < m_captureToCompounds.size( ) - 1; ++i1 ) {
        GRIN_captureToCompound const *GRIN_captureToCompound1 = m_captureToCompounds[i1];

        int index = GRIN_captureToCompound1->sampleCaptureLevel( a_protare, a_energy, a_rng, true );
        if( index > -1 ) return( index );
    }

    return( m_captureToCompounds.back( )->sampleCaptureLevel( a_protare, a_energy, a_rng, false ) );
}

/* *********************************************************************************************************//**
 * This method adds a null product to *a_products*. When running in multi-group mode, a sampled reaction may be rejected if the threshold 
 * is in the multi-group that the projectile is in. If this happens, only null products should be returned. This type of behavior was need
 * in MCAPM but is probably not needed for MCGIDI.
 *
 * @param a_protare                 [in]    The Protare this Reaction belongs to.
 * @param a_projectileEnergy        [in]    The energy of the projectile.
 * @param a_input                   [in]    Sample options requested by user.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_products                [in]    The object to add all sampled products to.
 ***********************************************************************************************************/

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::Reaction::sampleNullProducts( Protare const &a_protare, double a_projectileEnergy, Sampling::Input &a_input, 
        RNG && a_rng, PUSHBACK && a_push_back, Sampling::ProductHandler &a_products ) {

    a_input.m_GRIN_intermediateResidual = -1;
    a_input.m_sampledType = Sampling::SampledType::uncorrelatedBody;
    a_input.m_dataInTargetFrame = false;
    a_input.m_frame = GIDI::Frame::lab;
    a_input.m_delayedNeutronIndex = -1;
    a_input.m_delayedNeutronDecayRate = 0.0;

    a_input.m_energyOut1 = a_projectileEnergy;
    a_input.m_mu = 1.0;
    a_input.m_phi = 0.0;

    a_products.add( a_projectileEnergy, a_protare.projectileIntid( ), a_protare.projectileIndex( ), a_protare.projectileUserIndex( ), a_protare.projectileMass( ), a_input, a_rng, a_push_back, false );
}

/* *********************************************************************************************************//**
 * Returns the weight for a project with energy *a_energy_in* to cause this reaction to emitted a particle of index
 * *a_index* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy. 
 *
 * @param a_index                   [in]    The index of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The cumulative multiplicity.
 * @param a_checkOrphanProducts     [in]    If true, associated orphan products are also sampled.
 *
 * @return                                  The weight that the particle is emitted into mu *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Reaction::angleBiasing( int a_index, double a_temperature, double a_energy_in, double a_mu_lab, double &a_energy_out, 
                RNG && a_rng, double *a_cumulative_weight, bool a_checkOrphanProducts ) const {

    double cumulative_weight1 = 0.0;
    if( a_cumulative_weight == nullptr ) a_cumulative_weight = &cumulative_weight1;
    double weight1 = 0.0;

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->angleBiasing( this, a_index, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
#else

    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        (*productIter)->angleBiasing( this, a_index, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
    }

    if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_index == m_neutronIndex ) ) {
        for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
            DelayedNeutron const *delayedNeutron1 = m_delayedNeutrons[i1];
            Product const &product = delayedNeutron1->product( );

            product.angleBiasing( this, a_index, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
        }
    }
#endif

    if( a_checkOrphanProducts ) {
        for( auto productIter = m_associatedOrphanProducts.begin( ); productIter != m_associatedOrphanProducts.end( ); ++productIter ) {
            (*productIter)->angleBiasing( this, a_index, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, 
                    *a_cumulative_weight );
        }
    }

    return( weight1 );
}

/* *********************************************************************************************************//**
 * Returns the weight for a project with energy *a_energy_in* to cause this reaction to emitted a particle of intid
 * *a_intid* at angle *a_mu_lab* as seen in the lab frame. If a particle is emitted, *a_energy_out* is its sampled outgoing energy.
 *
 * @param a_intid                   [in]    The intid of the particle to emit.
 * @param a_temperature             [in]    Specifies the temperature of the material.
 * @param a_energy_in               [in]    The energy of the incident particle.
 * @param a_mu_lab                  [in]    The desired mu in the lab frame for the emitted particle.
 * @param a_energy_out              [in]    The energy of the emitted outgoing particle.
 * @param a_rng                     [in]    The random number generator function that returns a double in the range [0, 1.0).
 * @param a_cumulative_weight       [in]    The cumulative multiplicity.
 * @param a_checkOrphanProducts     [in]    If true, associated orphan products are also sampled.
 *
 * @return                                  The weight that the particle is emitted into mu *a_mu_lab*.
 ***********************************************************************************************************/

template <typename RNG>
LUPI_HOST_DEVICE double MCGIDI::Reaction::angleBiasingViaIntid( int a_intid, double a_temperature, double a_energy_in, double a_mu_lab, double &a_energy_out,
                RNG && a_rng, double *a_cumulative_weight, bool a_checkOrphanProducts ) const {

    double cumulative_weight1 = 0.0;
    if( a_cumulative_weight == nullptr ) a_cumulative_weight = &cumulative_weight1;
    double weight1 = 0.0;

#ifdef MCGIDI_USE_OUTPUT_CHANNEL
    m_outputChannel->angleBiasingViaIntid( this, a_intid, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
#else

    for( auto productIter = m_products.begin( ); productIter != m_products.end( ); ++productIter ) {
        (*productIter)->angleBiasingViaIntid( this, a_intid, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
    }

    if( ( m_totalDelayedNeutronMultiplicity != nullptr ) && ( a_intid == PoPI::Intids::neutron ) ) {
        for( std::size_t i1 = 0; i1 < (std::size_t) m_delayedNeutrons.size( ); ++i1 ) {
            DelayedNeutron const *delayedNeutron1 = m_delayedNeutrons[i1];
            Product const &product = delayedNeutron1->product( );

            product.angleBiasingViaIntid( this, a_intid, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng, *a_cumulative_weight );
        }
    }
#endif

    if( a_checkOrphanProducts ) {
        for( auto productIter = m_associatedOrphanProducts.begin( ); productIter != m_associatedOrphanProducts.end( ); ++productIter ) {
            (*productIter)->angleBiasingViaIntid( this, a_intid, a_temperature, a_energy_in, a_mu_lab, weight1, a_energy_out, a_rng,
                    *a_cumulative_weight );
        }
    }

    return( weight1 );
}

// From file: MCGIDI_sampling.cpp

template <typename RNG, typename PUSHBACK>
LUPI_HOST_DEVICE void MCGIDI::Sampling::ProductHandler::add( double a_projectileEnergy, int a_productIntid, int a_productIndex, int a_userProductIndex, 
                double a_productMass, Input &a_input, RNG && a_rng, PUSHBACK && a_push_back, bool a_isPhoton ) {

    Product product;

    if( a_isPhoton && ( a_input.m_sampledType != SampledType::unspecified ) ) a_input.m_sampledType = SampledType::photon;

    product.m_sampledType = a_input.m_sampledType;
    product.m_isVelocity = a_input.wantVelocity( );
    product.m_productIntid = a_productIntid;
    product.m_productIndex = a_productIndex;
    product.m_userProductIndex = a_userProductIndex;
    product.m_numberOfDBRC_rejections = a_input.m_numberOfDBRC_rejections;
    product.m_productMass = a_productMass;

    product.m_delayedNeutronIndex = a_input.m_delayedNeutronIndex;
    product.m_delayedNeutronDecayRate = a_input.m_delayedNeutronDecayRate;
    product.m_birthTimeSec = 0.;
    if( product.m_delayedNeutronDecayRate > 0. ) {
        product.m_birthTimeSec = -log( a_rng( ) ) / product.m_delayedNeutronDecayRate;
    }

    if( a_input.m_sampledType == SampledType::unspecified ) {
        product.m_kineticEnergy = 0.0;
        product.m_px_vx = 0.0;
        product.m_py_vy = 0.0;
        product.m_pz_vz = 0.0; }
    else if( a_input.m_sampledType == SampledType::uncorrelatedBody ) {
        if( a_input.m_frame == GIDI::Frame::centerOfMass ) {
            a_input.m_frame = GIDI::Frame::lab;

            double massRatio = a_input.m_projectileMass + a_input.m_targetMass;
            massRatio = a_input.m_projectileMass * a_productMass / ( massRatio * massRatio );
            double modifiedProjectileEnergy = massRatio * a_projectileEnergy;

            double sqrtModifiedProjectileEnergy = sqrt( modifiedProjectileEnergy );
            double sqrtEnergyOut_com = a_input.m_mu * sqrt( a_input.m_energyOut1 );

            a_input.m_energyOut1 += modifiedProjectileEnergy + 2. * sqrtModifiedProjectileEnergy * sqrtEnergyOut_com;
            if( a_input.m_energyOut1 != 0 ) a_input.m_mu = ( sqrtModifiedProjectileEnergy + sqrtEnergyOut_com ) / sqrt( a_input.m_energyOut1 );
        }

        product.m_kineticEnergy = a_input.m_energyOut1;

        double p_v = sqrt( a_input.m_energyOut1 * ( a_input.m_energyOut1 + 2. * a_productMass ) );
        if( product.m_isVelocity ) p_v *= MCGIDI_speedOfLight_cm_sec / ( a_input.m_energyOut1 + a_productMass );

        product.m_pz_vz = p_v * a_input.m_mu;
        p_v *= sqrt( 1. - a_input.m_mu * a_input.m_mu );
        product.m_px_vx = p_v * sin( a_input.m_phi );
        product.m_py_vy = p_v * cos( a_input.m_phi ); }
    else if( a_input.m_sampledType == SampledType::firstTwoBody ) {
        product.m_kineticEnergy = a_input.m_energyOut1;
        product.m_px_vx = a_input.m_px_vx1;
        product.m_py_vy = a_input.m_py_vy1;
        product.m_pz_vz = a_input.m_pz_vz1;
        a_input.m_sampledType = SampledType::secondTwoBody; }
    else if( a_input.m_sampledType == SampledType::secondTwoBody ) {
        product.m_kineticEnergy = a_input.m_energyOut2;
        product.m_px_vx = a_input.m_px_vx2;
        product.m_py_vy = a_input.m_py_vy2;
        product.m_pz_vz = a_input.m_pz_vz2; }
    else if( a_input.m_sampledType == SampledType::photon ) {
        product.m_kineticEnergy = a_input.m_energyOut1;

        double pz_vz_factor = a_input.m_energyOut1;
        if( product.m_isVelocity ) pz_vz_factor = MCGIDI_speedOfLight_cm_sec;
        product.m_pz_vz = a_input.m_mu * pz_vz_factor;

        double v_perp = sqrt( 1.0 - a_input.m_mu * a_input.m_mu ) * pz_vz_factor;
        product.m_px_vx = cos( a_input.m_phi ) * v_perp;
        product.m_py_vy = sin( a_input.m_phi ) * v_perp; }
    else {
        product.m_kineticEnergy = a_input.m_energyOut2;
        product.m_px_vx = a_input.m_px_vx2;
        product.m_py_vy = a_input.m_py_vy2;
        product.m_pz_vz = a_input.m_pz_vz2;
    }

    if( a_input.m_dataInTargetFrame && ( a_input.m_sampledType != SampledType::photon ) ) upScatterModelABoostParticle( a_input, a_rng, product );

    a_push_back( product );
}

#endif      // End of MCGIDI_headerSource_hpp_included
