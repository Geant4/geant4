/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_sampling_hpp_included
#define MCGIDI_sampling_hpp_included 1

#include <LUPI_declareMacro.hpp>
#include <MCGIDI_vector.hpp>
#include <MCGIDI_string.hpp>

namespace MCGIDI {

/*
============================================================
======================= DomainHash =========================
============================================================
*/
class DomainHash {

    private:
        std::size_t m_bins;                                                 /**< The number of bins for the hash. */
        double m_domainMin;                                                 /**< The minimum domain value for the hash. */
        double m_domainMax;                                                 /**< The maximum domain value for the hash. */
        double m_u_domainMin;                                               /**< The log of m_domainMin ). */
        double m_u_domainMax;                                               /**< The log of m_domainMax ). */
        double m_inverse_du;                                                /**< The value *m_bins* / ( *m_u_domainMax* - *m_u_domainMin* ). */

    public:
        LUPI_HOST_DEVICE DomainHash( );
        LUPI_HOST_DEVICE DomainHash( std::size_t a_bins, double a_domainMin, double a_domainMax );
        LUPI_HOST_DEVICE DomainHash( DomainHash const &a_domainHash );

        LUPI_HOST_DEVICE std::size_t bins( ) const { return( m_bins ); }                     /**< Returns the value of the **m_bins**. */
        LUPI_HOST_DEVICE double domainMin( ) const { return( m_domainMin ); }        /**< Returns the value of the **m_domainMax**. */
        LUPI_HOST_DEVICE double domainMax( ) const { return( m_domainMax ); }        /**< Returns the value of the **m_domainMax**. */
        LUPI_HOST_DEVICE double u_domainMin( ) const { return( m_u_domainMin ); }    /**< Returns the value of the **m_u_domainMin**. */
        LUPI_HOST_DEVICE double u_domainMax( ) const { return( m_u_domainMax ); }    /**< Returns the value of the **m_u_domainMax**. */
        LUPI_HOST_DEVICE double inverse_du( ) const { return( m_inverse_du ); }      /**< Returns the value of the **m_inverse_du**. */

        LUPI_HOST_DEVICE std::size_t index( double a_domain ) const ;
        LUPI_HOST_DEVICE Vector<std::size_t > map( Vector<double> const &a_domainValues ) const ;

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        LUPI_HOST void print( bool a_printValues ) const ;
};

namespace Sampling {

enum class SampledType { firstTwoBody, secondTwoBody, uncorrelatedBody, unspecified, photon };
class ProductHandler;

LUPI_HOST_DEVICE std::size_t evaluationForHashIndex( std::size_t a_hashIndex, Vector<std::size_t> const &a_hashIndices, double a_energy, 
                Vector<double> const &a_energies, double *a_energyFraction );

namespace Upscatter {

    enum class Model { none, A, B, BSnLimits, DBRC };

/*
============================================================
===================== ModelDBRC_data =======================
============================================================
*/

class ModelDBRC_data {

    public:
        double m_neutronMass;                   /**< The mass of the neutron. */
        double m_targetMass;                    /**< The mass of the target. */
        Vector<double> m_energies;              /**< The energy grid for the cross section. */
        Vector<double> m_crossSections;         /**< The cross sections corresponding to the energy grid. */
        Vector<std::size_t> m_hashIndices;      /**< The indicies for the energy hash function. */
        MCGIDI::DomainHash m_domainHash;        /**< The hash "function". */

    public:
        LUPI_HOST_DEVICE ModelDBRC_data( );
        LUPI_HOST ModelDBRC_data( double a_neutronMass, double a_targetMass, Vector<double> const &a_energies, Vector<double> const &a_crossSections,
                DomainHash const &a_domainHash );
        LUPI_HOST_DEVICE ~ModelDBRC_data( );

        LUPI_HOST_DEVICE double evaluate( double a_energy );
        LUPI_HOST_DEVICE double targetThermalSpeed( double a_temperature );
        LUPI_HOST_DEVICE double crossSectionMax( double a_energy, double a_targetThermalSpeed );

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

LUPI_HOST_DEVICE ModelDBRC_data *serializeModelDBRC_data( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode, ModelDBRC_data *a_modelDBRC_data );

}           // End of namespace Upscatter.

/*
============================================================
================ ClientRandomNumberGenerator ===============
============================================================
*/
class ClientRandomNumberGenerator {
    private:
        double (*m_generator)( void * );                    /**< User supplied generator. */
        void *m_state;                                      /**< User supplied state. */

    public:
        LUPI_HOST_DEVICE ClientRandomNumberGenerator( double (*a_generator)( void * ), void *a_state );

        LUPI_HOST_DEVICE double (*generator( ))( void * ) { return( m_generator ); }
        LUPI_HOST_DEVICE void *state( ) { return( m_state ); }
        LUPI_HOST_DEVICE double Double( ) { return( m_generator( m_state ) ); }

// The following are deprecated.
        LUPI_HOST_DEVICE double (*rng( ))( void * ) { return( generator( ) ); }
        LUPI_HOST_DEVICE void *rngState( ) { return( state( ) ); }
        LUPI_HOST_DEVICE double dRng( ) { return( Double( ) ); }
};

/*
============================================================
=================== Client Code RNG Data ===================
============================================================
*/
class ClientCodeRNGData : public ClientRandomNumberGenerator {

    public:
        LUPI_HOST_DEVICE ClientCodeRNGData( double (*a_generator)( void * ), void *a_state );
};

/*
============================================================
=========================== Input ==========================
============================================================
*/

class Input {

    friend ProtareSingle;
    friend Reaction;
    friend MCGIDI::Sampling::ProductHandler;
    friend GRIN_capture;
    friend GRIN_inelastic;

    private:
        bool m_wantVelocity = true ;                        /**< See member m_isVelocity in class Product for meaning. This is user input. */

        bool m_dataInTargetFrame = false;                   /**< **True** if the data are in the target's frame and **false** otherwise. */
        double m_modelTemperature = 0.0;                    /**< The temperature used when sampling product data. For example, in upscatter model A the projectile is boosted into a sampled target's frame and the modelled temperature is 0.0. */
        double m_modelEnergy = 0.0;                         /**< The projectile energy used when sampling product data (see comment for member **m_modelTemperature**. */

        SampledType m_sampledType = SampledType::uncorrelatedBody;  /**< For internal use only. Set by distributions and used in the method **MCGIDI::Sampling::ProductHandler::add**. */

                                                    // The next 2 members are set by the user via the setTemperatureAndEnergy method.
        double m_temperature = 0.0;                         /**< The temperature of the material. This member is set by the user. */
        double m_energy = 0.0;                              /**< The energy of the projectile. This member is set by the user. */

    public:

        Upscatter::Model m_upscatterModel = Upscatter::Model::none; /**< The upscatter model to use when sampling a target's velocity. */


                                                    // The rest of the members are set by MCGIDI methods.
                                                    // These five are used for upscatter model A and the last 4 also used by model B.
        double m_projectileBeta = 0.0;                      /**< The beta = speed / c of the projectile. */
        double m_muLab = 0.0;                               /**< The cosine of the angle between the projectile's and the sampled target's velocities. */
        double m_targetBeta = 0.0;                          /**< The beta = speed / c of the target. */
        double m_relativeBeta = 0.0;                        /**< The beta = speed / c of the relative speed between the projectile and the target. */

        Reaction const *m_reaction = nullptr;               /**< The current reaction whose products are being sampled. */

        double m_projectileMass = 0.0;                      /**< The mass of the projectile. */
        double m_targetMass = 0.0;                          /**< The mass of the target. */

        GIDI::Frame m_frame = GIDI::Frame::lab;             /**< The frame the product data are returned in. */
        int m_numberOfDBRC_rejections = 0;                  /**< For the DBRC upscattering model, this is the number of rejections + 1 per product sample. */

        double m_mu = 0.0;                                  /**< The sampled mu = cos( theta ) for the product. */
        double m_phi = 0.0;                                 /**< The sampled phi for the product. */

        double m_energyOut1 = 0.0;                          /**< The sampled energy of the product. */
        double m_px_vx1 = 0.0;                              /**< Variable used for two-body sampling. */
        double m_py_vy1 = 0.0;                              /**< Variable used for two-body sampling. */
        double m_pz_vz1 = 0.0;                              /**< Variable used for two-body sampling. */

        double m_energyOut2 = 0.0;                          /**< The sampled energy of the second product for a two-body interaction. */
        double m_px_vx2 = 0.0;                              /**< Variable used for two-body sampling. */
        double m_py_vy2 = 0.0;                              /**< Variable used for two-body sampling. */
        double m_pz_vz2 = 0.0;                              /**< Variable used for two-body sampling. */

        int m_delayedNeutronIndex = -1;                     /**< If the product is a delayed neutron, this is its index. */
        double m_delayedNeutronDecayRate = 0.0;             /**< If the product is a delayed neutron, this is its decay rate. */

        int m_GRIN_intermediateResidual = -1;               /**< For special GRIN product sampling, this is the GNDS intid of the intermediate residual. */

        LUPI_HOST_DEVICE Input( bool a_wantVelocity, Upscatter::Model a_upscatterModel );

        LUPI_HOST_DEVICE bool wantVelocity( ) const { return( m_wantVelocity ); }                       /**< Returns the value of the *m_wantVelocity* member. */
        LUPI_HOST_DEVICE double temperature( ) const { return( m_temperature ); }                       /**< Returns the value of the *m_temperature* member. */
        LUPI_HOST_DEVICE double energy( ) const { return( m_energy ); }                                 /**< Returns the value of the *m_energy* member. */
        LUPI_HOST_DEVICE void setTemperatureAndEnergy( double a_temperature, double a_energy );

        LUPI_HOST_DEVICE bool dataInTargetFrame( ) const { return( m_dataInTargetFrame ); }             /**< Returns the value of the *m_dataInTargetFrame*. */
        LUPI_HOST_DEVICE double modelTemperature( ) const { return( m_modelTemperature ); }             /**< Returns the value of the *m_dataInTargetFrame* member. */
        LUPI_HOST_DEVICE double modelEnergy( ) const { return( m_modelEnergy ); }                       /**< Returns the value of the *m_modelEnergy* member. */

        SampledType sampledType( ) const { return( m_sampledType ); }                                   /**< Returns the value of the *m_sampledType* member. */
        LUPI_HOST_DEVICE void setSampledType( SampledType a_sampledType ) { m_sampledType = a_sampledType; }             /**< Sets the member *m_sampledType* to *a_sampledType*. */
};

/*
============================================================
========================== Product =========================
============================================================
*/
class Product {

    public:
        SampledType m_sampledType;
        bool m_isVelocity;                      /**< If true, m_px_vx, m_py_vy and m_pz_vz are velocities otherwise momenta. */
        int m_productIntid;                     /**< The intid of the sampled product. */
        int m_productIndex;                     /**< The index of the sampled product. */
        int m_userProductIndex;                 /**< The user particle index of the sampled product. */
        int m_numberOfDBRC_rejections;          /**< For the DBRC upscattering model, this is the number of rejections + 1 per product sample. */
        double m_productMass;                   /**< The mass of the sampled product. */
        double m_kineticEnergy;                 /**< The kinetic energy of the sampled product. */
        double m_px_vx;                         /**< The velocity or momentum along the x-axis of the sampled product. */
        double m_py_vy;                         /**< The velocity or momentum along the y-axis of the sampled product. */
        double m_pz_vz;                         /**< The velocity or momentum along the z-axis of the sampled product. The z-axis is along the direction of the projectile's velolcity. */
        int m_delayedNeutronIndex;              /**< If the product is a delayed neutron, this is its index. */
        double m_delayedNeutronDecayRate;       /**< If the product is a delayed neutron, this is its decay rate. */
        double m_birthTimeSec;                  /**< Some products, like delayed fission neutrons, are to appear (be born) later. This is the time in seconds that such a particle should be born since the interaction. */
};

/*
============================================================
====================== ProductHandler ======================
============================================================
*/
class ProductHandler {

    public:
        LUPI_HOST_DEVICE ProductHandler( ) {}
        LUPI_HOST_DEVICE ~ProductHandler( ) {}

        template <typename RNG, typename PUSHBACK>
        LUPI_HOST_DEVICE void add( double a_projectileEnergy, int a_productIntid, int a_productIndex, int a_userProductIndex, double a_productMass, Input &a_input, 
                RNG && a_rng, PUSHBACK && push_back, bool isPhoton );
};

/*
============================================================
================ StdVectorProductHandler ===================
============================================================
*/
#ifdef __CUDACC__

#define MCGIDI_CUDACC_numberOfProducts 1000

class StdVectorProductHandler : public ProductHandler {

    private:
        std::size_t m_size;
        Product m_products[1024];

    public:
        LUPI_HOST_DEVICE StdVectorProductHandler( ) : m_size( 0 ) { }
        LUPI_HOST_DEVICE ~StdVectorProductHandler( ) { }

        LUPI_HOST_DEVICE std::size_t size( ) { return( m_size ); }
        LUPI_HOST_DEVICE Product &operator[]( long a_index ) { return( m_products[a_index] ); }
        LUPI_HOST_DEVICE void push_back( Product &a_product ) {
            if( m_size < MCGIDI_CUDACC_numberOfProducts ) {
                m_products[m_size] = a_product;
                ++m_size;
            }
        }
        LUPI_HOST_DEVICE void clear( ) { m_size = 0; }
};

#else
class StdVectorProductHandler : public ProductHandler {

    private:
        std::vector<Product> m_products;            /**< The list of products sampled. */

    public:
        LUPI_HOST_DEVICE StdVectorProductHandler( ) : m_products( ) { }
        LUPI_HOST_DEVICE ~StdVectorProductHandler( ) { }

        LUPI_HOST_DEVICE std::size_t size( ) { return( m_products.size( ) ); }
        LUPI_HOST_DEVICE Product &operator[]( std::size_t a_index ) { return( m_products[a_index] ); }
        LUPI_HOST_DEVICE std::vector<Product> &products( ) { return( m_products ); }
        LUPI_HOST_DEVICE void push_back( Product &a_product ) { m_products.push_back( a_product ); }
        LUPI_HOST_DEVICE void clear( ) { m_products.clear( ); }
};
#endif

/*
============================================================
============== MCGIDIVectorProductHandler ==================
============================================================
*/
class MCGIDIVectorProductHandler : public ProductHandler {

    private:
        Vector<Product> m_products;             /**< The list of products sampled. */

    public:
        LUPI_HOST_DEVICE MCGIDIVectorProductHandler( std::size_t a_size = 20 ) :
                m_products( ) {

            m_products.reserve( a_size );
        }
        LUPI_HOST_DEVICE ~MCGIDIVectorProductHandler( ) {}

        LUPI_HOST_DEVICE std::size_t size( ) { return( m_products.size( ) ); }
        LUPI_HOST_DEVICE Product const &operator[]( std::size_t a_index ) const { return( m_products[a_index] ); }
        LUPI_HOST_DEVICE Vector<Product> const &products( ) const { return( m_products ); }
        LUPI_HOST_DEVICE void push_back( Product &a_product ) { m_products.push_back( a_product ); }
        LUPI_HOST_DEVICE void clear( ) { m_products.clear( ); }
};

}       // End of namespace Sampling.

}       // End of namespace MCGIDI.

#endif      // End of MCGIDI_sampling_hpp_included
