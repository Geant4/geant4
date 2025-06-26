/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef MCGIDI_distributions_hpp_included
#define MCGIDI_distributions_hpp_included 1

#include <LUPI_declareMacro.hpp>

namespace MCGIDI {

namespace Distributions {

enum class Type { none, unspecified, angularTwoBody, KalbachMann, uncorrelated, branching3d, energyAngularMC, angularEnergyMC, 
        coherentPhotoAtomicScattering, incoherentPhotoAtomicScattering, incoherentPhotoAtomicScatteringElectron, incoherentBoundToFreePhotoAtomicScattering, pairProductionGamma,
        coherentElasticTNSL, incoherentElasticTNSL };

/*
============================================================
======================= Distribution =======================
============================================================
*/
class Distribution {

    private:
        Type m_type;                                    /**< Specifies the Type of the distribution. */
        GIDI::Frame m_productFrame;                     /**< Specifies the frame the product data are given in. */
        double m_projectileMass;                        /**< The mass of the projectile. */
        double m_targetMass;                            /**< The mass of the target. */
        double m_productMass;                           /**< The mass of the first product. */

    public:
        LUPI_HOST_DEVICE Distribution( );
        LUPI_HOST Distribution( Type a_type, GIDI::Distributions::Distribution const &a_distribution, SetupInfo &a_setupInfo );
        LUPI_HOST Distribution( Type a_type, GIDI::Frame a_productFrame, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION ~Distribution( ) MCGIDI_TRUE_VIRTUAL;

        LUPI_HOST_DEVICE Type type( ) const { return( m_type ); }                            /**< Returns the value of the **m_type**. */
        LUPI_HOST_DEVICE GIDI::Frame productFrame( ) const { return( m_productFrame ); }     /**< Returns the value of the **m_productFrame**. */

        LUPI_HOST_DEVICE double projectileMass( ) const { return( m_projectileMass ); }                  /**< Returns the value of the **m_projectileMass**. */
        LUPI_HOST_DEVICE double targetMass( ) const { return( m_targetMass ); }                          /**< Returns the value of the **m_targetMass**. */
        LUPI_HOST_DEVICE double productMass( ) const { return( m_productMass ); }                        /**< Returns the value of the **m_productMass**. */

        LUPI_HOST void setModelDBRC_data( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data );

        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const MCGIDI_TRUE_VIRTUAL;
        template <typename RNG>
        LUPI_HOST_DEVICE MCGIDI_VIRTUAL_FUNCTION double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const MCGIDI_TRUE_VIRTUAL;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== AngularTwoBody ======================
============================================================
*/
class AngularTwoBody : public Distribution {

    private:
        double m_residualMass;                                          /**< The mass of the second product (often the  residual). */
        double m_Q;                                                     /**< FIX ME. */
        double m_twoBodyThreshold;                                      /**< This is the T_1 value needed to do two-body kinematics (i.e., in the equation (K_{com,3_4} = m_2 * (K_1 - T_1) / (m_1 + m_2)). */
        bool m_Upscatter;                                               /**< Set to true if reaction is elastic which is the only reaction upscatter Model B is applied to. */
        Probabilities::ProbabilityBase2d_d1 *m_angular;                 /**< The 2d angular probability. */
        Sampling::Upscatter::ModelDBRC_data *m_modelDBRC_data;          /**< The cross section and other data needed for neutron elastic upscatter model DBRC. */

        template <typename RNG>
        LUPI_HOST_DEVICE bool upscatterModelB( double a_kineticLab, Sampling::Input &a_input, RNG && a_rng ) const ;

    public:
        LUPI_HOST_DEVICE AngularTwoBody( );
        LUPI_HOST AngularTwoBody( GIDI::Distributions::AngularTwoBody const &a_angularTwoBody, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~AngularTwoBody( );

        LUPI_HOST_DEVICE double residualMass( ) const { return( m_residualMass ); }                      /**< Returns the value of the **m_residualMass**. */
        LUPI_HOST_DEVICE double Q( ) const { return( m_Q ); }                                            /**< Returns the value of the **m_Q**. */
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *angular( ) const { return( m_angular ); }  /**< Returns the value of the **m_angular**. */
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
        LUPI_HOST_DEVICE bool Upscatter( ) const { return( m_Upscatter ); }                              /**< Returns the value of the **m_Upscatter**. */
        LUPI_HOST void setModelDBRC_data2( Sampling::Upscatter::ModelDBRC_data *a_modelDBRC_data );
};

/*
============================================================
======================= Uncorrelated =======================
============================================================
*/
class Uncorrelated : public Distribution {

    private:
        Probabilities::ProbabilityBase2d_d1 *m_angular;         /**< The angular probability P(mu|E). */
        Probabilities::ProbabilityBase2d *m_energy;             /**< The energy probability P(E'|E). */
        
    public:
        LUPI_HOST_DEVICE Uncorrelated( );
        LUPI_HOST Uncorrelated( GIDI::Distributions::Uncorrelated const &a_uncorrelated, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~Uncorrelated( );

        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *angular( ) const { return( m_angular ); }  /**< Returns the value of the **m_angular**. */
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d *energy( ) const { return( m_energy ); }       /**< Returns the value of the **m_energy**. */
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== Branching3d =======================
============================================================
*/
class Branching3d : public Distribution {

    private:
        int m_initialStateIndex;

    public:
        LUPI_HOST_DEVICE Branching3d( );
        LUPI_HOST Branching3d( GIDI::Distributions::Branching3d const &a_branching3d, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~Branching3d( );

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab,
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== EnergyAngularMC =====================
============================================================
*/
class EnergyAngularMC : public Distribution {

    private:
        Probabilities::ProbabilityBase2d_d1 *m_energy;              /**< The energy probability P(E'|E). */
        Probabilities::ProbabilityBase3d *m_angularGivenEnergy;     /**< The angular probability given E', P(mu|E,E'). */
        
    public:
        LUPI_HOST_DEVICE EnergyAngularMC( );
        LUPI_HOST EnergyAngularMC( GIDI::Distributions::EnergyAngularMC const &a_energyAngularMC, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~EnergyAngularMC( );

        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *energy( ) const { return( m_energy ); }       /**< Returns the value of the **m_energy**. */
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase3d *angularGivenEnergy( ) const { return( m_angularGivenEnergy ); }   /**< Returns the value of the **m_angularGivenEnergy**. */
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
====================== AngularEnergyMC =====================
============================================================
*/
class AngularEnergyMC : public Distribution {

    private:
        Probabilities::ProbabilityBase2d_d1 *m_angular;             /**< The angular probability P(mu|E). */
        Probabilities::ProbabilityBase3d *m_energyGivenAngular;     /**< The energy probability P(E'|E,mu). */
        
    public:
        LUPI_HOST_DEVICE AngularEnergyMC( );
        LUPI_HOST AngularEnergyMC( GIDI::Distributions::AngularEnergyMC const &a_angularEnergyMC, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~AngularEnergyMC( );

        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *angular( ) const { return( m_angular ); }     /**< Returns the value of the **m_angular**. */
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase3d *energyGivenAngular( ) const { return( m_energyGivenAngular ); }   /**< Returns the value of the **m_energyGivenAngular**. */
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
======================== KalbachMann =======================
============================================================
*/
class KalbachMann : public Distribution {

    private:
        double m_energyToMeVFactor;                                 /**< The factor that converts energies to MeV. */
        double m_eb_massFactor;                                     /**< FIX ME */
        Probabilities::ProbabilityBase2d_d1 *m_f;                   /**< The energy probability P(E'|E). */
        Functions::Function2d *m_r;                                 /**< The Kalbach-Mann r(E,E') function. */
        Functions::Function2d *m_a;                                 /**< The Kalbach-Mann a(E,E') function. */

    public:
        LUPI_HOST_DEVICE KalbachMann( );
        LUPI_HOST KalbachMann( GIDI::Distributions::KalbachMann const &a_KalbachMann, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~KalbachMann( );

        LUPI_HOST_DEVICE double energyToMeVFactor( ) const { return( m_energyToMeVFactor ); }    /**< Returns the value of the **m_energyToMeVFactor**. */
        LUPI_HOST_DEVICE double eb_massFactor( ) const { return( m_eb_massFactor ); }            /**< Returns the value of the **m_eb_massFactor**. */
        LUPI_HOST_DEVICE Probabilities::ProbabilityBase2d_d1 *f( ) const { return( m_f ); }      /**< Returns the value of the **m_f**. */
        LUPI_HOST_DEVICE Functions::Function2d *r( ) const { return( m_r ); }                    /**< Returns the value of the **m_r**. */
        LUPI_HOST_DEVICE Functions::Function2d *a( ) const { return( m_a ); }                    /**< Returns the value of the **m_a**. */
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        LUPI_HOST_DEVICE double evaluate( double E_in_lab, double E_out, double mu );
};

/*
============================================================
=============== CoherentPhotoAtomicScattering ==============
============================================================
*/
class CoherentPhotoAtomicScattering : public Distribution {

    private:
        bool m_anomalousDataPresent;                                /**< FIX ME */
        Vector<double> m_energies;                                  /**< FIX ME */
        Vector<double> m_formFactor;                                /**< FIX ME */
        Vector<double> m_a;                                         /**< FIX ME */
        Vector<double> m_integratedFormFactor;                      /**< FIX ME */
        Vector<double> m_integratedFormFactorSquared;               /**< FIX ME */
        Vector<double> m_probabilityNorm1_1;                        /**< FIX ME */
        Vector<double> m_probabilityNorm1_3;                        /**< FIX ME */
        Vector<double> m_probabilityNorm1_5;                        /**< FIX ME */
        Vector<double> m_probabilityNorm2_1;                        /**< FIX ME */
        Vector<double> m_probabilityNorm2_3;                        /**< FIX ME */
        Vector<double> m_probabilityNorm2_5;                        /**< FIX ME */
        Functions::Function1d_d1 *m_realAnomalousFactor;            /**< The real part of the anomalous scattering factor. */
        Functions::Function1d_d1 *m_imaginaryAnomalousFactor;       /**< The imaginary part of the anomalous scattering factor. */

        LUPI_HOST_DEVICE double Z_a( double a_Z, double a_a ) const ;

    public:
        LUPI_HOST_DEVICE CoherentPhotoAtomicScattering( );
        LUPI_HOST CoherentPhotoAtomicScattering( GIDI::Distributions::CoherentPhotoAtomicScattering const &a_coherentPhotoAtomicScattering, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~CoherentPhotoAtomicScattering( );

        LUPI_HOST_DEVICE double evaluate( double a_energyIn, double a_mu ) const ;
        LUPI_HOST_DEVICE double evaluateFormFactor( double a_energyIn, double a_mu ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;

        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
============== IncoherentPhotoAtomicScattering =============
============================================================
*/
class IncoherentPhotoAtomicScattering : public Distribution {

    private:
        Vector<double> m_energies;                                  /**< FIX ME */
        Vector<double> m_scatteringFactor;                          /**< FIX ME */
        Vector<double> m_a;                                         /**< FIX ME */

    public:
        LUPI_HOST_DEVICE IncoherentPhotoAtomicScattering( );
        LUPI_HOST IncoherentPhotoAtomicScattering( GIDI::Distributions::IncoherentPhotoAtomicScattering const &a_incoherentPhotoAtomicScattering, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~IncoherentPhotoAtomicScattering( );

        LUPI_HOST_DEVICE double energyRatio( double a_energyIn, double a_mu ) const ;
        LUPI_HOST_DEVICE double evaluateKleinNishina( double a_energyIn, double a_mu ) const ;
        LUPI_HOST_DEVICE double evaluateScatteringFactor( double a_X ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
/*
        LUPI_HOST_DEVICE double evaluate( double E_in_lab, double mu );
*/
};

/*
=======================================================================
============== IncoherentBoundToFreePhotoAtomicScattering =============
=======================================================================
*/
class IncoherentBoundToFreePhotoAtomicScattering : public Distribution {

    private:
        //Vector<double> m_energies;
        //Vector<double> m_ComptonProfile;
        Vector<double> m_occupationNumber;
        //Vector<double> m_a;
        Vector<double> m_pz;
        double m_bindingEnergy;

    public:
        LUPI_HOST_DEVICE IncoherentBoundToFreePhotoAtomicScattering( );
        LUPI_HOST IncoherentBoundToFreePhotoAtomicScattering( GIDI::Distributions::IncoherentBoundToFreePhotoAtomicScattering const &a_incoherentPhotoAtomicScattering, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~IncoherentBoundToFreePhotoAtomicScattering( );
        LUPI_HOST_DEVICE double energyRatio( double a_energyIn, double a_mu ) const ;
        LUPI_HOST_DEVICE double evaluateKleinNishina( double a_energyIn, double a_mu ) const ;
        LUPI_HOST_DEVICE double evaluateOccupationNumber( double a_X, double a_mu ) const;
        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab,
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========== IncoherentPhotoAtomicScatteringElectron =========
============================================================
*/
class IncoherentPhotoAtomicScatteringElectron : public Distribution {

    public:
        LUPI_HOST_DEVICE IncoherentPhotoAtomicScatteringElectron( );
        LUPI_HOST IncoherentPhotoAtomicScatteringElectron( SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~IncoherentPhotoAtomicScatteringElectron( );

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_energy, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab,
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
==================== PairProductionGamma ===================
============================================================
*/
class PairProductionGamma : public Distribution {

    private:
        bool m_firstSampled;                                    /**< When sampling photons for pair production, the photons must be emitted back-to-back. The flag help do this. */

    public:
        LUPI_HOST_DEVICE PairProductionGamma( );
        LUPI_HOST PairProductionGamma( SetupInfo &a_setupInfo, bool a_firstSampled );
        LUPI_HOST_DEVICE ~PairProductionGamma( );

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
==================== CoherentElasticTNSL ===================
============================================================
*/
class CoherentElasticTNSL : public Distribution {

    private:
        Interpolation m_temperatureInterpolation;
        Vector<double> m_temperatures;
        Vector<double> m_energies;
        Vector<double> m_S_table;

    public:
        LUPI_HOST_DEVICE CoherentElasticTNSL( );
        LUPI_HOST CoherentElasticTNSL( GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::CoherentElastic const *a_coherentElasticTNSL, 
                SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~CoherentElasticTNSL( ) {}

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_energy, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
==================== IncoherentElasticTNSL ===================
============================================================
*/
class IncoherentElasticTNSL : public Distribution {

    private:
        double m_temperatureToMeV_K;
        Functions::Function1d_d1 *m_DebyeWallerIntegral;

    public:
        LUPI_HOST_DEVICE IncoherentElasticTNSL( );
        LUPI_HOST IncoherentElasticTNSL( GIDI::DoubleDifferentialCrossSection::n_ThermalNeutronScatteringLaw::IncoherentElastic const *a_incoherentElasticTNSL, 
                SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~IncoherentElasticTNSL( ) {}

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_energy, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab,
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );

        Functions::Function1d       *DebyeWallerIntegral( )       { return( m_DebyeWallerIntegral ); }
        Functions::Function1d const *DebyeWallerIntegral( ) const { return( m_DebyeWallerIntegral ); }
};

/*
============================================================
======================= Unspecified ========================
============================================================
*/
class Unspecified : public Distribution {

    public:
        LUPI_HOST_DEVICE Unspecified( );
        LUPI_HOST Unspecified( GIDI::Distributions::Distribution const &a_distribution, SetupInfo &a_setupInfo );
        LUPI_HOST_DEVICE ~Unspecified( );

        template <typename RNG>
        LUPI_HOST_DEVICE void sample( double a_X, Sampling::Input &a_input, RNG && a_rng ) const ;
        template <typename RNG>
        LUPI_HOST_DEVICE double angleBiasing( Reaction const *a_reaction, double a_temperature, double a_energy_in, double a_mu_lab, 
                RNG && a_rng, double &a_energy_out ) const ;
        LUPI_HOST_DEVICE void serialize( LUPI::DataBuffer &a_buffer, LUPI::DataBuffer::Mode a_mode );
};

/*
============================================================
========================== Others ==========================
============================================================
*/
LUPI_HOST Distribution *parseGIDI( GIDI::Suite const &a_distribution, SetupInfo &a_setupInfo, Transporting::MC const &a_settings );
LUPI_HOST_DEVICE Type DistributionType( Distribution const *a_distribution );

}

}

#endif      // End of MCGIDI_distributions_hpp_included
