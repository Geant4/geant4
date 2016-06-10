/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef MCGIDI_h_included
#define MCGIDI_h_included

#define MCGIDI_VERSION_MAJOR 1
#define MCGIDI_VERSION_MINOR 0
#define MCGIDI_VERSION_PATCHLEVEL 0

#ifdef WIN32
#define M_PI 3.141592653589793238463
#endif

#include <GIDI_settings.hh>
#include <map>
#include <vector>

#include <statusMessageReporting.h>
#include <ptwXY.h>
#include <xDataTOM.h>

#include "MCGIDI_mass.h"
#include "MCGIDI_map.h"

/* Disable Effective C++ warnings in GIDI code. */
#if __INTEL_COMPILER > 1399
#pragma warning( disable:2021 )
#pragma warning( disable:593 )
#pragma warning( disable:111 )
#elif __INTEL_COMPILER > 1199
#pragma warning( disable:2304 )
#endif

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

typedef struct MCGIDI_GammaBranching_s MCGIDI_GammaBranching;
typedef struct MCGIDI_POP_s MCGIDI_POP;
typedef struct MCGIDI_POPs_s MCGIDI_POPs;
typedef struct MCGIDI_particle_s MCGIDI_particle;
typedef struct MCGIDI_target_s MCGIDI_target;
typedef struct MCGIDI_target_heated_info_s MCGIDI_target_heated_info;
typedef struct MCGIDI_target_heated_sorted_s MCGIDI_target_heated_sorted;
typedef struct MCGIDI_target_heated_s MCGIDI_target_heated;
typedef struct MCGIDI_reaction_s MCGIDI_reaction;
typedef struct MCGIDI_outputChannel_s MCGIDI_outputChannel;
typedef struct MCGIDI_product_s MCGIDI_product;
typedef struct MCGIDI_distribution_s MCGIDI_distribution;
typedef struct MCGIDI_KalbachMann_s MCGIDI_KalbachMann;
typedef struct MCGIDI_KalbachMann_ras_s MCGIDI_KalbachMann_ras;
typedef struct MCGIDI_pdfOfX_s MCGIDI_pdfOfX;
typedef struct MCGIDI_pdfsOfXGivenW_s MCGIDI_pdfsOfXGivenW;
typedef struct MCGIDI_pdfsOfXGivenW_sampled_s MCGIDI_pdfsOfXGivenW_sampled;
typedef struct MCGIDI_angular_s MCGIDI_angular;
typedef struct MCGIDI_energyWeightedFunctional_s MCGIDI_energyWeightedFunctional;
typedef struct MCGIDI_energyWeightedFunctionals_s MCGIDI_energyWeightedFunctionals;
typedef struct MCGIDI_energyNBodyPhaseSpace_s MCGIDI_energyNBodyPhaseSpace;
typedef struct MCGIDI_energy_s MCGIDI_energy;
typedef struct MCGIDI_energyAngular_s MCGIDI_energyAngular;
typedef struct MCGIDI_angularEnergy_s MCGIDI_angularEnergy;

typedef struct MCGIDI_decaySamplingInfo_s MCGIDI_decaySamplingInfo;
typedef struct MCGIDI_productsInfo_s MCGIDI_productsInfo;
typedef struct MCGIDI_productInfo_s MCGIDI_productInfo;
typedef struct MCGIDI_sampledProductsData_s MCGIDI_sampledProductsData;
typedef struct MCGIDI_sampledProductsDatas_s MCGIDI_sampledProductsDatas;

#if defined __cplusplus
    }
    }
#endif

enum MCGIDI_quantityLookupMode {      
    MCGIDI_quantityLookupMode_pointwise     /**< Pointwise data are used to determine a quantity's value an energy E. */,
    MCGIDI_quantityLookupMode_grouped       /**< Grouped data are used to determine a quantity's value an energy E. */
};

class MCGIDI_quantitiesLookupModes {

    private:
        int mProjectilesPOPID;
        double mProjectileEnergy;
        int mGroupIndex;
        double mProjectileEnergyForGroupIndex;
        double mTemperature;
        enum MCGIDI_quantityLookupMode mCrossSectionMode;
        enum MCGIDI_quantityLookupMode mMultiplicityMode;

    public:
        MCGIDI_quantitiesLookupModes( int projectilesPOPID );
        ~MCGIDI_quantitiesLookupModes( );

        inline double getProjectileEnergy( void ) const { return( mProjectileEnergy ); }
        void setProjectileEnergy( double e_in ) { mProjectileEnergy = e_in; }

        inline int getGroupIndex( void ) const { return( mGroupIndex ); }
        int setGroupIndex( GIDI_settings const &settings, bool encloseOutOfRange );

        inline double getTemperature( void ) const { return( mTemperature ); }
        void setTemperature( double temperature ) { mTemperature = temperature; }

        enum MCGIDI_quantityLookupMode getMode( std::string const &quantity ) const;
        enum MCGIDI_quantityLookupMode getCrossSectionMode( void ) const { return( mCrossSectionMode ); };
        std::vector<std::string> getListOfLookupQuanities( ) const;
        void setMode( std::string const &quantity, enum MCGIDI_quantityLookupMode mode );
        void setCrossSectionMode( enum MCGIDI_quantityLookupMode mode ) { mCrossSectionMode = mode; };
        void setModeAll( enum MCGIDI_quantityLookupMode mode );
};

typedef struct MCGIDI_samplingMultiplicityBias_s MCGIDI_samplingMultiplicityBias;

struct MCGIDI_samplingMultiplicityBias_s {
    int PoPID;
    double multiplicityFactor;
};

class MCGIDI_samplingMethods {

    public:
        MCGIDI_samplingMethods( );
        ~MCGIDI_samplingMethods( );
};

class MCGIDI_samplingSettings {

    private:                    // This is user input.
        enum GIDI::xDataTOM_frame mWantFrame;
        bool mWantVelocities;
        double (*mRng)( void * );
        void *mRngState;
        std::vector<struct MCGIDI_samplingMultiplicityBias_s> mSamplingMultiplicityBiases;

    public:                     // Temporary variables used in MCGIDI sampling routines.
        enum GIDI::xDataTOM_frame mGotFrame;
        GIDI::MCGIDI_POP *mPoP;
        double mMu;
        double mEp;

    public:
        MCGIDI_samplingSettings( enum GIDI::xDataTOM_frame frame, bool wantVelocities, double (*rng)( void * ), void *rngState );
        ~MCGIDI_samplingSettings( );

        inline double getProductMultiplicityBias( int PoPID ) const {
                for( int i1 = 0; i1 < (int) mSamplingMultiplicityBiases.size( ); ++i1 ) {
                    if( PoPID == mSamplingMultiplicityBiases[i1].PoPID ) return( mSamplingMultiplicityBiases[i1].multiplicityFactor );
                }
                return( 1. ); }
        int setProductMultiplicityBias( GIDI::statusMessageReporting *smr, int PoPID, double fractor );
};

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

enum MCGIDI_transportability {                          /**< This enum is used to give the transportability status for a particle in a reaction or target. */
    MCGIDI_transportability_unknown,                    /**< Particle is not a product of this reaction or target. */
    MCGIDI_transportability_none,                       /**< Particle is a product but has not distribution data. */
    MCGIDI_transportability_partial,                    /**< Particle is a product and has some distribution data. */
    MCGIDI_transportability_full };                     /**< Particle is a product and all needed distribution data. */

#if defined __cplusplus
    }
    }
#endif

typedef std::map<int, enum GIDI::MCGIDI_transportability> transportabilitiesMap;

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

#define MCGIDI_crossSectionType_grouped 1
#define MCGIDI_crossSectionType_pointwise 2

#define MCGIDI_nullReaction -10001

#define MCGIDI_speedOfLight_cm_sec 2.99792458e10
#define MCGIDI_AMU2MeV 931.494028

enum MCGIDI_reactionType { 
    MCGIDI_reactionType_unknown_e,                      /* This should never happen. */
    MCGIDI_reactionType_null_e,                         /* Only occurs when sampling with from grouped cross sections and the projectile is below threshold. */
    MCGIDI_reactionType_elastic_e,                      /* A nuclear elastic reaction. */
    MCGIDI_reactionType_scattering_e,                   /* A nuclear reaction where the projectile and target are products as well as gammas, 
                                                            excluding reactions that are MCGIDI_reactionType_elastic_e and 
                                                            MCGIDI_reactionType_nuclearLevelTransition_e. */
    MCGIDI_reactionType_nuclearIsomerTransmutation_e,   /* A nuclear that changes N or Z and is not one of the others.*/
    MCGIDI_reactionType_nuclearLevelTransition_e,       /* Reaction in which the residual is the same isotope as the target but in a
                                                            different nuclear level. Mainly for meta-stables. */
    MCGIDI_reactionType_capture_e,                      /* A nuclear capture reaction. */
    MCGIDI_reactionType_fission_e,                      /* A nuclear fission reaction. */
    MCGIDI_reactionType_sumOfRemainingOutputChannels_e, /* ENDF MT 5 reactions. */
    MCGIDI_reactionType_atomic_e                        
};

enum MCGIDI_channelGenre { MCGIDI_channelGenre_undefined_e, MCGIDI_channelGenre_twoBody_e, MCGIDI_channelGenre_uncorrelated_e, 
    MCGIDI_channelGenre_sumOfRemaining_e, MCGIDI_channelGenre_twoBodyDecay_e, MCGIDI_channelGenre_uncorrelatedDecay_e };

enum MCGIDI_productMultiplicityType { MCGIDI_productMultiplicityType_invalid_e, MCGIDI_productMultiplicityType_unknown_e, MCGIDI_productMultiplicityType_integer_e,
    MCGIDI_productMultiplicityType_energyDependent_e, MCGIDI_productMultiplicityType_gammaBranching_e, MCGIDI_productMultiplicityType_mixed_e };

enum MCGIDI_distributionType { MCGIDI_distributionType_none_e, MCGIDI_distributionType_unknown_e, MCGIDI_distributionType_angular_e, 
    MCGIDI_distributionType_KalbachMann_e, MCGIDI_distributionType_uncorrelated_e, MCGIDI_distributionType_energyAngular_e, 
    MCGIDI_distributionType_angularEnergy_e };

enum MCGIDI_angularType { MCGIDI_angularType_isotropic, MCGIDI_angularType_recoil, MCGIDI_angularType_linear };

enum MCGIDI_energyType { MCGIDI_energyType_unknown, MCGIDI_energyType_primaryGamma, MCGIDI_energyType_discreteGamma, 
    MCGIDI_energyType_linear, MCGIDI_energyType_generalEvaporation, MCGIDI_energyType_simpleMaxwellianFission, MCGIDI_energyType_evaporation, 
        MCGIDI_energyType_Watt, MCGIDI_energyType_MadlandNix, MCGIDI_energyType_NBodyPhaseSpace, MCGIDI_energyType_weightedFunctional };

extern const char *MCGIDI_productGenre_unknown, *MCGIDI_productGenre_twoBody_angular, *MCGIDI_productGenre_twoBody_formFactor,
    *MCGIDI_productGenre_NBody_angular_energy, *MCGIDI_productGenre_NBody_pairProduction;

#define MCGIDI_particleLevel_continuum -1
#define MCGIDI_particleLevel_sum -2

struct MCGIDI_GammaBranching_s {
    MCGIDI_POP *finalLevel;
    double probability;
};

struct MCGIDI_POP_s {
    MCGIDI_POP *next;
    MCGIDI_POP *parent;
    char *name;
    int globalPoPsIndex;        /* Index of particle in the PoPs library if particle can be return to packages using */
    int Z, A, level, m;         /* this library. Otherwise, -1. */
    double mass_MeV;
    double level_MeV;
    int numberOfGammaBranchs;
    MCGIDI_GammaBranching *gammas;
};

struct MCGIDI_POPs_s {
    int numberOfPOPs, size, increment;
    MCGIDI_POP *first, *last, **sorted;
};

struct MCGIDI_particle_s {
    MCGIDI_particle *prior;
    MCGIDI_particle *next;
    int ordinal;
    int Z, A, m;
    double mass_MeV;
    char *name;
};

struct MCGIDI_decaySamplingInfo_s {
    enum xDataTOM_frame frame;                  /* The frame the product data are in. */
    int isVelocity;                             /* See struct MCGIDI_sampledProductsData_s for meaning. This is user input. */
    double (*rng)( void * );                    /* User supplied rng. */
    void *rngState;                             /* User supplied rng state. */
    MCGIDI_POP *pop;                            /* pop for the sampled product. */
    double mu;                                  /* mu = cos( theta ) for the sampled product. Frame is given by frame member. */
    double Ep;                                  /* Energy of the product. Frame is given by frame member. */
};

struct MCGIDI_productInfo_s {
    int globalPoPsIndex;
    enum MCGIDI_productMultiplicityType productMultiplicityType;
    int multiplicity;
    int transportable;
};

struct MCGIDI_productsInfo_s {
    int numberOfProducts;
    int numberOfAllocatedProducts;
    MCGIDI_productInfo *productInfo;
};

struct MCGIDI_sampledProductsData_s {
    int isVelocity;             /* If true, px_vx, py_vy and pz_vz are velocities otherwise momenta. */
    MCGIDI_POP *pop;
    double kineticEnergy;
    double px_vx;
    double py_vy;
    double pz_vz;
    int delayedNeutronIndex;
    double delayedNeutronRate;
    double birthTimeSec;        /* Some products, like delayed fission neutrons, are to appear (be born) later. */
};

struct MCGIDI_sampledProductsDatas_s {
    int numberOfProducts;
    int numberAllocated;
    int incrementSize;
    MCGIDI_sampledProductsData *products;
};

struct MCGIDI_pdfOfX_s {
    int numberOfXs;
    double *Xs;
    double *pdf;
    double *cdf;
};

struct MCGIDI_pdfsOfXGivenW_s {
    int numberOfWs;
    ptwXY_interpolation interpolationWY, interpolationXY;
    double *Ws;
    MCGIDI_pdfOfX *dist;
};

struct MCGIDI_pdfsOfXGivenW_sampled_s {
    statusMessageReporting *smr;
    ptwXY_interpolation interpolationWY, interpolationXY;
    int iW, iX1, iX2;
    double x, w, frac;
};

struct MCGIDI_angular_s {
    enum xDataTOM_frame frame;
    enum MCGIDI_angularType type;
    MCGIDI_angular *recoilProduct;
    MCGIDI_pdfsOfXGivenW dists;
    double projectileMass_MeV, targetMass_MeV, productMass_MeV, residualMass_MeV;
};

struct MCGIDI_energyWeightedFunctional_s {
    ptwXYPoints *weight;
    MCGIDI_energy *energy;
};

struct MCGIDI_energyWeightedFunctionals_s {
    int numberOfWeights;
    MCGIDI_energyWeightedFunctional weightedFunctional[4];      /* ??????????? Hardwired for no good reason. Will handle up to a (z,4n) reaction. */
};

struct MCGIDI_energyNBodyPhaseSpace_s {
    int numberOfProducts;
    double mass, massFactor, e_inCOMFactor, Q_MeV;
};

struct MCGIDI_energy_s {
    enum xDataTOM_frame frame;
    enum MCGIDI_energyType type;
    double gammaEnergy_MeV;
    double primaryGammaMassFactor;
    double e_inCOMFactor;
    MCGIDI_pdfsOfXGivenW dists;
    double U;
    ptwXYPoints *theta, *Watt_a, *Watt_b;
    ptwXY_interpolation gInterpolation;
    MCGIDI_pdfOfX g;
    MCGIDI_energyWeightedFunctionals weightedFunctionals;
    MCGIDI_energyNBodyPhaseSpace NBodyPhaseSpace;
};

struct MCGIDI_energyAngular_s {
    enum xDataTOM_frame frame;
    MCGIDI_pdfsOfXGivenW pdfOfEpGivenE;
    MCGIDI_pdfsOfXGivenW *pdfOfMuGivenEAndEp;   /* The number of MCGIDI_pdfsOfXGivenW allocated is given by pdfOfEpGivenE.numberOfWs. */
};

struct MCGIDI_angularEnergy_s {
    enum xDataTOM_frame frame;
    MCGIDI_pdfsOfXGivenW pdfOfMuGivenE;
    MCGIDI_pdfsOfXGivenW *pdfOfEpGivenEAndMu;   /* The number of MCGIDI_pdfsOfXGivenW allocated is given by pdfOfMuGivenE.numberOfWs. */
};

struct MCGIDI_KalbachMann_ras_s {
    double *rs;
    double *as;
};

struct MCGIDI_KalbachMann_s {
    enum xDataTOM_frame frame;
    double energyToMeVFactor, massFactor, Sa, Sb, Ma, mb;           /* Needed if a(E,E') is caluclated from the formula. */
    MCGIDI_pdfsOfXGivenW dists;                             /* Sa currently not used. */
    MCGIDI_KalbachMann_ras *ras;
};

struct MCGIDI_distribution_s {
    MCGIDI_product *product;
    enum MCGIDI_distributionType type;
    MCGIDI_angular *angular;                /* All distribution forms must have a frame member. */
    MCGIDI_energy *energy;
    MCGIDI_energyAngular *energyAngular;
    MCGIDI_angularEnergy *angularEnergy;
    MCGIDI_KalbachMann *KalbachMann;
};

struct MCGIDI_outputChannel_s {
    enum MCGIDI_channelGenre genre;
    MCGIDI_reaction *reaction;              /* This is only used for output channels. */
    MCGIDI_product *parent;                 /* This is only used for decay channels. */
    int QIsFloat;
    double Q;
    int numberOfProducts;
    MCGIDI_product *products;
};

struct MCGIDI_product_s {
    MCGIDI_POP *pop;
    char *label;
    MCGIDI_outputChannel *outputChannel;
    int multiplicity;                                       /* If 0, the multiplicity is either 'energyDependent' or 'partialProduction'. */
    int delayedNeutronIndex;
    double delayedNeutronRate;
    ptwXYPoints *multiplicityVsEnergy;
    ptwXYPoints *norms;
    int numberOfPiecewiseMultiplicities;
    ptwXYPoints **piecewiseMultiplicities;
    MCGIDI_distribution distribution;
    MCGIDI_outputChannel decayChannel;
};

struct MCGIDI_reaction_s {
    MCGIDI_target_heated *target;
    int ENDF_MT, ENDL_C, ENDL_S;
    enum MCGIDI_reactionType reactionType;
    char const *outputChannelStr;
    xDataTOM_attributionList attributes;            /* Do not free, owned by attributes. */
    int domainValuesPresent;                        /* True if cross section data defined so EMin and EMax are value. */
    int thresholdGroupIndex;                        /* For grouped data, the group index where threshold starts. */
    double thresholdGroupDomain;                    /* This is groupEnergy[thresholdGroupIndex+1] - EMin. */
    double thresholdGroupedDeltaCrossSection;       /* The adjusted group cross section in group thresholdGroupIndex. */
    double EMin, EMax, finalQ;                      /* BRB, EMin is used as threshold. However, some reactions, especially charged particle */
    ptwXYPoints *crossSection;                      /* have effective thresholds much higher than EMin, may need to handle these differently??????? */
    ptwXPoints *crossSectionGrouped;
    MCGIDI_outputChannel outputChannel;
    MCGIDI_productsInfo productsInfo;               /* See MCGIDI_reaction_ParseDetermineReactionProducts for description. */
    transportabilitiesMap *transportabilities;
};

struct MCGIDI_target_heated_s {
    int ordinal;
    char *path;            /* Partial path of input file. */
    char *absPath;         /* Full absolute path of input file. */
    MCGIDI_POPs pops;
    MCGIDI_POP *projectilePOP;
    MCGIDI_POP *targetPOP;
    xDataTOM_attributionList attributes;
    char *contents;
    double temperature_MeV;
    double EMin, EMax;
    ptwXYPoints *crossSection;
    ptwXPoints *crossSectionGrouped;
    ptwXPoints *crossSectionGroupedForSampling;
    int numberOfReactions;
    MCGIDI_reaction *reactions;
    transportabilitiesMap *transportabilities;
};

struct MCGIDI_target_heated_info_s {
    int ordinal;
    double temperature;
    char *path;                 /* Full path of input file. */
    char *contents;
    MCGIDI_target_heated *heatedTarget;
};

struct MCGIDI_target_s {
    char *path;                 /* Full path of input file. */
    char *absPath;              /* Full absolute path of input file. */
    MCGIDI_POP *projectilePOP;
    MCGIDI_POP *targetPOP;
    xDataTOM_attributionList attributes;
    int nHeatedTargets, nReadHeatedTargets;
    MCGIDI_target_heated *baseHeatedTarget;   /* The lowest temperature whose contents is "all" data, (e.g, not just "crossSection"). */
    MCGIDI_target_heated_info *heatedTargets;         /* List of heated targets in order by temperature. */
    MCGIDI_target_heated_info **readHeatedTargets;    /* List of "read in" heated targets in order by temperature. */
};

char const *MCGIDI_version( void );
int MCGIDI_versionMajor( void );
int MCGIDI_versionMinor( void );
int MCGIDI_versionPatchLevel( void );

/*
* Routines in MCGIDI_target.c
*/
MCGIDI_target *MCGIDI_target_new( statusMessageReporting *smr );
int MCGIDI_target_initialize( statusMessageReporting *smr, MCGIDI_target *target );
MCGIDI_target *MCGIDI_target_newRead( statusMessageReporting *smr, const char *fileName );
int MCGIDI_target_readFromMapViaPoPIDs( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_map *map, const char *evaluation,
        int projectile_PoPID, int target_PoPID );
int MCGIDI_target_readFromMap( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_map *map, const char *evaluation, const char *projectileName, 
    const char *targetName );
MCGIDI_target *MCGIDI_target_newReadFromMapViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation,
        int projectile_PoPID, int target_PoPID );
MCGIDI_target *MCGIDI_target_newReadFromMap( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, const char *projectileName, 
    const char *targetName );
MCGIDI_target *MCGIDI_target_free( statusMessageReporting *smr, MCGIDI_target *target );
int MCGIDI_target_release( statusMessageReporting *smr, MCGIDI_target *target );
int MCGIDI_target_read( statusMessageReporting *smr, MCGIDI_target *target, const char *fileName );
char const *MCGIDI_target_getAttributesValue( statusMessageReporting *smr, MCGIDI_target *target, char const *name );
int MCGIDI_target_getTemperatures( statusMessageReporting *smr, MCGIDI_target *target, double *temperatures );
int MCGIDI_target_readHeatedTarget( statusMessageReporting *smr, MCGIDI_target *target, int index );
MCGIDI_target_heated *MCGIDI_target_getHeatedTargetAtIndex_ReadIfNeeded( statusMessageReporting *smr, MCGIDI_target *target, int index );
MCGIDI_target_heated *MCGIDI_target_getHeatedTargetAtTIndex( statusMessageReporting *smr, MCGIDI_target *target, int index );

int MCGIDI_target_numberOfReactions( statusMessageReporting *smr, MCGIDI_target *target );
enum MCGIDI_reactionType MCGIDI_target_getReactionTypeAtIndex( statusMessageReporting *smr, MCGIDI_target *target, int index );
MCGIDI_reaction *MCGIDI_target_getReactionAtIndex( MCGIDI_target *target, int index );
MCGIDI_reaction *MCGIDI_target_getReactionAtIndex_smr( statusMessageReporting *smr, MCGIDI_target *target, int index );
int MCGIDI_target_numberOfProductionReactions( statusMessageReporting *smr, MCGIDI_target *target );

transportabilitiesMap const *MCGIDI_target_getUniqueProducts( statusMessageReporting *smr, MCGIDI_target *target );
int MCGIDI_target_recast( statusMessageReporting *smr, MCGIDI_target *target, GIDI_settings &settings );

int MCGIDI_target_getDomain( statusMessageReporting *smr, MCGIDI_target *target, double *EMin, double *EMax );
double MCGIDI_target_getTotalCrossSectionAtTAndE( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_quantitiesLookupModes &modes,
        bool sampling );
double MCGIDI_target_getIndexReactionCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target *target, int index, MCGIDI_quantitiesLookupModes &modes,
        bool sampling );
int MCGIDI_target_sampleReaction( statusMessageReporting *smr, MCGIDI_target *target, MCGIDI_quantitiesLookupModes &modes, double totalXSec, 
        double (*userrng)( void * ), void *rngState );
int MCGIDI_target_sampleNullReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target *target,
    MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productDatas );
int MCGIDI_target_sampleIndexReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target *target, int index, 
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productData );
double MCGIDI_target_getIndexReactionFinalQ( statusMessageReporting *smr, MCGIDI_target *target, int index, MCGIDI_quantitiesLookupModes &modes );

/*
* Routines in MCGIDI_target_heated.c
*/
MCGIDI_target_heated *MCGIDI_target_heated_new( statusMessageReporting *smr );
int MCGIDI_target_heated_initialize( statusMessageReporting *smr, MCGIDI_target_heated *target );
MCGIDI_target_heated *MCGIDI_target_heated_newRead( statusMessageReporting *smr, const char *fileName );
MCGIDI_target_heated *MCGIDI_target_heated_free( statusMessageReporting *smr, MCGIDI_target_heated *target );
int MCGIDI_target_heated_release( statusMessageReporting *smr, MCGIDI_target_heated *target );
int MCGIDI_target_heated_read( statusMessageReporting *smr, MCGIDI_target_heated *target, const char *fileName );
int MCGIDI_target_heated_numberOfReactions( statusMessageReporting *smr, MCGIDI_target_heated *target );
int MCGIDI_target_heated_numberOfProductionReactions( statusMessageReporting *smr, MCGIDI_target_heated *target );
MCGIDI_reaction *MCGIDI_target_heated_getReactionAtIndex( MCGIDI_target_heated *target, int index );
MCGIDI_reaction *MCGIDI_target_heated_getReactionAtIndex_smr( statusMessageReporting *smr, MCGIDI_target_heated *target, int index );
#if 0
MCGIDI_reaction *MCGIDI_target_heated_getProductionReactionAtIndex( MCGIDI_target_heated *target, int index );
#endif
MCGIDI_POP *MCGIDI_target_heated_getPOPForProjectile( statusMessageReporting *smr, MCGIDI_target_heated *target );
MCGIDI_POP *MCGIDI_target_heated_getPOPForTarget( statusMessageReporting *smr, MCGIDI_target_heated *target );
double MCGIDI_target_heated_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_target_heated *target );
double MCGIDI_target_heated_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_target_heated *target );
int MCGIDI_target_heated_getEnergyGrid( statusMessageReporting *smr, MCGIDI_target_heated *target, double **energyGrid );
double MCGIDI_target_heated_getTotalCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, MCGIDI_quantitiesLookupModes &modes,
        bool sampling );
double MCGIDI_target_heated_getIndexReactionCrossSectionAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes, bool sampling );
int MCGIDI_target_heated_sampleIndexReactionProductsAtE( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productData );
double MCGIDI_target_heated_getReactionsThreshold( statusMessageReporting *smr, MCGIDI_target_heated *target, int index );
int MCGIDI_target_heated_getReactionsDomain( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, double *EMin, double *EMax );
double MCGIDI_target_heated_getIndexReactionFinalQ( statusMessageReporting *smr, MCGIDI_target_heated *target, int index, 
        MCGIDI_quantitiesLookupModes &modes );

transportabilitiesMap const *MCGIDI_target_heated_getUniqueProducts( statusMessageReporting *smr, MCGIDI_target_heated *target );
int MCGIDI_target_heated_recast( statusMessageReporting *smr, MCGIDI_target_heated *target, GIDI_settings &settings );

/*
* Routines in MCGIDI_reaction.c
*/
MCGIDI_reaction *MCGIDI_reaction_new( statusMessageReporting *smr );
int MCGIDI_reaction_initialize( statusMessageReporting *smr, MCGIDI_reaction *reaction );
MCGIDI_reaction *MCGIDI_reaction_free( statusMessageReporting *smr, MCGIDI_reaction *reaction );
int MCGIDI_reaction_release( statusMessageReporting *smr, MCGIDI_reaction *reaction );
int MCGIDI_reaction_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_target_heated *target, 
    MCGIDI_POPs *pops, MCGIDI_reaction *reaction );
enum MCGIDI_reactionType MCGIDI_reaction_getReactionType( statusMessageReporting *smr, MCGIDI_reaction *reaction );
MCGIDI_target_heated *MCGIDI_reaction_getTargetHeated( statusMessageReporting *smr, MCGIDI_reaction *reaction );
double MCGIDI_reaction_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_reaction *reaction );
double MCGIDI_reaction_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_reaction *reaction );
int MCGIDI_reaction_getDomain( statusMessageReporting *smr, MCGIDI_reaction *reaction, double *EMin, double *EMax );
int MCGIDI_reaction_fixDomains( statusMessageReporting *smr, MCGIDI_reaction *reaction, double EMin, double EMax, nfu_status *status );
double MCGIDI_reaction_getCrossSectionAtE( statusMessageReporting *smr, MCGIDI_reaction *reaction, MCGIDI_quantitiesLookupModes &modes, bool sampling );
double MCGIDI_reaction_getFinalQ( statusMessageReporting *smr, MCGIDI_reaction *reaction, MCGIDI_quantitiesLookupModes &modes );
int MCGIDI_reaction_getENDF_MTNumber( MCGIDI_reaction *reaction );
int MCGIDI_reaction_getENDL_CSNumbers( MCGIDI_reaction *reaction, int *S );
int MCGIDI_reaction_recast( statusMessageReporting *smr, MCGIDI_reaction *reaction, GIDI_settings &settings, 
    GIDI_settings_particle const *projectileSettings, double temperature_MeV, ptwXPoints *totalGroupedCrossSection );

MCGIDI_productsInfo *MCGIDI_reaction_getProductsInfo( MCGIDI_reaction *reaction );
int MCGIDI_productsInfo_getNumberOfUniqueProducts( MCGIDI_productsInfo *productsInfo );
int MCGIDI_productsInfo_getPoPsIndexAtIndex( MCGIDI_productsInfo *productsInfo, int index );
enum MCGIDI_productMultiplicityType MCGIDI_productsInfo_getMultiplicityTypeAtIndex( MCGIDI_productsInfo *productsInfo, int index );
int MCGIDI_productsInfo_getIntegerMultiplicityAtIndex( MCGIDI_productsInfo *productsInfo, int index );
int MCGIDI_productsInfo_getTransportableAtIndex( MCGIDI_productsInfo *productsInfo, int index );

/*
* Routines in MCGIDI_pop.c
*/
MCGIDI_POPs *MCGIDI_POPs_new( statusMessageReporting *smr, int size );
int MCGIDI_POPs_initial( statusMessageReporting *smr, MCGIDI_POPs *pops, int size );
void *MCGIDI_POPs_free( MCGIDI_POPs *pops );
int MCGIDI_POPs_release( MCGIDI_POPs *pops );
MCGIDI_POP *MCGIDI_POPs_addParticleIfNeeded( statusMessageReporting *smr, MCGIDI_POPs *pops, char const *name, double mass_MeV, 
    double level_MeV, MCGIDI_POP *parent, int globalParticle );
int MCGIDI_POPs_findParticleIndex( MCGIDI_POPs *pops, char const *name );
MCGIDI_POP *MCGIDI_POPs_findParticle( MCGIDI_POPs *pops, char const *name );
void MCGIDI_POPs_writeSortedList( MCGIDI_POPs *pops, FILE *f );
void MCGIDI_POPs_printSortedList( MCGIDI_POPs *pops );

MCGIDI_POP *MCGIDI_POP_new( statusMessageReporting *smr, char const *name, double mass_MeV, double level_MeV, MCGIDI_POP *parent );
MCGIDI_POP *MCGIDI_POP_free( MCGIDI_POP *pop );
MCGIDI_POP *MCGIDI_POP_release( MCGIDI_POP *pop );
double MCGIDI_POP_getMass_MeV( MCGIDI_POP *pop );

/*
* Routines in MCGIDI_particle.c
*/
MCGIDI_particle *MCGIDI_particle_new( statusMessageReporting *smr );
int MCGIDI_particle_initialize( statusMessageReporting *smr, MCGIDI_particle *particle );
MCGIDI_particle *MCGIDI_particle_free( statusMessageReporting *smr, MCGIDI_particle *particle );
int MCGIDI_particle_release( statusMessageReporting *smr, MCGIDI_particle *particle );
int MCGIDI_particle_freeInternalList( statusMessageReporting *smr );
MCGIDI_particle *MCGIDI_particle_getInternalID( statusMessageReporting *smr, const char * const name, MCGIDI_POPs *pops );
int MCGIDI_particle_printInternalSortedList( statusMessageReporting *smr );

/*
* Routines in MCGIDI_outputChannel.c
*/
MCGIDI_outputChannel *MCGIDI_outputChannel_new( statusMessageReporting *smr );
int MCGIDI_outputChannel_initialize( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
MCGIDI_outputChannel *MCGIDI_outputChannel_free( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
int MCGIDI_outputChannel_release( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
int MCGIDI_outputChannel_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_POPs *pops, MCGIDI_outputChannel *outputChannel,
    MCGIDI_reaction *reaction, MCGIDI_product *parent );
int MCGIDI_outputChannel_numberOfProducts( MCGIDI_outputChannel *outputChannel );
MCGIDI_product *MCGIDI_outputChannel_getProductAtIndex( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, int i );
int MCGIDI_outputChannel_getDomain( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, double *EMin, double *EMax );
MCGIDI_target_heated *MCGIDI_outputChannel_getTargetHeated( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
double MCGIDI_outputChannel_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
double MCGIDI_outputChannel_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel );
double MCGIDI_outputChannel_getQ_MeV( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, double e_in );
double MCGIDI_outputChannel_getFinalQ( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, double e_in );
int MCGIDI_outputChannel_sampleProductsAtE( statusMessageReporting *smr, MCGIDI_outputChannel *outputChannel, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo, MCGIDI_sampledProductsDatas *productDatas, double *masses );

/*
* Routines in MCGIDI_product.c
*/
MCGIDI_product *MCGIDI_product_new( statusMessageReporting *smr );
int MCGIDI_product_initialize( statusMessageReporting *smr, MCGIDI_product *product );
MCGIDI_product *MCGIDI_product_free( statusMessageReporting *smr, MCGIDI_product *product );
int MCGIDI_product_release( statusMessageReporting *smr, MCGIDI_product *product );
int MCGIDI_product_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_outputChannel *outputChannel,
        MCGIDI_POPs *pops, MCGIDI_product *product, int *delayedNeutronIndex );
int MCGIDI_product_getDomain( statusMessageReporting *smr, MCGIDI_product *product, double *EMin, double *EMax );
int MCGIDI_product_setTwoBodyMasses( statusMessageReporting *smr, MCGIDI_product *product, double projectileMass_MeV, double targetMass_MeV,
    double productMass_MeV, double residualMass_MeV );
double MCGIDI_product_getMass_MeV( statusMessageReporting *smr, MCGIDI_product *product );
MCGIDI_target_heated *MCGIDI_product_getTargetHeated( statusMessageReporting *smr, MCGIDI_product *product );
double MCGIDI_product_getProjectileMass_MeV( statusMessageReporting *smr, MCGIDI_product *product );
double MCGIDI_product_getTargetMass_MeV( statusMessageReporting *smr, MCGIDI_product *product );
int MCGIDI_product_sampleMultiplicity( statusMessageReporting *smr, MCGIDI_product *product, double e_in, double r );
int MCGIDI_product_sampleMu( statusMessageReporting *smr, MCGIDI_product *product, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

int MCGIDI_sampledProducts_initialize( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas, int incrementSize );
int MCGIDI_sampledProducts_release( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas );
int MCGIDI_sampledProducts_remalloc( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas );
int MCGIDI_sampledProducts_addProduct( statusMessageReporting *smr, MCGIDI_sampledProductsDatas *sampledProductsDatas, 
    MCGIDI_sampledProductsData *sampledProductsData );
int MCGIDI_sampledProducts_number( MCGIDI_sampledProductsDatas *sampledProductsDatas );
MCGIDI_sampledProductsData *MCGIDI_sampledProducts_getProductAtIndex( MCGIDI_sampledProductsDatas *sampledProductsDatas, int index );

/*
* Routines in MCGIDI_distribution.c
*/
MCGIDI_distribution *MCGIDI_distribution_new( statusMessageReporting *smr );
int MCGIDI_distribution_initialize( statusMessageReporting *smr, MCGIDI_distribution *distribution );
MCGIDI_distribution *MCGIDI_distribution_free( statusMessageReporting *smr, MCGIDI_distribution *distribution );
int MCGIDI_distribution_release( statusMessageReporting *smr, MCGIDI_distribution *distribution );
int MCGIDI_distribution_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_product *product, MCGIDI_POPs *pops, ptwXYPoints *norms );

/*
* Routines in MCGIDI_angular.c
*/
MCGIDI_angular *MCGIDI_angular_new( statusMessageReporting *smr );
int MCGIDI_angular_initialize( statusMessageReporting *smr, MCGIDI_angular *angular );
MCGIDI_angular *MCGIDI_angular_free( statusMessageReporting *smr, MCGIDI_angular *angular );
int MCGIDI_angular_release( statusMessageReporting *smr, MCGIDI_angular *angular );
int MCGIDI_angular_setTwoBodyMasses( statusMessageReporting *smr, MCGIDI_angular *angular, double projectileMass_MeV, double targetMass_MeV,
    double productMass_MeV, double residualMass_MeV );
int MCGIDI_angular_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms );
int MCGIDI_angular_sampleMu( statusMessageReporting *smr, MCGIDI_angular *angular, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_energy.c
*/
MCGIDI_energy *MCGIDI_energy_new( statusMessageReporting *smr );
int MCGIDI_energy_initialize( statusMessageReporting *smr, MCGIDI_energy *energy );
MCGIDI_energy *MCGIDI_energy_free( statusMessageReporting *smr, MCGIDI_energy *energy );
int MCGIDI_energy_release( statusMessageReporting *smr, MCGIDI_energy *energy );
int MCGIDI_energy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms,
    enum MCGIDI_energyType energyType, double gammaEnergy_MeV );
int MCGIDI_energy_sampleEnergy( statusMessageReporting *smr, MCGIDI_energy *energy, MCGIDI_quantitiesLookupModes &modes, 
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_energyAngular.c
*/
int MCGIDI_energyAngular_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );
MCGIDI_energyAngular *MCGIDI_energyAngular_new( statusMessageReporting *smr );
int MCGIDI_energyAngular_initialize( statusMessageReporting *smr, MCGIDI_energyAngular *energyAngular );
MCGIDI_energyAngular *MCGIDI_energyAngular_free( statusMessageReporting *smr, MCGIDI_energyAngular *energyAngular );
int MCGIDI_energyAngular_release( statusMessageReporting *smr, MCGIDI_energyAngular *energyAngular );
int MCGIDI_energyAngular_sampleDistribution( statusMessageReporting *smr, MCGIDI_distribution *distribution, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_angularEnergy.c
*/
MCGIDI_angularEnergy *MCGIDI_angularEnergy_new( statusMessageReporting *smr );
int MCGIDI_angularEnergy_initialize( statusMessageReporting *smr, MCGIDI_angularEnergy *energyAngular );
MCGIDI_angularEnergy *MCGIDI_angularEnergy_free( statusMessageReporting *smr, MCGIDI_angularEnergy *energyAngular );
int MCGIDI_angularEnergy_release( statusMessageReporting *smr, MCGIDI_angularEnergy *energyAngular );
int MCGIDI_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );
int MCGIDI_angularEnergy_sampleDistribution( statusMessageReporting *smr, MCGIDI_angularEnergy *angularEnergy, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_KalbachMann.c
*/
MCGIDI_KalbachMann *MCGIDI_KalbachMann_new( statusMessageReporting *smr, ptwXY_interpolation interpolationWY, ptwXY_interpolation interpolationXY );
int MCGIDI_KalbachMann_initialize( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann, ptwXY_interpolation interpolationWY, ptwXY_interpolation interpolationXY );
MCGIDI_KalbachMann *MCGIDI_KalbachMann_free( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann );
int MCGIDI_KalbachMann_release( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann );
int MCGIDI_KalbachMann_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );
int MCGIDI_KalbachMann_sampleEp( statusMessageReporting *smr, MCGIDI_KalbachMann *KalbachMann, MCGIDI_quantitiesLookupModes &modes, 
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_uncorrelated.c
*/
int MCGIDI_uncorrelated_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution, ptwXYPoints *norms,
    enum MCGIDI_energyType energyType, double gammaEnergy_MeV );
int MCGIDI_uncorrelated_sampleDistribution( statusMessageReporting *smr, MCGIDI_distribution *distribution, MCGIDI_quantitiesLookupModes &modes,
    MCGIDI_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in MCGIDI_LLNLAngular_angularEnergy.c
*/
int MCGIDI_LLNLAngular_angularEnergy_parseFromTOM( statusMessageReporting *smr, xDataTOM_element *element, MCGIDI_distribution *distribution );

/*
* Routines in MCGIDI_kinetics.c
*/
int MCGIDI_kinetics_2BodyReaction( statusMessageReporting *smr, MCGIDI_angular *angular, double K, double mu, double phi,
        MCGIDI_sampledProductsData *outgoingData );
int MCGIDI_kinetics_COMKineticEnergy2LabEnergyAndMomentum( statusMessageReporting *smr, double beta, double e_kinetic_com, double mu, double phi,
        double m3cc, double m4cc, MCGIDI_sampledProductsData *outgoingData );
int MCGIDI_kinetics_COM2Lab( statusMessageReporting *smr, MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo, double masses[3] );

/*
* Routines in MCGIDI_sampling.c
*/
int MCGIDI_sampling_pdfsOfXGivenW_initialize( statusMessageReporting *smr, MCGIDI_pdfsOfXGivenW *dists );
int MCGIDI_sampling_pdfsOfXGivenW_release( statusMessageReporting *smr, MCGIDI_pdfsOfXGivenW *dists );
int MCGIDI_sampling_pdfsOfX_release( statusMessageReporting *smr, MCGIDI_pdfOfX *dist );
int MCGIDI_sampling_sampleX_from_pdfsOfXGivenW( MCGIDI_pdfsOfXGivenW *dists, MCGIDI_pdfsOfXGivenW_sampled *sampled, double r );
int MCGIDI_sampling_sampleX_from_pdfOfX( MCGIDI_pdfOfX *dist, MCGIDI_pdfsOfXGivenW_sampled *sampled, double r );
int MCGIDI_sampling_doubleDistribution( statusMessageReporting *smr, MCGIDI_pdfsOfXGivenW *pdfOfWGivenV, MCGIDI_pdfsOfXGivenW *pdfOfXGivenVAndW,  
        MCGIDI_quantitiesLookupModes &modes, MCGIDI_decaySamplingInfo *decaySamplingInfo );
int MCGIDI_sampling_interpolationValues( statusMessageReporting *smr, ptwXY_interpolation interpolation, double *ws, double y1, double y2, double *y );
double MCGIDI_sampling_ptwXY_getValueAtX( ptwXYPoints *ptwXY, double x1 );

/*
* Routines in MCGIDI_misc.c
*/
int MCGIDI_misc_NumberOfZSymbols( void );
const char *MCGIDI_misc_ZToSymbol( int iZ );
int MCGIDI_misc_symbolToZ( const char *Z );
int MCGIDI_miscNameToZAm( statusMessageReporting *smr, const char *name, int *Z, int *A, int *m, int *level );
xDataTOM_Int MCGIDI_misc_binarySearch( xDataTOM_Int n, double *ds, double d );
int MCGIDI_misc_PQUStringToDouble( statusMessageReporting *smr, char const *str, char const *unit, double conversion, double *value );
int MCGIDI_misc_PQUStringToDoubleInUnitOf( statusMessageReporting *smr, char const *str, char const *toUnit, double *value );
void MCGIDI_misc_updateTransportabilitiesMap( transportabilitiesMap *transportabilities, int PoPID, enum MCGIDI_transportability transportability );
void MCGIDI_misc_updateTransportabilitiesMap2( transportabilitiesMap *transportabilities, int PoPID, int transportable );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of MCGIDI_h_included. */
