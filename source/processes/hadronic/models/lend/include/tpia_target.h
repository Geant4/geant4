/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#ifndef tpia_target_h_included
#define tpia_target_h_included

#if defined __cplusplus
    extern "C" {
#endif

#include <xData.h>
#include <tpia_defs.h>
#include <tpia_map.h>
#include <tpi_IDs.h>

#if defined __cplusplus
    namespace GIDI {
#endif

typedef struct tpia_particle_s tpia_particle;
typedef struct tpia_target_s tpia_target;
typedef struct tpia_target_heated_info_s tpia_target_heated_info;
typedef struct tpia_target_heated_sorted_s tpia_target_heated_sorted;
typedef struct tpia_target_heated_s tpia_target_heated;
typedef struct tpia_decaySamplingInfo_s tpia_decaySamplingInfo;
typedef struct tpia_productOutgoingData_s tpia_productOutgoingData;
typedef struct tpia_multiplicity_s tpia_multiplicity;
typedef struct tpia_samplingMethods_s tpia_samplingMethods;
typedef struct tpia_EqualProbableBinSpectra_s tpia_EqualProbableBinSpectra;
typedef struct tpia_EqualProbableBinSpectrum_s tpia_EqualProbableBinSpectrum;
typedef struct tpia_LegendreBin_s tpia_LegendreBin;
typedef struct tpia_angularEnergyBin_s tpia_angularEnergyBin;
typedef struct tpia_angular_s tpia_angular;
typedef struct tpia_Legendre_s tpia_Legendre;
typedef struct tpia_angularEnergy_s tpia_angularEnergy;
typedef struct tpia_decayChannel_s tpia_decayChannel;
typedef struct tpia_product_s tpia_product;
typedef struct tpia_1dData_s tpia_1dData;
typedef struct tpia_2dData_s tpia_2dData;
typedef struct tpia_channel_s tpia_channel;
typedef struct tpia_data_frame_s tpia_data_frame;

//#include <xData.h>
#include <tpia_mass.h>
//#include <tpia_map.h>
//#include <tpi_IDs.h>

#define tpia_mode_MonteCarlo 1
#define tpia_mode_Pn 2
#define tpia_maxNumberOfFrames 8
#define tpia_referenceFrame_None 0
#define tpia_referenceFrame_COM 1
#define tpia_referenceFrame_lab 2
#define tpia_referenceFrame_Max tpia_referenceFrame_lab

#define tpia_crossSectionType_grouped 1
#define tpia_crossSectionType_pointwise 2

#define tpia_m_depositionEnergy 1
#define tpia_m_multiplicity ( 1 << 1 )
#define tpia_m_decayChannel ( 1 << 2 )
#define tpia_m_commonShift 3
#define tpia_m_angular ( 1 << 3 )
#define tpia_m_formFactor ( 1 << 4 )
#define tpia_m_Legendre ( 1 << 5 )
#define tpia_m_angular_energy ( 1 << 6 )

#define tpia_speedOfLight_cm_sec 2.99792458e10
#define tpia_AMU2MeV 931.494028

extern DLL_LEND const char *tpia_productGenre_unknown, *tpia_productGenre_twoBody_angular, *tpia_productGenre_twoBody_formFactor,
    *tpia_productGenre_NBody_Legendre, *tpia_productGenre_NBody_angular_energy, *tpia_productGenre_NBody_uncorrelate_Legendre,
    *tpia_productGenre_NBody_pairProduction;

extern DLL_LEND const char *tpia_samplingMethods_constant, *tpia_samplingMethods_linear;

//extern const char *tpia_samplingMethods_constant = "constant";
//extern const char *tpia_samplingMethods_linear = "linear";

//#define tpia_samplingMethods_isConstant( method ) ( method == tpia_samplingMethods_constant )
//#define tpia_samplingMethods_isLinear( method ) ( method == tpia_samplingMethods_linear )

//tpia_samplingMethods_constant and tpia_samplingMethods_linear are defined in src/tpia_samplingMethods.cc in namespace GIDI.
//So TK directory write the values in follwoing macros. 110527
//#define tpia_samplingMethods_isConstant( method ) ( method == "constant" )
//#define tpia_samplingMethods_isLinear( method ) ( method == "linear" )
//Fix above bug. 110602
#define tpia_samplingMethods_isConstant( method ) ( strcmp ( method , "constant" ) == 0 )
#define tpia_samplingMethods_isLinear( method ) ( strcmp ( method , "linear" ) == 0 )


struct tpia_particle_s {
    tpia_particle *prior;
    tpia_particle *next;
    int ordinal;
    int Z, A, m;
    double mass;
    double fullMass_MeV;
    char *name;
    tpi_spectralID *spectralID;
};

struct tpia_samplingMethods_s {
    const char *angular_equalProbableBinMethod;
};

struct tpia_data_frame_s {
    unsigned int frames;
};

struct tpia_1dData_s {
    xData_Int start, end, length;
    double *data;
};

struct tpia_decaySamplingInfo_s {
    double e_in;
    tpia_samplingMethods *samplingMethods;
    int isVelocity;
    tpia_data_frame frame;
    double (*rng)( void * );
    void *rngState;
    tpia_particle *productID;
    char const *genre;
    double mu;
    double Ep;
};

struct tpia_productOutgoingData_s {
    char const *genre;
    int isVelocity;
    tpia_data_frame frame;
    tpia_particle *productID;
    double kineticEnergy;
    double px_vx;
    double py_vy;
    double pz_vz;
    tpia_decayChannel *decayChannel;
};

struct tpia_multiplicity_s {
    tpia_multiplicity *next;
    tpia_data_frame frame;
    double timeScale;                           /* Only for delay neutrons from fission, otherwise -1. */
    xData_Int numberOfPointwise;
    tpia_1dData grouped;
    double *pointwise;
};
/*
*  This data is used for discrete 2-body data where mu is give in equal probable bins vs E and
*  for isotropic scattering where the outgoing E' is give in equal probable bins vs E.
*/
struct tpia_EqualProbableBinSpectrum_s {
    int nBins, index;
    double value;
    double *bins;
};

struct tpia_EqualProbableBinSpectra_s {
    int iValue;
    int nBins;
    int numberOfEs;
    double dValue;
    tpia_EqualProbableBinSpectrum *energies;
};

struct tpia_LegendreBin_s {
    int numberOfLs;
    tpia_EqualProbableBinSpectra *ls;
};

struct tpia_angularEnergyBin_s {
    int nBins, numberOfEs;
    tpia_EqualProbableBinSpectra *energies;
};

struct tpia_decayChannel_s {
    int numberOfProducts;
    tpia_product *products;
    double m1_fullMass_MeV;     /* Currently, only used for 2-body decay. Either the mass of the projectile and target, respectively; */
    double m2_fullMass_MeV;     /* or, the mass of the decaying particle and 0., respectively. */
};

struct tpia_angular_s {
    tpia_data_frame frame;
    tpia_EqualProbableBinSpectra binned;
};

struct tpia_Legendre_s {
    tpia_data_frame frame;
    tpia_LegendreBin binned;
};

struct tpia_angularEnergy_s {
    tpia_data_frame frame;
    tpia_angularEnergyBin binned;
};

struct tpia_product_s {
    tpia_product *next;
    tpia_channel *channel;
    tpia_product *parentProduct;                            /* If this product is the result of a particle decaying, then this the the parent particle. */
    xData_attributionList attributes;
    tpia_particle *productID;
    char const *genre;
    long b_dataPresent;
    long b_dataRequired;
    tpia_1dData depositionEnergyGrouped;
    int multiplicity;                                       /* If 0, the multiplicity is either 'energyDependent' or 'partialProduction'. */
    tpia_multiplicity *multiplicityVsEnergy;
    tpia_multiplicity *delayedNeutronMultiplicityVsEnergy;   /* Fission delayed neutron multiplicities. */
    tpia_angular angular;
    tpia_Legendre Legendre;
    tpia_angularEnergy angularEnergy;
    tpia_decayChannel decayChannel;
};

struct tpia_channel_s {
    tpia_target_heated *target;
    char *outputChannel;
    char *genre;
    xData_attributionList attributes;
    int ENDL_C, ENDF_MT;
    char *QString;
    char *fission;
    int QIsFloat;
    double Q;
    tpia_data_frame crossSectionFrame;
    tpia_1dData crossSectionPointwise;
    tpia_1dData crossSectionGrouped;
    tpia_1dData availableEnergyGrouped;
    tpia_decayChannel decayChannel;
};

struct tpia_target_heated_s {
    int ordinal;
    char *path;            /* Partial path of input file. */
    char *absPath;         /* Full absolute path of input file. */
    tpia_particle *projectileID;
    tpia_particle *targetID;
    int nGroups;
    xData_Int energyGridLength;
    double *energyGrid;
    tpia_1dData totalCrossSectionPointwise;
    tpia_1dData totalCrossSectionGrouped;
    xData_attributionList attributes;
    char *contents;
    int nChannels, nProductionChannels;
    tpia_channel **channels;
    tpia_channel **productionChannels;
    double *kerma;
};

struct tpia_target_heated_info_s {
    int ordinal;
    double temperature;
    char *path;                 /* Full path of input file. */
    char *contents;
    tpia_target_heated *heatedTarget;
};

struct tpia_target_s {
    char *path;                 /* Full path of input file. */
    char *absPath;              /* Full absolute path of input file. */
    tpia_particle *projectileID;
    tpia_particle *targetID;
    xData_attributionList attributes;
    tpia_samplingMethods samplingMethods;
    int nHeatedTargets, nReadHeatedTargets;
    tpia_target_heated *baseHeatedTarget;   /* The lowest temperature whose contents is "all" data, (e.g, not just "crossSection"). */
    tpia_target_heated_info *heatedTargets;         /* List of heated targets in order by temperature. */
    tpia_target_heated_info **readHeatedTargets;    /* List of "read in" heated targets in order by temperature. */
};

/*
* Routines in tpia_target.c
*/
tpia_target *tpia_target_create( statusMessageReporting *smr );
int tpia_target_initialize( statusMessageReporting *smr, tpia_target *target );
tpia_target *tpia_target_createRead( statusMessageReporting *smr, const char *fileName );
int tpia_target_readFromMap( statusMessageReporting *smr, tpia_target *target, tpia_map *map, const char *evaluation, const char *projectileName, 
    const char *targetName );
tpia_target *tpia_target_createReadFromMap( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectileName, 
    const char *targetName );
tpia_target *tpia_target_free( statusMessageReporting *smr, tpia_target *target );
int tpia_target_release( statusMessageReporting *smr, tpia_target *target );
int tpia_target_read( statusMessageReporting *smr, tpia_target *target, const char *fileName );
char *tpia_target_getAttributesValue( statusMessageReporting *smr, tpia_target *target, char const *name );
int tpia_target_getTemperatures( statusMessageReporting *smr, tpia_target *target, double *temperatures );
int tpia_target_readHeatedTarget( statusMessageReporting *smr, tpia_target *target, int index, int checkElememtsForAccess );
tpia_target_heated *tpia_target_getHeatedTargetAtIndex_ReadIfNeeded( statusMessageReporting *smr, tpia_target *target, int index );
int tpia_target_numberOfChannels( statusMessageReporting *smr, tpia_target *target );
int tpia_target_numberOfProductionChannels( statusMessageReporting *smr, tpia_target *target );
xData_Int tpia_target_getEnergyGridAtTIndex( statusMessageReporting *smr, tpia_target *target, int index, double **energyGrid );
tpia_1dData *tpia_target_getTotalCrossSectionAtTIndex( statusMessageReporting *smr, tpia_target *target, int index, int crossSectionType );
double tpia_target_getTotalCrossSectionAtTAndE( statusMessageReporting *smr, tpia_target *target, double T, xData_Int iEg, double e_in,
        int crossSectionType );
double tpia_target_getIndexChannelCrossSectionAtE( statusMessageReporting *smr, tpia_target *target, int index, double T, xData_Int iEg, double e_in,
        int crossSectionType );
int tpia_target_sampleIndexChannelProductsAtE( statusMessageReporting *smr, tpia_target *target, int index, double T, 
        tpia_decaySamplingInfo *decaySamplingInfo, int nProductData, tpia_productOutgoingData *productData );

/*
* Routines in tpia_target_heated.c
*/
tpia_target_heated *tpia_target_heated_create( statusMessageReporting *smr );
int tpia_target_heated_initialize( statusMessageReporting *smr, tpia_target_heated *target );
tpia_target_heated *tpia_target_heated_createRead( statusMessageReporting *smr, const char *fileName, int checkElememtsForAccess );
tpia_target_heated *tpia_target_heated_free( statusMessageReporting *smr, tpia_target_heated *target );
int tpia_target_heated_release( statusMessageReporting *smr, tpia_target_heated *target );
int tpia_target_heated_read( statusMessageReporting *smr, tpia_target_heated *target, const char *fileName, int checkElememtsForAccess );
int tpia_target_heated_numberOfChannels( statusMessageReporting *smr, tpia_target_heated *target );
int tpia_target_heated_numberOfProductionChannels( statusMessageReporting *smr, tpia_target_heated *target );
tpia_channel *tpia_target_heated_getChannelAtIndex( tpia_target_heated *target, int index );
tpia_channel *tpia_target_heated_getChannelAtIndex_smr( statusMessageReporting *smr, tpia_target_heated *target, int index );
tpia_channel *tpia_target_heated_getProductionChannelAtIndex( tpia_target_heated *target, int index );
xData_Int tpia_target_heated_getEnergyGrid( statusMessageReporting *smr, tpia_target_heated *target, double **energyGrid );
xData_Int tpia_target_heated_getEIndex( tpia_target_heated *target, double e_in );
double tpia_target_heated_getTotalCrossSectionAtE( statusMessageReporting *smr, tpia_target_heated *target, xData_Int gE, double e_in, 
        int crossSectionType );
double tpia_target_heated_getIndexChannelCrossSectionAtE( statusMessageReporting *smr, tpia_target_heated *target, int index, xData_Int iEg, double e_in,
        int crossSectionType );
int tpia_target_heated_sampleIndexChannelProductsAtE( statusMessageReporting *smr, tpia_target_heated *target, int index, 
        tpia_decaySamplingInfo *decaySamplingInfo, int nProductData, tpia_productOutgoingData *productData );

/*
* Routines in tpia_channel.c
*/
tpia_channel *tpia_channel_create( statusMessageReporting *smr );
int tpia_channel_initialize( statusMessageReporting *smr, tpia_channel *channel );
tpia_channel *tpia_channel_createGetFromElement( statusMessageReporting *smr, tpia_target_heated *target, xData_element *channelElement,
    int pointwiseRequired );
tpia_channel *tpia_channel_free( statusMessageReporting *smr, tpia_channel *channel );
int tpia_channel_release( statusMessageReporting *smr, tpia_channel *channel );
int tpia_channel_getFromElement( statusMessageReporting *smr, tpia_target_heated *target, xData_element *channelElement, tpia_channel *channel,
    int pointwiseRequired );
tpia_product *tpia_channel_getFirstProduct( tpia_channel *channel );
tpia_product *tpia_channel_getProductByIndex( statusMessageReporting *smr, tpia_channel *channel, int index );
int tpia_channel_numberOfProducts( statusMessageReporting *smr, tpia_channel *channel );
int tpia_channel_isProduction( statusMessageReporting *smr, tpia_channel *channel );
double tpia_channel_getCrossSectionAtE( statusMessageReporting *smr, tpia_channel *channel, xData_Int iEg, double e_in, int crossSectionType );

/*
* Routines in tpia_particle.c
*/
tpia_particle *tpia_particle_create( statusMessageReporting *smr );
int tpia_particle_initialize( statusMessageReporting *smr, tpia_particle *particle );
tpia_particle *tpia_particle_free( statusMessageReporting *smr, tpia_particle *particle );
int tpia_particle_release( statusMessageReporting *smr, tpia_particle *particle );
int tpia_particle_freeInternalList( statusMessageReporting *smr );
tpia_particle *tpia_particle_getInternalID( statusMessageReporting *smr, const char * const name );
int tpia_particle_printInternalSortedList( statusMessageReporting *smr );

/*
* Routines in tpia_product.c
*/
tpia_product *tpia_product_create( statusMessageReporting *smr );
int tpia_product_initialize( statusMessageReporting *smr, tpia_product *product );
tpia_product *tpia_product_createGetFromElement( statusMessageReporting *smr, tpia_channel *channel, tpia_product *parentProduct, 
    xData_element *productElement );
tpia_product *tpia_product_free( statusMessageReporting *smr, tpia_product *product );
int tpia_product_release( statusMessageReporting *smr, tpia_product *product);
int tpia_product_getFromElement( statusMessageReporting *smr, tpia_channel *channel, tpia_product *parentProduct, xData_element *productElement, 
    tpia_product *product );
int tpia_product_getDecayChannelFromElement( statusMessageReporting *smr, xData_element *parentElement, tpia_channel *channel, tpia_product *parentProduct,
    tpia_product **priorProductNext );
long tpia_product_dataRequired( statusMessageReporting *smr, tpia_product *product );
tpia_product *tpia_product_getFirstProduct( tpia_product *product );
tpia_product *tpia_product_getProductByIndex( statusMessageReporting *smr, tpia_product *product, int index );
int tpia_product_doesDecay( statusMessageReporting *smr, tpia_product *product );
int tpia_product_numberOfProducts( statusMessageReporting *smr, tpia_product *product );
int tpia_product_isDataPresent( statusMessageReporting *smr, tpia_product *product, int b_data );
int tpia_product_sampleMultiplicity( statusMessageReporting *smr, tpia_product *product, double e_in, double r );

/*
* Routines in tpia_decayChannel.c
*/
tpia_product *tpia_decayChannel_getFirstProduct( tpia_decayChannel *decayChannel );
tpia_product *tpia_decayChannel_getNextProduct( tpia_product *product );
int tpia_decayChannel_sampleProductsAtE( statusMessageReporting *smr, tpia_decayChannel *decayChannel, tpia_decaySamplingInfo *decaySamplingInfo,
        int nProductData, tpia_productOutgoingData *productData );

/*
* Routines in tpia_kinetics.c
*/
int tpia_kinetics_2BodyReaction( statusMessageReporting *smr, tpia_decayChannel *decayChannel, double K, double mu, double phi, 
        tpia_productOutgoingData *outgoingData );
int tpia_kinetics_COMKineticEnergy2LabEnergyAndMomentum( statusMessageReporting *smr, double beta, double e_kinetic_com, double mu, double phi,
        double m3cc, double m4cc, tpia_productOutgoingData *outgoingData );

/*
* Routines in tpia_frame.c
*/
int tpia_frame_clear( statusMessageReporting *smr, tpia_data_frame *frame );
int tpia_frame_setFromElement( statusMessageReporting *smr, xData_element *element, int dimension, tpia_data_frame *frame );
int tpia_frame_setFromString( statusMessageReporting *smr, const char *forItem, const char *value, int dimension, tpia_data_frame *frame );
int tpia_frame_getDimensions( statusMessageReporting *smr, tpia_data_frame *frame );
char *tpia_frame_toString( statusMessageReporting *smr, const char *fromItem, tpia_data_frame *frame );
int tpia_frame_setColumns( statusMessageReporting *smr, tpia_data_frame *frame, int nColumns, int *values );
int tpia_frame_setColumn( statusMessageReporting *smr, tpia_data_frame *frame, int column, int value );
int tpia_frame_getColumn( statusMessageReporting *smr, tpia_data_frame *frame, int column );

/*
* Routines in tpia_multiplicity.c
*/
tpia_multiplicity *tpia_multiplicity_create( statusMessageReporting *smr );
int tpia_multiplicity_initialize( statusMessageReporting *smr, tpia_multiplicity *multiplicity );
tpia_multiplicity *tpia_multiplicity_free( statusMessageReporting *smr, tpia_multiplicity *multiplicity );
int tpia_multiplicity_release( statusMessageReporting *smr, tpia_multiplicity *multiplicity );
tpia_multiplicity *tpia_multiplicity_createGetFromElement( statusMessageReporting *smr, xData_element *multiplicityElement, int nGroups );
int tpia_multiplicity_getFromElement( statusMessageReporting *smr, xData_element *multiplicityElement, tpia_multiplicity *multiplicity, int nGroups );
int tpia_multiplicity_getTimeScaleFromElement( statusMessageReporting *smr, xData_element *element, const char **timeScale, int *isDelayedNeutrons,
    double *dTimeScale );

/*
* Routines in tpia_angular.c
*/
int tpia_angular_initialize( statusMessageReporting *smr, tpia_angular *angular );
int tpia_angular_release( statusMessageReporting *smr, tpia_angular *angular );
int tpia_angular_getFromElement( statusMessageReporting *smr, xData_element *angularElement, tpia_angular *angular );
int tpia_angular_SampleMu( statusMessageReporting *smr, tpia_angular *angular, tpia_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in tpia_Legendre.c
*/
int tpia_Legendre_initialize( statusMessageReporting *smr, tpia_Legendre *Legendre );
int tpia_Legendre_release( statusMessageReporting *smr, tpia_Legendre *Legendre );
int tpia_Legendre_getFromElement( statusMessageReporting *smr, xData_element *LegendreElement, tpia_Legendre *Legendre );
int tpia_Legendre_SampleEp( statusMessageReporting *smr, tpia_Legendre *Legendre, int sampleMu, tpia_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in tpia_angularEnergy.c
*/
int tpia_angularEnergy_initialize( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy );
int tpia_angularEnergy_release( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy );
int tpia_angularEnergy_getFromElement( statusMessageReporting *smr, xData_element *angularEnergyElement, tpia_angularEnergy *angularEnergy );
int tpia_angularEnergy_SampleEp( statusMessageReporting *smr, tpia_angularEnergy *angularEnergy, tpia_decaySamplingInfo *decaySamplingInfo );

/*
* Routines in tpia_samplingMethods.c
*/
int tpia_samplingMethods_initialize( statusMessageReporting *smr, tpia_samplingMethods *samplingMethods );

/*
* Routines in tpia_misc.c
*/
int tpia_misc_NumberOfZSymbols( void );
const char *tpia_misc_ZToSymbol( int iZ );
int tpia_misc_symbolToZ( const char *Z );
int tpia_miscNameToZAm( statusMessageReporting *smr, const char *name, int *Z, int *A, int *m );
char *tpia_misc_pointerToAttributeIfAllOk( statusMessageReporting *smr, xData_element *element, const char *path, int required, 
        xData_attributionList *attributes, const char *name, const char *file, int line );
double *tpia_misc_get2dx_y_data( statusMessageReporting *smr, xData_element *element, xData_Int *length );
double *tpia_misc_get2dxindex_y_data( statusMessageReporting *smr, xData_element *element, xData_Int *start, xData_Int *end, double *xValues );
double *tpia_misc_get2d_xShared_yHistogram_data( statusMessageReporting *smr, xData_element *element, xData_Int *start, xData_Int *end, xData_Int *length );
int tpia_misc_get2d_xShared_yHistogram_data_Grouped( statusMessageReporting *smr, xData_element *element, tpia_1dData *group );
double tpia_misc_getPointwiseCrossSectionAtE( statusMessageReporting *smr, tpia_1dData *crossSection, double *energyGrid, xData_Int index, double e_in );
double tpia_misc_drng( double (*rng)( void * ), void *rngState );
int tpia_misc_sampleEqualProbableBin( statusMessageReporting *smr, tpia_decaySamplingInfo *decaySamplingInfo, double e_in, int nBins,
        tpia_EqualProbableBinSpectra *binned, double *value_ );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of tpia_target_h_included. */
