/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef GIDI_hpp_included
#define GIDI_hpp_included 1

#include <string>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <iostream>
#include <stdexcept>
#include <regex>

#include <GUPI.hpp>
#include <RISI.hpp>
#include <PoPI.hpp>

#include <nf_utilities.h>
#include <nf_buffer.h>
#include <ptwXY.h>

#include "GIDI_data.hpp"

namespace GIDI {

class SetupInfo;
class Form;
class Suite;
class FissionFragmentData;
class OutputChannel;
class Protare;
class ProtareSingle;
class ParticleInfo;
class MultiGroupCalulationInformation;
namespace  GRIN {
class GRIN_continuumGammas;
}

typedef std::set<int> ExcludeReactionsSet;

namespace Functions {
    class XYs1d;
    class Xs_pdf_cdf1d;
    class Branching1d;
    class Function2dForm;
}                   // End namespace Functions.

namespace Map {
    class ProtareBase;
    class TNSL;
    class Map;
}                   // End of namespace Map.

typedef bool (*MapWalkCallBack)( Map::ProtareBase const *a_protareEntry, std::string const &a_library, void *a_userData, int a_level );

namespace Construction {
    class Settings;
}                   // End of namespace Construction.

namespace Documentation_1_10 {
    class Suite;
}                   // End of namespace Documentation_1_10.

namespace ExternalFiles {
    class Suite;
}                   // End of namespace ExternalFiles

namespace Table {
    class Column;
}                   // End of namespace Table

namespace Styles {
    class Suite;
    class MultiGroup;
    class HeatedMultiGroup;
}                   // End of namespace Styles.

typedef Form *(*parseSuite)( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
        PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );

enum class GNDS_FileType { uninitialized, unknown, pops, protare, covarianceSuite, map };

class GNDS_FileTypeInfo {

    private:
        GNDS_FileType m_GNDS_fileType;
        std::string m_projectileID;
        std::string m_targetID;
        std::string m_evaluation;
        std::string m_interaction;

    public:
        GNDS_FileTypeInfo( );
        GNDS_FileTypeInfo( GNDS_FileType a_GNDS_fileType, std::string a_projectileID = "", std::string a_targetID = "", std::string a_evaluation = "",
                        std::string a_interaction = "" );
        GNDS_FileTypeInfo( GNDS_FileTypeInfo const &a_GNDS_fileTypeInfo );
        GNDS_FileTypeInfo &operator=( GNDS_FileTypeInfo const &a_rhs );

        GNDS_FileType GNDS_fileType( ) const { return( m_GNDS_fileType ); }
        void setGNDS_fileType( GNDS_FileType a_GNDS_fileType ) { m_GNDS_fileType = a_GNDS_fileType; }
        std::string const &projectileID( ) const { return( m_projectileID ); }
        std::string const &targetID( ) const { return( m_targetID ); }
        std::string const &evaluation( ) const { return( m_evaluation ); }
        std::string const &interaction( ) const { return( m_interaction ); }
};

enum class ProtareType { single, composite, TNSL };

enum class FormType { generic, lazyParsingHelperForm, group, groups, transportable, flux, fluxes, externalFile, style,
                reaction, product, delayedNeutron, fissionFragmentData, rate,
                physicalQuantity, axisDomain, axis, grid, axes,
                flattenedArrayData, array3d, array,
                    // 1d functions.
                constant1d, XYs1d, Ys1d, polynomial1d, Legendre1d, gridded1d, reference1d, xs_pdf_cdf1d, regions1d, 
                resonancesWithBackground1d, resonanceBackground1d, resonanceBackgroundRegion1d, URR_probabilityTables1d,
                fissionEnergyRelease1d, branching1d, branching1dPids, thermalNeutronScatteringLaw1d, unspecified1d,
                    // 2d functions.
                XYs2d, recoil2d, isotropic2d, discreteGamma2d, primaryGamma2d, regions2d, gridded2d,
                generalEvaporation2d, simpleMaxwellianFission2d, evaporation2d, Watt2d, MadlandNix2d, 
                weighted_function2d, weightedFunctionals2d, NBodyPhaseSpace2d,
                    // 3d functions.
                XYs3d, regions3d, gridded3d,
                    // distributions.
                angularTwoBody, KalbachMann, uncorrelated, unspecified, reference3d, multiGroup3d, 
                energyAngular, energyAngularMC, angularEnergy, angularEnergyMC, LLNL_angularEnergy,
                coherentPhotonScattering, incoherentPhotonScattering, incoherentBoundToFreePhotonScattering, thermalNeutronScatteringLaw, branching3d,
                coherentElastic, incoherentElastic, incoherentInelastic, CoulombPlusNuclearElastic3d, LLNLLegendre,
                    // Sums stuff.
                crossSectionSum, multiplicitySum, summands,
                    // ACE style URR stuff currently in the applicationData node.
                ACE_URR_probabilityTable, ACE_URR_incidentEnergy,
                    // Table stuff.
                table, columnHeaders, column,
                    // Non-GNDS compliant GRIN forms.
                GRIN_inelasticIncidentEnergy, GRIN_captureLevelProbability };

enum class Frame { lab, centerOfMass };
enum class TransportCorrectionType { None, Pendlebury, LLNL, Ferguson };
enum class FileType { XML, HDF };

#define GIDI_emptyFileNameChars ""

#define GIDI_mapFormatVersion_0_1Chars "0.1"
#define GIDI_mapFormatVersion_0_2Chars "0.2"

#define GIDI_LLNL_Chars "LLNL"
#define GIDI_LLNL_multiGroupReactions_Chars "LLNL::multiGroupReactions"
#define GIDI_LLNL_multiGroupDelayedNeutrons_Chars "LLNL::multiGroupDelayedNeutrons"
#define GIDI_LLNL_URR_probability_tables_Chars "LLNL::URR_probability_tables"
#define GIDI_LLNL_pointwiseAverageProductEnergies "LLNL::pointwiseAverageProductEnergies"
#define GIDI_LLNL_GRIN_continuumGammas "LLNL::GRIN_continuumGammas"

#define GIDI_mapChars "map"
#define GIDI_importChars "import"
#define GIDI_protareChars "protare"
#define GIDI_TNSLChars "TNSL"

#define GIDI_formatChars "format"

#define GIDI_topLevelChars "reactionSuite"
#define GIDI_covarianceSuiteChars "covarianceSuite"

#define GIDI_externalFilesChars "externalFiles"
#define GIDI_externalFileChars "externalFile"

#define GIDI_documentationChars "documentation"
#define GIDI_documentations_1_10_Chars "documentations"
#define GIDI_stylesChars "styles"
#define GIDI_PoPsChars "PoPs"
#define GIDI_reactionsChars "reactions"
#define GIDI_reactionChars "reaction"
#define GIDI_orphanProductsChars "orphanProducts"
#define GIDI_orphanProductChars "orphanProduct"
#define GIDI_incompleteReactionsChars "incompleteReactions"
#define GIDI_fissionComponentsChars "fissionComponents"
#define GIDI_fissionComponentChars "fissionComponent"
#define GIDI_ACE_URR_probabilityTablesChars "probabilityTables"
#define GIDI_ACE_URR_probabilityTableChars "probabilityTable"
#define GIDI_LLNL_photoAtomicIncoherentDoppler_Chars "LLNL::photoAtomicIncoherentDoppler"

#define GIDI_tableChars "table"
#define GIDI_rowsChars "rows"
#define GIDI_columnsChars "columns"
#define GIDI_columnHeadersChars "columnHeaders"
#define GIDI_columnChars "column"
#define GIDI_nameChars "name"
#define GIDI_typesChars "types"
#define GIDI_dataChars "data"
#define GIDI_sepChars "sep"

#define GIDI_applicationDataChars "applicationData"
#define GIDI_institutionChars "institution"
#define GIDI_nuclearPlusCoulombInterferenceChars "nuclearPlusCoulombInterference"

#define GIDI_sumsChars "sums"
#define GIDI_sumsCrossSectionsChars "crossSections"
#define GIDI_sumsMultiplicitiesChars "multiplicities"
#define GIDI_sumsAddChars "add"
#define GIDI_sumsSummandsChars "summands"
#define GIDI_crossSectionSumsChars "crossSectionSums"
#define GIDI_crossSectionSumChars "crossSectionSum"
#define GIDI_multiplicitySumsChars "multiplicitySums"
#define GIDI_multiplicitySumChars "multiplicitySum"

#define GIDI_doubleDifferentialCrossSectionChars "doubleDifferentialCrossSection"
#define GIDI_crossSectionChars "crossSection"
#define GIDI_availableEnergyChars "availableEnergy"
#define GIDI_availableMomentumChars "availableMomentum"

#define GIDI_QChars "Q"
#define GIDI_productsChars "products"
#define GIDI_productChars "product"

#define GIDI_multiplicityChars "multiplicity"
#define GIDI_distributionChars "distribution"
#define GIDI_averageEnergyChars "averageProductEnergy"
#define GIDI_averageMomentumChars "averageProductMomentum"
#define GIDI_outputChannelChars "outputChannel"

#define GIDI_fissionFragmentDataChars "fissionFragmentData"
#define GIDI_delayedNeutronsChars "delayedNeutrons"
#define GIDI_delayedNeutronChars "delayedNeutron"
#define GIDI_fissionEnergyReleasesChars "fissionEnergyReleases"
#define GIDI_fissionEnergyReleaseChars "fissionEnergyRelease"
#define GIDI_rateChars "rate"

#define GIDI_groupsChars "groups"
#define GIDI_groupChars "group"
#define GIDI_fluxesChars "fluxes"

#define GIDI_evaluatedStyleChars "evaluated"
#define GIDI_crossSectionReconstructedStyleChars "crossSectionReconstructed"
#define GIDI_angularDistributionReconstructedStyleChars "angularDistributionReconstructed"
#define GIDI_CoulombPlusNuclearElasticMuCutoffStyleChars "CoulombPlusNuclearElasticMuCutoff"
#define GIDI_averageProductDataStyleChars "averageProductData"
#define GIDI_MonteCarlo_cdfStyleChars "MonteCarlo_cdf"
#define GIDI_multiGroupStyleChars "multiGroup"
#define GIDI_transportablesChars "transportables"
#define GIDI_transportableChars "transportable"
#define GIDI_realizationChars "realization"
#define GIDI_heatedStyleChars "heated"
#define GIDI_griddedCrossSectionStyleChars "griddedCrossSection"
#define GIDI_URR_probabilityTablesStyleChars "URR_probabilityTables"
#define GIDI_heatedMultiGroupStyleChars "heatedMultiGroup"
#define GIDI_SnElasticUpScatterStyleChars "SnElasticUpScatter"
#define GIDI_projectileEnergyDomainChars "projectileEnergyDomain"

// array monikers and allowed values.
#define GIDI_arrayChars "array"
#define GIDI_noneChars "none"
#define GIDI_valuesChars "values"
#define GIDI_shapeChars "shape"

#define GIDI_compressionChars "compression"
#define GIDI_diagonalChars "diagonal"
#define GIDI_startingIndices "startingIndices"
#define GIDI_flattenedChars "flattened"
#define GIDI_startsChars "starts"
#define GIDI_lengthsChars "lengths"
#define GIDI_embeddedChars "embedded"

#define GIDI_symmetryChars "symmetry"
#define GIDI_lowerChars "lower"
#define GIDI_upperChars "upper"

#define GIDI_permutationChars "permutation"
#define GIDI_plusOneChars "+1"
#define GIDI_minusOneChars "-1"

#define GIDI_storageOrderChars "storageOrder"
#define GIDI_rowMajorChars "rowMajor"
#define GIDI_columnMajorChars "columnMajor"

#define GIDI_offsetChars "offset"
#define GIDI_startIndexChars "startIndex"

// 1d Function monikers.
#define GIDI_constant1dChars "constant1d"
#define GIDI_XYs1dChars "XYs1d"
#define GIDI_Ys1dChars "Ys1d"
#define GIDI_polynomial1dChars "polynomial1d"
#define GIDI_LegendreChars "Legendre"
#define GIDI_regions1dChars "regions1d"
#define GIDI_gridded1dChars "gridded1d"
#define GIDI_referenceChars "reference"
#define GIDI_xs_pdf_cdf1dChars "xs_pdf_cdf1d"
#define GIDI_branching1dChars "branching1d"
#define GIDI_TNSL1dChars "thermalNeutronScatteringLaw1d"

// 2d Function monikers.
#define GIDI_XYs2dChars "XYs2d"
#define GIDI_recoilChars "recoil"
#define GIDI_isotropic2dChars "isotropic2d"
#define GIDI_discreteGammaChars "discreteGamma"
#define GIDI_primaryGammaChars "primaryGamma"
#define GIDI_generalEvaporationChars "generalEvaporation"
#define GIDI_simpleMaxwellianFissionChars "simpleMaxwellianFission"
#define GIDI_evaporationChars "evaporation"
#define GIDI_WattChars "Watt"
#define GIDI_MadlandNixChars "MadlandNix"
#define GIDI_weightedFunctionalsChars "weightedFunctionals"
#define GIDI_NBodyPhaseSpaceChars "NBodyPhaseSpace"
#define GIDI_regions2dChars "regions2d"
#define GIDI_gridded2dChars "gridded2d"

// 3d Function monikers.
#define GIDI_XYs3dChars "XYs3d"
#define GIDI_gridded3dChars "gridded3d"

// Double differentials
#define GIDI_optionsChars "options"
#define GIDI_S_alpha_betaChars "S_alpha_beta"
#define GIDI_S_tableChars "S_table"
#define GIDI_formFactorChars "formFactor"
#define GIDI_realAnomalousFactorChars "realAnomalousFactor"
#define GIDI_imaginaryAnomalousFactorChars "imaginaryAnomalousFactor"
#define GIDI_scatteringFactorChars "scatteringFactor"
#define GIDI_ComptonProfileChars "ComptonProfile"
#define GIDI_boundAtomCrossSectionChars "boundAtomCrossSection"
#define GIDI_characteristicCrossSectionChars "characteristicCrossSection"
#define GIDI_DebyeWallerIntegralChars "DebyeWallerIntegral"
#define GIDI_DebyeWallerChars "DebyeWaller"
#define GIDI_massChars "mass"
#define GIDI_freeAtomCrossSectionChars "freeAtomCrossSection"
#define GIDI_e_criticalChars "e_critical"
#define GIDI_e_maxChars "e_max"
#define GIDI_T_effectiveChars "T_effective"
#define GIDI_UChars "U"
#define GIDI_thetaChars "theta"
#define GIDI_gChars "g"

// Distribution forms.
#define GIDI_multiGroup3dChars "multiGroup3d"
#define GIDI_angularTwoBodyChars "angularTwoBody"
#define GIDI_uncorrelatedChars "uncorrelated"
#define GIDI_angularChars "angular"
#define GIDI_energyChars "energy"
#define GIDI_KalbachMannChars "KalbachMann"
#define GIDI_energyAngularChars "energyAngular"
#define GIDI_energyAngularMCChars "energyAngularMC"
#define GIDI_angularEnergyChars "angularEnergy"
#define GIDI_angularEnergyMCChars "angularEnergyMC"
#define GIDI_LLNLAngularEnergyChars "LLNLAngularEnergy"
#define GIDI_LLNLAngularOfAngularEnergyChars "LLNLAngularOfAngularEnergy"
#define GIDI_LLNLAngularEnergyOfAngularEnergyChars "LLNLAngularEnergyOfAngularEnergy"
#define GIDI_coherentPhotonScatteringChars "coherentPhotonScattering"
#define GIDI_incoherentPhotonScatteringChars "incoherentPhotonScattering"
#define GIDI_incoherentBoundToFreePhotonScatteringChars "incoherentBoundToFreePhotonScattering"
#define GIDI_TNSL_coherentElasticChars "thermalNeutronScatteringLaw_coherentElastic"
#define GIDI_TNSL_incoherentElasticChars "thermalNeutronScatteringLaw_incoherentElastic"
#define GIDI_TNSL_incoherentInelasticChars "thermalNeutronScatteringLaw_incoherentInelastic"
#define GIDI_thermalNeutronScatteringLawChars "thermalNeutronScatteringLaw"
#define GIDI_branching3dChars "branching3d"
#define GIDI_unspecifiedChars "unspecified"

#define GIDI_scatteringAtomsChars "scatteringAtoms"
#define GIDI_scatteringAtomChars "scatteringAtom"

#define GIDI_resonancesWithBackgroundChars "resonancesWithBackground"
#define GIDI_resonancesChars "resonances"
#define GIDI_resonanceBackground1dChars   "background"
#define GIDI_resolvedRegionChars "resolvedRegion"
#define GIDI_unresolvedRegionChars "unresolvedRegion"
#define GIDI_fastRegionChars "fastRegion"

#define GIDI_CoulombPlusNuclearElasticChars "CoulombPlusNuclearElastic"
#define GIDI_RutherfordScatteringChars "RutherfordScattering"
#define GIDI_nuclearPlusInterferenceChars "nuclearPlusInterference"

#define GIDI_URR_probabilityTables1dChars "URR_probabilityTables1d"
#define GIDI_LLNLLegendreChars "LLNLLegendre"

#define GIDI_axesChars "axes"
#define GIDI_axisChars "axis"
#define GIDI_gridChars "grid"
#define GIDI_fluxNodeChars "flux"

#define GIDI_function1dsChars "function1ds"
#define GIDI_function2dsChars "function2ds"
#define GIDI_uncertaintyChars "uncertainty"
#define GIDI_fChars "f"
#define GIDI_rChars "r"
#define GIDI_aChars "a"
#define GIDI_bChars "b"
#define GIDI_EFL_Chars "EFL"
#define GIDI_EFH_Chars "EFH"
#define GIDI_T_M_Chars "T_M"
#define GIDI_weightedChars "weighted"

#define GIDI_promptProductKEChars "promptProductKE"
#define GIDI_promptNeutronKEChars "promptNeutronKE"
#define GIDI_delayedNeutronKEChars  "delayedNeutronKE"
#define GIDI_promptGammaEnergyChars "promptGammaEnergy"
#define GIDI_delayedGammaEnergyChars "delayedGammaEnergy"
#define GIDI_delayedBetaEnergyChars "delayedBetaEnergy"
#define GIDI_neutrinoEnergyChars "neutrinoEnergy"
#define GIDI_nonNeutrinoEnergyChars "nonNeutrinoEnergy"
#define GIDI_totalEnergyChars "totalEnergy"

#define GIDI_trueChars "true"
#define GIDI_fissionGenreChars "fissionGenre"
#define GIDI_libraryChars "library"
#define GIDI_startChars "start"
#define GIDI_projectileChars "projectile"
#define GIDI_targetChars "target"
#define GIDI_evaluationChars "evaluation"
#define GIDI_interactionChars "interaction"
#define GIDI_standardTargetChars "standardTarget"
#define GIDI_standardEvaluationChars "standardEvaluation"
#define GIDI_projectileFrameChars "projectileFrame"
#define GIDI_ENDF_MT_Chars "ENDF_MT"
#define GIDI_dateChars "date"
#define GIDI_derivedFromChars "derivedFrom"
#define GIDI_versionChars "version"
#define GIDI_temperatureChars "temperature"
#define GIDI_muCutoffChars "muCutoff"
#define GIDI_lMaxChars "lMax"
#define GIDI_parametersChars "parameters"
#define GIDI_upperCalculatedGroupChars "upperCalculatedGroup"
#define GIDI_calculatedAtThermalChars "calculatedAtThermal"
#define GIDI_asymmetricChars "asymmetric"
#define GIDI_valueTypeChars "valueType"

#define GIDI_productFrameChars "productFrame"
#define GIDI_interpolationChars "interpolation"
#define GIDI_interpolationQualifierChars "interpolationQualifier"
#define GIDI_outerDomainValueChars "outerDomainValue"
#define GIDI_indexChars "index"
#define GIDI_labelChars "label"
#define GIDI_unitChars "unit"
#define GIDI_hrefChars "href"
#define GIDI_initialChars "initial"
#define GIDI_finalChars "final"
#define GIDI_minChars "min"
#define GIDI_maxChars "max"
#define GIDI_valueChars "value"
#define GIDI_domainMinChars "domainMin"
#define GIDI_domainMaxChars "domainMax"
#define GIDI_finalStateChars "finalState"
#define GIDI_numberOfProductsChars "numberOfProducts"
#define GIDI_pathChars "path"
#define GIDI_styleChars "style"
#define GIDI_genreChars "genre"
#define GIDI_processChars "process"
#define GIDI_pidChars "pid"
#define GIDI_countChars "count"

#define GIDI_inverseSpeedChars "inverseSpeed"

#define GIDI_centerOfMassChars "centerOfMass"
#define GIDI_labChars "lab"
#define GIDI_twoBodyChars "twoBody"
#define GIDI_NBodyChars "NBody"

// Allowed values for the 'conserve' attribute
#define GIDI_conserveNumberChars "number"
#define GIDI_conserveEnergyOutChars "energyOut"

// GRIN.
#define GIDI_GRIN_continuumGammasChars "GRIN_continuumGammas"
#define GIDI_captureNeutronSeparationEnergyChars "captureNeutronSeparationEnergy"
#define GIDI_maximumIncidentEnergyChars "maximumIncidentEnergy"
#define GIDI_inelasticIncidentEnergiesChars "inelasticIncidentEnergies"
#define GIDI_inelasticIncidentEnergyChars "inelasticIncidentEnergy"
#define GIDI_captureLevelProbabilitiesChars "captureLevelProbabilities"
#define GIDI_captureLevelProbabilityChars "captureLevelProbability"
#define GIDI_probabilityChars "probability"
#define GIDI_spinUnitChars "spinUnit"
#define GIDI_capturePrimaryToContinuaChars "capturePrimaryToContinua"

typedef std::pair<std::string, double> stringAndDoublePair;
typedef std::vector<stringAndDoublePair> stringAndDoublePairs;
typedef std::map<std::string, ParticleInfo> ParticleSubstitution;

#ifdef _WIN32
#define GIDI_FILE_SEPARATOR   "\\"
#else
#define GIDI_FILE_SEPARATOR   "/"
#endif

std::vector<std::string> vectorOfStrings( std::string const &a_string );

/*
============================================================
========================= Exception ========================
============================================================
*/
class Exception : public std::runtime_error {

    public :
        explicit Exception( std::string const &a_message );
};

namespace Construction {

/* *********************************************************************************************************
 * This enum allows a user to limit the data read in by various constructors. Limiting the data speeds up the reading
 * and parsing, and uses less memory.
 ***********************************************************************************************************/

enum class ParseMode : int { all,                             /**< Read and parse all data. */
                             multiGroupOnly,                  /**< Only read and parse data needed for multi-group transport. */
                             MonteCarloContinuousEnergy,      /**< Only read and parse data needed for continuous energy Monte Carlo. */
                             excludeProductMatrices,          /**< Read and parse all data but multi-group product matrices. */
                             readOnly,                        /**< Only read and parse all the data but do no calculations. Useful for reading an incomplete GNDS file. */
                             outline,                         /**< Does parse any component data (e.g., cross section, multiplicity, distribution). */
                             noParsing                        /**< Only timing development, must be used with caution. The mode allows one to test the io time without any parsing of the **HAPI::File**. */ };

enum class PhotoMode : int { nuclearAndAtomic,                /**< Instructs method Map::protare to create a Protare with both photo-nuclear and photo-atomic data when the projectile is photon. */
                             nuclearOnly,                     /**< Instructs method Map::protare to create a Protare with only photo-nuclear data when the projectile is photon. */
                             atomicOnly                       /**< Instructs method Map::protare to create a Protare with only photo-atomic data when the projectile is photon. */ };

/* *********************************************************************************************************
 * This enum specifies what fission redisual products will be added to the list of products produced in a fission reaction.
 ***********************************************************************************************************/

enum class FissionResiduals : int {
        none,               /**< No additional product is produced to the fission reaction. */
        ENDL99120,          /**< A LLNL ENDL 99120 fission product will be produced with multiplicity 2. */
        ENDL99125           /**< A LLNL ENDL 99125 fission product will be produced with multiplicity 2. */ };

/*
============================================================
========================= Settings =========================
============================================================
*/
class Settings {

    private:
        ParseMode m_parseMode;                                      /**< Parameter used by various constructors to limit data read into. */
        PhotoMode m_photoMode;                                      /**< Determines whether photo-nuclear and/or photo-atomic are included a Protare when the projectile is photon. */
        int m_useSystem_strtod;                                     /**< Flag passed to the function nfu_stringToListOfDoubles of the numericalFunctions library. */
        bool m_lazyParsing;                                         /**< It **true**, **Component** suites are lazy parsed. */
        bool m_decayPositronium;                                    /**< If **true**, whenever a positron is created, it is assumed to immediately form positronium and decay into 2 511 KeV photons. Ergo, the photons are produced in the reaction and not a positron. */
        bool m_usePhotoAtomicIncoherentDoppler;
        FissionResiduals m_fissionResiduals;                        /**< This member specifies what fission redisual products will be added to the list of products produced in a fission reaction. */
        bool m_GRIN_continuumGammas;                                /**< If true and institution/LLNL::GRIN_continuumGammas are loaded and used in MCGIDI. */

    public:
        Settings( ParseMode a_parseMode, PhotoMode a_photoMode );
        Settings( Settings const &a_settings );

        ParseMode parseMode( ) const { return( m_parseMode ); }             /**< Returns the value of the *m_parseMode* member. */

        PhotoMode photoMode( ) const { return( m_photoMode ); }             /**< Returns the value of the *m_photoMode* member. */
        void setPhotoMode( PhotoMode a_photoMode ) { m_photoMode = a_photoMode; }
                                                                            /**< Set the *m_photoMode* member to *a_photoMode*. */

        bool lazyParsing( ) const { return( m_lazyParsing ); }              /**< Returns the value of the *m_lazyParsing* member. */
        void setLazyParsing( bool a_lazyParsing ) { m_lazyParsing = a_lazyParsing; }
                                                                            /**< Set the *m_lazyParsing* member to *a_lazyParsing*. */

        bool decayPositronium( ) const { return( m_decayPositronium ); }    /**< Returns the value of the *m_decayPositronium* member. */
        void setDecayPositronium( bool a_decayPositronium ) { m_decayPositronium = a_decayPositronium; }
                                                                            /**< Set the *m_decayPositronium* member to *a_decayPositronium*. */

        FissionResiduals fissionResiduals( ) const { return( m_fissionResiduals ); }    /**< Returns the value of the *m_fissionResiduals* member. */
        void setFissionResiduals( FissionResiduals a_fissionResiduals ) { m_fissionResiduals = a_fissionResiduals ; }
                                                                            /**< Set the *m_fissionResiduals* member to *a_fissionResiduals*. */

        bool GRIN_continuumGammas( void ) const { return( m_GRIN_continuumGammas ); }   /**< Returns the value of the *m_GRIN_continuumGammas* member. */
        void setGRIN_continuumGammas( bool a_GRIN_continuumGammas ) { m_GRIN_continuumGammas = a_GRIN_continuumGammas; }
                                                                            /**< Set the *m_GRIN_continuumGammas* member to *a_GRIN_continuumGammas*. */

        int useSystem_strtod( ) const { return( m_useSystem_strtod ); }     /**< Returns the value of the *m_useSystem_strtod* member. */
        void setUseSystem_strtod( bool a_useSystem_strtod ) { m_useSystem_strtod = a_useSystem_strtod ? 1 : 0; }
                                                                            /**< Set the *m_useSystem_strtod* member to *a_useSystem_strtod*. */

        bool usePhotoAtomicIncoherentDoppler( ) const { return( m_usePhotoAtomicIncoherentDoppler ); }
        void setUsePhotoAtomicIncoherentDoppler( bool a_usePhotoAtomicIncoherentDoppler ) { m_usePhotoAtomicIncoherentDoppler = a_usePhotoAtomicIncoherentDoppler; }
};

}               // End namespace Construction.

/*
============================================================
========================= SetupInfo ========================
============================================================
*/
class SetupInfo {

    public:
        ProtareSingle *m_protare;
        ParticleSubstitution *m_particleSubstitution;
        LUPI::FormatVersion m_formatVersion;
        Styles::MultiGroup *m_multiGroup;
        Styles::HeatedMultiGroup *m_heatedMultiGroup;
        bool m_isENDL_C_9;
        int m_outputChannelLevel;
        std::string m_initialState;

        SetupInfo( ProtareSingle *a_protare ) :
                m_protare( a_protare ),
                m_particleSubstitution( nullptr ),
                m_formatVersion( ),
                m_multiGroup( nullptr ),
                m_heatedMultiGroup( nullptr ),
                m_isENDL_C_9( false ),
                m_outputChannelLevel( 0 ),
                m_initialState( "" ) {

        }

        SetupInfo( SetupInfo const &a_setupInfo ) :
                m_protare( a_setupInfo.m_protare ),
                m_particleSubstitution( a_setupInfo.m_particleSubstitution ),
                m_formatVersion( a_setupInfo.m_formatVersion ),
                m_multiGroup( a_setupInfo.m_multiGroup ),
                m_heatedMultiGroup( a_setupInfo.m_heatedMultiGroup ),
                m_isENDL_C_9( a_setupInfo.m_isENDL_C_9 ),
                m_outputChannelLevel( a_setupInfo.m_outputChannelLevel ),
                m_initialState( a_setupInfo.m_initialState ) {

        }
        ~SetupInfo( ) { }
};

/*
============================================================
=========================== Form ===========================
============================================================
*/
class Form : public GUPI::Ancestry {

    friend class Table::Column;

    private:
        Suite *m_parent;                        /**< The parent for the form. */
        FormType m_type;                        /**< The type of the form. */
        std::string m_keyName;                  /**< The key name used when the **Form** resides in a **Suite**. */
        mutable std::string m_keyValue;         /**< The value associated with the key. */
        std::string m_label;                    /**< The label for the form. */

    public:
        Form( FormType a_type );
        Form( std::string const &a_moniker, FormType a_type, std::string const &a_label );
        Form( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_suite = nullptr );
        Form( Form const &a_form );
        virtual ~Form( );
        Form &operator=( Form const &a_rhs );

        Suite *parent( ) const { return( m_parent ); }                                          /**< Returns the value of the *m_parent* member. */

        std::string const &label( ) const { return( m_label ); }                                /**< Returns the value of the *m_label* member. */
        void setLabel( std::string const &a_label );
        virtual std::string actualMoniker( ) const { return( moniker( ) ); }                    /**< Returns the value of the moniker. */

        std::string const &keyName( ) const ;
        void setKeyName( std::string const &a_keyName );
        std::string const &keyValue( ) const ;
        virtual void setKeyValue( std::string const &a_keyName ) const ;

        FormType type( ) const { return( m_type ); }                                            /**< Returns the value of the *m_type* member. */
        Form const *sibling( std::string a_label ) const ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }
        std::string xlinkItemKey( ) const {

            if( m_label == "" ) return( "" );
            return( buildXLinkItemKey( GIDI_labelChars, m_label ) );
        }                                                                                       /**< Returns the value of *this*'s key. */
};

/*
============================================================
================== LazyParsingHelperForm ===================
============================================================
*/
class LazyParsingHelperForm : public Form {

    private:
        Construction::Settings m_construction;
        HAPI::Node const m_node;
        SetupInfo m_setupInfo;
        PoPI::Database const *m_pops;
        PoPI::Database const *m_internalPoPs;
        std::string m_name;
        Styles::Suite const *m_styles;
        parseSuite m_parser;

    public:
        LazyParsingHelperForm( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node,
                SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name,
                Styles::Suite const *a_styles, parseSuite a_parser );
        ~LazyParsingHelperForm( );

        std::string actualMoniker( ) const { return( m_name ); }        /**< Returns the value of the *m_name* member which is the moniker of the actual form. */

        Form *parse( );
};

/*
============================================================
===================== PhysicalQuantity =====================
============================================================
*/
class PhysicalQuantity : public Form {

    private:
        double m_value;                                                 /**< The value for the physical quantity. */
        std::string m_unit;                                             /**< The unit for the physical quantity. */

    public:
        PhysicalQuantity( HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        PhysicalQuantity( double a_value, std::string a_unit );
        PhysicalQuantity( PhysicalQuantity const &a_physicalQuantity ) : 
                Form( FormType::physicalQuantity ),
                m_value( a_physicalQuantity.value( ) ),
                m_unit( a_physicalQuantity.unit( ) ) { }
        ~PhysicalQuantity( );
        PhysicalQuantity &operator=( PhysicalQuantity const &a_rhs );

        double value( ) const { return( m_value ); }                    /**< Returns the value of the *m_value* member. */
        std::string const &unit( ) const { return( m_unit ); }          /**< Returns the value of the *m_unit* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;

        friend std::ostream &operator<<( std::ostream& a_os, PhysicalQuantity const &a_physicalQuantity );
};

/*
============================================================
======================= ParticleInfo =======================
============================================================
*/
class ParticleInfo {

    public:
        static std::string const IDPortion( std::string const &a_id );
        static std::string const qualifierPortion( std::string const &a_id );

    private:
        std::string m_id;                                       /**< The particle's PoPs id. */
        std::string m_qualifier;                                /**< The particle's qualifier. For example the "1s1/2" in "Th{1s1/2}". */
        std::string m_pid;                                      /**< The same as *m_id* unless particle is an alias, then the final particle's id. */
        PhysicalQuantity m_mass;                                /**< The mass of the particle including nuclear excitation energy. */
        PhysicalQuantity m_excitationEnergy;                    /**< If the particle is a PoPI::Nuclide or PoPI::Nucleus, this is it nuclear excitation energy. Otherwise, it is 0. */

    public:
        ParticleInfo( std::string const &a_id, std::string const &a_pid, double a_mass, double a_excitationEnergy = 0.0 );
        ParticleInfo( std::string const &a_id, PoPI::Database const &a_globalPoPs, PoPI::Database const &a_internalPoPs, bool a_requiredInGlobalPoPs );
        ParticleInfo( ParticleInfo const &a_particleInfo );
        ParticleInfo &operator=( ParticleInfo const &a_rhs );

        std::string const &ID( ) const { return( m_id  ); }                     /**< Returns a const reference to *m_id* member. */
        std::string const &qualifier( ) const { return( m_qualifier ); }        /**< Returns a const reference to *m_qualifier**. */
        std::string const &pid( ) const { return( m_pid  ); }                   /**< Returns a const reference to *m_pid* member. */
        bool isAlias( ) const { return( m_pid != "" ); }                        /**< Returns true if particle id is an alias and false otherwise. */

        PhysicalQuantity const &mass( ) const { return( m_mass ); }             /**< Returns a const reference to *m_mass* member. */
        PhysicalQuantity const &excitationEnergy( ) const { return( m_excitationEnergy ); }     /**< Returns a const reference to *m_excitationEnergy* member. */
        double mass( std::string const &a_unit ) const ;
};

/*
============================================================
======================== AxisDomain ========================
============================================================
*/
class AxisDomain : public Form {

    private:
        double m_minimum;                                           /**< The minimum value for the domain. */
        double m_maximum;                                           /**< The maximum value for the domain. */
        std::string m_unit;                                         /**< The unit for the domain. */

    public:
        AxisDomain( HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        AxisDomain( double m_minimum, double m_maximum, std::string const &a_unit );
        ~AxisDomain( );

        double minimum( ) const { return( m_minimum ); }            /**< Returns the value of the *m_minimum* member. */
        double maximum( ) const { return( m_maximum ); }            /**< Returns the value of the *m_maximum* member. */
        std::string const &unit( ) const { return( m_unit ); }      /**< Returns the value of the *m_unit* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
=========================== Axis ===========================
============================================================
*/
class Axis : public Form {

    private:
        int m_index;                                                            /**< The index for the axis. */
        std::string m_unit;                                                     /**< The unit for the axis. */
        std::string m_href;                                                     /**< The **GNDS**'s href if instance points to another Axis or Grid instance. */

    public:
        Axis( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type = FormType::axis );
        Axis( int a_index, std::string a_label, std::string a_unit, FormType a_type = FormType::axis );
        Axis( Axis const &a_axis );
        virtual ~Axis( );

        int index( ) const { return( m_index ); }                               /**< Returns the value of the *m_index* member. */
        std::string const &unit( ) const { return( m_unit ); }                  /**< Returns the value of the *m_unit* member. */

        std::string const &href( ) const { return( m_href ); }                  /**< Returns the value of the *m_href* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Grid ===========================
============================================================
*/
class Grid : public Axis {

    private:
        std::string m_style;                                                        /**< The **GNDS grid**'s style. */
        std::string m_keyName;                                                      /**< **FIXME**. */
        std::string m_keyValue;                                                     /**< **FIXME**. */
        std::string m_valueType;                                                    /**< The type of data in m_values. Can be "Integer32". */
        std::string m_interpolation;
        nf_Buffer<double> m_values;                                                 /**< The **GNDS grid**'s values. */

    public:
        Grid( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod );
        Grid( Grid const &a_grid );

        std::size_t size( ) const { return( m_values.size( ) ); }                   /**< Returns the number of values in the *m_values* member. */
        inline double &operator[]( std::size_t a_index ) noexcept { return( m_values[a_index] ); }  /**< Returns the value at m_values[a_index]. */

        std::string const &style( ) const { return( m_style ); }                    /**< Returns the value of the *m_style* member. */
        std::string keyName( ) const { return( m_keyName ); }                       /**< Returns the value of the *m_keyName* member. */
        std::string keyValue( ) const { return( m_keyValue ); }                     /**< Returns the value of the *m_keyValue* member. */
        std::string valueType( ) const { return( m_valueType ); }                     /**< Returns the value of the *m_valueType* member. */

        std::string const &interpolation( ) const { return( m_interpolation ); }

        nf_Buffer<double> const &values( ) const { return( m_values ); }          /**< Returns the value of the *m_values* member. */
        nf_Buffer<double> const &data( ) const { return( m_values ); }            /**< Returns the value of the *m_values* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Axes ===========================
============================================================
*/
class Axes : public Form {

    private:
        std::vector<Axis *> m_axes;                                                 /**< Stores the list of Axis nodes. */

    public:
        Axes( );
        Axes( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod );
        Axes( Axes const &a_axes );
        ~Axes( );
        Axes &operator=( Axes const &a_rhs );

        std::size_t size( ) const { return( m_axes.size( ) ); }                     /**< Returns the number of *Axis* instances in *this*. */
        Axis const *operator[]( std::size_t a_index ) const { return( (m_axes[a_index]) ); }    /**< Returns m_axes[a_index]. */
        std::size_t dimension( ) const { return( m_axes.size( ) - 1 ); }            /**< Returns the dimension of the instance. */

        void append( Axis *a_axis ) { m_axes.push_back( a_axis ); }                 /**< Appends *a_axis* to the* list of *Axis* nodes. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;

        static Axes makeAxes( std::vector<std::pair<std::string, std::string>> const &a_labelsAndUnits );
};

namespace Array {

/*
============================================================
========================= FullArray ========================
============================================================
*/

class FullArray {

    public:
        FullArray( std::vector<int> a_shape );
        FullArray( std::vector<int> a_shape, std::vector<double> a_flattenedValues );
        ~FullArray( ) {}

        std::vector<int> m_shape;                           /**< The shape of the array. */
        std::vector<double> m_flattenedValues;              /**< A *std::vector<double>* representing the flattened arrary. */

        std::size_t size( ) const { return( m_flattenedValues.size( ) ); }
};

/*
============================================================
=========================== Array ==========================
============================================================
*/

class Array : public Form {

    private:
        std::vector<int> m_shape;               /**< The shape of the array. */
        std::string m_compression;              /**< The compression of the array. Allowed values are *none*, *diagonal*, *flattened* or *embedded*. */
        std::string m_symmetry;                 /**< The symmetry of the array. Allowed values are *none*, *lower* or *upper*. */
        std::string m_permutation;              /**< The permutation of the array. Allowed values are *none*, *-1*, and *1*. */
        std::string m_storageOrder;             /**< The storage order of the array. Allowed values are *row-major* or *colunn-major*. */
        nf_Buffer<double> m_values;             /**< The list of *values* of the array. */
        nf_Buffer<int> m_starts;                /**< If *compression* is *flattened*, this is the *starts* node. If *compression* is *diagonal* this is the *startingIndices* node. Otherwise, this member is not used. */
        nf_Buffer<int> m_length;                /**< If *compression* is *flattened*, this is the *lengths* node, otherwise, not used. */
        nf_Buffer<int> m_offset;                /**< The offset of the array. Only used if *this* array is embedded into another array. */
        std::vector<Array *> m_array;           /**< The list of embedded arrays. This is only used if *compression* is *embedded*. */

    public:
        Array( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod );
        ~Array( );

        std::size_t dimension( ) const { return( m_shape.size( ) ); }                   /**< Returns the dimension of the array. */
        std::size_t size( ) const ;
        std::vector<int> const &shape( ) const { return( m_shape ); }                   /**< Returns a const reference to member *m_shape*. */
        FullArray constructArray( ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

}               // End namespace Array.

/*
============================================================
====================== FlattenedArrayData ==================
============================================================
*/
class FlattenedArrayData : public Form {

    public:
        std::vector<int> m_shape;                                               /**< The shape of the flattened array. */
        std::size_t m_numberOfStarts;                                           /**< The number of start values. */
        std::size_t m_numberOfLengths;                                          /**< The number of length values. */
        nf_Buffer<int> m_starts;                                                /**< The start values. */
        nf_Buffer<int> m_lengths;                                               /**< The length values. */
        nf_Buffer<double> m_dValues;                                            /**< The given array data. */
//        int32_t *m_starts;                                                      /**< The start values. */
//        int32_t *m_lengths;                                                     /**< The length values. */
//        std::vector<double> m_dValues;                                          /**< The given array data. */

        FlattenedArrayData( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_dimensions, int a_useSystem_strtod );
        ~FlattenedArrayData( );

        std::vector<int> const &shape( ) const { return( m_shape ); }
        void setToValueInFlatRange( int a_start, int a_end, double a_value );
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================= Array3d ==========================
============================================================
*/
class Array3d : public Form {

    private:
        FlattenedArrayData m_array;                                             /**< The 3d array as a FlattenedArrayData instance. */

    public:
        Array3d( HAPI::Node const &a_node, SetupInfo &a_setupInfo, int a_useSystem_strtod );
        ~Array3d( );

        std::size_t size( ) const { return( m_array.m_shape.back( ) ); }        /**< The length of the 3d diminsion. */

        Matrix matrix( std::size_t a_index ) const ;

        void modifiedMultiGroupElasticForTNSL( int maxTNSL_index );
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const { m_array.toXMLList( a_writeInfo, a_indent ); }
};

namespace Functions {

/*
============================================================
====================== FunctionForm ======================
============================================================
*/
class FunctionForm : public Form {

    private:
        int m_dimension;                                    /**< The dimension of the function (i.e., the number of independent axes. */
        Axes m_axes;                                        /**< The axes node for the function. */
        ptwXY_interpolation m_interpolation;                /**< The interpolation for functions highest independent axis and its dependent axis. */
        std::string m_interpolationString;                  /**< The interpolation for functions highest independent axis and its dependent axis. */
        int m_index;                                        /**< Currently not used. */
        double m_outerDomainValue;                          /**< If function is part of a higher dimensional function, this is the next higher dimensions domain value. */

    public:
        FunctionForm( std::string const &a_moniker, FormType a_type, int a_dimension, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        FunctionForm( std::string const &a_moniker, FormType a_type, int a_dimension, Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        FunctionForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, int a_dimension, Suite *a_suite = nullptr );
        FunctionForm( FunctionForm const &a_form );
        ~FunctionForm( );
        FunctionForm &operator=( FunctionForm const &a_rhs );

        int dimension( ) const { return( m_dimension ); }                                       /**< Returns the value of the *m_dimension* member. */

        int index( ) const { return( m_index ); }                                               /**< Returns the value of the *m_index* member. */
        double outerDomainValue( ) const { return( m_outerDomainValue ); }                      /**< Returns the value of the *m_outerDomainValue* member. */
        void setOuterDomainValue( double a_outerDomainValue ) { m_outerDomainValue = a_outerDomainValue; }
        Axes const &axes( ) const { return( m_axes ); }                                         /**< Returns a const reference to the *m_axes* member. */
        Axes &axes( ) { return( m_axes ); }                                                     /**< Returns a reference to the *m_axes* member. */

        ptwXY_interpolation interpolation( ) const { return( m_interpolation ); }               /**< Returns the value of the *m_interpolation* member. */
        void setInterpolation( ptwXY_interpolation a_interpolation );
        std::string interpolationString( ) const { return( m_interpolationString ); }           /**< Returns the value of the *m_interpolationString* member. */

        virtual double domainMin( ) const = 0;
        virtual double domainMax( ) const = 0;

        virtual void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const { toXMLList_func( a_writeInfo, a_indent, false, false ); }
};

/*
============================================================
======================= Function1dForm =====================
============================================================
*/
class Function1dForm : public FunctionForm {

    public:
        Function1dForm( std::string const &a_moniker, FormType a_type, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function1dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function1dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_suite = nullptr );
        Function1dForm( Function1dForm const &a_form );
        ~Function1dForm( );
        Function1dForm &operator=( Function1dForm const &a_rhs );

        virtual double evaluate( double a_x1 ) const = 0;
        virtual void mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const ;
        virtual XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        virtual void write( FILE *a_file, std::string const &a_format ) const ;
        void print( std::string const &a_format ) const ;
};

/*
============================================================
========================= Constant1d =======================
============================================================
*/
class Constant1d : public Function1dForm {

    private:
        double m_value;                                                 /**< The constant value of the function. */
        double m_domainMin;                                             /**< The minimum domain value the function is valid. */
        double m_domainMax;                                             /**< The maximum domain value the function is valid. */

    public:
        Constant1d( Axes const &a_axes, double value, double a_domainMin, double a_domainMax, int a_index = 0, double a_outerDomainValue = 0.0 );
        Constant1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Constant1d( );

        double value( ) const { return( m_value ); }                    /**< Returns the value of the *m_value* member. */
        double domainMin( ) const { return( m_domainMin ); }            /**< Returns the value of the *m_domainMin* member. */
        double domainMax( ) const { return( m_domainMax ); }            /**< Returns the value of the *m_domainMax* member. */

        double evaluate( double a_x1 ) const ;
        void mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== XYs1d ===========================
============================================================
*/
class XYs1d : public Function1dForm {

    private:
        mutable ptwXYPoints *m_ptwXY;                                               /**< The ptwXYPoints instance that stores points and is used to do calculations. */

    public:
        XYs1d( );
        XYs1d( Axes const &a_axes, ptwXY_interpolation m_interpolation, int a_index = 0, double a_outerDomainValue = 0.0 );
        XYs1d( Axes const &a_axes, ptwXY_interpolation m_interpolation, std::vector<double> const &a_values, int a_index = 0, 
                double a_outerDomainValue = 0.0 );
        XYs1d( Axes const &a_axes, ptwXY_interpolation m_interpolation, std::vector<double> const &a_xs, 
                std::vector<double> const &a_ys, int a_index = 0, double a_outerDomainValue = 0.0 );
        XYs1d( Axes const &a_axes, ptwXYPoints *a_ptwXY, int a_index = 0, double a_outerDomainValue = 0.0 );
        XYs1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        XYs1d( XYs1d const &a_XYs1d );
        ~XYs1d( );
        XYs1d &operator=( XYs1d const &a_rhs );

        std::size_t size( ) const { return( ptwXY_length( nullptr, m_ptwXY ) ); }   /**< Returns the number of points (i.e., x,y pairs) in this. */
        ptwXYPoints const *ptwXY( ) const { return( m_ptwXY ); }                    /**< Returns the value of the *m_ptwXY* member. */
        ptwXYPoints *ptwXY( ) { return( m_ptwXY ); }                                /**< Returns the value of the *m_ptwXY* member. */

        std::pair<double, double> operator[]( std::size_t a_index ) const ;
        XYs1d operator+( XYs1d const &a_XYs1d ) const ;
        XYs1d &operator+=( XYs1d const &a_XYs1d );
        XYs1d operator-( XYs1d const &a_XYs1d ) const ;
        XYs1d &operator-=( XYs1d const &a_XYs1d );
        XYs1d operator*( double a_value ) const ;
        XYs1d operator*( XYs1d const &a_XYs1d ) const ;
        XYs1d &operator*=( double a_value );
        XYs1d &operator*=( XYs1d const &a_XYs1d );

        double domainMin( ) const { return( (*this)[0].first ); }                   /**< Returns first x1 value of this. */
        double domainMax( ) const { return( (*this)[size( )-1].first ); }           /**< Returns last x1 value of this. */
        std::vector<double> xs( ) const ;
        std::vector<double> ys( ) const ;
        std::vector<double> ysMappedToXs( std::vector<double> const &a_xs, std::size_t *a_offset ) const ;
        XYs1d domainSlice( double a_domainMin, double a_domainMax, bool a_fill ) const ;
        XYs1d domainSliceMax( double a_domainMax ) const ;

        double evaluate( double a_x1 ) const ;
        void mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        double integrate( double a_dommainMin, double a_dommainMax );
        double normalize( );
        Xs_pdf_cdf1d toXs_pdf_cdf1d( );

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
        void write( FILE *a_file, std::string const &a_format ) const ;

        static XYs1d *makeConstantXYs1d( Axes const &a_axes, double a_domainMin, double a_domainMax, double a_value );
};

/*
============================================================
=========================== Ys1d ===========================
============================================================
*/
class Ys1d : public Function1dForm {

    private:
        std::size_t m_start;                                /**< The index in the grid for the x1 value for the first y value in *m_Ys*. */
        std::vector<double> m_Ys;                           /**< This list of y values. */

    public:
        Ys1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index = 0, double a_outerDomainValue = 0.0 );
        Ys1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, std::size_t a_start, std::vector<double> const &a_Ys, int a_index = 0, double a_outerDomainValue = 0.0 );
        Ys1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        Ys1d( Ys1d const &a_Ys1d );
        ~Ys1d( );

        std::size_t size( ) const { return( m_Ys.size( ) ); }                           /**< Returns the number of values in *m_Ys*. */

        double operator[]( std::size_t a_index ) const { return( m_Ys[a_index] ); }     /**< Returns the y value at *m_Ys*[a_index]. */
        void push_back( double a_y ) { m_Ys.push_back( a_y ); }
        Ys1d operator+( Ys1d const &a_Ys1d ) const ;
        Ys1d &operator+=( Ys1d const &a_Ys1d );

        double domainMin( ) const ;
        double domainMax( ) const ;
        std::size_t start( ) const { return( m_start ); }                               /**< Returns the value of the *m_start* member. */
        void setStart( std::size_t a_start ) { m_start = a_start; }                     /**< Sets the *m_start* member to *a_start. */
        std::size_t length( ) const { return( m_start + m_Ys.size( ) ); }               /**< Returns the sum of m_start and size( ). */
        std::vector<double> const &Ys( ) const { return( m_Ys ); }                      /**< Returns a reference to the list of y-values. */
        std::vector<double> &Ys( ) { return( m_Ys ); }                                  /**< Returns a reference to the list of y-values. */

        double evaluate( double a_x1 ) const ;
        void set( std::size_t a_index, double a_value ) { m_Ys[a_index] = a_value; }    /**< Set the value at *m_Ys*[a_index] to a_value. */

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
        void write( FILE *a_file, std::string const &a_format ) const ;
};

/*
============================================================
======================= Polynomial1d =======================
============================================================
*/
class Polynomial1d : public Function1dForm {

    private:
        double m_domainMin;                                             /**< The minimum domain value the function is valid. */
        double m_domainMax;                                             /**< The maximum domain value the function is valid. */
        std::vector<double> m_coefficients;                             /**< The coefficients of the polynomial. */

    public:
        Polynomial1d( Axes const &a_axes, double a_domainMin, double a_domainMax, std::vector<double> const &a_coefficients, int a_index = 0, double a_outerDomainValue = 0.0 );
        Polynomial1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        Polynomial1d( Polynomial1d const &a_polynomial1d );
        ~Polynomial1d( );

        double domainMin( ) const { return( m_domainMin ); }                                /**< Returns the value of the *m_domainMin* member. */
        double domainMax( ) const { return( m_domainMax ); }                                /**< Returns the value of the *m_domainMax* member. */

        std::vector<double> const &coefficients( ) const { return( m_coefficients ); }      /**< Returns the value of the *m_coefficients* member. */

        double evaluate( double a_x1 ) const ;
        void mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================= Legendre1d =======================
============================================================
*/
class Legendre1d : public Function1dForm {

    private:
        std::vector<double> m_coefficients;                                             /**< The Legendre coefficients. */

    public:
        Legendre1d( Axes const &a_axes, int a_index = 0, double a_outerDomainValue = 0.0 );
        Legendre1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        Legendre1d( Legendre1d const &a_Legendre1d );
        ~Legendre1d( );

        double domainMin( ) const { return( -1.0 ); }                                   /**< Returns the value of the *domainMin* which is always -1.0. */
        double domainMax( ) const { return( 1.0 ); }                                    /**< Returns the value of the *domainMax* which is always 1.0. */

        std::vector<double> const &coefficients( ) const { return( m_coefficients ); }  /**< Returns the value of the *m_coefficients* member. */
        std::vector<double> &coefficients( ) { return( m_coefficients ); }              /**< Returns the value of the *m_coefficients* member. */

        double evaluate( double a_x1 ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== Gridded1d =======================
============================================================
*/
class Gridded1d : public Function1dForm {

    private:
        Vector m_grid;                                                          /**< The grid for the gridded 1d function. Can be a link. */
        Vector m_data;                                                          /**< The value of the function on the grid. */
// BRB should have <array compression="flattened"> ... instead of m_data.

    public:
        Gridded1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Gridded1d( );

        double domainMin( ) const { return( m_grid[0] ); }                      /**< Returns the value of the *domainMin*. */
        double domainMax( ) const { return( m_grid[m_grid.size( )-1] ); }       /**< Returns the value of the *domainMax*. */

        Vector const &grid( ) const { return( m_grid ); }                       /**< Returns the value of the *m_grid* member. */
        Vector const &data( ) const { return( m_data ); }                       /**< Returns the value of the *m_data* member. */
        void setData( Vector const &a_data ) { m_data = a_data; }               /**< Sets the *m_data* member to *a_data*. */

        void modifiedMultiGroupElasticForTNSL( int a_maxTNSL_index );
        double evaluate( double a_x1 ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
        void write( FILE *a_file, std::string const &a_format ) const ;
};

/*
============================================================
======================== Reference1d =======================
============================================================
*/
class Reference1d : public Function1dForm {

    private:
        std::string m_xlink;                                                    /**< Link to the other function. */

    public:
        Reference1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Reference1d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        std::string const &xlink( ) const { return( m_xlink ); }                /**< Returns the value of the *m_xlink* member. */
        double evaluate( double a_x1 ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= Xs_pdf_cdf1d =======================
============================================================
*/
class Xs_pdf_cdf1d : public Function1dForm {

    private:
        std::vector<double> m_xs;                                               /**< List of x1 values. */
        std::vector<double> m_pdf;                                              /**< The pdf evaluated at the x1 values. */
        std::vector<double> m_cdf;                                              /**< The cdf evaluated at the x1 values. */
// BRB m_xs, m_pdf and m_cdf need to be a class like ListOfDoubles.

    public:
        Xs_pdf_cdf1d( );
        Xs_pdf_cdf1d( Axes const &a_axes, ptwXY_interpolation a_interpolation, std::vector<double> const &a_Xs, 
                std::vector<double> const &a_pdf, std::vector<double> const &a_cdf, int a_index = 0, double a_outerDomainValue = 0.0 );
        Xs_pdf_cdf1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Xs_pdf_cdf1d( );
        Xs_pdf_cdf1d &operator=( Xs_pdf_cdf1d const &a_rhs );

        double domainMin( ) const { return( m_xs[0] ); }                        /**< Returns the value of the *domainMin*. */
        double domainMax( ) const { return( m_xs[m_xs.size( )-1] ); }           /**< Returns the value of the *domainMax*. */

        std::vector<double> const &Xs( ) const { return( m_xs ); }              /**< Returns the value of the *m_xs* member. */
        std::vector<double> const &pdf( ) const { return( m_pdf ); }            /**< Returns the value of the *m_pdf* member. */
        std::vector<double> const &cdf( ) const { return( m_cdf ); }            /**< Returns the value of the *m_cdf* member. */
        double evaluate( double a_x1 ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== Regions1d =======================
============================================================
*/
class Regions1d : public Function1dForm {

    private:
        std::vector<double> m_Xs;                                               /**< List of *x1* domain values that bounds each region. */
        std::vector<Function1dForm *> m_function1ds;                            /**< List of regions. */

    public:
        Regions1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Regions1d( );

        std::size_t size( ) const { return( m_function1ds.size( ) ); }                          /**< Returns number of regions. */
        Function1dForm const *operator[]( std::size_t a_index ) const { return( m_function1ds[a_index] ); } /**< Returns the region at index *a_index* - 1. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        void append( Function1dForm *a_function );
        double evaluate( double a_x1 ) const ;
        void mapToXsAndAdd( int a_offset, std::vector<double> const &a_Xs, std::vector<double> &a_results, double a_scaleFactor ) const ;
        XYs1d *asXYs1d( bool a_asLinlin, double a_accuray, double a_lowerEps, double a_upperEps ) const ;

        std::vector<double> const &Xs( ) const { return( m_Xs ); }                              /**< Returns the value of the *m_Xs* member. */
        std::vector<Function1dForm *> const &function1ds( ) const { return( m_function1ds ); }  /**< Returns the value of the *m_function1ds* member. */
        std::vector<Function1dForm *> &function1ds( ) { return( m_function1ds ); }              /**< Returns the value of the *m_function1ds* member. */

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
        void write( FILE *a_file, std::string const &a_format ) const ;
};

/*
============================================================
========================= Branching1d ======================
============================================================
*/
class Branching1d : public Function1dForm {

    private:
        std::string m_initialState;                                         /**< The nuclide level that decays, emitting a photon. */
        double m_multiplicity;                                              /**< The average number of photons emitted when transitioning from the initial to the final state. */

    public:
        Branching1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Branching1d( );

        std::string const &initialState( ) const { return( m_initialState ); }      /**< Returns the value of the *m_initialState* member. */

        double multiplicity( ) const { return( m_multiplicity ); }                  /**< Returns the value of the *m_multiplicity* member. */
        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;
};

/*
============================================================
================ ResonanceBackgroundRegion1d ===============
============================================================
*/

class ResonanceBackgroundRegion1d : public Function1dForm {

    private:
        Function1dForm *m_function1d;                                           /**< The 1-d function representing *this*. */

    public:
        ResonanceBackgroundRegion1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~ResonanceBackgroundRegion1d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
=================== ResonanceBackground1d ==================
============================================================
*/
class ResonanceBackground1d : public Function1dForm {

    private:
        ResonanceBackgroundRegion1d *m_resolvedRegion;                      /**< The 1-d function for the resolved region. */
        ResonanceBackgroundRegion1d *m_unresolvedRegion;                    /**< The 1-d function for the unresolved region. */
        ResonanceBackgroundRegion1d *m_fastRegion;                          /**< The 1-d function for the fast region. */

    public:
        ResonanceBackground1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~ResonanceBackground1d( );

        Function1dForm const *resolvedRegion( ) const { return( m_resolvedRegion ); }       /**< Returns the value of the *m_resolvedRegion* member. */
        Function1dForm const *unresolvedRegion( ) const { return( m_unresolvedRegion ); }   /**< Returns the value of the *m_unresolvedRegion* member. */
        Function1dForm const *fastRegion( ) const { return( m_fastRegion ); }               /**< Returns the value of the *m_fastRegion* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
================= ResonancesWithBackground1d ===============
============================================================
*/
class ResonancesWithBackground1d : public Function1dForm {

    private:
        std::string m_resonances;                                                           /**< The reference to the resonance data for *this*.*/
        ResonanceBackground1d m_background;                                                 /**< The background .*/

    public:
        ResonancesWithBackground1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~ResonancesWithBackground1d( );

        double domainMin( ) const { return( m_background.domainMin( ) ); }                  /**< Returns *this* function's domain mimimun value. */
        double domainMax( ) const { return( m_background.domainMax( ) ); }                  /**< Returns *this* function's domain mimimun value. */

        double evaluate( double a_x1 ) const { return( m_background.evaluate( a_x1 ) ); }   /**< Returns the value *this* evaluated at *a_x1*. */

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
================= URR_probabilityTables1d ==================
============================================================
*/
class URR_probabilityTables1d : public Function1dForm {

    private:
        Function2dForm *m_function2d;                                                       /**< The URR probability tables. */

    public:
        URR_probabilityTables1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~URR_probabilityTables1d( );

        Function2dForm const *function2d( ) const { return( m_function2d ); }               /**< Returns the pointer to the *m_function2d* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
=============== ThermalNeutronScatteringLaw1d ================
============================================================
*/
class ThermalNeutronScatteringLaw1d : public Function1dForm {

    private:
        std::string m_href;                                                 /**< xlink to the IncoherentPhotoAtomicScattering instance under the *m_doubleDifferentialCrossSection* node. */

    public:
        ThermalNeutronScatteringLaw1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~ThermalNeutronScatteringLaw1d( );

        std::string const &href( ) const { return( m_href ); }              /**< Returns the value of the *m_href* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;
};

/*
============================================================
====================== Unspecified1d =======================
============================================================
*/
class Unspecified1d : public Function1dForm {

    public:
        Unspecified1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Unspecified1d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x1 ) const ;

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
======================= Function2dForm =====================
============================================================
*/
class Function2dForm : public FunctionForm {

    public:
        Function2dForm( std::string const &a_moniker, FormType a_type, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function2dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function2dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_suite = nullptr );
        Function2dForm( Function2dForm const &a_form );
        ~Function2dForm( );

        virtual double evaluate( double a_x2, double a_x1 ) const = 0;
};

/*
============================================================
=========================== XYs2d ==========================
============================================================
*/
class XYs2d : public Function2dForm {

    private:
        std::string m_interpolationQualifier;                       /**< The interpolation qualifier for *this*. */
        std::vector<double> m_Xs;                                   /**< The list of *x2* values for each function. */
        std::vector<Function1dForm *> m_function1ds;                /**< The list of 1d functions. */

    public:
        XYs2d( Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index = 0, double a_outerDomainValue = 0.0 );
        XYs2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~XYs2d( );

        std::string interpolationQualifier( ) const { return( m_interpolationQualifier ); }         /**< Returns the value of the *m_interpolationQualifier* member. */
        void setInterpolationQualifier( std::string a_interpolationQualifier ) { m_interpolationQualifier = a_interpolationQualifier; }
                                                                                                    /**< Sets the *m_interpolationQualifier* member to *a_interpolationQualifier*. */

        double domainMin( ) const ;
        double domainMax( ) const ;
        double evaluate( double a_x2, double a_x1 ) const ;

        std::vector<double> const &Xs( ) const { return( m_Xs ); }                                  /**< Returns the value of the *m_Xs* member. */
        std::vector<Function1dForm *> const &function1ds( ) const { return( m_function1ds ); }      /**< Returns the value of the *m_function1ds* member. */
        std::vector<Function1dForm *>       &function1ds( )       { return( m_function1ds ); }      /**< Returns the value of the *m_function1ds* member. */
        void append( Function1dForm *a_function1d );

        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================= Recoil2d =========================
============================================================
*/
class Recoil2d : public Function2dForm {

    private:
        std::string m_xlink;                                        /**< Link to the recoil product. */

    public:
        Recoil2d( std::string const &a_label, std::string const &a_href );
        Recoil2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Recoil2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        std::string const &xlink( ) const { return( m_xlink ); }    /**< Returns the value of the *m_xlink* member. */
        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================= Isotropic2d ======================
============================================================
*/
class Isotropic2d : public Function2dForm {

    public:
        Isotropic2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Isotropic2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, LUPI_maybeUnused bool a_embedded, LUPI_maybeUnused bool a_inRegions ) const { 
                a_writeInfo.addNodeStarterEnder( a_indent, moniker( ) ); }
};

/*
============================================================
======================= DiscreteGamma2d ====================
============================================================
*/
class DiscreteGamma2d : public Function2dForm {

    private:
        double m_domainMin;                                             /**< The minimum domain value the function is valid. */
        double m_domainMax;                                             /**< The maximum domain value the function is valid. */
        double m_value;                                                 /**< The energy of the discrete gamma. */

    public:
        DiscreteGamma2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~DiscreteGamma2d( );

        double domainMin( ) const { return( m_domainMin ); }            /**< Returns the value of the *m_domainMin* member. */
        double domainMax( ) const { return( m_domainMax ); }            /**< Returns the value of the *m_domainMax* member. */
        double value( ) const { return( m_value ); }                    /**< Returns the value of the *m_value* member. */

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
======================= PrimaryGamma2d ====================
============================================================
*/
class PrimaryGamma2d : public Function2dForm {

    private:
        double m_domainMin;                                             /**< The minimum domain value the function is valid. */
        double m_domainMax;                                             /**< The maximum domain value the function is valid. */
        double m_value;                                                 /**< The binding energy needed to calculate the energy of the primary gamma. */
        std::string m_finalState;                                       /**< The nuclear state the compound is in after the primary photon (gamma) is emitted. */

    public:
        PrimaryGamma2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~PrimaryGamma2d( );

        double domainMin( ) const { return( m_domainMin ); }            /**< Returns the value of the *m_domainMin* member. */
        double domainMax( ) const { return( m_domainMax ); }            /**< Returns the value of the *m_domainMax* member. */
        double value( ) const { return( m_value ); }                    /**< Returns the value of the *m_value* member. */
        std::string const &finalState( ) const { return( m_finalState ); }  /**< Returns a reference to the *m_finalState* member. */

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
=================== GeneralEvaporation2d ===================
============================================================
*/
class GeneralEvaporation2d : public Function2dForm {

    private:
        PhysicalQuantity m_U;                                           /**< The *U* value for the general evaporation function. */
        Function1dForm *m_theta;                                        /**< The *theta* function for the general evaporation function. */
        Function1dForm *m_g;                                            /**< The *g* function for the general evaporation function. */

    public:
        GeneralEvaporation2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~GeneralEvaporation2d( );

        double U( ) const { return( m_U.value( ) ); }                   /**< Returns the GNDS *U* value for *this*. */
        Function1dForm const *theta( ) const { return( m_theta ); }     /**< Returns the value of the *m_theta* member. */
        Function1dForm const *g( ) const { return( m_g ); }             /**< Returns the value of the *m_g* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
================= SimpleMaxwellianFission2d ================
============================================================
*/
class SimpleMaxwellianFission2d : public Function2dForm {

    private:
        PhysicalQuantity m_U;                                           /**< The *U* value for the simple Maxwellian function. */
        Function1dForm *m_theta;                                        /**< The *theta* function for the simple Maxwellian function. */

    public:
        SimpleMaxwellianFission2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~SimpleMaxwellianFission2d( );

        double U( ) const { return( m_U.value( ) ); }                   /**< Returns the GNDS *U* value for *this*. */
        Function1dForm const *theta( ) const { return( m_theta ); }     /**< Returns the value of the *m_theta* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
====================== Evaporation2d =======================
============================================================
*/
class Evaporation2d : public Function2dForm {

    private:
        PhysicalQuantity m_U;                                           /**< The *U* value for the evaporation function. */
        Function1dForm *m_theta;                                        /**< The *theta* function for the evaporation function. */

    public:
        Evaporation2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Evaporation2d( );

        double U( ) const { return( m_U.value( ) ); }                   /**< Returns the *m_U* value for *this*. */
        Function1dForm const *theta( ) const { return( m_theta ); }     /**< Returns the value of the *m_theta* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== Watt2d ==========================
============================================================
*/
class Watt2d : public Function2dForm {

    private:
        PhysicalQuantity m_U;                                           /**< The *U* value for the Watt function. */
        Function1dForm *m_a;                                            /**< The *a* function for the Watt function. */
        Function1dForm *m_b;                                            /**< The *b* function for the Watt function. */

    public:
        Watt2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Watt2d( );

        double U( ) const { return( m_U.value( ) ); }                   /**< Returns the GNDS *U* value for *this*. */
        Function1dForm const *a( ) const { return( m_a ); }             /**< Returns the value of the *m_a* member. */
        Function1dForm const *b( ) const { return( m_b ); }             /**< Returns the value of the *m_b* member. */

        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
======================= MadlandNix2d =======================
============================================================
*/
class MadlandNix2d : public Function2dForm {

    private:
        PhysicalQuantity m_EFL;                                         /**< The *EFL* value for the Madland/Nix function. */
        PhysicalQuantity m_EFH;                                         /**< The *EFH* value for the Madland/Nix function. */
        Function1dForm *m_T_M;                                          /**< The *T_M* function for the Madland/Nix function. */

    public:
        MadlandNix2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~MadlandNix2d( );

        double EFL( ) const { return( m_EFL.value( ) ); }               /**< Returns the GNDS *EFL* value for *this*. */
        double EFH( ) const { return( m_EFH.value( ) ); }               /**< Returns the GNDS *EFH* value for *this*. */
        Function1dForm const *T_M( ) const { return( m_T_M ); }         /**< Returns the value of the *m_T_M* member. */


        double domainMin( ) const ;
        double domainMax( ) const ;

        double evaluate( double a_x2, double a_x1 ) const ;
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
=================== Weighted_function2d ====================
============================================================
*/
class Weighted_function2d : public Function2dForm {

    private:
        Function1dForm *m_weight;                                       /**< The weight for this function. */
        Function2dForm *m_energy;                                       /**< The energy functional. */

    public:
        Weighted_function2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Weighted_function2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        Function1dForm const *weight( ) const { return( m_weight ); }       /**< Returns the value of the *m_weight* member. */
        Function2dForm const *energy( ) const { return( m_energy ); }       /**< Returns the value of the *m_energy* member. */
        double evaluate( double a_x2, double a_x1 ) const ;
};

/*
============================================================
================== WeightedFunctionals2d ===================
============================================================
*/
class WeightedFunctionals2d : public Function2dForm {

    private:
        std::vector<Weighted_function2d *> m_weighted_function2d;       /**< The list of Weighted_function2d. */

    public:
        WeightedFunctionals2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~WeightedFunctionals2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        std::vector<Weighted_function2d *> const &weighted_function2d( ) const { return( m_weighted_function2d ); } /**< Returns the value of the *m_weighted_function2d* member. */
        double evaluate( double a_x2, double a_x1 ) const ;
};

/*
============================================================
==================== NBodyPhaseSpace2d =====================
============================================================
*/
class NBodyPhaseSpace2d : public Function2dForm {

    private:
        int m_numberOfProducts;                                         /**< The number of products for the NBodyPhaseSpace function. */
        PhysicalQuantity m_mass;                                        /**< The mass for the NBodyPhaseSpace function. */

    public:
        NBodyPhaseSpace2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~NBodyPhaseSpace2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        int numberOfProducts( ) const { return( m_numberOfProducts ); }     /**< Returns the value of the *m_numberOfProducts* member. */
        PhysicalQuantity const &mass( ) const { return( m_mass ); }         /**< Returns the value of the *m_mass* member. */

        double evaluate( double a_x2, double a_x1 ) const ;
};

/*
============================================================
========================== Regions2d =======================
============================================================
*/
class Regions2d : public Function2dForm {

    private:
        std::vector<double> m_Xs;                                           /**< List of *x2* domain values that bounds each region. */
        std::vector<Function2dForm *> m_function2ds;                        /**< List of 2d regions. */

    public:
        Regions2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Regions2d( );

        double domainMin( ) const ;
        double domainMax( ) const ;

        void append( Function2dForm *a_function );
        double evaluate( double a_x2, double a_x1 ) const ;

        std::vector<double> const &Xs( ) const { return( m_Xs ); }                                  /**< Returns the value of the *m_Xs* member. */
        std::vector<Function2dForm *> const &function2ds( ) const { return( m_function2ds ); }      /**< Returns the value of the *m_function2ds* member. */
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== Gridded2d =======================
============================================================
*/
class Gridded2d : public Function2dForm {

    private:
        Array::Array m_array;                                                               /**< The multi-group transfer matrix. */

    public:
        Gridded2d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Gridded2d( );

        double domainMin( ) const { return( 0.0 ); }                                        /**< Not properly implemented. */
        double domainMax( ) const { return( 0.0 ); }                                        /**< Not properly implemented. */
        double evaluate( LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const { return( 0.0 ); }                /**< Not properly implemented. */

        Array::Array const &array( ) const { return( m_array ); }                           /**< Returns the value of the *m_array* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= Function3dForm =====================
============================================================
*/
class Function3dForm : public FunctionForm {

    public:
        Function3dForm( std::string const &a_moniker, FormType a_type, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function3dForm( std::string const &a_moniker, FormType a_type, Axes const &a_axes, ptwXY_interpolation a_interpolation, int a_index, double a_outerDomainValue );
        Function3dForm( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_suite = nullptr );
        Function3dForm( Function3dForm const &a_form );
        ~Function3dForm( );

        virtual double evaluate( double a_x3, double a_x2, double a_x1 ) const = 0;
};

/*
============================================================
=========================== XYs3d ==========================
============================================================
*/
class XYs3d : public Function3dForm {

    private:
        std::string m_interpolationQualifier;                       /**< The interpolation qualifier for *this*. */
        std::vector<double> m_Xs;                                   /**< The list of *x3* values for each function. */
        std::vector<Function2dForm *> m_function2ds;                /**< The list of 2d functions. */

    public:
        XYs3d( Axes const &a_axes, ptwXY_interpolation a_interpolation = ptwXY_interpolationLinLin, int a_index = 0, double a_outerDomainValue = 0.0 );
        XYs3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~XYs3d( );

        std::string interpolationQualifier( ) const { return( m_interpolationQualifier ); }         /**< Returns the value of the *m_interpolationQualifier* member. */
        void setInterpolationQualifier( std::string a_interpolationQualifier ) { m_interpolationQualifier = a_interpolationQualifier; }
                                                                                                    /**< Sets the *m_interpolationQualifier* member to *a_interpolationQualifier*. */

        double domainMin( ) const ;
        double domainMax( ) const ;
        double evaluate( double a_x3, double a_x2, double a_x1 ) const ;

        std::vector<double> const &Xs( ) const { return( m_Xs ); }                                  /**< Returns the value of the *m_Xs* member. */
        std::vector<Function2dForm *> const &function2ds( ) const { return( m_function2ds ); }      /**< Returns a const reference to the *m_function2ds* member. */

        void append( Function2dForm *a_function2d );                                                /**< Appends the 2d function *a_function2d* to the end the *this*. */
        void toXMLList_func( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, bool a_embedded, bool a_inRegions ) const ;
};

/*
============================================================
========================== Gridded3d =======================
============================================================
*/
class Gridded3d : public Function3dForm {

    private:
        Array3d m_data;                                                                 /**< The multi-group transfer matrix. */

    public:
        Gridded3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Gridded3d( );

        double domainMin( ) const { return( 0.0 ); }                                        /**< Not properly implemented. */
        double domainMax( ) const { return( 0.0 ); }                                        /**< Not properly implemented. */
        double evaluate( LUPI_maybeUnused double a_x3, LUPI_maybeUnused double a_x2, LUPI_maybeUnused double a_x1 ) const { return( 0.0 ); }   /**< Not properly implemented. */

        Array3d const &data( ) const { return( m_data ); }                                  /**< Returns the value of the *m_data* member. */

        void modifiedMultiGroupElasticForTNSL( int maxTNSL_index );
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

}               // End namespace Functions.

/*
============================================================
=========== DoubleDifferentialCrossSection stuff ===========
============================================================
*/
namespace DoubleDifferentialCrossSection {

/*
============================================================
============================= Base =========================
============================================================
*/
class Base : public Form {

    public:
        Base( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_parent );
};

/*
============================================================
================ CoherentPhotoAtomicScattering =============
============================================================
*/
class CoherentPhotoAtomicScattering : public Base {

    private:
        Functions::Function1dForm *m_formFactor;                                   /**< The form factor for coherent photo-atomic scattering. */
        Functions::Function1dForm *m_realAnomalousFactor;                          /**< The real anomalous factor of coherent photo-atomic scattering. */
        Functions::Function1dForm *m_imaginaryAnomalousFactor;                     /**< The imaginary anomalous factor of coherent photo-atomic scattering. */

    public:
        CoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~CoherentPhotoAtomicScattering( );

        Functions::Function1dForm *formFactor( ) { return( m_formFactor ); }                                            /**< Returns the value of the *m_formFactor* member. */
        Functions::Function1dForm const *formFactor( ) const { return( m_formFactor ); }                                /**< Returns the value of the *m_formFactor* member. */
        Functions::Function1dForm *realAnomalousFactor( ) { return( m_realAnomalousFactor ); }                          /**< Returns the value of the *m_realAnomalousFactor* member. */
        Functions::Function1dForm const *realAnomalousFactor( ) const { return( m_realAnomalousFactor ); }              /**< Returns the value of the *m_realAnomalousFactor* member. */
        Functions::Function1dForm *imaginaryAnomalousFactor( ) { return( m_imaginaryAnomalousFactor ); }                /**< Returns the value of the *m_imaginaryAnomalousFactor* member. */
        Functions::Function1dForm const *imaginaryAnomalousFactor( ) const { return( m_imaginaryAnomalousFactor ); }    /**< Returns the value of the *m_imaginaryAnomalousFactor* member. */
};

/*
============================================================
============== IncoherentPhotoAtomicScattering =============
============================================================
*/
class IncoherentPhotoAtomicScattering : public Base {

    private:
        Functions::Function1dForm *m_scatteringFactor;                          /**< The scattering factor for incoherent photo-atomic scattering. */

    public:
        IncoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, 
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~IncoherentPhotoAtomicScattering( );

        Functions::Function1dForm const *scatteringFactor( ) const { return( m_scatteringFactor); }     /**< Returns the value of the *m_scatteringFactor* member. */
};

/*
=======================================================================
============== IncoherentBoundToFreePhotoAtomicScattering =============
=======================================================================
*/
class IncoherentBoundToFreePhotoAtomicScattering : public Base {

    private:
        Functions::Function1dForm *m_ComptonProfile;

    public:
        IncoherentBoundToFreePhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~IncoherentBoundToFreePhotoAtomicScattering( );

        Functions::Function1dForm *ComptonProfile( ) { return( m_ComptonProfile); }
        Functions::Function1dForm const *ComptonProfile( ) const { return( m_ComptonProfile); }
};

namespace n_ThermalNeutronScatteringLaw {

/*
============================================================
========================== S_table =========================
============================================================
*/
class S_table : public Form {


    private:
        Functions::Function2dForm *m_function2d;           /**< The cumulative scattering factor \f$S(T,E)\f$. */

    public:
        S_table( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~S_table( );

        Functions::Function2dForm *function2d( ) { return( m_function2d ); }                /**< Returns the value of the *m_function2d* member. */
        Functions::Function2dForm const *function2d( ) const { return( m_function2d ); }    /**< Returns the value of the *m_function2d* member. */
};

/*
============================================================
====================== CoherentElastic =====================
============================================================
*/
class CoherentElastic : public Base {
    
    private:
        S_table m_S_table;                                                  /**< The cumulative scattering factor \f$S(T,E)\f$. */

    public:
        CoherentElastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~CoherentElastic( );

        S_table const &s_table( ) const { return( m_S_table ); }            /**< Returns the value of the *m_S_table* member. */
};

/*
============================================================
=================== DebyeWallerIntegral ====================
============================================================
*/
class DebyeWallerIntegral : public Form {


    private:
        Functions::Function1dForm *m_function1d;                        /**< The 1-d function representing the Debye-Waller integral function \f$W(T)\f$. */

    public:
        DebyeWallerIntegral( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~DebyeWallerIntegral( );

        Functions::Function1dForm       *function1d( )       { return( m_function1d ); }       /**< Returns the value of the *m_function1d* member. */
        Functions::Function1dForm const *function1d( ) const { return( m_function1d ); }       /**< Returns the value of the *m_function1d* member. */
};

/*
============================================================
====================== IncoherentElastic =====================
============================================================
*/
class IncoherentElastic : public Base {

    private:
        PhysicalQuantity m_boundAtomCrossSection;                       /**< The characteristic bound cross section. */
        DebyeWallerIntegral m_DebyeWallerIntegral;                      /**< The Debye-Waller integral function \f$W(T)\f$. */

    public:
        IncoherentElastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~IncoherentElastic( );

        PhysicalQuantity const &boundAtomCrossSection( ) { return( m_boundAtomCrossSection ); }   /**< Returns the value of the *m_boundAtomCrossSection* member. */
        DebyeWallerIntegral const &debyeWallerIntegral( ) const { return( m_DebyeWallerIntegral ); }  /**< Returns the value of the *m_DebyeWallerIntegral* member. */
};

/*
============================================================
========================== Options =========================
============================================================
*/

class Options : public Form {

    private:
        bool m_calculatedAtThermal;                                     /**< If *true* calculate at 0.0253 eV/k. */
        bool m_asymmetric;                                              /**< If *true* S(alpha,beta) is asymmetric, otherwise it is symmetric. */

    public:
        Options( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Options( );

        bool calculatedAtThermal( ) { return( m_calculatedAtThermal ); }    /**< Returns the value of the *m_calculatedAtThermal* member. */
        bool asymmetric( ) { return( m_asymmetric ); }                      /**< Returns the value of the *m_asymmetric* member. */
};

/*
============================================================
======================== T_effective =======================
============================================================
*/
class T_effective : public Form {

    private:
        Functions::Function1dForm *m_function1d;                               /**< The 1-d function representing effective temperature \f$T_{\rm eff}(T)\f$. */

    public:
        T_effective( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~T_effective( );

        Functions::Function1dForm const *function1d( ) const { return( m_function1d ); }   /**< Returns the value of the *m_function1d* member. */
};

/*
============================================================
====================== ScatteringAtom ======================
============================================================
*/
class ScatteringAtom : public Form {

    private:
        PhysicalQuantity m_mass;                                /**< The mass of the atom. */
        PhysicalQuantity m_freeAtomCrossSection;                /**< The free atom scattering cross section. */
        PhysicalQuantity m_e_critical;                          /**< The energy value above which the static model of elastic scattering is adequate. */
        PhysicalQuantity m_e_max;                               /**< The upper energy limit for the constant. */
        T_effective m_T_effective;                              /**< The effective temperatures for the shortcollision-time approximation given as a function of moderator temperature for the atom. */

    public:
        ScatteringAtom( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~ScatteringAtom( );

        PhysicalQuantity const &mass( ) const { return( m_mass ); }                                     /**< Returns the value of the *m_mass* member. */
        PhysicalQuantity const &freeAtomCrossSection( ) const { return( m_freeAtomCrossSection ); }     /**< Returns the value of the *m_freeAtomCrossSection* member. */
        PhysicalQuantity const &e_critical( ) const { return( m_e_critical ); }                         /**< Returns the value of the *m_e_critical* member. */
        PhysicalQuantity const &e_max( ) const { return( m_e_max ); }                                   /**< Returns the value of the *m_e_max* member. */
        T_effective const &t_effective( ) const { return( m_T_effective ); }                            /**< Returns the value of the *m_T_effective* member. */
};

/*
============================================================
======================= S_alpha_beta =======================
============================================================
*/
class S_alpha_beta : public Form {

    private:
        Functions::Function3dForm *m_function3d;                           /**< The \f$S(T,\alpha,\beta)\f$ function. */

    public:
        S_alpha_beta( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~S_alpha_beta( );

        Functions::Function3dForm *function3d( ) { return( m_function3d ); }                   /**< Returns the value of the *m_function3d* member. */
};

}               // End namespace n_ThermalNeutronScatteringLaw.

}               // End namespace DoubleDifferentialCrossSection.

namespace Distributions {

/*
============================================================
========================= Distribution =====================
============================================================
*/
class Distribution : public Form {

    private:
        Frame m_productFrame;                                                   /**< The product frame for the distribution form. */

    public:
        Distribution( std::string const &a_moniker, FormType a_type, std::string const &a_label, Frame a_productFrame );
        Distribution( HAPI::Node const &a_node, SetupInfo &a_setupInfo, FormType a_type, Suite *a_parent );

        Frame productFrame( ) const { return( m_productFrame ); }               /**< Returns the value of the *m_productFrame* member. */
        void toXMLNodeStarter( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
======================= AngularTwoBody =====================
============================================================
*/
class AngularTwoBody : public Distribution {

    private:
        Functions::Function2dForm *m_angular;                                              /**< The P(mu|E) distribution as a Function2dForm. */

    public:
        AngularTwoBody( std::string const &a_label, Frame a_productFrame, Functions::Function2dForm *a_angular = nullptr );
        AngularTwoBody( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~AngularTwoBody( );

        Functions::Function2dForm const *angular( ) const { return( m_angular ); }         /**< Returns the value of the *m_angular* member as a const pointer. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
========================= KalbachMann ======================
============================================================
*/
class KalbachMann : public Distribution {

    private:
        Functions::Function2dForm *m_f;                                                    /**< The P(E'|E) distribution as a Function2dForm. */
        Functions::Function2dForm *m_r;                                                    /**< The Kalbach/Mann r(E,E') function as a Function2dForm. */
        Functions::Function2dForm *m_a;                                                    /**< The Kalbach/Mann a(E,E') function as a Function2dForm. */

    public:
        KalbachMann( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~KalbachMann( );

        Functions::Function2dForm const *f( ) const { return( m_f ); }                     /**< Returns the value of the *m_f* member. */
        Functions::Function2dForm const *r( ) const { return( m_r ); }                     /**< Returns the value of the *m_r* member. */
        Functions::Function2dForm const *a( ) const { return( m_a ); }                     /**< Returns the value of the *m_a* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================== EnergyAngular =====================
============================================================
*/
class EnergyAngular : public Distribution {

    private:
        Functions::Function3dForm *m_energyAngular;                                                /**< The P(E',mu|E) distribution as a Function3dForm. */

    public:
        EnergyAngular( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~EnergyAngular( );

        Functions::Function3dForm const *energyAngular( ) const { return( m_energyAngular ); }     /**< Returns the value of the *m_energyAngular* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= EnergyAngularMC ====================
============================================================
*/
class EnergyAngularMC : public Distribution {

    private:
        Functions::Function2dForm *m_energy;                                               /**< The P(E'|E) distribution as a Function2dForm. */
        Functions::Function3dForm *m_energyAngular;                                        /**< The P(mu|E,E') distribution as a Function3dForm. */

    public:
        EnergyAngularMC( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~EnergyAngularMC( );

        Functions::Function2dForm const *energy( ) const { return( m_energy ); }                   /**< Returns the value of the *m_energy* member. */
        Functions::Function3dForm const *energyAngular( ) const { return( m_energyAngular ); }     /**< Returns the value of the *m_energyAngular* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================== AngularEnergy =====================
============================================================
*/
class AngularEnergy : public Distribution {

    private:
        Functions::Function3dForm *m_angularEnergy;                                                /**< The P(mu,E'|E) distribution as a Function3dForm. */

    public:
        AngularEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~AngularEnergy( );

        Functions::Function3dForm const *angularEnergy( ) const { return( m_angularEnergy ); }     /**< Returns the value of the *m_angularEnergy* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= AngularEnergyMC ====================
============================================================
*/
class AngularEnergyMC : public Distribution {

    private:
        Functions::Function2dForm *m_angular;                                                      /**< The P(mu|E) distribution as a Function2dForm. */
        Functions::Function3dForm *m_angularEnergy;                                                /**< The P(E'|E,mu) distribution as a Function3dForm. */

    public:
        AngularEnergyMC( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~AngularEnergyMC( );

        Functions::Function2dForm const *angular( ) const { return( m_angular ); }                 /**< Returns the value of the *m_angular* member. */
        Functions::Function3dForm const *angularEnergy( ) const { return( m_angularEnergy ); }     /**< Returns the value of the *m_angularEnergy* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================= Uncorrelated =====================
============================================================
*/
class Uncorrelated : public Distribution {

    private:
        Functions::Function2dForm *m_angular;                                              /**< The P(mu|E) distribution as a Function2dForm. */
        Functions::Function2dForm *m_energy;                                               /**< The P(E'|E) distribution as a Function2dForm. */

    public:
        Uncorrelated( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Uncorrelated( );

        Functions::Function2dForm const *angular( ) const { return( m_angular ); }         /**< Returns the value of the *m_angular* member. */
        Functions::Function2dForm const *energy( ) const { return( m_energy ); }           /**< Returns the value of the *m_energy* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================= MultiGroup3d =====================
============================================================
*/
class MultiGroup3d : public Distribution {

    private:
        Functions::Gridded3d m_gridded3d;                                              /**< The multi-group Legendre distribution as a Gridded3d instance. */

    public:
        MultiGroup3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        Functions::Gridded3d const &data( ) const { return( m_gridded3d ); }           /**< Returns the value of the *m_gridded3d* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
====================== LLNLAngularEnergy ===================
============================================================
*/
class LLNLAngularEnergy : public Distribution {

    private:
        Functions::Function2dForm *m_angular;                                          /**< The P(mu|E) distribution as a Function2dForm. */
        Functions::Function3dForm *m_angularEnergy;                                    /**< The P(E'|E,mu) distribution as a Function3dForm. */

    public:
        LLNLAngularEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~LLNLAngularEnergy( );

        Functions::Function2dForm const *angular( ) const { return( m_angular ); }                 /**< Returns the value of the *m_angular* member. */
        Functions::Function3dForm const *angularEnergy( ) const { return( m_angularEnergy ); }     /**< Returns the value of the *m_angularEnergy* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
============== CoherentPhotoAtomicScattering ===============
============================================================
*/
class CoherentPhotoAtomicScattering : public Distribution {

    private:
        std::string m_href;                                                 /**< xlink to the IncoherentPhotoAtomicScattering instance under the *m_doubleDifferentialCrossSection* node. */

    public:
        CoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                          /**< Returns the value of the *m_href* member. */
};

/*
============================================================
============== IncoherentPhotoAtomicScattering =============
============================================================
*/
class IncoherentPhotoAtomicScattering : public Distribution {

    private:
        std::string m_href;                                                 /**< xlink to the IncoherentPhotoAtomicScattering instance under the *m_doubleDifferentialCrossSection* node. */

    public:
        IncoherentPhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                          /**< Returns the value of the *m_href* member. */
};

/*
============================================================
======== IncoherentBoundToFreePhotoAtomicScattering ========
============================================================
*/
class IncoherentBoundToFreePhotoAtomicScattering : public Distribution {

    private:
        std::string m_href;                                                 /**< xlink to the IncoherentPhotoAtomicScattering instance under the *m_doubleDifferentialCrossSection* node. */

    public:
        IncoherentBoundToFreePhotoAtomicScattering( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                          /**< Returns the value of the *m_href* member. */
};

/*
============================================================
=============== ThermalNeutronScatteringLaw ================
============================================================
*/
class ThermalNeutronScatteringLaw : public Distribution {

    private:
        std::string m_href;                                                 /**< xlink to the IncoherentPhotoAtomicScattering instance under the *m_doubleDifferentialCrossSection* node. */

    public:
        ThermalNeutronScatteringLaw( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                          /**< Returns the value of the *m_href* member. */
};

/*
============================================================
======================== Branching3d =======================
============================================================
*/
class Branching3d : public Distribution {

    private:
        std::string m_initialState;                                         /**< The nuclide level that decays, emitting a photon. */

    public:
        Branching3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &initialState( ) const { return( m_initialState ); }        /**< Returns the value of the *m_initialState* member. */
};

/*
============================================================
======================= Reference3d ========================
============================================================
*/
class Reference3d : public Distribution {

    private:
        std::string m_href;                                                     /**< Link to the other function. */

    public:
        Reference3d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                  /**< Returns the value of the *m_xlink* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
================ CoulombPlusNuclearElastic =================
============================================================
*/

class CoulombPlusNuclearElastic : public Distribution {

    private:
        std::string m_href;                                                     /**< Link to the other function. */

    public:
        CoulombPlusNuclearElastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        std::string const &href( ) const { return( m_href ); }                  /**< Returns the value of the *m_xlink* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= LLNLLegendre =======================
============================================================
*/

class LLNLLegendre : public Distribution {
//
// This class is woefully inadequate but some form is needed by the method Product::isCompleteParticle.
//

    public:
        LLNLLegendre( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================= Unspecified ======================
============================================================
*/
class Unspecified : public Distribution {

    public:
        Unspecified( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

}                     // End of namespace Distributions.

/*
============================================================
=========================== Suite ==========================
============================================================
*/
class Suite : public GUPI::Ancestry {

    public:
        typedef std::vector<Form *> Forms;                              /**< The typedef the the *m_forms* member. */

    private:
        std::string m_keyName;                                          /**< The name of the key used to look up items in the suite. */
        mutable Forms m_forms;                                          /**< The list of nodes stored within *this*. */
        std::map<std::string,int> m_map;                                /**< A map of *this* node labels to their index in *m_forms*. */
        Styles::Suite const *m_styles;                                  /**< The Styles::Suite for the Protare that *this* resides in. */
        bool m_allowsLazyParsing;                                       /**< If **true**, the suite allows its elements to be lazy parsed. */
        std::string m_href;                                             /**< xlink to the to a Suite that has the elements for this Suite. */

            // FIXME should we make public or private copy constructor?

    public:
        Suite( std::string const &a_keyName = GIDI_labelChars );
        Suite( std::string const &a_moniker, std::string const &a_keyName );
        Suite( Construction::Settings const &a_construction, std::string const &a_moniker, std::string const &a_keyName, HAPI::Node const &a_node, 
                        SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, parseSuite a_parseSuite, 
                        Styles::Suite const *a_styles, bool a_allowsLazyParsing = false );
        ~Suite( );

        std::string const &keyName( ) const { return( m_keyName ); }                        /**< Returns a *const* reference to the *m_keyName* member. */
        std::size_t size( ) const { return( m_forms.size( ) ); }                            /**< Returns the number of node contained by *this*. */
        typedef Forms::iterator iterator;
        typedef Forms::const_iterator const_iterator;
        iterator begin( ) { return m_forms.begin( ); }                                      /**< The C++ *begin iterator* for *this*. */
        const_iterator begin( ) const { return m_forms.begin( ); }                          /**< The C++ const *begin iterator* for *this*. */
        iterator end( ) { return m_forms.end( ); }                                          /**< The C++ *end iterator* for *this*. */
        const_iterator end( ) const { return m_forms.end( ); }                              /**< The C++ const *end iterator* for *this*. */
        int operator[]( std::string const &a_label ) const ;
        template<typename T> T       *get( std::size_t a_Index );
        template<typename T> T const *get( std::size_t a_Index ) const ;
        template<typename T> T       *get( std::string const &a_label );
        template<typename T> T const *get( std::string const &a_label ) const ;
        template<typename T> T *getViaLineage( std::string const &a_label );
        template<typename T> T       *pop( std::size_t a_Index );
        template<typename T> T       *pop( std::string const &a_label );

        Styles::Suite const *styles( ) { return( m_styles ); }                              /**< Returns the value of the *m_styles* member. */
        std::string const &href( ) const { return( m_href ); }                              /**< Returns a reference to the *m_ref* member. */

        void parse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, 
                        parseSuite a_parseSuite, Styles::Suite const *a_styles );
        void add( Form *a_form );
        iterator find( std::string const &a_label, bool a_convertLazyParsingHelperForm = false );
        const_iterator find( std::string const &a_label, bool a_convertLazyParsingHelperForm = false ) const ;
        bool has( std::string const &a_label ) const { return( find( a_label ) != m_forms.end( ) ); }
        Form *checkLazyParsingHelperForm( std::size_t a_index );
        Form *checkLazyParsingHelperForm( std::size_t a_index ) const ;
        iterator checkLazyParsingHelperFormIterator( iterator a_iter ) ;
        const_iterator checkLazyParsingHelperFormIterator( const_iterator a_iter ) const ;

        void modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex );
        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;
        std::vector<iterator> findAllOfMoniker( std::string const &a_moniker ) ;
        std::vector<const_iterator> findAllOfMoniker( std::string const &a_moniker ) const ;
        Form const *findInstanceOfTypeInLineage( std::string const &_label, std::string const &a_moniker ) const ;
        Form       *findInstanceOfTypeInLineage( Styles::Suite const &a_styles, std::string const &_label, std::string const &a_moniker ) ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
        void printFormLabels( std::string const &a_header ) const ;
};

/* *********************************************************************************************************//**
 * Returns the node at index *a_index*.
 *
 * @param a_index               [in]    The index of the node to return.
 *
 * @return                              The node at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T *Suite::get( std::size_t a_index ) {

    Form *__form = checkLazyParsingHelperForm( a_index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::get( std::size_t ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node at index *a_index*.
 *
 * @param a_index               [in]    The index of the node to return.
 *
 * @return                              The node at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T const *Suite::get( std::size_t a_index ) const {

    Form *__form = checkLazyParsingHelperForm( a_index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::get( std::size_t ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node with label *a_label*.
 *
 * @param a_label               [in]    The label of the node to return.
 *
 * @return                              The node with label *a_label*.
 ***********************************************************************************************************/

template<typename T> T *Suite::get( std::string const &a_label ) {

    int index = (*this)[a_label];
    Form *__form = checkLazyParsingHelperForm( index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::get( std::string const & ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Returns the node with label *a_label*.
 *
 * @param a_label               [in]    The label of the node to return.
 *
 * @return                              The node with label *a_label*.
 ***********************************************************************************************************/

template<typename T> T const *Suite::get( std::string const &a_label ) const {

    int index = (*this)[a_label];
    Form *__form = checkLazyParsingHelperForm( index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::get( std::string const & ): invalid cast" );

    return( object );
}

/* *********************************************************************************************************//**
 * Removes the form at index *a_index* and returns it. It is up to the calling function to delete the form,
 * otherwise there will be memory leak.
 *
 * @param a_index               [in]    The index of the node to return.
 *
 * @return                              The node at index *a_index*.
 ***********************************************************************************************************/

template<typename T> T *Suite::pop( std::size_t a_index ) {

    Form *__form = checkLazyParsingHelperForm( a_index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::pop( std::size_t ): invalid cast" );

    for( std::size_t index = a_index + 1; index < m_forms.size( ); ++index ) {
        m_forms[index-1] = m_forms[index];
        m_map[m_forms[index-1]->label( )] = index - 1;
    }
    m_forms.resize( m_forms.size( ) - 1 );

    return( object );
}

/* *********************************************************************************************************//**
 * Removes the form with label *a_label* and returns it. It is up to the calling function to delete the form,
 * otherwise there will be memory leak.
 *
 * @param a_label               [in]    The label of the node to return.
 *
 * @return                              The node at index *a_label*.
 ***********************************************************************************************************/

template<typename T> T *Suite::pop( std::string const &a_label ) {

    int index = (*this)[a_label];                           // This will throw an exception if *a_label* is not in *this*.
    Form *__form = checkLazyParsingHelperForm( index );
    T *object = dynamic_cast<T *>( __form );

    if( object == nullptr ) throw Exception( "GIDI::Suite::pop( std::size_t ): invalid cast" );

    for( std::size_t index2 = index + 1; index2 < m_forms.size( ); ++index2 ) {
        m_forms[index2-1] = m_forms[index2];
        m_map[m_forms[index2-1]->label( )] = index2 - 1;
    }
    m_forms.resize( m_forms.size( ) - 1 );

    return( object );

}

/*
============================================================
======================== Component =========================
============================================================
*/
class Component : public Suite {

    public:
        Component( Construction::Settings const &a_construction, std::string const &a_moniker, std::string const &a_keyName, 
                        HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, 
                        parseSuite a_parseSuite, Styles::Suite const *a_styles );
        Component( std::string const &a_moniker, std::string const &a_keyName = GIDI_labelChars );
};

namespace Table {

/*
============================================================
========================== Column ==========================
============================================================
*/

class Column : public Form {

    private:
        std::string m_index;                                                    /**< The index of the column. */
        std::string m_name;                                                     /**< The name of the column. */
        std::string m_unit;                                                     /**< The unit of the data in the column. */
        std::string m_types;                                                    /**< The types of the data in the column. */

    public:
        Column( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~Column( );

        std::string const &index( ) const { return( m_index ); }                /**< Returns a *const* reference of the *m_index* member. */
        std::string const &name( ) const { return( m_name ); }                  /**< Returns a *const* reference of the *m_name* member. */
        std::string const &unit( ) const { return( m_unit ); }                  /**< Returns a *const* reference of the *m_unit* member. */
        std::string const &types( ) const { return( m_types ); }                /**< Returns a *const* reference of the *m_types* member. */

        void setKeyValue( std::string const &a_keyName ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Data ===========================
============================================================
*/

class Data : public GUPI::Ancestry {

    private:
        std::string m_sep;
        std::string m_body;

    public:
        Data( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Data( );

        std::string const &sep( ) const { return( m_sep ); }                      /**< Returns a *const* reference of the *m_sep* member. */
        std::string const &body( ) const { return( m_body ); }                    /**< Returns a *const* reference of the *m_body* member. */

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Table ==========================
============================================================
*/

class Table : public Form {

    private:
        int m_rows;                                     /**< The number of rows in the table. */
        int m_coluns;                                   /**< The number of columns in the table. */
        std::string m_storageOrder;                     /**< The storageOrder for the data in the table. */
        Suite m_columnHeaders;                          /**< The column header for the table. */
        Data m_data;                                    /**< The data for the table. */

    public:
        Table( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Table( );

        int rows( ) const { return( m_rows ); }                                     /**< Returns the value of the *m_rows* member. */
        int columns( ) const { return( m_coluns ); }                                /**< Returns the value of the *m_coluns* member. */
        std::string const &storageOrder( ) const { return( m_storageOrder ); }      /**< Returns the value of the *m_storageOrder* member. */
        Suite const &columnHeaders( ) const { return( m_columnHeaders ); }          /**< Returns the value of the *m_columnHeaders* member. */
        Data const &data( ) const { return( m_data ); }                             /**< Returns the value of the *m_data* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

}               // End Table namespace.

/*
============================================================
=========================== Flux ===========================
============================================================
*/
class Flux : public Form {

    private:
        Functions::Function2dForm *m_flux;                                          /**< The flux f(E,mu). */

    public:
        Flux( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Flux( );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================== Group ===========================
============================================================
*/
class Group : public Form {

    private:
        Grid m_grid;                                                /**< Multi-group boundaries for this Group. */

    public:
        Group( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops );
        Group( Group const &a_group );

        std::size_t size( ) const { return( m_grid.size( ) ); }                             /**< Returns the number of multi-group boundaries. */
        inline double &operator[]( std::size_t a_index ) { return( m_grid[a_index] ); }     /**< Returns the multi-group boundary at index *a_index*. */
        std::vector<double> data( ) const { return( m_grid.data().vector() ); }              /**< Returns the multi-group boundaries. */
        Grid const &grid( ) const { return( m_grid ); }                                     /**< Returns the value of the *m_grid* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
====================== Transportable =======================
============================================================
*/
class Transportable : public Form {

    private:
        std::string m_conserve;                                     /**< Conservation flag for the transfer matrices for this particle. Currently, only "*number*" is allowed. */
        Group m_group;                                              /**< Multi-group boundaries for this Transportable. */

    public:
        Transportable( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, Suite *a_parent );
        Transportable( Transportable const &a_transportable );

        std::string pid( ) const { return( label( ) ); }                                    /**< Returns the value of the particle id for the *Transportable*. */
        std::string const &conserve( ) const { return( m_conserve ); }                      /**< Returns a const reference to member *m_conserve*. */
        Group const &group( ) const { return( m_group ); }                                  /**< Returns the value of the *m_group* member. */
        std::vector<double> groupBoundaries( ) const { return( m_group.data( ) ); }         /**< Returns the multi-group boundaries for this transportable particle. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= ExternalFile =======================
============================================================
*/

class ExternalFile : public Form {

    private:
        std::string m_path;                         /**< The path to the external file. */

    public:
        ExternalFile( std::string const &a_label, std::string const &a_path );
        ExternalFile( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~ExternalFile( );

        std::string const &path( ) const { return( m_path ); }

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
==================== Documentation_1_10 ====================
============================================================
*/

namespace Documentation_1_10 {

class Documentation : public Form {

    private:
        std::string m_label;                /**< The label for the documentation. */
        std::string m_text;                 /**< The documentation text. */

    public:
        Documentation( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~Documentation( ) { }

        std::string const &label( ) const { return m_label; }       /**< Returns a const reference to the *m_label* member. */
        std::string const &text( ) const { return m_text; }         /**< Returns a const reference to the *m_text* member. */

};

/*
============================================================
=========================== Suite ==========================
============================================================
*/
class Suite : public GIDI::Suite {

    public:
        Suite( );

        void parse( HAPI::Node const &a_node, SetupInfo &a_setupInfo );
};

}                     // End of namespace Documentation_1_10.

/*
============================================================
===================== ExternalFiles stuff ==================
============================================================
*/

namespace ExternalFiles {

/*
============================================================
========================== Suite ===========================
============================================================
*/
class Suite : public GIDI::Suite {

    public:
        void registerBinaryFiles(std::string a_parentDir, SetupInfo &a_setupInfo);

};

}                     // End of namespace ExternalFiles.

/*
============================================================
========================= Styles stuff =====================
============================================================
*/

namespace Styles {

/*
============================================================
========================== Base ============================
============================================================
*/
class Base : public Form {

    private:
        std::string m_date;                     /**< The GNDS <**date**> attribute. */
        std::string m_label;                    /**< The GNDS <**label**> attribute. */
        std::string m_derivedStyle;             /**< The GNDS <**derivedFrom**> attribute. */
        GUPI::Documentation *m_documentation;

    public:
        Base( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~Base( );

        std::string const &date( ) const { return( m_date ); }                      /**< Returns the value of the *m_date* member. */
        std::string const &label( ) const { return( m_label ); }                    /**< Returns the value of the *m_label* member. */
        std::string const &derivedStyle( ) const { return( m_derivedStyle ); }      /**< Returns the value of the *m_derivedStyle* member. */
        bool hasDocumentation( ) { return ( m_documentation != nullptr ); }
        GUPI::Documentation *documentation( ) { return ( m_documentation ); }
        virtual PhysicalQuantity const &temperature( ) const = 0;
        Base const *getDerivedStyle( ) const ;
        Base const *getDerivedStyle( std::string const &a_moniker ) const ;

        std::string baseXMLAttributes( GUPI::WriteInfo &a_writeInfo ) const ;
};

/*
============================================================
======================== Evaluated =========================
============================================================
*/
class Evaluated : public Base {

    private:
        std::string m_library;                      /**< The GNDS <**library**> attribute. */
        std::string m_version;                      /**< The GNDS <**version**> attribute. */
        PhysicalQuantity m_temperature;             /**< The GNDS <**temperature**> node data. */
        AxisDomain m_projectileEnergyDomain;        /**< The GNDS <**projectileEnergyDomain**> node data. */

    public:
        Evaluated( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );

        PhysicalQuantity const &temperature( ) const { return( m_temperature ); }   /**< Returns the value of the *m_temperature* member. */
        AxisDomain const &projectileEnergyDomain( ) const { return( m_projectileEnergyDomain ); }
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
================ CrossSectionReconstructed =================
============================================================
*/
class CrossSectionReconstructed : public Base {

    private:
        PhysicalQuantity *m_temperature;            /**< The GNDS <**temperature**> node data. */

    public:
        CrossSectionReconstructed( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~CrossSectionReconstructed( );

        PhysicalQuantity const &temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
============= AngularDistributionReconstructed =============
============================================================
*/
class AngularDistributionReconstructed : public Base {

    private:
        PhysicalQuantity *m_temperature;            /**< The GNDS <**temperature**> node data. */

    public:
        AngularDistributionReconstructed( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~AngularDistributionReconstructed( );

        PhysicalQuantity const &temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
============= CoulombPlusNuclearElasticMuCutoff ============
============================================================
*/
class CoulombPlusNuclearElasticMuCutoff : public Base {

    private:
        double m_muCutoff;                      /**< The GNDS <**muCutoff**> attribute. */

    public:
        CoulombPlusNuclearElasticMuCutoff( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );

        PhysicalQuantity const &temperature( ) const ;
        double muCutoff( ) const { return( m_muCutoff ); }          /**< Returns the value of the *m_muCutoff* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================= Realization ========================
============================================================
*/
class Realization : public Base {

    public:
        Realization( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );

        PhysicalQuantity const & temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
=================== AverageProductData =====================
============================================================
*/
class AverageProductData : public Base {

    private:
        PhysicalQuantity *m_temperature;            /**< The GNDS <**temperature**> node data. */

    public:
        AverageProductData( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        ~AverageProductData( );

        PhysicalQuantity const &temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
===================== MonteCarlo_cdf =======================
============================================================
*/
class MonteCarlo_cdf : public Base {

    public:
        MonteCarlo_cdf( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );

        PhysicalQuantity const &temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
======================== MultiGroup ========================
============================================================
*/
class MultiGroup : public Base {

    private:
        int m_maximumLegendreOrder;         /**< The GNDS <**lMax**> attribute. */
        GIDI::Suite m_transportables;       /**< The GNDS <**transportables**> node. */

    public:
        MultiGroup( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, GIDI::Suite *a_parent );
        ~MultiGroup( );

        int maximumLegendreOrder( ) const { return( m_maximumLegendreOrder ); }     /**< Returns the value of the *m_maximumLegendreOrder* member. */
        PhysicalQuantity const &temperature( ) const ;

        std::vector<double> groupBoundaries( std::string const &a_productID ) const ;
        GIDI::Suite const &transportables( ) const { return( m_transportables ); }  /**< Returns the value of the *m_transportables* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================= Heated ===========================
============================================================
*/
class Heated : public Base {

    private:
        PhysicalQuantity m_temperature;                                 /**< The GNDS <**temperature**> node data. */

    public:
        Heated( HAPI::Node const &a_node, SetupInfo &a_setupInfo, GIDI::Suite *a_parent );
        PhysicalQuantity const & temperature( ) const { return( m_temperature ); }  /**< Returns the value of the *m_temperature* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
===================== HeatedMultiGroup =====================
============================================================
*/
class HeatedMultiGroup : public Base {

    private:
        GIDI::Suite m_transportables;                   /**< The GNDS <**transportables**> node. For GNDS 2.0 and above. */
        Flux m_flux;                                    /**< The GNDS <**flux**> node. */
        Functions::Gridded1d m_inverseSpeed;            /**< The GNDS <**inverseSpeed**> node data. */
        std::string m_parameters;                       /**< The GNDS <**parameters**> attribute. Only used for GNDS 1.10. */

    public:
        HeatedMultiGroup( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, GIDI::Suite *a_parent );
        ~HeatedMultiGroup( );

        PhysicalQuantity const &temperature( ) const ;

        GIDI::Suite const &transportables( ) const { return( m_transportables ); }   /**< Returns a const reference to *m_transportables*. */
        Transportable const &transportable( std::string const &a_ID ) const ;
        std::vector<double> groupBoundaries( std::string const &a_ID ) const ;
        Flux const &flux( ) const { return( m_flux ); }                         /**< Returns a const reference to member *m_flux*. */
        std::string const &parameters( ) const { return( m_parameters ); }      /**< Returns a const reference to member *m_parameters*. Only used for GNDS 1.10. */

        Vector inverseSpeedData( ) const { return( m_inverseSpeed.data( ) ); }      /**< Returns the value of the *m_inverseSpeed* data. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
==================== SnElasticUpScatter ====================
============================================================
*/
class SnElasticUpScatter : public Base {

    private:
        int m_upperCalculatedGroup;             /**< The GNDS <**upperCalculatedGroup**> attribute. */

    public:
        SnElasticUpScatter( HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, GIDI::Suite *a_parent );
        ~SnElasticUpScatter( );

        PhysicalQuantity const &temperature( ) const ;
        int upperCalculatedGroup( ) const { return( m_upperCalculatedGroup ); }     /**< Returns the value of the *m_upperCalculatedGroup* data. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
==================== GriddedCrossSection ===================
============================================================
*/
class GriddedCrossSection : public Base {

    private:
        Grid m_grid;                        /**< The GNDS <**grid**> node. */

    public:
        GriddedCrossSection( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, GIDI::Suite *a_parent );
        ~GriddedCrossSection( );

        PhysicalQuantity const &temperature( ) const ;
        Grid const &grid( ) const { return( m_grid ); }     /**< Returns the value of the *m_grid*. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
=================== URR_probabilityTables ==================
============================================================
*/
class URR_probabilityTables : public Base {

    public:
        URR_probabilityTables( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, GIDI::Suite *a_parent );
        ~URR_probabilityTables( );

        PhysicalQuantity const &temperature( ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
========================== Suite ===========================
============================================================
*/
class Suite : public GIDI::Suite {

    public:
        Suite( );

        std::string const *findLabelInLineage( GIDI::Suite const &a_suite, std::string const &a_label ) const ;
};

/*
============================================================
===================== TemperatureInfo ======================
============================================================
*/

class TemperatureInfo {

    private:
        PhysicalQuantity m_temperature;                 /**< The temperature for this TemperatureInfo. */
        std::string m_heatedCrossSection;               /**< The label for the *heatedCrossSection* data for this temperature. */
        std::string m_griddedCrossSection;              /**< The label for the *griddedCrossSection* data for this temperature. */
        std::string m_URR_probabilityTables;            /**< The label for the *URR_probabilityTables* data for this temperature. */
        std::string m_heatedMultiGroup;                 /**< The label for the *heatedMultiGroup* data for this temperature. */
        std::string m_SnElasticUpScatter;               /**< The label for the *SnElasticUpScatter* data for this temperature. */

    public:
        TemperatureInfo( );
        TemperatureInfo( PhysicalQuantity const &a_temperature, std::string const &a_heatedCrossSection, std::string const &a_griddedCrossSection,
                std::string const &a_URR_probabilityTables, std::string const &a_heatedMultiGroup, std::string const &a_SnElasticUpScatter );

        PhysicalQuantity const &temperature( ) const { return( m_temperature ); }                   /**< Returns the value of the *m_temperature*. */
        std::string const &heatedCrossSection( ) const { return( m_heatedCrossSection ); }          /**< Returns the value of the *m_heatedCrossSection*. */
        std::string const &griddedCrossSection( ) const { return( m_griddedCrossSection ); }        /**< Returns the value of the *m_griddedCrossSection*. */
        std::string const &URR_probabilityTables( ) const { return( m_URR_probabilityTables ); }    /**< Returns the value of the *m_URR_probabilityTables*. */
        std::string const &heatedMultiGroup( ) const { return( m_heatedMultiGroup ); }              /**< Returns the value of the *m_heatedMultiGroup*. */
        std::string const &SnElasticUpScatter( ) const { return( m_SnElasticUpScatter ); }          /**< Returns the value of the *m_SnElasticUpScatter*. */

        void print( ) const ;
};

typedef std::vector<Styles::TemperatureInfo> TemperatureInfos;

}               // End of namespace Styles.

/*
=========================================================
*/
template<typename T> T *Suite::getViaLineage( std::string const &a_label ) {

    std::string const *label = m_styles->findLabelInLineage( (Styles::Suite &) *this, a_label );

    return( get<T>( *label ) );
}

/*
============================================================
==================== Transporting stuff ====================
============================================================
*/

namespace Transporting {

class ProcessedFlux;

enum class Mode { multiGroup, multiGroupWithSnElasticUpScatter, MonteCarloContinuousEnergy };
enum class DelayedNeutrons { off, on };
enum class Conserve { number, energyOut };

/*
============================================================
========================== MultiGroup ======================
============================================================
*/
class MultiGroup {

    private:
        std::string m_label;                                        /**< The label for the multi-group. */
        std::vector<double> m_boundaries;                           /**< The list of boundaries for the multi-group. */

    public:
        MultiGroup( );
        MultiGroup( std::string const &a_label, int a_length, double const *a_values );
        MultiGroup( std::string const &a_label, std::vector<double> const &a_boundaries );
        MultiGroup( Group const &a_group );
        MultiGroup( MultiGroup const &a_multiGroup );
        ~MultiGroup( );
        MultiGroup &operator=( MultiGroup const &a_rhs );

        double operator[]( int const a_index ) const { return( m_boundaries[a_index] ); }           /**< Returns the multi-group boundary at index *a_index*. */
        std::size_t size( ) const { return( m_boundaries.size( ) ); }                               /**< Returns the number of multi-group boundaries. */
        int numberOfGroups( ) const { return( (int) ( m_boundaries.size( ) - 1 ) ); }               /**< Returns the number of multi-group groups. */
        std::vector<double> const &boundaries( ) const { return( m_boundaries ); }                  /**< Returns the value of the *m_boundaries* member. */
        double const *pointer( ) const { return( &(m_boundaries[0]) ); }                            /**< Returns a pointer to the beginning of the multi-group boundaries. */

        void set( std::string const &a_label, std::vector<double> const &a_boundaries );
        std::string const &label( ) const { return( m_label ); }                                    /**< Returns the value of the *m_label* member. */
        int multiGroupIndexFromEnergy( double a_energy, bool a_encloseOutOfRange ) const ;
        void print( std::string const &a_indent, bool a_outline = false, int a_valuesPerLine = 10 ) const ;
};

/*
============================================================
==================== Groups_from_bdfls =====================
============================================================
*/
class Groups_from_bdfls {

    private:
        std::vector<MultiGroup> m_multiGroups;                                          /**< List of MultiGroup's read in from the bdfls file. */

    public:
        Groups_from_bdfls( std::string const &a_fileName );
        Groups_from_bdfls( char const *a_fileName );
        ~Groups_from_bdfls( );

        MultiGroup viaLabel( std::string const &a_label ) const ;
        MultiGroup getViaGID( int a_gid ) const;
        std::vector<std::string> labels( ) const;
        std::vector<int> GIDs( ) const;
        void print( bool a_outline = true, int a_valuesPerLine = 10 ) const;

    private:
        void initialize( char const *a_fileName );
};

/*
============================================================
========================= Flux_order =======================
============================================================
*/
class Flux_order {

    private:
        int m_order;                        /**< The Legendre order of the flux. */
        std::vector<double> m_energies;     /**< List of flux energies. */
        std::vector<double> m_fluxes;       /**< List of flux values - one for each element of m_energies. */

    public:
        Flux_order( int a_order, int a_length, double const *a_energies, double const *a_fluxes );
        Flux_order( int a_order, std::vector<double> const &a_energies, std::vector<double> const &a_fluxes );
        Flux_order( Flux_order const  &a_fluxOrder  );
        ~Flux_order( );

        int order( ) const { return( m_order ); }                                   /**< Returns the value of the *m_order* member. */
        int size( ) const { return( (int) m_energies.size( ) ); }                   /**< Returns the number of energy, flux pairs. */
        double const *energies( ) const { return( &(m_energies[0]) ); }             /**< Returns a pointer to the beginning of the energy data. */
        std::vector<double> const &v_energies( ) const { return( m_energies ); }    /**< Returns the value of the *m_energies* member. */
        double const *fluxes( ) const { return( &(m_fluxes[0]) ); }                 /**< Returns a pointer to the beginning of the flux data. */
        std::vector<double> const &v_fluxes( ) const { return( m_fluxes ); }        /**< Returns the value of the *m_fluxes* member. */
        void print( int a_valuesPerLine = 10 ) const;
};

/*
============================================================
============================ Flux ==========================
============================================================
*/
class Flux {

    private:
        std::string m_label;                        /**< Label for the flux. */
        double m_temperature;                       /**< Temperature of the material that produced this flux. */
        std::vector<Flux_order> m_fluxOrders;       /**< List of fluxes for each Legendre order, *l*, sorted by Legendre order starting with *l* = 0. */

    public:
        Flux( std::string const &a_label, double a_temperature_MeV );
        Flux( char const *a_label, double a_temperature_MeV );
        Flux( Flux const &a_flux );
        ~Flux( );

        Flux_order const &operator[]( int a_order ) const { return( m_fluxOrders[a_order] ); }  /**< Returns the Flux_order for Legendre order *a_order*. */
        int maxOrder( ) const { return( (int) m_fluxOrders.size( ) - 1 ); }                     /**< Returns the maximum number of Legendre orders for *this*. */
        int size( ) const { return( (int) m_fluxOrders.size( ) ); }                             /**< Returns the number of stored Legendre orders. */

        std::string const &label( ) const { return( m_label ); }                                /**< Returns the value of the *m_label* member. */
        double temperature( ) const { return( m_temperature ); }                                /**< Returns the value of the *m_temperature* member. */
        void addFluxOrder( Flux_order const &a_fluxOrder );
        ProcessedFlux process( std::vector<double> const &a_multiGroup ) const ;
        void print( std::string const &a_indent, bool a_outline = true, int a_valuesPerLine = 10 ) const ;
};

/*
============================================================
===================== Fluxes_from_bdfls ====================
============================================================
*/
class Fluxes_from_bdfls {

    private:
        std::vector<Flux> m_fluxes;                     /**< The list of Flux read in from the *bdfls* file. */

    public:
        Fluxes_from_bdfls( std::string const &a_fileName, double a_temperature_MeV );
        Fluxes_from_bdfls( char const *a_fileName, double a_temperature_MeV );
        ~Fluxes_from_bdfls( );

        Flux getViaFID( int a_fid ) const ;
        Functions::XYs3d *get3dViaFID( int a_fid ) const ;
        std::vector<std::string> labels( ) const ;
        std::vector<int> FIDs( ) const ;
        void print( bool a_outline = true, int a_valuesPerLine = 10 ) const ;

    private:
        void initialize( char const *a_fileName, double a_temperature_MeV );
};

/*
============================================================
======================= ProcessedFlux ======================
============================================================
*/
class ProcessedFlux {

    private:
        double m_temperature;                                           /**< The temperature of the material that produced the flux. */
        std::vector<double> m_multiGroupFlux;                           /**< The Legendre order = 0 multi-grouped flux. */

    public:
        ProcessedFlux( double a_temperature, std::vector<double> const &a_multiGroupFlux );
        ProcessedFlux( ProcessedFlux const &a_processedFlux );
        ~ProcessedFlux( );

        double temperature( ) const { return( m_temperature ); }                            /**< Returns the value of the *m_temperature* member. */
        std::vector<double> const &multiGroupFlux( ) const { return( m_multiGroupFlux ); }  /**< Returns the value of the *m_multiGroupFlux* member. */
};

/*
============================================================
========================= Particle =========================
============================================================
*/
class Particle {

    private:
        std::string m_pid;                                                  /**< The PoPs id for the particle. */
        Transporting::Mode m_mode;                                          /**< Indicates the type of transport the user is likely, but not guaranteed, to do. */
        Transporting::Conserve m_conserve;                                  /**< Indicates the conservation option for this transportable. */
        MultiGroup m_multiGroup;                                            /**< Coarse multi-group to collapse to. */
        MultiGroup m_fineMultiGroup;                                        /**< Fine multi-group to collapse from. For internal use only. */
        std::vector<int> m_collapseIndices;                                 /**< Indices for collapsing to m_multiGroup. */
        std::vector<Flux> m_fluxes;                                         /**< One flux for each temperature. */
        std::vector<ProcessedFlux> m_processedFluxes;                       /**< One processed flux for each temperature. */

    public:
        Particle( std::string const &a_pid, MultiGroup const &a_multiGroup, Functions::Function3dForm const &a_fluxes, 
                        Transporting::Mode a_mode = Transporting::Mode::multiGroup );
        Particle( std::string const &a_pid, Transporting::Mode a_mode = Transporting::Mode::multiGroup );
        Particle( std::string const &a_pid, MultiGroup const &a_multiGroup, Transporting::Mode a_mode = Transporting::Mode::multiGroup );
        Particle( Particle const &a_particle );
        ~Particle( );

        std::string const &pid( ) const { return( m_pid ); }                                /**< Returns the value of the *m_pid* member. */
        Transporting::Mode mode( ) const { return( m_mode ); }                              /**< Returns the value of the *m_mode* member. */
        Transporting::Conserve conserve( ) const { return( m_conserve ); }                  /**< Returns the value of the *m_conserve* member. */
        int multiGroupIndexFromEnergy( double a_e_in, bool a_encloseOutOfRange ) const { return( m_multiGroup.multiGroupIndexFromEnergy( a_e_in, a_encloseOutOfRange ) ); }
                                                                                            /**< Returns the coarse multi-group index corresponding to energy *a_e_in*. See MultiGroup::multiGroupIndexFromEnergy. */
        int numberOfGroups( ) const { return( m_multiGroup.numberOfGroups( ) ); }           /**< Returns the number of coarse multi-group groups. */
        MultiGroup multiGroup( ) const { return( m_multiGroup ); }                          /**< Returns the value of the *m_multiGroup* member. */
        MultiGroup fineMultiGroup( ) const { return( m_fineMultiGroup ); }                  /**< Returns the value of the *m_fineMultiGroup* member. */
        int appendFlux( Flux const &a_flux );
        ProcessedFlux const *nearestProcessedFluxToTemperature( double a_temperature ) const;
        std::vector<int> const &collapseIndices( ) const { return( m_collapseIndices ); }   /**< Returns the value of the *m_collapseIndices* member. */

        void process( Transportable const &a_transportable, double a_epsilon = 1e-6 );
        void print( std::string const &a_indent ) const ;
};

/*
============================================================
======================== Particles =========================
============================================================
*/
class Particles {

    private:
        std::map<std::string, Particle> m_particles;

    public:
        Particles( );
        ~Particles( );

        std::map<std::string, Particle> &particles( ) { return( m_particles ); }                /**< Returns the value of the *m_particles* member. */
        std::map<std::string, Particle> const &particles( ) const { return( m_particles ); }    /**< Returns the value of the *m_particles* member. */
        Particle const *particle( std::string const &a_particleID ) const;
        bool add( Particle const &a_particle );
        bool remove( std::string const &a_particleID );
        void clear( ) { m_particles.clear( ); }
        bool hasParticle( std::string const &a_id ) const ;

        void process( Protare const &a_protare, std::string const &a_label );

        std::vector<std::string> sortedIDs( bool a_orderIsAscending = true ) const ;

        void print( ) const ;
};

/*
============================================================
========================= Settings =========================
============================================================
*/
class Settings {

    private:
        std::string m_projectileID;                                 /**< The PoPs id of the projectile. */
        DelayedNeutrons m_delayedNeutrons;                          /**< If **true**, include delayed neutrons when returning or setting up data. */
        bool m_nuclearPlusCoulombInterferenceOnly;                  /**< If **true**, for charge particle as projectile and elastic scattering, the Rutherford term is excluded from the elastic reaction. */
        bool m_throwOnError;                                        /**< For methods that have an argument of type *LUPI:StatusMessageReporting**, if this member is true, an error will cause a thorw; otherwise, the error will be ignored and reported to the *LUPI:StatusMessageReporting** instance. */
        bool m_zeroDepositionIfAllProductsTracked;                  /**< For a reaction, if **true* and all products are tracked, then the deposition energy will be set to zero, independent of that the data may yield. Otherwise, the data results are returned. */

    public:
        Settings( std::string const &a_projectileID, DelayedNeutrons a_delayedNeutrons );
        ~Settings( );

        std::string const &projectileID( ) const { return( m_projectileID ); }                                      /**< Returns the value of the *m_projectileID* member. */

        DelayedNeutrons delayedNeutrons( ) const { return( m_delayedNeutrons ); }                                   /**< Returns the value of the *m_delayedNeutrons* member. */
        void setDelayedNeutrons( DelayedNeutrons a_delayedNeutrons ) { m_delayedNeutrons = a_delayedNeutrons; }     /**< Sets the *m_delayedNeutrons* member to *a_delayedNeutrons*. */

        bool nuclearPlusCoulombInterferenceOnly( ) const { return( m_nuclearPlusCoulombInterferenceOnly ); }        /**< Returns the value of the *m_nuclearPlusCoulombInterferenceOnly* member. */
        void setNuclearPlusCoulombInterferenceOnly( bool a_nuclearPlusCoulombInterferenceOnly )
            { m_nuclearPlusCoulombInterferenceOnly = a_nuclearPlusCoulombInterferenceOnly; }                        /**< Sets the *m_nuclearPlusCoulombInterferenceOnly* to *a_nuclearPlusCoulombInterferenceOnly*. */

        bool zeroDepositionIfAllProductsTracked( ) const { return( m_zeroDepositionIfAllProductsTracked ); }        /**< Returns the value of the *m_zeroDepositionIfAllProductsTracked* member. */
        void setZeroDepositionIfAllProductsTracked( bool a_zeroDepositionIfAllProductsTracked ) 
            { m_zeroDepositionIfAllProductsTracked = a_zeroDepositionIfAllProductsTracked; }                        /**< Sets the *m_zeroDepositionIfAllProductsTracked* to *a_zeroDepositionIfAllProductsTracked*. */

        bool throwOnError( ) const { return( m_throwOnError ); }
        void setThrowOnError( bool a_throwOnError ) { m_throwOnError = a_throwOnError; }

        Vector multiGroupZeroVector( Particles const &a_particles, bool a_collapse = true ) const ;
        Matrix multiGroupZeroMatrix( Particles const &a_particles, std::string const &a_particleID, bool a_collapse = true ) const ;

//        void print( ) const ;
};

/*
============================================================
============================ MG ============================
============================================================
*/
class MG : public Settings {

    private:
        Mode m_mode;                                    /**< Specifies the type of data to use or retrieve for transport codes. */
        bool m_useMultiGroupSummedData;                 /**< If **true** and multi-grouped summed data available in protare, use it instead of summing data over reactions. */

    public:
        MG( std::string const &a_projectileID, Mode a_mode, DelayedNeutrons a_delayedNeutrons );

        Mode mode( ) const { return( m_mode ); }                /**< Returns the value of the *m_mode* member. */
        void setMode( Mode a_mode ) { m_mode = a_mode; }        /**< Sets the *m_mode* member to *a_mode*. */

        bool useMultiGroupSummedData( ) const { return( m_useMultiGroupSummedData ); }  /**< Returns the value of the *m_useMultiGroupSummedData* member. */
        void setUseMultiGroupSummedData( bool a_useMultiGroupSummedData ) { m_useMultiGroupSummedData = a_useMultiGroupSummedData; }
                                                                /**< Sets the *m_useMultiGroupSummedData* member to *a_useMultiGroupSummedData*. */

        Form const *form( LUPI::StatusMessageReporting &a_smr, GIDI::Suite const &a_suite, Styles::TemperatureInfo const &a_temperatureInfo,
                std::string a_dataType ) const ;
};

}           // End of namespace Transporting.

namespace GRIN {

/*
============================================================
================= InelasticIncidentEnergy ==================
============================================================
*/

class InelasticIncidentEnergy : public Form {

    private:
        double m_energy;
        std::string m_unit;
        Table::Table m_table;

    public:
        InelasticIncidentEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~InelasticIncidentEnergy( );

        double energy( ) const { return( m_energy ); }
        std::string const &unit( ) const { return( m_unit ); }
        Table::Table const &table( ) const { return( m_table ); }
};

/*
============================================================
================= CaptureLevelProbability ==================
============================================================
*/

class CaptureLevelProbability : public Form {

    private:
        double m_probabilty;
        double m_spin;
        std::string m_spinUnit;
        int m_parity;
        std::string m_capturePrimaryToContinua;
        Table::Table m_table;

    public:
        CaptureLevelProbability( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~CaptureLevelProbability( );

        double probabilty( ) const { return( m_probabilty ); }
        double spin( ) const { return( m_spin ); }
        std::string const &spinUnit( ) const { return( m_spinUnit ); }
        int parity( ) const { return( m_parity ); }
        std::string const &capturePrimaryToContinua( ) const { return( m_capturePrimaryToContinua ); }
        Table::Table const &table( ) const { return( m_table ); }
};

/*
============================================================
=================== GRIN_continuumGammas ===================
============================================================
*/

class GRIN_continuumGammas : public GUPI::Ancestry {

    private:
        PhysicalQuantity m_captureNeutronSeparationEnergy;
        PhysicalQuantity m_maximumCaptureIncidentEnergy;
        PoPI::Database m_pops;
        Suite m_inelasticIncidentEnergies;
        Suite m_captureLevelProbabilities;
        std::string m_captureResidualId;                                                /**< The GNDS PoPs' id of the heavy capture residual particle. */
        int m_captureResidualIntid;                                                     /**< The intid of the heavy capture residual particle. */
        int m_captureResidualIndex;                                                     /**< The PoPI index of the heavy capture residual particle. */
        double m_captureResidualMass;                                                   /**< The mass if the heavy capture residual particle. */

    public:
        GRIN_continuumGammas( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, 
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, ProtareSingle const &a_protare, Styles::Suite const *a_styles );
        ~GRIN_continuumGammas( );

        PhysicalQuantity const &captureNeutronSeparationEnergy( ) const { return( m_captureNeutronSeparationEnergy ); }
        PhysicalQuantity const &maximumCaptureIncidentEnergy( ) const { return( m_maximumCaptureIncidentEnergy ); }
        PoPI::Database const &pops( ) const { return( m_pops ); }
        Suite const &inelasticIncidentEnergies( ) const { return( m_inelasticIncidentEnergies ); }
        Suite const &captureLevelProbabilities( ) const { return( m_captureLevelProbabilities ); }
        std::string captureResidualId( ) const { return( m_captureResidualId ); }       /**< Returns the value of the *m_captureResidualId* member. */
        int captureResidualIntid( ) const { return( m_captureResidualIntid ); }         /**< Returns the value of the *m_captureResidualIntid* member. */
        int captureResidualIndex( ) const { return( m_captureResidualIndex ); }         /**< Returns the value of the *m_captureResidualIndex* member. */
        double captureResidualMass( ) const { return( m_captureResidualMass ); }        /**< Returns the value of the *m_captureResidualMass* member. */

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;
};

}               // End of namespace GRIN.

/*
============================================================
========================= Product ==========================
============================================================
*/
class Product : public Form {

    private:
        ParticleInfo m_particle;                    /**< The products *ParticleInfo* data. */
        ParticleInfo m_GNDS_particle;               /**< The products *ParticleInfo* data. This is the product's equivalent of the Protare::m_GNDS_target member. */

        int m_productMultiplicity;                  /**< Product integer multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer. */
        bool m_treatProductAsIfInfinityMass;        /**< If **true**, the product is photo-atomic or TNSL target and should be handled as if it has infinite mass. Ths is, energy and momentum data are returned with 0 value. */
        Component m_multiplicity;                   /**< The GNDS <**multiplicity**> node. */
        Component m_distribution;                   /**< The GNDS <**distribution**> node. */
        Component m_averageEnergy;                  /**< The GNDS <**averageEnergy**> node. */
        Component m_averageMomentum;                /**< The GNDS <**averageMomentum**> node. */
        OutputChannel *m_outputChannel;             /**< The GNDS <**outputChannel**> node if present. */

    public:
        Product( PoPI::Database const &a_pops, std::string const &a_productID, std::string const &a_label );
        Product( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, 
                PoPI::Database const &a_internalPoPs, Suite *a_parent, Styles::Suite const *a_styles );
        ~Product( );

        ParticleInfo const &particle( ) const { return( m_particle ); }                     /**< Returns the value of the *m_particle* member. */
        void setParticle( ParticleInfo const &a_particle ) { m_particle = a_particle; }     /**< Sets *m_particle* to *a_particle*. */
        std::string pid( ) const { return( m_particle.ID( ) ); }
        ParticleInfo const &GNDS_particle( ) const { return( m_GNDS_particle ); }           /**< Returns a const reference to the *m_GNDS_particle* member. */
        ParticleInfo &GNDS_particle( ) { return( m_GNDS_particle ); }                       /**< Returns the value of the *m_GNDS_particle* member. */
        int depth( ) const ;

        Component &multiplicity( ) { return( m_multiplicity ); }                            /**< Returns a reference to the *m_multiplicity* member. */
        Component const &multiplicity( ) const { return( m_multiplicity ); }                /**< Returns a const reference to the *m_multiplicity* member. */
        Component &distribution( ) { return( m_distribution ); }                            /**< Returns a reference to the *m_distribution* member. */
        Component const &distribution( ) const { return( m_distribution ); }                /**< Returns a reference to the *m_distribution* member. */
        Component &averageEnergy( ) { return( m_averageEnergy ); }                          /**< Returns a reference to the *m_averageEnergy* member. */
        Component const &averageEnergy( ) const { return( m_averageEnergy ); }              /**< Returns a const reference to the *m_averageEnergy* member. */
        Component &averageMomentum( ) { return( m_averageMomentum ); }                      /**< Returns a reference to the *m_averageMomentum* member. */
        Component const &averageMomentum( ) const { return( m_averageMomentum ); }          /**< Returns a const reference to the *m_averageMomentum* member. */
        OutputChannel *outputChannel( ) const { return( m_outputChannel ); }                /**< Returns a reference to the *m_outputChannel* member. */

        void modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex );

        bool hasFission( ) const ;
        bool isDelayedFissionNeutronComplete( bool a_isDelayedNeutron ) const ;
        bool areAllProductsTracked( Transporting::Particles const &a_particles ) const ;

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int productMultiplicity( std::string const &a_productID ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        bool a_final ) const ;
        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;
        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, std::string const &a_productID, 
                        int a_order ) const ;

        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        void continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const ;
        void mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const ;

        bool isCompleteParticle( ) const ;
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
        void calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
====================== DelayedNeutron ======================
============================================================
*/
class DelayedNeutron : public Form {

    private:
        int m_delayedNeutronIndex;                  /**< If this is a delayed fission neutron, this is its index. */
        Suite m_rate;                               /**< The GNDS <**rate**> node. */
        Product m_product;                          /**< The GNDS <**product**> node. */

    public:
        DelayedNeutron( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent, Styles::Suite const *a_styles );
        ~DelayedNeutron( );

        int delayedNeutronIndex( ) const { return( m_delayedNeutronIndex ); }
        void setDelayedNeutronIndex( int a_delayedNeutronIndex ) { m_delayedNeutronIndex = a_delayedNeutronIndex; }
        Suite &rate( ) { return( m_rate ); }
        Suite const &rate( ) const { return( m_rate ); }
        Product &product( ) { return( m_product ); }
        Product const &product( ) const { return( m_product ); }

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        void productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int productMultiplicity( std::string const &a_productID ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;
        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;

        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
        void continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const ;
        void mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const ;
        void calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=================== DelayedNeutronProduct ==================
============================================================
*/
class DelayedNeutronProduct {

    private:
        int m_delayedNeutronIndex;                  /**< If this is a delayed fission neutron, this is its index. */
        PhysicalQuantity m_rate;
        Product const *m_product;

    public:
        DelayedNeutronProduct( int a_delayedNeutronIndex, PhysicalQuantity a_rate, Product const *a_product ) : 
                m_delayedNeutronIndex( a_delayedNeutronIndex ),
                m_rate( a_rate ),
                m_product( a_product ) {
        }
        DelayedNeutronProduct( DelayedNeutronProduct const &a_delayedNeutronProduct ) :
                m_delayedNeutronIndex( a_delayedNeutronProduct.delayedNeutronIndex( ) ),
                m_rate( a_delayedNeutronProduct.rate( ) ),
                m_product( a_delayedNeutronProduct.product( ) ) {
        }
        ~DelayedNeutronProduct( ) {}

        int delayedNeutronIndex( ) const { return( m_delayedNeutronIndex ); }
        PhysicalQuantity rate( ) const { return( m_rate ); }
        Product const *product( ) const { return( m_product ); }
};

typedef std::vector<DelayedNeutronProduct> DelayedNeutronProducts;

/*
============================================================
==================== FissionFragmentData ===================
============================================================
*/
class FissionFragmentData : public GUPI::Ancestry {

    private:
        Suite m_delayedNeutrons;                            /**< The GNDS <**delayedNeutrons**> node. This members stores a list of DelayedNeutron instances. */
        Component m_fissionEnergyReleases;                      /**< The GNDS <**fissionEnergyReleases**> node. */

    public:
        FissionFragmentData( );
        FissionFragmentData( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Styles::Suite const *a_styles );
        ~FissionFragmentData( );

        Suite &delayedNeutrons( ) { return( m_delayedNeutrons ); }
        Suite const &delayedNeutrons( ) const { return( m_delayedNeutrons ); }
        Component &fissionEnergyReleases( ) { return( m_fissionEnergyReleases ); }
        Component const &fissionEnergyReleases( ) const { return( m_fissionEnergyReleases ); }

        bool isDelayedFissionNeutronComplete( ) const ;

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        void productIDs( std::set<std::string> &a_indices, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int productMultiplicity( std::string const &a_productID ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        bool a_final ) const ;
        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;

        void delayedNeutronProducts( DelayedNeutronProducts &a_delayedNeutronProducts ) const ;
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
        void continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const ;
        void mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const ;
        void calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
======================= OutputChannel ======================
============================================================
*/
class OutputChannel : public GUPI::Ancestry {

    private:
        bool m_twoBody;                                     /**< true if the output channel is two-body and false otherwise. */
        bool m_fissions;                                    /**< true if the output channel is a fission channel and false otherwise. */
        std::string m_process;                              /**< The GNDS *process* attribute for the channel. */

        Component m_Q;                                      /**< The GNDS <**Q**> node. */
        Suite m_products;                                   /**< The GNDS <**products**> node. */
        FissionFragmentData m_fissionFragmentData;          /**< The GNDS <**fissionFragmentData**> node. */
        Construction::FissionResiduals m_fissionResiduals;  /**< This member specifies what fission redisual products will be added to the list of products produced in a fission reaction. */

    public:
        OutputChannel( bool a_twoBody, bool a_fissions, std::string a_process );
        OutputChannel( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                PoPI::Database const &a_internalPoPs, Styles::Suite const *a_styles, bool a_isFission, bool a_addFissionResiduals );
        ~OutputChannel( );

        bool twoBody( ) const { return( m_twoBody ); }                              /**< Returns the value of the *m_twoBody* member. */
        std::string process( ) const { return( m_process ); }                       /**< Returns the value of the *m_process* member. */
        int depth( ) const ;

        Component &Q( ) { return( m_Q ); }                                          /**< Returns a reference to the *m_Q* member. */
        Component const &Q( ) const { return( m_Q ); }                              /**< Returns a reference to the *m_Q* member. */
        Suite &products( ) { return( m_products ); }                                /**< Returns a reference to the *m_products* member. */
        Suite const &products( ) const { return( m_products ); }                    /**< Returns a reference to the *m_products* member. */
        FissionFragmentData &fissionFragmentData( ) { return( m_fissionFragmentData ); }
        FissionFragmentData const &fissionFragmentData( ) const { return( m_fissionFragmentData ); }

        Construction::FissionResiduals fissionResiduals( ) const { return( m_fissionResiduals ); }  /**< Returns the value of the *m_fissionResiduals* member. */

        void modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex );
        bool areAllProductsTracked( Transporting::Particles const &a_particles ) const ;

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        bool isFission( ) const { return( m_fissions ); }                           /**< Returns true if the output channel is a fission output channel. */
        bool hasFission( ) const ;
        bool isDelayedFissionNeutronComplete( ) const ;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int productMultiplicity( std::string const &a_productID ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        bool a_final ) const ;
        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;

        void delayedNeutronProducts( DelayedNeutronProducts &a_delayedNeutronProducts ) const { m_fissionFragmentData.delayedNeutronProducts( a_delayedNeutronProducts ); }
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
        void continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const ;
        void mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const ;
        void calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, 
                Functions::XYs1d const &a_crossSectionXYs1d );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

namespace DoubleDifferentialCrossSection {

namespace n_ThermalNeutronScatteringLaw {

/*
============================================================
==================== IncoherentInelastic ===================
============================================================
*/
class IncoherentInelastic : public Base {
    
    private:
        Options m_options;                              /**< Options for *this*. */
        Suite m_scatteringAtoms;                        /**< The list of atoms and there information. */
        S_alpha_beta m_S_alpha_beta;                    /**< The S(alpha,beta,T) function. */
    
    public:
        IncoherentInelastic( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Suite *a_parent );
        ~IncoherentInelastic( );
        
        Options &options( ) { return( m_options ); }                                    /**< Returns the value of the *m_options* */
        Suite &scatteringAtoms( ) { return( m_scatteringAtoms ); }                      /**< Returns the value of the *m_scatteringAtoms* */
        S_alpha_beta const &s_alpha_beta( ) const { return( m_S_alpha_beta ); }         /**< Returns the value of the *m_S_alpha_beta* */
};

}               // End namespace n_ThermalNeutronScatteringLaw.

}               // End namespace DoubleDifferentialCrossSection.

namespace ACE_URR {

/*
============================================================
====================== IncidentEnergy ======================
============================================================
*/
class IncidentEnergy: public Form {

    private:
        double m_value;
        std::string m_unit;
        Table::Table m_table;

    public:
        IncidentEnergy( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~IncidentEnergy( );

        double value( ) const { return( m_value ); }      /**< Returns the value of the *m_value* member. */
        std::string const &unit( ) const { return( m_unit ); }  /**< Returns a *const* reference to the *m_unit* member. */
        Table::Table const &table( ) const { return( m_table ); }            /**< Returns a *const* reference to the *m_table* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

/*
============================================================
===================== ProbabilityTable =====================
============================================================
*/
class ProbabilityTable : public Form {

    public:
        typedef std::vector<IncidentEnergy *> Forms;                              /**< The typedef the the *m_forms* member. */

    private:
        mutable Forms m_forms;                                          /**< The list of nodes stored within *this*. */

    public:
        ProbabilityTable( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~ProbabilityTable( );

        Forms &forms( ) { return( m_forms ); }                          /**< Returns a *const* reference to the *m_forms* member. */
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent ) const ;
};

}               // End namespace ACE_URR.

/*
============================================================
========================= Reaction =========================
============================================================
*/
class Reaction : public Form {

    friend class ProtareSingle;

    private:
        mutable int m_reactionIndex;                    /**< The index of the reaction in the ProtareSingle. */
        bool m_active;                                  /**< If true, this reaction is used for calcualtion (e.g., its cross section is added to the total for its protare), otherwise, this reaction is ignored. */
        int m_ENDF_MT;                                  /**< The ENDF MT value for the reaction. */
        int m_ENDL_C;                                   /**< The ENDL C value for the reaction. */
        int m_ENDL_S;                                   /**< The ENDL S value for the reaction. */
        std::string m_fissionGenre;                     /**< If the reaction is fission, this is its genre. */
        double m_QThreshold;                            /**< Threshold value calculated from the Q and the protare's m_thresholdFactor. */
        double m_crossSectionThreshold;                 /**< Threshold value derived from cross section data via *evaluated* or *griddedCrossSection*. */
        double m_twoBodyThreshold;                      /**< This is the T_1 value needed by MCGIDI to do two-body kinematics (i.e., in the equation (K_{com,3_4} = m_2 * (K_1 - T_1) / (m_1 + m_2)). */
        bool m_isPairProduction;                        /**< Kludge! Currently needed because GNDS specification unclear about how to specify photo-atomic pair production reaction. */
        bool m_isPhotoAtomicIncoherentScattering;       /**< **true** if the reaction is photo-atomic incoherent scattering and **false** otherwise. Helpful for MCGIDI. */
        bool m_RutherfordScatteringPresent;             /**> For charged particle elastic scattering, this member is *true* if Rutherford scattering is present and *false* otherwise. */
        bool m_onlyRutherfordScatteringPresent;         /**> For charged particle elastic scattering, this member is *true* if only Rutherford scattering is present and *false* otherwise. */
        bool m_nuclearPlusInterferencePresent;          /**> For charged particle elastic scattering, this member is *true* if nuclear plus interference is present and *false* otherwise. */
        bool m_decayPositronium;                        /**< If **true**, whenever a positron is created, it is assumed to immediately form positronium and decay into 2 511 KeV photons. Ergo, the photons are produced in the reaction and not a positron. */

        Component m_doubleDifferentialCrossSection;     /**< The GNDS <**doubleDifferentialCrossSection**> node. */
        Component m_crossSection;                       /**< The GNDS <**crossSection**> node. */
        Component m_availableEnergy;                    /**< The GNDS <**availableEnergy**> node. */
        Component m_availableMomentum;                  /**< The GNDS <**availableMomentum**> node. */
        OutputChannel *m_outputChannel;                 /**< The reaction's output channel. */
        void setReactionIndex( int a_reactionIndex ) const 
                { m_reactionIndex = a_reactionIndex ; } /**< Sets *m_reactionIndex* to *a_reactionIndex*. */

    public:
        Reaction( int a_ENDF_MT, std::string a_fissionGenre );
        Reaction( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, Protare const &a_protare,
                        Styles::Suite const *a_styles );
        ~Reaction( );

        bool active( ) const { return( m_active ); }                                    /**< Returns the value of the *m_active* member. */
        void setActive( bool a_active ) { m_active = a_active; }                        /**< Sets *m_active* to *a_active*. */
        int reactionIndex( ) const { return( m_reactionIndex ); }                       /**< Returns the value of the *m_reactionIndex* member. */
        int depth( ) const { return( m_outputChannel->depth( ) ); }                     /**< Returns the maximum product depth for this reaction. */
        int ENDF_MT( ) const { return( m_ENDF_MT ); }                                   /**< Returns the value of the *m_ENDF_MT* member. */
        int ENDL_C( ) const { return( m_ENDL_C ); }                                     /**< Returns the value of the *m_ENDL_C* member. */
        int ENDL_S( ) const { return( m_ENDL_S ); }                                     /**< Returns the value of the *m_ENDL_S* member. */
        std::string const &fissionGenre( ) const { return( m_fissionGenre ); }
        bool isPairProduction( ) const { return( m_isPairProduction ); }                /**< Returns the value of the *m_isPairProduction* member. */
        bool isPhotoAtomicIncoherentScattering( ) const { return( m_isPhotoAtomicIncoherentScattering ); }                /**< Returns the value of the *m_isPhotoAtomicIncoherentScattering* member. */
        bool RutherfordScatteringPresent( ) const { return( m_RutherfordScatteringPresent ); }
                                                                                        /**< Returns the value of *m_RutherfordScatteringPresent* member. */
        bool onlyRutherfordScatteringPresent( ) const { return( m_onlyRutherfordScatteringPresent ); }
                                                                                        /**< Returns the value of *m_onlyRutherfordScatteringPresent* member. */
        bool nuclearPlusInterferencePresent( ) const { return( m_nuclearPlusInterferencePresent ); }
                                                                                        /**< Returns the value of *m_nuclearPlusInterferencePresent* member. */

        Component &doubleDifferentialCrossSection( ) { return( m_doubleDifferentialCrossSection ); }    /**< Returns a reference to the *m_doubleDifferentialCrossSection* member. */
        Component const &doubleDifferentialCrossSection( ) const { return( m_doubleDifferentialCrossSection ); }    /**< Returns a reference to the *m_doubleDifferentialCrossSection* member. */
        Component &crossSection( ) { return( m_crossSection ); }                            /**< Returns a reference to the *m_crossSection* member. */
        Component const &crossSection( ) const { return( m_crossSection ); }                            /**< Returns a reference to the *m_crossSection* member. */

        Component &availableEnergy( ) { return( m_availableEnergy ); }                      /**< Returns a reference to the *m_availableEnergy* member. */
        Component const &availableEnergy( ) const { return( m_availableEnergy ); }          /**< Returns a reference to the *m_availableEnergy* member. */
        Component &availableMomentum( ) { return( m_availableMomentum ); }                  /**< Returns a reference to the *m_availableMomentum* member. */
        Component const &availableMomentum( ) const { return( m_availableMomentum ); }      /**< Returns a reference to the *m_availableMomentum* member. */

        OutputChannel *outputChannel( ) const { return( m_outputChannel ); }            /**< Returns a reference to the *m_outputChannel* member. */
        void setOutputChannel( OutputChannel *a_outputChannel );

        void modifiedMultiGroupElasticForTNSL( std::map<std::string,std::size_t> a_maximumTNSL_MultiGroupIndex );

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;
        std::string xlinkItemKey( ) const { return( GUPI::Ancestry::buildXLinkItemKey( GIDI_labelChars, label( ) ) ); }   /**< Returns the result of calling "GUPI::Ancestry::buildXLinkItemKey( GIDI_labelChars, label() )". */

        bool hasFission( ) const ;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int productMultiplicity( std::string const &a_productID ) const {
                return( m_outputChannel->productMultiplicity( a_productID ) ); }               /**< Returns the product multiplicity (e.g., 0, 1, 2, ...) or -1 if energy dependent or not an integer. */
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        double threshold( ) const { return( m_QThreshold ); }                           /**< Returns the value of the *m_QThreshold* member. */
        double crossSectionThreshold( ) const { return( m_crossSectionThreshold ); }    /**< Returns the value of the *m_crossSectionThreshold* member. */
        double twoBodyThreshold( ) const { return( m_twoBodyThreshold ); }              /**< Returns the value of the *m_twoBodyThreshold* member. */

        bool areAllProductsTracked( Transporting::Particles const &a_particles ) const ;

        Vector multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo ) 
                        const ;
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        bool a_final ) const ;
        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;

        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles, std::string const &a_productID, int a_order ) const ;
        Matrix multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles, int a_order ) const ;

        Vector multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo ) 
                        const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles ) const ;

        Vector multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo ) 
                        const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        Transporting::Particles const &a_particles ) const ;

        Vector multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID, std::string const &a_projectileID ) const ;

        void delayedNeutronProducts( DelayedNeutronProducts &a_delayedNeutronProducts ) const ;
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
        void continuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, double a_energy, 
                double &a_productEnergy, double &a_productMomentum, double &a_productGain, bool a_ignoreIncompleteParticles ) const ;
        void mapContinuousEnergyProductData( Transporting::Settings const &a_settings, std::string const &a_particleID, 
                std::vector<double> const &a_energies, int a_offset, std::vector<double> &a_productEnergies, std::vector<double> &a_productMomenta, 
                std::vector<double> &a_productGains, bool a_ignoreIncompleteParticles ) const ;

        bool modifyCrossSection( Functions::XYs1d const *a_offset, Functions::XYs1d const *a_slope, bool a_updateMultiGroup = false );
        bool modifiedCrossSection( Functions::XYs1d const *a_offset, Functions::XYs1d const *a_slope );
        void recalculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo );
        void calculateMultiGroupData( ProtareSingle const *a_protare, Styles::TemperatureInfo const &a_temperatureInfo, 
                std::string const &a_heatedMultiGroupLabel, MultiGroupCalulationInformation const &a_multiGroupCalulationInformation );

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

namespace Sums {

namespace Summand {

/*
============================================================
=========================== Base ===========================
============================================================
*/
class Base : public GUPI::Ancestry {

    private:
        std::string m_href;                                                     /**< xlink for the summand. */

    public:
        Base( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Base( );

        std::string const &href( ) const { return( m_href ); }                  /**< Returns the value of the *m_href* member. */
        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Add ============================
============================================================
*/
class Add : public Base {

    public:
        Add( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
};

}           // End of namespace Summand.

/*
============================================================
========================= Summands =========================
============================================================
*/
class Summands : public Form {

    private:
        std::vector<Summand::Base *> m_summands;                            /**< List of summand for *this*. */

    public:
        Summands( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo );
        ~Summands( );

        std::size_t size( ) const { return( m_summands.size( ) ); }         /**< Returns the number of summands in *this*. */
        Summand::Base const *operator[]( std::size_t a_index ) const { return( m_summands[a_index] ); } /**< Returns the summand at index *a_index*. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Base ===========================
============================================================
*/
class Base : public Form {

    private:
        int m_ENDF_MT;                                                      /**< ENDF MT value for the sum. */
        Summands m_summands;                                                /**< List of Summands for *this*. */

    public:
        Base( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                FormType a_type );

        int ENDF_MT( ) const { return( m_ENDF_MT ); }                       /**< Returns the value of the *m_ENDF_MT* member. */
        Summands const &summands( ) const { return( m_summands ); }         /**< Returns the value of the *m_summands* member. */
};

/*          
============================================================
====================== CrossSectionSum =====================
============================================================
*/
class CrossSectionSum : public Base {

    private:
        Component m_Q;                                                          /**< The GNDS <**Q**> node. */
        Component m_crossSection;                                               /**< The GNDS <**crossSection**> node. */

    public:
        CrossSectionSum( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs );
        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        Component &Q( ) { return( m_Q ); }                                      /**< Returns a reference to the *m_Q* member. */
        Component &crossSection( ) { return( m_crossSection ); }                /**< Returns a reference to the *m_crossSection* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*          
============================================================
====================== MultiplicitySum =====================
============================================================
*/
class MultiplicitySum : public Base {

    private:
        Suite m_multiplicity;                                               /**< The GNDS <**multiplicity**> node. */

    public:
        MultiplicitySum( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs );

        Suite &multiplicity( ) { return( m_multiplicity ); }                /**< Returns a reference to the *m_multiplicity* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Sums ===========================
============================================================
*/
class Sums : public GUPI::Ancestry {

    private:
        Suite m_crossSectionSums;                                               /**< The GNDS <**crossSectionSums**> node. */
        Suite m_multiplicitySums;                                               /**< The GNDS <**multiplicitySums**> node. */

    public:
        Sums( );
        ~Sums( );

        Suite &crossSectionSums( ) { return( m_crossSectionSums ); }                /**< Returns the value of the *m_crossSectionSums* member. */
        Suite const &crossSectionSums( ) const { return( m_crossSectionSums ); }    /**< Returns the value of the *m_crossSectionSums* member. */
        Suite &multiplicitySums( ) { return( m_multiplicitySums ); }                /**< Returns the value of the *m_multiplicitySums* member. */
        Suite const &multiplicitySums( ) const { return( m_multiplicitySums ); }    /**< Returns the value of the *m_multiplicitySums* member. */

        void parse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs );
        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

}           // End of namespace Sums.

/*
============================================================
========================== Protare =========================
============================================================
*/
class Protare : public GUPI::Ancestry {

    private:
        ParticleInfo m_projectile;              /**< Information about the projectile. */
        ParticleInfo m_target;                  /**< Information about the target. */
        ParticleInfo m_GNDS_target;             /**< Information about the target as specified in the GNDS file. For example, for requested target 'H1' for a photo-atomic GNDS file, the GNDS target will be 'H'. */

    protected:
        void initialize( HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, bool a_targetRequiredInGlobalPoPs, 
                        bool a_requiredInPoPs = true );

    public:
        Protare( );
        ~Protare( );

        ParticleInfo const &projectile( ) const { return( m_projectile ); }         /**< Returns the value of the *m_projectile* member. */
        void setProjectile( ParticleInfo const &a_projectile ) { m_projectile = a_projectile; }    /**< Sets *m_projectile* to *a_projectile*. */
        ParticleInfo const &target( ) const { return( m_target ); }                 /**< Returns the value of the *m_target* member. */
        void setTarget( ParticleInfo const &a_target ) { 
                m_target = a_target;
                if( m_GNDS_target.ID( ) == "" ) m_GNDS_target = a_target; }         /**< Sets *m_target* to *a_target* and m_GNDS_target if it is an empty string. */
        ParticleInfo const &GNDS_target( ) const { return( m_GNDS_target ); }       /**< Returns the value of the *m_GNDS_target* member. */

        virtual ProtareType protareType( ) const = 0;                               /**< Returns the type of the protare. */
        virtual bool isTNSL_ProtareSingle( ) const { return( false ); }             /**< Returns *true* if the instance is a ProtareSingle instance with only TNSL data and *false* otherwise. */
        virtual std::size_t numberOfProtares( ) const = 0;                          /**< Returns the number of protares contained in *this*. */
        virtual ProtareSingle *protare( std::size_t a_index ) = 0;                  /**< Returns the *a_index* - 1 Protare contained in *this*. */
        virtual ProtareSingle const *protare( std::size_t a_index ) const = 0;      /**< Returns the *a_index* - 1 Protare contained in *this*. */

        virtual LUPI::FormatVersion const &formatVersion( std::size_t a_index = 0 ) const = 0;
        virtual std::string const &fileName( std::size_t a_index = 0 ) const = 0;
        virtual std::string const &realFileName( std::size_t a_index = 0 ) const = 0;

        virtual std::vector<std::string> libraries( std::size_t a_index = 0 ) const = 0;
        virtual std::string const &evaluation( std::size_t a_index = 0 ) const = 0;
        virtual Frame projectileFrame( std::size_t a_index = 0 ) const = 0;
        virtual int numberOfLazyParsingHelperForms( ) const = 0;
        virtual int numberOfLazyParsingHelperFormsReplaced( ) const = 0;
        virtual double thresholdFactor( ) const = 0;

        virtual Documentation_1_10::Suite &documentations( ) = 0;

        virtual Styles::Base &style( std::string const a_label ) = 0;
        virtual Styles::Suite &styles( ) = 0;
        virtual Styles::Suite const &styles( ) const = 0;

        virtual int intid( std::string const &a_id ) const = 0;
        virtual void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const = 0;
        virtual int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const = 0;

        virtual Styles::TemperatureInfos temperatures( ) const  = 0;

        virtual std::size_t numberOfReactions( ) const = 0;
        virtual Reaction *reaction( std::size_t a_index ) = 0;
        virtual Reaction const *reaction( std::size_t a_index ) const = 0;
        virtual Reaction const *reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual std::size_t numberOfOrphanProducts( ) const = 0;
        virtual Reaction *orphanProduct( std::size_t a_index ) = 0;
        virtual Reaction const *orphanProduct( std::size_t a_index ) const = 0;
        virtual void updateReactionIndices( int a_offset ) const = 0;

        virtual bool hasFission( ) const = 0;
        virtual bool isDelayedFissionNeutronComplete( ) const = 0;

        virtual GUPI::Ancestry *findInAncestry3( std::string const &a_item ) = 0;
        virtual GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const = 0;

        virtual std::vector<double> groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const = 0;
        virtual Vector multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo ) const = 0;

        virtual Vector multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic = true,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                        std::string const &a_productID, int a_order, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Matrix multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order, 
                        TransportCorrectionType a_transportCorrectionType, double a_temperature,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual Vector multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual Vector multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;
        virtual Vector multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual Vector multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const = 0;

        virtual void TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::XYs1d &a_crossSectionSum );
        virtual void TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::Ys1d &a_crossSectionSum );
        virtual void TNSL_crossSectionSumCorrection( std::string const &a_label, Vector &a_crossSectionSum );

        virtual stringAndDoublePairs muCutoffForCoulombPlusNuclearElastic( ) const = 0;
        virtual DelayedNeutronProducts delayedNeutronProducts( ) const = 0;
        virtual void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const = 0;
        ExcludeReactionsSet reactionIndicesMatchingENDLCValues( std::set<int> const &a_CValues, bool a_checkActiveState = true );
};

/*
============================================================
====================== ProtareSingle =======================
============================================================
*/
class ProtareSingle : public Protare {

    private:
        HAPI::File *m_doc;                      /**< If data read from file, this member is a pointer to the opened **HAPI::File** instance. */
        HAPI::DataManager *m_dataManager;       /**< If data read from hybrid file, this member is a pointer to the **HAPI::DataManager** instance. */
        int m_numberOfLazyParsingHelperForms;  /**< This counts the number of LazyParsingHelperForms instantiated. */
        int m_numberOfLazyParsingHelperFormsReplaced;   /**< This counts the number of LazyParsingHelperForms replaced with the appropriate form. */
        LUPI::FormatVersion m_formatVersion;    /**< Store the GNDS format version. */
        PoPI::Database m_internalPoPs;          /**< The *PoPs* specified under the protare (e.g., reactionSuite) node. */

        std::vector<std::string> m_libraries;   /**< The list of libraries *this* was found in. */
        std::string m_evaluation;               /**< The protare's evaluation string. */
        std::string m_interaction;              /**< The protare's interaction string. */
        std::string m_fileName;                 /**< The path to the protare's file. May be relative. */
        std::string m_realFileName;             /**< The real path to the protare's file. Equivalent to the value returned by the C-function *realpath( m_fileName )* on Unix systems. */
        Frame m_projectileFrame;                /**< The frame the projectile data are given in. */
        double m_projectileEnergyMin;           /**< The projectile's minimum energy for which data are complete as specified in the evaluated style. */
        double m_projectileEnergyMax;           /**< The projectile's maximum energy for which data are complete as specified in the evaluated style. */
        bool m_isTNSL_ProtareSingle;            /**< If *this* is a ProtareSingle instance with TNSL data *true* and otherwise *false*. */
        bool m_isPhotoAtomic;                   /**< true if photo-atomic protare and false otherwise. */
        bool m_decayPositronium;                /**< If **true**, whenever a positron is created, it is assumed to immediately form positronium and decay into 2 511 KeV photons. Ergo, the photons are produced in the reaction and not a positron. */

        double m_thresholdFactor;               /**< The non-relativistic factor that converts a Q-value into a threshold. */

        PoPI::NuclideGammaBranchStateInfos m_nuclideGammaBranchStateInfos;  /**< Simplified list of gamma branching data from nuclide level decays derived from the internal PoPI::Database. */

        ExternalFiles::Suite m_externalFiles;   /**< The GNDS <**externalFiles**> node. */
        Styles::Suite m_styles;                 /**< The GNDS <**styles**> node. */
        Documentation_1_10::Suite m_documentations;  /**< The GNDS <**documentations**> node. */
        Suite m_reactions;                      /**< The GNDS <**reactions**> node. */
        Suite m_orphanProducts;                 /**< The GNDS <**orphanProducts**> node. */
        Suite m_incompleteReactions;            /**< The GNDS <**incompleteReactions**> node. */

        Sums::Sums m_sums;                      /**< The GNDS <**sums**> node. */
        Suite m_fissionComponents;              /**< The GNDS <**fissionComponents**> node. */

        bool m_RutherfordScatteringPresent;     /**> For charged particle elastic scattering, this member is *true* if Rutherford scattering is present and *false* otherwise. */
        bool m_onlyRutherfordScatteringPresent; /**> For charged particle elastic scattering, this member is *true* if only Rutherford scattering is present and *false* otherwise. */

//  The following are non-GNDS 2.0 data types that are stored in the applicationData node.
        Reaction *m_nuclearPlusCoulombInterferenceOnlyReaction;     /**< The nuclear + interference (ENDL C=9) reaction in the applicationData node. */
        Reaction *m_multiGroupSummedReaction;                       /**< This reaction contains the sum multi-group data from all other reactions. */
        OutputChannel *m_multiGroupSummedDelayedNeutrons;           /**< This reaction contains the sum multi-group data from delayed neutrons. */
        Suite m_ACE_URR_probabilityTables;                          /**< This suite stores ACE style URR probability tables. */
        Suite m_photoAtomicIncoherentDoppler;                       /**< This suite stores the data for the impulse approximation photon doppler broadening reaction (MT 1534-1572) */
        Component m_pointwiseAverageProductEnergy;                  /**< This suite stores upscatter model B pointwise energy deposition data for the outgoing neutron. */
        GRIN::GRIN_continuumGammas *m_GRIN_continuumGammas;         /**< This stores continuum gamma information from the GRIN project. */

        void initialize( );
        void initialize( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                bool a_targetRequiredInGlobalPoPs, bool a_requiredInPoPs = true );

        bool useMultiGroupSummedData( Transporting::MG const &a_settings, ExcludeReactionsSet const &a_reactionsToExclude ) const ;
        bool useMultiGroupSummedDelayedNeutronsData( Transporting::MG const &a_settings ) const ;

    public:
        ProtareSingle( PoPI::Database const &a_pops, std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_evaluation,
                std::string const &a_interaction, std::string const &a_formatVersion = GNDS_formatVersion_1_10Chars );
        ProtareSingle( Construction::Settings const &a_construction, std::string const &a_fileName, FileType a_fileType, PoPI::Database const &a_pops, 
                ParticleSubstitution const &a_particleSubstitution, std::vector<std::string> const &a_libraries, std::string const &a_interaction,
                bool a_targetRequiredInGlobalPoPs = true, bool a_requiredInPoPs = true );
        ProtareSingle( Construction::Settings const &a_construction, HAPI::Node const &a_protare, PoPI::Database const &a_pops, 
                ParticleSubstitution const &a_particleSubstitution, std::vector<std::string> const &a_libraries, std::string const &a_interaction,
                bool a_targetRequiredInGlobalPoPs = true, bool a_requiredInPoPs = true );
        ~ProtareSingle( );

        PoPI::NuclideGammaBranchStateInfos const &nuclideGammaBranchStateInfos( ) const { return( m_nuclideGammaBranchStateInfos ); }
                                                                                    /**< Returns the value of the *m_nuclideGammaBranchStateInfos* member. */

        HAPI::DataManager *dataManager( ) { return( m_dataManager ); }              /**< Returns the value of the *m_dataManager* member. */
        void setDataManager( HAPI::DataManager *a_dataManager ) { m_dataManager = a_dataManager; }
                                                                                    /**< Sets the member *m_dataManager* to *a_dataManager*. */
        void incrementNumberOfLazyParsingHelperForms( ) { ++m_numberOfLazyParsingHelperForms; }
                                                /**> Increments the *m_numberOfLazyParsingHelperForms* member of this by 1. */
        void incrementNumberOfLazyParsingHelperFormsReplaced( ) { ++m_numberOfLazyParsingHelperFormsReplaced; }
                                                /**> Increments the *m_numberOfLazyParsingHelperFormsReplaced* member of this by 1. */

        double projectileEnergyMin( ) const { return( m_projectileEnergyMin ); }
        double projectileEnergyMax( ) const { return( m_projectileEnergyMax ); }
        bool isTNSL_ProtareSingle( ) const { return( m_isTNSL_ProtareSingle ); }    /**< Returns *true* if the instance is a ProtareSingle instance with only TNSL data and *false* otherwise. */
        bool isPhotoAtomic( ) const { return( m_isPhotoAtomic ); }                  /**< Returns the value of the *m_isPhotoAtomic* member. */

        Suite &reactions( ) { return( m_reactions ); }                              /**< Returns a reference to the *m_reactions* member. */
        Suite const &reactions( ) const { return( m_reactions ); }                  /**< Returns a *const* reference to the *m_reactions* member. */
        Suite &orphanProducts( ) { return( m_orphanProducts ); }                    /**< Returns a reference to the *m_orphanProducts* member. */
        Suite const &orphanProducts( ) const { return( m_orphanProducts ); }        /**< Returns a *const* reference to the *m_orphanProducts* member. */
        Suite &incompleteReactions( ) { return( m_incompleteReactions ); }          /**< Returns a reference to the *m_incompleteReactions* member. */
        Suite const &incompleteReactions( ) const { return( m_incompleteReactions ); }  /**< Returns a *const* reference to the *m_incompleteReactions* member. */

        Sums::Sums &sums( ) { return( m_sums ); }                                   /**< Returns a reference to the *m_sums* member. */
        Sums::Sums const &sums( ) const { return( m_sums ); }                       /**< Returns a reference to the *m_sums* member. */
        Suite &fissionComponents( ) { return( m_fissionComponents ); }              /**< Returns a reference to the *m_fissionComponents* member. */

        bool RutherfordScatteringPresent( ) const { return( m_RutherfordScatteringPresent ); }
                                                                                    /**< Returns the value of *m_RutherfordScatteringPresent*. */
        bool onlyRutherfordScatteringPresent( ) const { return( m_onlyRutherfordScatteringPresent ); }
                                                                                    /**< Returns the value of *m_onlyRutherfordScatteringPresent*. */
        Reaction const *nuclearPlusCoulombInterferenceOnlyReaction( ) const { return( m_nuclearPlusCoulombInterferenceOnlyReaction ); }
                                                                                    /**< Returns the *m_nuclearPlusCoulombInterferenceOnlyReaction* member which is a pointer. */
        Reaction const *checkIf_nuclearPlusCoulombInterferenceWanted( Transporting::MG const &a_settings, Reaction const *a_reaction ) const ;
        Reaction const *reactionToMultiGroup( Transporting::MG const &a_settings, std::size_t a_index,
                ExcludeReactionsSet const &a_reactionsToExclude ) const ;
        Reaction const *multiGroupSummedReaction( ) const { return( m_multiGroupSummedReaction ); }     /**< Returns the *m_multiGroupSummedReaction* member which is a pointer. */
        OutputChannel const *multiGroupSummedDelayedNeutrons( ) const { return( m_multiGroupSummedDelayedNeutrons ); }  /**< Returns the *m_multiGroupSummedReaction* member which is a pointer. */
        Suite const &ACE_URR_probabilityTables( ) const { return( m_ACE_URR_probabilityTables ); }      /**< Returns a *const* reference to the *m_ACE_URR_probabilityTables* member. */
        Suite const &photoAtomicIncoherentDoppler( ) const { return( m_photoAtomicIncoherentDoppler ); }
        GRIN::GRIN_continuumGammas const *GRIN_continuumGammas2( ) const { return( m_GRIN_continuumGammas ); }    /**< Returns a *const* pointer to the *m_GRIN_continuumGammas* member. */

// The rest are virtual methods defined in the Protare class.

        ProtareType protareType( ) const { return( ProtareType::single ); }                                 /**< Returns the type of the protare. */
        std::size_t numberOfProtares( ) const { return( 1 ); }                                              /**< Returns 1. */
        ProtareSingle *protare( std::size_t a_index );
        ProtareSingle const *protare( std::size_t a_index ) const ;

        LUPI::FormatVersion const &formatVersion( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_formatVersion ); }  /**< Returns the value of the *m_formatVersion* member. */
        std::string const &fileName( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_fileName ); }              /**< Returns the value of the *m_fileName* member. */
        std::string const &realFileName( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_realFileName ); }      /**< Returns the value of the *m_realFileName* member. */

        std::vector<std::string> libraries( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_libraries ); }      /**< Returns the libraries that *this* resided in. */
        std::string const &evaluation( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_evaluation ); }          /**< Returns the value of the *m_evaluation* member. */
        std::string const &interaction( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_interaction ); }        /**< Returns the value of the *m_interaction* member. */
        Frame projectileFrame( LUPI_maybeUnused std::size_t a_index = 0 ) const { return( m_projectileFrame ); }             /**< Returns the value of the *m_projectileFrame* member. */
        int numberOfLazyParsingHelperForms( ) const { return( m_numberOfLazyParsingHelperForms ); }
                                                                                    /**< Returns the value of the *m_numberOfLazyParsingHelperForms* member. */
        int numberOfLazyParsingHelperFormsReplaced( ) const { return( m_numberOfLazyParsingHelperFormsReplaced ); }
                                                                                    /**< Returns the value of the *m_numberOfLazyParsingHelperFormsReplaced* member. */
        double thresholdFactor( ) const { return( m_thresholdFactor ); }                                    /**< Returns the value of the *m_thresholdFactor* member. */

        Documentation_1_10::Suite &documentations( ) { return( m_documentations ); }                        /**< Returns the value of the *m_documentations* member. */

        ExternalFile const &externalFile( std::string const a_label ) const { return( *m_externalFiles.get<ExternalFile>( a_label ) ); }      /**< Returns the external file with label *a_label*. */
        ExternalFiles::Suite const &externalFiles( ) const { return( m_externalFiles ); }                   /**< Returns the value of the *m_externalFiles* member. */

        Styles::Base &style( std::string const a_label ) { return( *m_styles.get<Styles::Base>( a_label ) ); }              /**< Returns the style with label *a_label*. */
        Styles::Base const &style( std::string const a_label ) const { return( *m_styles.get<Styles::Base const>( a_label ) ); }  /**< Returns the const style with label *a_label*. */
        Styles::Suite &styles( ) { return( m_styles ); }                                                    /**< Returns the value of the *m_styles* member. */
        Styles::Suite const &styles( ) const { return( m_styles ); }                                        /**< Returns a *const* reference to the *m_styles* member. */

        PoPI::Database const &internalPoPs( ) const { return( m_internalPoPs ); }                           /**< Returns a *const* reference to the *m_internalPoPs* member. */
        PoPI::Database &internalPoPs( ) { return( m_internalPoPs ); }                                       /**< Returns a reference to the *m_internalPoPs* member. */

        int intid( std::string const &a_id ) const;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        Styles::TemperatureInfos temperatures( ) const ;

        std::size_t numberOfReactions( ) const { return( m_reactions.size( ) ); }                                   /**< Returns the number of reactions in the **Protare**. */
        Reaction *reaction( std::size_t a_index ) { return( m_reactions.get<Reaction>( a_index ) ); }               /**< Returns the *a_index* - 1 reaction. */
        Reaction const *reaction( std::size_t a_index ) const { return( m_reactions.get<Reaction>( a_index ) ); }   /**< Returns the *a_index* - 1 reaction. */
        Reaction const *reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const;
        std::size_t numberOfInactiveReactions( ) const ;

        std::size_t numberOfOrphanProducts( ) const { return( m_orphanProducts.size( ) ); }                         /**< Returns the number of orphan product reactions in the **Protare**. */
        Reaction *orphanProduct( std::size_t a_index ) { return( m_orphanProducts.get<Reaction>( a_index ) ); }     /**< Returns the *a_index* - 1 orphan product reaction. */
        Reaction const *orphanProduct( std::size_t a_index ) const { return( m_orphanProducts.get<Reaction>( a_index ) ); }     /**< Returns the *a_index* - 1 orphan product reaction. */

        std::size_t numberOfIncompleteReactions( ) const { return( m_incompleteReactions.size( ) ); }                                   /**< Returns the number of incomplete reactions in the **Protare**. */
        Reaction *incompleteReaction( std::size_t a_index ) { return( m_incompleteReactions.get<Reaction>( a_index ) ); }               /**< Returns the *a_index* - 1 reaction. */
        Reaction const *incompleteReaction( std::size_t a_index ) const { return( m_incompleteReactions.get<Reaction>( a_index ) ); }   /**< Returns the *a_index* - 1 reaction. */
        void updateReactionIndices( int a_offset ) const;

        bool hasFission( ) const ;
        bool isDelayedFissionNeutronComplete( ) const ;

        GUPI::Ancestry *findInAncestry3( std::string const &a_item );
        GUPI::Ancestry const *findInAncestry3( std::string const &a_item ) const ;

        std::vector<double> groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo ) const ;

        Vector multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic = true,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                        std::string const &a_productID, int a_order, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Matrix multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order, 
                        TransportCorrectionType a_transportCorrectionType, double a_temperature,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        stringAndDoublePairs muCutoffForCoulombPlusNuclearElastic( ) const ;
        DelayedNeutronProducts delayedNeutronProducts( ) const ;
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;

        void saveAs( std::string const &a_fileName ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
===================== ProtareComposite =====================
============================================================
*/
class ProtareComposite : public Protare {

    private:
        std::vector<Protare *> m_protares;                                      /**< List of protares added to *this* instance. */

    public:
        ProtareComposite( Construction::Settings const &a_construction );
        ~ProtareComposite( );

        std::vector<Protare *> &protares( ) { return( m_protares ); }           /**< Returns the value of the *m_protares* member. */
        void append( Protare *a_protare );

// The rest are virtual methods defined in the Protare class.

        ProtareType protareType( ) const { return( ProtareType::composite ); }  /**< Returns the type of the protare. */
        std::size_t numberOfProtares( ) const ;
        ProtareSingle *protare( std::size_t a_index );
        ProtareSingle const *protare( std::size_t a_index ) const ;

        LUPI::FormatVersion const &formatVersion( std::size_t a_index = 0 ) const ;
        std::string const &fileName( std::size_t a_index = 0 ) const ;
        std::string const &realFileName( std::size_t a_index = 0 ) const ;

        std::vector<std::string> libraries( std::size_t a_index = 0 ) const ;
        std::string const &evaluation( std::size_t a_index = 0 ) const ;
        Frame projectileFrame( std::size_t a_index = 0 ) const ;
        int numberOfLazyParsingHelperForms( ) const ;
        int numberOfLazyParsingHelperFormsReplaced( ) const ;
        double thresholdFactor( ) const ;

        Documentation_1_10::Suite &documentations( );

        Styles::Base &style( std::string const a_label );
        Styles::Suite &styles( );
        Styles::Suite const &styles( ) const ;

        int intid( std::string const &a_id ) const;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        Styles::TemperatureInfos temperatures( ) const ;

        std::size_t numberOfReactions( ) const ;
        Reaction *reaction( std::size_t a_index );
        Reaction const *reaction( std::size_t a_index ) const ;
        Reaction const *reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const;
        std::size_t numberOfOrphanProducts( ) const ;
        Reaction *orphanProduct( std::size_t a_index );
        Reaction const *orphanProduct( std::size_t a_index ) const ;
        void updateReactionIndices( int a_offset ) const;

        bool hasFission( ) const ;
        bool isDelayedFissionNeutronComplete( ) const ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }  /**< Always returns *nullptr*. */
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }  /**< Always returns *nullptr*. */

        std::vector<double> groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo ) const ;

        Vector multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, 
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic = true,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                        std::string const &a_productID, int a_order, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Matrix multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order, 
                        TransportCorrectionType a_transportCorrectionType, double a_temperature,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        stringAndDoublePairs muCutoffForCoulombPlusNuclearElastic( ) const ;
        DelayedNeutronProducts delayedNeutronProducts( ) const ;
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
};

/*
============================================================
======================= ProtareTNSL ========================
============================================================
*/
class ProtareTNSL : public Protare {

    private:
        ProtareSingle *m_protare;                                           /**< Protare with non thermal neutron scattering law data. */
        ProtareSingle *m_TNSL;                                              /**< Protare with thermal neutron scattering law data. */
        Reaction *m_elasticReaction;                                        /**< The elastic reaction from the non TNSL protare. */
        std::map<std::string,std::size_t> m_maximumTNSL_MultiGroupIndex;    /**< For each neutron multi-group data, this the number of valid groups for the TNSL data. */

    public:
        ProtareTNSL( Construction::Settings const &a_construction, ProtareSingle *a_protare, ProtareSingle *a_TNSL );
        ~ProtareTNSL( );

        ProtareSingle *TNSL( ) { return( m_TNSL ); }                        /**< Returns the *m_TNSL* member. */
        ProtareSingle const *TNSL( ) const { return( m_TNSL ); }            /**< Returns the *m_TNSL* member. */
        Reaction *elasticReaction( ) { return( m_elasticReaction ); }       /**< Returns the *m_elasticReaction* member. */
        std::size_t maximumTNSL_MultiGroupIndex( Styles::TemperatureInfo const &a_temperatureInfo ) const ;
        void combineVectors( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, Vector &a_vector, 
                        Vector const &a_vectorElastic, Vector const &a_vectorTNSL ) const ;
        void combineMatrices( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, Matrix &a_matrix, 
                        Matrix const &a_matrixElastic, Matrix const &a_matrixTNSL ) const ;

// The rest are virtual methods defined in the Protare class.

        ProtareType protareType( ) const { return( ProtareType::TNSL ); }   /**< Returns the type of the protare. */
        std::size_t numberOfProtares( ) const { return( 2 ); }              /**< Always returns 2. */
        ProtareSingle *protare( std::size_t a_index = 0 );
        ProtareSingle const *protare( std::size_t a_index = 0 ) const ;

        LUPI::FormatVersion const &formatVersion( std::size_t a_index = 0 ) const ;
        std::string const &fileName( std::size_t a_index = 0 ) const ;
        std::string const &realFileName( std::size_t a_index = 0 ) const ;

        std::vector<std::string> libraries( std::size_t a_index = 0 ) const ;
        std::string const &evaluation( std::size_t a_index = 0 ) const ;
        Frame projectileFrame( std::size_t a_index = 0 ) const ;
        int numberOfLazyParsingHelperForms( ) const ;
        int numberOfLazyParsingHelperFormsReplaced( ) const ;
        double thresholdFactor( ) const ;

        Documentation_1_10::Suite &documentations( );

        Styles::Base &style( std::string const a_label );
        Styles::Suite &styles( );
        Styles::Suite const &styles( ) const ;

        int intid( std::string const &a_id ) const;
        void productIDs( std::set<std::string> &a_ids, Transporting::Particles const &a_particles, bool a_transportablesOnly ) const ;
        int maximumLegendreOrder( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID ) const ;

        Styles::TemperatureInfos temperatures( ) const ;

        std::size_t numberOfReactions( ) const ;
        Reaction *reaction( std::size_t a_index );
        Reaction const *reaction( std::size_t a_index ) const ;
        Reaction const *reaction( std::size_t a_index, Transporting::MG const &a_settings, 
                    ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const;
        std::size_t numberOfOrphanProducts( ) const ;
        Reaction *orphanProduct( std::size_t a_index );
        Reaction const *orphanProduct( std::size_t a_index ) const ;
        void updateReactionIndices( int a_offset ) const;

        bool hasFission( ) const ;
        bool isDelayedFissionNeutronComplete( ) const ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }                      /**< Always returns *nullptr*. */
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }          /**< Always returns *nullptr*. */

        std::vector<double> groupBoundaries( Transporting::MG const &a_settings, Styles::TemperatureInfo const &a_temperatureInfo, 
                        std::string const &a_productID ) const ;
        Vector multiGroupInverseSpeed( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo ) const ;

        Vector multiGroupCrossSection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, bool a_final, bool a_effectivePhotoAtomic = true,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionNeutronMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupFissionGammaMultiplicity( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Matrix multiGroupProductMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, 
                        std::string const &a_productID, int a_order, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Matrix multiGroupFissionMatrix( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupTransportCorrection( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles, int a_order, 
                        TransportCorrectionType a_transportCorrectionType, double a_temperature,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionEnergy( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupAvailableMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupAverageMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;
        Vector multiGroupDepositionMomentum( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, Transporting::Particles const &a_particles,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        Vector multiGroupGain( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo, std::string const &a_productID,
                        ExcludeReactionsSet const &a_reactionsToExclude = ExcludeReactionsSet {} ) const ;

        void TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::XYs1d &a_crossSectionSum );
        void TNSL_crossSectionSumCorrection( std::string const &a_label, Functions::Ys1d &a_crossSectionSum );
        void TNSL_crossSectionSumCorrection( std::string const &a_label, Vector &a_crossSectionSum ) {
            return( Protare::TNSL_crossSectionSumCorrection( a_label, a_crossSectionSum ) );
        }

        stringAndDoublePairs muCutoffForCoulombPlusNuclearElastic( ) const ;
        DelayedNeutronProducts delayedNeutronProducts( ) const { return( m_protare->delayedNeutronProducts( ) ); }
        void incompleteParticles( Transporting::Settings const &a_settings, std::set<std::string> &a_incompleteParticles ) const ;
};

namespace Map {

enum class EntryType { import, protare, TNSL };
#define GIDI_MapInteractionNuclearChars "nuclear"
#define GIDI_MapInteractionAtomicChars "atomic"
#define GIDI_MapInteractionTNSLChars "thermalNeutronScatteringLaw"

typedef std::vector<ProtareBase const *> FindProtareEntries;

/*
============================================================
========================= BaseEntry ========================
============================================================
*/
class BaseEntry : public GUPI::Ancestry {

    public:
        enum class PathForm { entered, cumulative, real };

    private:
        std::string m_name;                                 /**< Designates the entry as either a protare or a map. */
        Map const *m_parent;                                /**< Pointer to map containing *this*. */
        std::string m_path;                                 /**< Absolute or relative (to map file) path of the protare or map file. */
        std::string m_cumulativePath;                       /**< Currently not used. */

    public:
        BaseEntry( HAPI::Node const &a_node, std::string const &a_basePath, Map const *a_parent );
        virtual ~BaseEntry( ) = 0;

        std::string const &name( ) const { return( m_name ); }              /**< Returns the value of the *m_name* member. */
        Map const *parent( ) const { return( m_parent ); }                  /**< Returns the value of the *m_parent* member. */
        std::string path( PathForm a_form = PathForm::real ) const ;

        virtual EntryType entryType( ) const = 0;

        void libraries( std::vector<std::string> &a_libraries ) const ;
        virtual ProtareBase const *findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID, 
                        std::string const &a_library = "", std::string const &a_evaluation = "" ) const = 0 ;
        virtual void findProtareEntries( FindProtareEntries &a_protareEntries, std::regex const &a_projectileID,
                        std::regex const &a_targetID, std::regex const &a_library = std::regex( ".*" ), 
                        std::regex const &a_evaluation = std::regex( ".*" ) ) const = 0 ;

        virtual void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const = 0;
};

/*
============================================================
========================== Import ==========================
============================================================
*/
class Import : public BaseEntry {

    private:
        Map *m_map;                                         /**< Map instance for this Import. */

    public:
        Import( HAPI::Node const &a_node, PoPI::Database const &a_pops, std::string const &a_basePath, Map const *a_parent );
        ~Import( );

        EntryType entryType( ) const { return( EntryType::import ); }   /**< Returns EntryType::import. */

        Map const *map( ) const { return( m_map ); }                    /**< Returns the value of the *m_map* member. */

        ProtareBase const *findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID,
                std::string const &a_library = "", std::string const &a_evaluation = "" ) const ;
        void findProtareEntries( FindProtareEntries &a_protareEntries, std::regex const &a_projectileID,
                std::regex const &a_targetID, std::regex const &a_library = std::regex( ".*" ), std::regex const &a_evaluation = std::regex( ".*" ) ) const ;
        std::string protareFilename( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library = "",
                std::string const &a_evaluation = "", PathForm a_form = PathForm::real ) const ;
        bool isProtareAvailable( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library = "",
                std::string const &a_evaluation = "" ) const {
            return( protareFilename( a_projectileID, a_targetID, a_library, a_evaluation ) != GIDI_emptyFileNameChars ); }
                                                                        /**< Returns the value of the *m_map* member. */
        std::vector<std::string> availableEvaluations( std::string const &a_projectileID, std::string const &a_targetID ) const ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }                  /**< Always returns *nullptr*. */
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }      /**< Always returns *nullptr*. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
======================= ProtareBase ========================
============================================================
*/
class ProtareBase : public BaseEntry {

    private:
        std::string m_projectileID;             /**< Projectile id for protare. */
        std::string m_targetID;                 /**< Target id for protare. */
        std::string m_evaluation;               /**< Evaluation string for protare. */
        std::string m_interaction;              /**< The interaction type for the protare. */

    public:
        ProtareBase( HAPI::Node const &a_node, std::string const &a_basePath, Map const *const a_map );
        ~ProtareBase( );

        std::string const &projectileID( ) const { return( m_projectileID ); }      /**< Returns the value of the *m_projectileID* member. */
        std::string const &targetID( ) const { return( m_targetID ); }              /**< Returns the value of the *m_targetID* member. */
        std::string const &evaluation( ) const { return( m_evaluation ); }          /**< Returns the value of the *m_evaluation* member. */
        std::string const &interaction( ) const { return( m_interaction ); }        /**< Returns the value of the *m_interaction* member. */
        void setInteraction( std::string const &a_interaction ) { m_interaction = a_interaction; }  /**< Set the *m_interaction* member to *a_interaction*. */

        bool isMatch( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_evaluation = "" ) const ;
        std::string const &library( ) const ;
        std::string const &resolvedLibrary( ) const ;

        ProtareBase const *findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID,
                std::string const &a_library = "", std::string const &a_evaluation = "" ) const ;
        void findProtareEntries( FindProtareEntries &a_protareEntries, std::regex const &a_projectileID,
                std::regex const &a_targetID, std::regex const &a_library = std::regex( ".*" ), std::regex const &a_evaluation = std::regex( ".*" ) ) const ;
        virtual GIDI::Protare *protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, 
                        ParticleSubstitution const &a_particleSubstitution ) const = 0 ;
        virtual GIDI::ProtareSingle *protareSingle( Construction::Settings const &a_construction, PoPI::Database const &a_pops, 
                        ParticleSubstitution const &a_particleSubstitution ) const = 0 ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }                  /**< Always returns *nullptr*. */
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }      /**< Always returns *nullptr*. */
};

/*
============================================================
========================= Protare ==========================
============================================================
*/
class Protare : public ProtareBase {

    private:
        bool m_isPhotoAtomic;                   /**< true if photo-atomic protare and false otherwise. */

    public:
        Protare( HAPI::Node const &a_node, PoPI::Database const &a_pops, std::string const &a_basePath, Map const *const a_parent );
        ~Protare( );

        EntryType entryType( ) const { return( EntryType::protare ); }              /**< Returns EntryType::protare. */

        bool isPhotoAtomic( ) const { return( m_isPhotoAtomic ); }                  /**< Returns the value of the *m_isPhotoAtomic* member. */
        GIDI::Protare *protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution ) const ;
        GIDI::ProtareSingle *protareSingle( Construction::Settings const &a_construction, PoPI::Database const &a_pops, 
                ParticleSubstitution const &a_particleSubstitution ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== TNSL ===========================
============================================================
*/
class TNSL : public ProtareBase {

    private:
        std::string m_standardTarget;                           /**< The non-TNSL target. */
        std::string m_standardEvaluation;                       /**< The non-TNSL evaluation. */

    public:
        TNSL( HAPI::Node const &a_node, PoPI::Database const &a_pops, std::string const &a_basePath, Map const *const a_parent );
        ~TNSL( );

        EntryType entryType( ) const { return( EntryType::TNSL ); }                 /**< Returns EntryType::TNSL. */

        std::string const &standardTarget( ) const { return( m_standardTarget ); }
        std::string const &standardEvaluation( ) const { return( m_standardEvaluation ); }
        GIDI::Protare *protare( Construction::Settings const &a_construction, PoPI::Database const &a_pops, ParticleSubstitution const &a_particleSubstitution ) const ;
        GIDI::ProtareSingle *protareSingle( Construction::Settings const &a_construction, PoPI::Database const &a_pops, 
                ParticleSubstitution const &a_particleSubstitution ) const ;

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

/*
============================================================
=========================== Map ============================
============================================================
*/
class Map : public GUPI::Ancestry {

    private:
        Map const *m_parent;                            /**< Pointer to map containing *this* if this is an imported map. */
        std::string m_fileName;                         /**< Specified path to Map file. */
        std::string m_realFileName;                     /**< Absolute, real path to Map file. */
        std::string m_library;                          /**< The name of the library. */
        std::vector<BaseEntry *> m_entries;             /**< List of Map entries. */
        RISI::Projectiles m_projectiles;                /**< **RISI::Projectiles** loaded when method RIS_load is called. */
        bool m_projectilesLoaded;                       /**< If **true** data for **m_projectiles** have been read in, otherwise they have not been read in. */

        void initialize( std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent );
        void initialize( HAPI::Node const &a_node, std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent );

    public:
        Map( std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent = nullptr );
        Map( HAPI::Node const &a_node, std::string const &a_fileName, PoPI::Database const &a_pops, Map const *a_parent = nullptr );
        ~Map( );

        Map const *parent( ) const { return( m_parent ); }                      /**< Returns the value of the *m_parent* member. */
        std::string const &fileName( ) const { return( m_fileName ); }          /**< Returns the value of the *m_fileName* member. */
        std::string const &realFileName( ) const { return( m_realFileName ); }  /**< Returns the value of the *m_realFileName* member. */

        std::string const &library( ) const { return( m_library ); }            /**< Returns the value of the *m_library* member. */
        std::string const &resolvedLibrary( ) const ;
        void libraries( std::vector<std::string> &a_libraries ) const ;

        std::size_t size( ) const { return( m_entries.size( ) ); }              /**< Returns the number of entries in *this*. Does not descend map entries. */
        BaseEntry const *operator[]( std::size_t a_index ) const { return( m_entries[a_index] ); }
                                                                                /**< Returns the map entry at index *a_index*. */

        ProtareBase const *findProtareEntry( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library = "",
                std::string const &a_evaluation = "" ) const ;
        void findProtareEntries( FindProtareEntries &a_protareEntries, std::regex const &a_projectileID,
                std::regex const &a_targetID, std::regex const &a_library = std::regex( ".*" ), std::regex const &a_evaluation = std::regex( ".*" ) ) const ;
        std::string protareFilename( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library = "",
                std::string const &a_evaluation = "", BaseEntry::PathForm a_form = BaseEntry::PathForm::real ) const ;

        bool isProtareAvailable( std::string const &a_projectileID, std::string const &a_targetID, std::string const &a_library = "",
                std::string const &a_evaluation = "" ) const {
            return( protareFilename( a_projectileID, a_targetID, a_library, a_evaluation, BaseEntry::PathForm::entered ) != GIDI_emptyFileNameChars ); }
                                /**< Returns true if the map contains a Protare matching *a_projectileID*, *a_targetID*, *a_library* and *a_evaluation*, and false otherwise. */
        bool isTNSL_target( std::string const &a_targetID ) const ;
        std::vector<std::string> availableEvaluations( std::string const &a_projectileID, std::string const &a_targetID ) const ;

        GIDI::Protare *protare(             Construction::Settings const &a_construction, PoPI::Database const &a_pops, std::string const &a_projectileID, 
                std::string const &a_targetID, std::string const &a_library = "", std::string const &a_evaluation = "", 
                bool a_targetRequiredInGlobalPoPs = true, bool a_requiredInPoPs = true ) const ;

        std::vector<ProtareBase const *> directory( std::string const &a_projectileID = "", std::string const &a_targetID = "", 
                std::string const &a_library = "", std::string const &a_evaluation = "" ) const ;
        bool walk( MapWalkCallBack a_mapWalkCallBack, void *a_userData, int a_level = 0 ) const ;

        GUPI::Ancestry *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) { return( nullptr ); }                  /**< Always returns *nullptr*. */
        GUPI::Ancestry const *findInAncestry3( LUPI_maybeUnused std::string const &a_item ) const { return( nullptr ); }      /**< Always returns *nullptr*. */

        void saveAs( std::string const &a_fileName ) const ;
        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;

        std::string RIS_fileName( );
        bool RIS_fileExist( );
        RISI::Projectiles const &RIS_load( std::string const &a_energyUnit );
        std::string replacementTarget( PoPI::Database const &a_pops, std::string const &a_projectile, std::string const &a_target );
};

}           // End of namespace Map.

namespace Functions {
/*
============================================================
================= FissionEnergyRelease ==================
============================================================
*/
class FissionEnergyRelease : public Function1dForm {

    private:
        Function1dForm *m_promptProductKE;                  /**< The **ENDF** prompt total product kinetic energy released. */
        Function1dForm *m_promptNeutronKE;                  /**< The **ENDF** prompt neutron kinetic energy released. */
        Function1dForm *m_delayedNeutronKE;                 /**< The **ENDF** delayed neutron kinetic energy released. */
        Function1dForm *m_promptGammaEnergy;                /**< The **ENDF** prompt gamma energy released. */
        Function1dForm *m_delayedGammaEnergy;               /**< The **ENDF** delayed gamma energy released. */
        Function1dForm *m_delayedBetaEnergy;                /**< The **ENDF** delayed beta kinetic energy released. */
        Function1dForm *m_neutrinoEnergy;                   /**< The **ENDF** neutrino energy released. */
        Function1dForm *m_nonNeutrinoEnergy;                /**< The **ENDF** non neutrino energy released. */
        Function1dForm *m_totalEnergy;                      /**< The **ENDF** total energy released. */

        void energyReleaseToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_moniker, std::string const &a_indent, Function1dForm *a_function1d ) const ;

    public:
        FissionEnergyRelease( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
        ~FissionEnergyRelease( );

        double domainMin( ) const { return( m_nonNeutrinoEnergy->domainMin( ) ); }                  /**< Returns the minimum domain value for the energy released. */
        double domainMax( ) const { return( m_nonNeutrinoEnergy->domainMax( ) ); }                  /**< Returns the maximum domain value for the energy released. */
        double evaluate( double a_x1 ) const { return( m_nonNeutrinoEnergy->evaluate( a_x1 ) ); }   /**< Returns the value of *m_nonNeutrinoEnergy* evaluated at *a_x1*. */
        Vector multiGroupQ( LUPI::StatusMessageReporting &a_smr, Transporting::MG const &a_settings, 
                        Styles::TemperatureInfo const &a_temperatureInfo ) const ;

        Function1dForm const *promptProductKE( ) const { return( m_promptProductKE ); }             /**< Returns the value of the *m_promptProductKE* member. */
        Function1dForm       *promptProductKE( )       { return( m_promptProductKE ); }             /**< Returns the value of the *m_promptProductKE* member. */
        Function1dForm const *promptNeutronKE( ) const { return( m_promptNeutronKE ); }             /**< Returns the value of the *m_promptNeutronKE* member. */
        Function1dForm       *promptNeutronKE( )       { return( m_promptNeutronKE ); }             /**< Returns the value of the *m_promptNeutronKE* member. */
        Function1dForm const *delayedNeutronKE( ) const { return( m_delayedNeutronKE ); }           /**< Returns the value of the *m_delayedNeutronKE* member. */
        Function1dForm       *delayedNeutronKE( )       { return( m_delayedNeutronKE ); }           /**< Returns the value of the *m_delayedNeutronKE* member. */
        Function1dForm const *promptGammaEnergy( ) const { return( m_promptGammaEnergy ); }         /**< Returns the value of the *m_promptGammaEnergy* member. */
        Function1dForm       *promptGammaEnergy( )       { return( m_promptGammaEnergy ); }         /**< Returns the value of the *m_promptGammaEnergy* member. */
        Function1dForm const *delayedGammaEnergy( ) const { return( m_delayedGammaEnergy ); }       /**< Returns the value of the *m_delayedGammaEnergy* member. */
        Function1dForm       *delayedGammaEnergy( )       { return( m_delayedGammaEnergy ); }       /**< Returns the value of the *m_delayedGammaEnergy* member. */
        Function1dForm const *delayedBetaEnergy( ) const { return( m_delayedBetaEnergy ); }         /**< Returns the value of the *m_delayedBetaEnergy* member. */
        Function1dForm       *delayedBetaEnergy( )       { return( m_delayedBetaEnergy ); }         /**< Returns the value of the *m_delayedBetaEnergy* member. */
        Function1dForm const *neutrinoEnergy( ) const { return( m_neutrinoEnergy ); }               /**< Returns the value of the *m_neutrinoEnergy* member. */
        Function1dForm       *neutrinoEnergy( )       { return( m_neutrinoEnergy ); }               /**< Returns the value of the *m_neutrinoEnergy* member. */
        Function1dForm const *nonNeutrinoEnergy( ) const { return( m_nonNeutrinoEnergy ); }         /**< Returns the value of the *m_neutrinoEnergy* member. */
        Function1dForm       *nonNeutrinoEnergy( )       { return( m_nonNeutrinoEnergy ); }         /**< Returns the value of the *m_neutrinoEnergy* member. */
        Function1dForm const *totalEnergy( ) const { return( m_totalEnergy ); }                     /**< Returns the value of the *m_totalEnergy* member. */
        Function1dForm       *totalEnergy( )       { return( m_totalEnergy ); }                     /**< Returns the value of the *m_totalEnergy* member. */

        void toXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent = "" ) const ;
};

}           // End of namespace Functions.

/*
============================================================
========================== Groups ==========================
============================================================
*/
class Groups : public Suite {

    public:
        Groups( );
        Groups( std::string const &a_fileName );

        void addFile( std::string const &a_fileName );
};

/*
============================================================
========================== Fluxes ==========================
============================================================
*/
class Fluxes : public Suite {

    public:
        Fluxes( );
        Fluxes( std::string const &a_fileName );

        void addFile( std::string const &a_fileName );
};

/*
============================================================
============= MultiGroupCalulationInformation ==============
============================================================
*/

class MultiGroupCalulationInformation {

    public:
        Transporting::MultiGroup const &m_multiGroup;       /**< The multi-group boundaries. */
        Transporting::Flux const &m_flux;                   /**< The flux weighting. */
        ptwXPoints *m_boundaries_xs;                        /**< This is an **ptwXPoints** representation of *m_heatedMultiGroupLabel* as needed by numerical functions. */
        ptwXYPoints *m_fluxes_xys;                          /**< This is an **ptwXYPoints** representation of *m_flux* as needed by numerical functions. */
        ptwXPoints *m_multiGroupFlux;                       /**< This is the grouped representation of *m_flux* as needed by numerical functions. */

        MultiGroupCalulationInformation( Transporting::MultiGroup const &a_multiGroup, Transporting::Flux const &a_flux );
        ~MultiGroupCalulationInformation( );
};

/*
============================================================
========================== others ==========================
============================================================
*/
Form *parseExternalFilesSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseStylesSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseTransportablesSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseReaction( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseOrphanProduct( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseFissionComponent( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );
Form *parseReactionType( std::string const &a_moniker, Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseSumsCrossSectionsSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );
Form *parseSumsMultiplicitiesSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops,
                PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );
Form *parseDoubleDifferentialCrossSectionSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseScatteringAtom( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseCrossSectionSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseDelayedNeutronsSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseFissionEnergyReleasesSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parsePhysicalQuantitySuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseAvailableSuite( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseQSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseProductSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs, 
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseMultiplicitySuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseDistributionSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseAverageEnergySuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseAverageMomentumSuite( Construction::Settings const &a_construction, Suite *parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo, PoPI::Database const &a_pop, PoPI::Database const &a_internalPoPs,
                std::string const &a_name, Styles::Suite const *a_styles );
Form *parseACE_URR_probabilityTables( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );
Form *parseColumnHeaders( Construction::Settings const &a_construction, Suite *a_parent, HAPI::Node const &a_node, SetupInfo &a_setupInfo,
                PoPI::Database const &a_pops, PoPI::Database const &a_internalPoPs, std::string const &a_name, Styles::Suite const *a_styles );
Functions::Function1dForm *data1dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *parent );
Functions::Function1dForm *data1dParseAllowEmpty( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *a_parent );
void data1dListParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, std::vector<Functions::Function1dForm *> &a_function1ds );
Functions::Function2dForm *data2dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *parent );
void data2dListParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, std::vector<Functions::Function2dForm *> &a_function2ds );
Functions::Function3dForm *data3dParse( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Suite *parent );
void checkOuterDomainValues1d( std::vector<Functions::Function1dForm *> &a_functions, std::vector<double> &a_Xs );
void checkOuterDomainValues2d( std::vector<Functions::Function2dForm *> &a_functions, std::vector<double> &a_Xs );
void checkSequentialDomainLimits1d( std::vector<Functions::Function1dForm *> &a_functions, std::vector<double> &a_Xs );
void checkSequentialDomainLimits2d( std::vector<Functions::Function2dForm *> &a_functions, std::vector<double> &a_Xs );

int parseFlattened1d( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, Vector &data );

Vector collapse( Vector const &a_vector, Transporting::Settings const &a_settings, Transporting::Particles const &a_particles, double a_temperature );
Matrix collapse( Matrix const &a_matrix, Transporting::Settings const &a_settings, Transporting::Particles const &a_particles, double a_temperature, 
                std::string const &a_productID );

Vector transportCorrect( Vector const &a_vector, Vector const &a_transportCorrection );
Matrix transportCorrect( Matrix const &a_matrix, Vector const &a_transportCorrection );

Vector multiGroupXYs1d( Transporting::MultiGroup const &a_boundaries, Functions::XYs1d const &a_function, Transporting::Flux const &a_flux );
Vector *multiGroupTwoXYs1ds( MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Functions::XYs1d const &a_function1,
                Functions::XYs1d const &a_function2 );
void calculate1dMultiGroupDataInComponent( ProtareSingle const *a_protare, std::string const &a_heatedMultiGroupLabel,
                MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Component &a_component, Functions::XYs1d const &a_crossSection );
void calculate1dMultiGroupFissionEnergyRelease( MultiGroupCalulationInformation const &a_multiGroupCalulationInformation, Functions::XYs1d const &a_weight,
                Functions::Function1dForm const *a_evaluated, Functions::Function1dForm *a_gridded1d );


int ENDL_CFromENDF_MT( int ENDF_MT, int *ENDL_C, int *ENDL_S );

GNDS_FileType GNDS_fileType( std::string const &a_fileName, GNDS_FileTypeInfo &a_GNDS_fileTypeInfo );

/*
*   The following are in the file GIDI_misc.cpp.
*/
long binarySearchVector( double a_x, std::vector<double> const &a_Xs );
void intsToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, std::vector<int> a_values, std::string const &a_attributes );
void parseValuesOfDoubles( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<double> &a_vector );
void parseValuesOfDoubles( HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<double> &a_vector, int a_useSystem_strtod );
void parseValuesOfInts( Construction::Settings const &a_construction, HAPI::Node const &a_node, SetupInfo &a_setupInfo, std::vector<int> &a_vector );
void parseValuesOfInts( HAPI::Node const &a_node, SetupInfo &a_setupInfo, nf_Buffer<int> &a_vector );
void doublesToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_indent, std::vector<double> a_values, std::size_t a_start = 0, bool a_newLine = true,
        std::string const &a_valueType = "" );
Frame parseFrame( HAPI::Node const &a_node, SetupInfo &a_setupInfo, std::string const &a_name );
std::string frameToString( Frame a_frame );
std::string intToString( int a_value );
std::string size_t_ToString( std::size_t a_value );
std::string nodeWithValuesToDoubles( GUPI::WriteInfo &a_writeInfo, std::string const &a_nodeName, std::vector<double> const &a_values );
void excludeReactionsSetAdjust( ExcludeReactionsSet a_excludeReactionsSet, Protare const &a_protare );

Functions::Ys1d gridded1d2GIDI_Ys1d( Functions::Function1dForm const &a_function1d );
Functions::Ys1d vector2GIDI_Ys1d( Axes const &a_axes, Vector const &a_vector );

std::string LLNL_gidToLabel( int a_gid );
std::string LLNL_fidToLabel( int a_fid );

std::vector<std::string> sortedListOfStrings( std::vector<std::string> const &a_strings, bool a_orderIsAscending = true );

void energy2dToXMLList( GUPI::WriteInfo &a_writeInfo, std::string const &a_moniker, std::string const &a_indent, Functions::Function1dForm *a_function );

std::vector<Transporting::Flux> settingsFluxesFromFunction3d( Functions::Function3dForm const &a_function3d );

}           // End of namespace GIDI.

#endif      // End of GIDI_hpp_included
