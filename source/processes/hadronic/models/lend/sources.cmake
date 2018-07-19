#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_lend
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_lend
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106463 2017-10-11 08:05:17Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${EXPAT_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/pre_equilibrium/exciton_model/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)


#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4had_lend
    HEADERS
        G4GIDI.hh
        G4GIDI_Misc.hh
        G4GIDI_map.hh
        G4GIDI_mass.hh
        G4GIDI_target.hh
        G4LENDCapture.hh
        G4LENDCaptureCrossSection.hh
        G4LENDCrossSection.hh
        G4LENDCombinedCrossSection.hh
        G4LENDCombinedModel.hh
        G4LENDElastic.hh
        G4LENDElasticCrossSection.hh
        G4LENDFission.hh
        G4LENDFissionCrossSection.hh
        G4LENDHeader.hh
        G4LENDInelastic.hh
        G4LENDInelasticCrossSection.hh
        G4LENDManager.hh
        G4LENDModel.hh
        G4LENDUsedTarget.hh
        GIDI_settings.hh
        MCGIDI.h
        MCGIDI_fromTOM.h
        MCGIDI_map.h
        MCGIDI_mass.h
        MCGIDI_misc.h
        MCGIDI_private.h
        PoPs.h
        PoPs_Bcast_private.h
        PoPs_data.h
        PoPs_mass.h
        PoPs_private.h
        nf_Legendre.h
        nf_integration.h
        nf_specialFunctions.h
        nf_utilities.h
        ptwX.h
        ptwXY.h
        statusMessageReporting.h
        xDataTOM.h
        xDataTOM_importXML_private.h
        xDataTOM_private.h
    SOURCES
        G4GIDI.cc
        G4GIDI_Misc.cc
        G4GIDI_map.cc
        G4GIDI_mass.cc
        G4GIDI_target.cc
        G4LENDCapture.cc
        G4LENDCaptureCrossSection.cc
        G4LENDCombinedCrossSection.cc
        G4LENDCombinedModel.cc
        G4LENDCrossSection.cc
        G4LENDElastic.cc
        G4LENDElasticCrossSection.cc
        G4LENDFission.cc
        G4LENDFissionCrossSection.cc
        G4LENDInelastic.cc
        G4LENDInelasticCrossSection.cc
        G4LENDManager.cc
        G4LENDModel.cc
        G4LENDUsedTarget.cc
        GIDI_settings.cc
        GIDI_settings_flux.cc
        GIDI_settings_group.cc
        GIDI_settings_particle.cc
        MCGIDI_KalbachMann.cc
        MCGIDI_LLNLAngular_angularEnergy.cc
        MCGIDI_angular.cc
        MCGIDI_angularEnergy.cc
        MCGIDI_distribution.cc
        MCGIDI_energy.cc
        MCGIDI_energyAngular.cc
        MCGIDI_fromTOM.cc
        MCGIDI_kinetics.cc
        MCGIDI_map.cc
        MCGIDI_mass.cc
        MCGIDI_misc.cc
        MCGIDI_outputChannel.cc
        MCGIDI_particle.cc
        MCGIDI_pop.cc
        MCGIDI_product.cc
        MCGIDI_quantitiesLookupMode.cc
        MCGIDI_reaction.cc
        MCGIDI_sampling.cc
        MCGIDI_samplingSettings.cc
        MCGIDI_target.cc
        MCGIDI_target_heated.cc
        MCGIDI_uncorrelated.cc
        MCGIDI_version.cc
        PoPs.cc
        PoPs_Bcast.cc
        PoPs_data.cc
        PoPs_mass.cc
        lPoPs.cc
        nf_GnG_adaptiveQuadrature.cc
        nf_Legendre.cc
        nf_Legendre_GaussianQuadrature.cc
        nf_angularMomentumCoupling.cc
        nf_exponentialIntegral.cc
        nf_gammaFunctions.cc
        nf_incompleteGammaFunctions.cc
        nf_polevl.cc
        nf_stringToDoubles.cc
        nf_stringToDoubles_main.cc
        nf_utilities.cc
        ptwXY_binaryOperators.cc
        ptwXY_convenient.cc
        ptwXY_core.cc
        ptwXY_functions.cc
        ptwXY_integration.cc
        ptwXY_interpolation.cc
        ptwXY_methods.cc
        ptwXY_misc.cc
        ptwXY_unitaryOperators.cc
        ptwX_core.cc
        ptwX_misc.cc
        statusMessageReporting.cc
        xDataTOM.cc
        xDataTOM_KalbachMann.cc
        xDataTOM_LegendreSeries.cc
        xDataTOM_Misc.cc
        xDataTOM_V_W_XYs.cc
        xDataTOM_V_W_XYs_LegendreSeries.cc
        xDataTOM_W_XYs.cc
        xDataTOM_W_XYs_LegendreSeries.cc
        xDataTOM_XYs.cc
        xDataTOM_axes.cc
        xDataTOM_importXML.cc
        xDataTOM_importXML_KalbachMann.cc
        xDataTOM_importXML_V_W_XYs.cc
        xDataTOM_importXML_V_W_XYs_LegendreSeries.cc
        xDataTOM_importXML_W_XYs.cc
        xDataTOM_importXML_W_XYs_LegendreSeries.cc
        xDataTOM_importXML_XYs.cc
        xDataTOM_importXML_axes.cc
        xDataTOM_importXML_polynomial.cc
        xDataTOM_importXML_regionsW_XYs_LegendreSeries.cc
        xDataTOM_importXML_regionsXYs.cc
        xDataTOM_interpolation.cc
        xDataTOM_polynomial.cc
        xDataTOM_regionsW_XYs_LegendreSeries.cc
        xDataTOM_regionsXYs.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_LE
        G4hadronic_deex_management
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
        ${EXPAT_LIBRARIES}
)

# List any source specific properties here
#
# Not this again... Need to add *another* compile definition to lend sources
# for DLLs....
# We do this quick and dirty, add the definition directly. I *don't want it 
# polluting compile of all other processes modules
# in global mode
#
set_source_files_properties(
    ${G4had_lend_SOURCES}
    PROPERTIES COMPILE_DEFINITIONS G4PROCESSES_EXPORT
)


