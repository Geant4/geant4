#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phys_lists
# Package: Geant4.src.G4physicslists.G4physlist_ctors.G4physlist_ctor_ions
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 10/01/2013
#
# $Id: sources.cmake 97179 2016-05-27 12:46:40Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/shortlived/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/builders/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/physics_lists/constructors/factory/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/transportation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/decay/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/stopping/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/binary_cascade/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/cascade/cascade/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/qmd/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/coherent_elastic/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fermi_breakup/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fission/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/handler/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/multifragmentation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/im_r_matrix/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/inclxx/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/inclxx/incl_physics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/inclxx/interface/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/lend/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/lepto_nuclear/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/particle_hp/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/parton_string/diffraction/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/parton_string/hadronization/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/parton_string/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/parton_string/qgsm/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/pre_equilibrium/exciton_model/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/radioactive_decay/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/theo_high_energy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4phys_ctor_ions
    HEADERS
        G4IonBinaryCascadePhysics.hh
        G4IonINCLXXPhysics.hh
        G4IonPhysics.hh
        G4IonPhysicsPHP.hh
        G4IonQMDPhysics.hh
    SOURCES
        G4IonBinaryCascadePhysics.cc
        G4IonINCLXXPhysics.cc
        G4IonPhysics.cc
        G4IonPhysicsPHP.cc
        G4IonQMDPhysics.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4decay
        G4digits
        G4emhighenergy
        G4emlowenergy
        G4emstandard
        G4emutils
        G4event
        G4geometrymng
        G4globman
        G4had_im_r_matrix
        G4had_lend
        G4had_mod_man
        G4had_mod_util
        G4had_lept_nuclear
        G4had_neu_hp
        G4had_preequ_exciton
        G4had_string_diff
        G4had_string_frag
        G4had_string_man
        G4had_theo_max
        G4hadronic_HE
        G4hadronic_LE
        G4hadronic_bert_cascade
        G4hadronic_binary
        G4hadronic_coherent_elastic
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_fission
        G4hadronic_deex_gem_evaporation
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_inclxx_interface
        G4hadronic_inclxx_physics
        G4hadronic_inclxx_utils
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_qgstring
        G4hadronic_radioactivedecay
        G4hadronic_stop
        G4hadronic_util
        G4hadronic_xsect
        G4hits
        G4intercoms
        G4ions
        G4leptons
        G4magneticfield
        G4materials
        G4mesons
        G4muons
        G4navigation
        G4optical
        G4partman
        G4phys_builders
        G4phys_ctor_factory
        G4procman
        G4run
        G4shortlived
        G4track
        G4tracking
        G4transportation
        G4volumes
        G4xrays
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

