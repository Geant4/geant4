#------------------------------------------------------------------------------
# sources.cmake
# Module : G4materials
# Package: Geant4.src.G4materials
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.3 2010-11-29 17:15:46 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4materials
    HEADERS
        G4Isotope.hh                   
        G4NistManager.hh
        G4AtomicShells.hh       
        G4IsotopeVector.hh             
        G4NistMaterialBuilder.hh
        G4DensityEffectData.hh  
        G4Material.hh                  
        G4NistMessenger.hh
        G4Element.hh            
        G4MaterialPropertiesTable.hh   
        G4OpticalSurface.hh
        G4ElementTable.hh       
        G4MaterialPropertiesTable.icc  
        G4SandiaTable.hh
        G4ElementVector.hh      
        G4MaterialPropertyVector.hh    
        G4StaticSandiaData.hh
        G4ExtDEDXTable.hh       
        G4MaterialPropertyVector.icc   
        G4SurfaceProperty.hh
        G4IonisParamElm.hh      
        G4MaterialTable.hh             
        G4VIonDEDXTable.hh
        G4IonisParamMat.hh      
        G4MPVEntry.hh
        G4IonStoppingData.hh    
        G4NistElementBuilder.hh
    SOURCES
        G4IonStoppingData.cc
        G4NistManager.cc
        G4AtomicShells.cc       
        G4Isotope.cc
        G4NistMaterialBuilder.cc
        G4DensityEffectData.cc  
        G4Material.cc              
        G4NistMessenger.cc
        G4Element.cc          
        G4MaterialPropertiesTable.cc  
        G4OpticalSurface.cc
        G4ExtDEDXTable.cc       
        G4MaterialPropertyVector.cc   
        G4SandiaTable.cc
        G4IonisParamElm.cc      
        G4MPVEntry.cc                 
        G4SurfaceProperty.cc
        G4IonisParamMat.cc      
        G4NistElementBuilder.cc       
        G4VIonDEDXTable.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

