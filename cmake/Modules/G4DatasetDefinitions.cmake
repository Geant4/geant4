# - Define datasets known and used by Geant4
# We keep this separate from the Geant4InstallData module for conveniance
# when updating and patching because datasets may change more rapidly.
# It allows us to decouple the dataset definitions from how they are
# checked/installed/configured
#

# - NDL
geant4_add_dataset(
  NAME      G4NDL
  VERSION   4.7.1
  FILENAME  G4NDL
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONHPDATA
  MD5SUM    54f0ed3995856f02433d42ec96d70bc6
  )

# - Low energy electromagnetics
geant4_add_dataset(
  NAME      G4EMLOW
  VERSION   8.6.1
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    9db67a37acc3eae9b0ffdace41a23b74
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
  VERSION   6.1
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    92d68b937cdad0fd49892a66878863de 
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
  VERSION   6.1.2
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    20d494f73d4bddabd7fab5c06a58895c
  )

# - Particle XS - replaces Neutron XS
geant4_add_dataset(
  NAME      G4PARTICLEXS
  VERSION   4.1
  FILENAME  G4PARTICLEXS
  EXTENSION tar.gz
  ENVVAR    G4PARTICLEXSDATA
  MD5SUM    878252a464ba6b38f085741840f053e6
  )

# - PII
geant4_add_dataset(
  NAME      G4PII
  VERSION   1.3
  FILENAME  G4PII
  EXTENSION tar.gz
  ENVVAR    G4PIIDATA
  MD5SUM    05f2471dbcdf1a2b17cbff84e8e83b37
  )

# - Optical Surfaces
geant4_add_dataset(
  NAME      RealSurface
  VERSION   2.2
  FILENAME  G4RealSurface
  EXTENSION tar.gz
  ENVVAR    G4REALSURFACEDATA
  MD5SUM    ea8f1cfa8d8aafd64b71fb30b3e8a6d9
  )

# - SAID
geant4_add_dataset(
  NAME      G4SAIDDATA
  VERSION   2.0
  FILENAME  G4SAIDDATA
  EXTENSION tar.gz
  ENVVAR    G4SAIDXSDATA
  MD5SUM    d5d4e9541120c274aeed038c621d39da
  )

# - ABLA
geant4_add_dataset(
  NAME      G4ABLA
  VERSION   3.3
  FILENAME  G4ABLA
  EXTENSION tar.gz
  ENVVAR    G4ABLADATA
  MD5SUM    b25d093339e1e4532e31038653580ca6
  )


# - INCL
geant4_add_dataset(
  NAME      G4INCL
  VERSION   1.2
  FILENAME  G4INCL
  EXTENSION tar.gz
  ENVVAR    G4INCLDATA
  MD5SUM    0a76df936839bb557dae7254117eb58e
  )

# - ENSDFSTATE
geant4_add_dataset(
  NAME      G4ENSDFSTATE
  VERSION   3.0
  FILENAME  G4ENSDFSTATE
  EXTENSION tar.gz
  ENVVAR    G4ENSDFSTATEDATA
  MD5SUM    c500728534ce3e9fb2fefa0112eb3a74
  )

# - Channeling
geant4_add_dataset(
  NAME      G4CHANNELING
  VERSION   1.0
  FILENAME  G4CHANNELING
  EXTENSION tar.gz
  ENVVAR    G4CHANNELINGDATA
  MD5SUM    b2f692ec7109418c6354ea1ecbc62da7
  )

# - TENDL
option(GEANT4_INSTALL_DATASETS_TENDL "Install optional TENDL dataset" OFF)
mark_as_advanced(GEANT4_INSTALL_DATASETS_TENDL)

if(GEANT4_INSTALL_DATASETS_TENDL)
  geant4_add_dataset(
  NAME      G4TENDL
  VERSION   1.4
  FILENAME  G4TENDL
  EXTENSION tar.gz
  ENVVAR    G4PARTICLEHPDATA
  MD5SUM    ab17a1eec5ad4c8fcfafcfbc0de30ed8
  )
endif()

# - NuDEX
option(GEANT4_INSTALL_DATASETS_NUDEXLIB "Install optional NuDEXLib dataset" OFF)
mark_as_advanced(GEANT4_INSTALL_DATASETS_NUDEXLIB)

if(GEANT4_INSTALL_DATASETS_NUDEXLIB)
  geant4_add_dataset(
  NAME      G4NUDEXLIB
  VERSION   1.0
  FILENAME  G4NUDEXLIB
  EXTENSION tar.gz
  ENVVAR    G4NUDEXLIBDATA
  MD5SUM    09a85f907d2282dbf234d1784f436db3
  )
endif()

# - G4URRPT
option(GEANT4_INSTALL_DATASETS_URRPT "Install optional URRPT dataset" OFF)
mark_as_advanced(GEANT4_INSTALL_DATASETS_URRPT)

if(GEANT4_INSTALL_DATASETS_URRPT)
  geant4_add_dataset(
    NAME      G4URRPT
    VERSION   1.1
    FILENAME  G4URRPT
    EXTENSION tar.gz
    ENVVAR    G4URRPTDATA
    MD5SUM    ec9a2acb0745c8bfb6365ca3434bd3b8
    )
endif()
