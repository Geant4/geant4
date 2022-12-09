# - Define datasets known and used by Geant4
# We keep this separate from the Geant4InstallData module for conveniance
# when updating and patching because datasets may change more rapidly.
# It allows us to decouple the dataset definitions from how they are
# checked/installed/configured
#

# - NDL
geant4_add_dataset(
  NAME      G4NDL
  VERSION   4.7
  FILENAME  G4NDL
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONHPDATA
  MD5SUM    b001a2091bf9392e6833830347672ea2
  )

# - Low energy electromagnetics
geant4_add_dataset(
  NAME      G4EMLOW
  VERSION   8.2
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    07773e57be3f6f2ebb744da5ed574f6d
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
  VERSION   5.7
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    81ff27deb23af4aa225423e6b3a06b39
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
  VERSION   5.6
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    acc1dbeb87b6b708b2874ced729a3a8f
  )

# - Particle XS - replaces Neutron XS
geant4_add_dataset(
  NAME      G4PARTICLEXS
  VERSION   4.0
  FILENAME  G4PARTICLEXS
  EXTENSION tar.gz
  ENVVAR    G4PARTICLEXSDATA
  MD5SUM    d82a4d171d50f55864e28b6cd6f433c0
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
  VERSION   3.1
  FILENAME  G4ABLA
  EXTENSION tar.gz
  ENVVAR    G4ABLADATA
  MD5SUM    180f1f5d937733b207f8d5677f76296e
  )


# - INCL
geant4_add_dataset(
  NAME      G4INCL
  VERSION   1.0
  FILENAME  G4INCL
  EXTENSION tar.gz
  ENVVAR    G4INCLDATA
  MD5SUM    85fe937b6df46d41814f07175d3f5b51
  )

# - ENSDFSTATE
geant4_add_dataset(
  NAME      G4ENSDFSTATE
  VERSION   2.3
  FILENAME  G4ENSDFSTATE
  EXTENSION tar.gz
  ENVVAR    G4ENSDFSTATEDATA
  MD5SUM    6f18fce8f217e7aaeaa3711be9b2c7bf
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

