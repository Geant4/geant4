# - Define datasets known and used by Geant4
# We keep this separate from the Geant4InstallData module for conveniance
# when updating and patching because datasets may change more rapidly.
# It allows us to decouple the dataset definitions from how they are
# checked/installed/configured
#

# - NDL
geant4_add_dataset(
  NAME      G4NDL
  VERSION   4.6
  FILENAME  G4NDL
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONHPDATA
  MD5SUM    d07e43499f607e01f2c1ce06d7a09f3e
  )

# - Low energy electromagnetics
geant4_add_dataset(
  NAME      G4EMLOW
  VERSION   7.9.1
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    1308dc5d73b5539557e2ec6b96f7e5ef
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
  VERSION   5.5
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    707514c864414089af9671db0f656e35
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
  VERSION   5.4
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    08abe2bcc0bcd1ac4bbe09f5ae69cdbe
  )

# - Particle XS - replaces Neutron XS
geant4_add_dataset(
  NAME      G4PARTICLEXS
  VERSION   2.1
  FILENAME  G4PARTICLEXS
  EXTENSION tar.gz
  ENVVAR    G4PARTICLEXSDATA
  MD5SUM    24a68bb627a95629e2edcd098131d6b3
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
  VERSION   2.1.1
  FILENAME  G4RealSurface
  EXTENSION tar.gz
  ENVVAR    G4REALSURFACEDATA
  MD5SUM    1d0fcc24c7082edae1e22a3d43fbb4d9
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
  VERSION   2.2
  FILENAME  G4ENSDFSTATE
  EXTENSION tar.gz
  ENVVAR    G4ENSDFSTATEDATA
  MD5SUM    495439cf600225753d7bd99825e5c6bc
  )

