# - Define datasets known and used by Geant4
# We keep this separate from the Geant4InstallData module for conveniance
# when updating and patching because datasets may change more rapidly.
# It allows us to decouple the dataset definitions from how they are
# checked/installed/configured
#

# - NDL
geant4_add_dataset(
  NAME      G4NDL
  VERSION   4.5
  FILENAME  G4NDL
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONHPDATA
  MD5SUM    fd29c45fe2de432f1f67232707b654c0
  )

# - Low energy electromagnetics
geant4_add_dataset(
  NAME      G4EMLOW
  VERSION   6.48
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    844064faa16a063a6a08406dc7895b68
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
  VERSION   3.2
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    01d5ba17f615d3def01f7c0c6b19bd69
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
  VERSION   4.3
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    9f1630a5d9f00b4ba1ffc5f7df174827
  )

# - Neutron XS
geant4_add_dataset(
  NAME      G4NEUTRONXS
  VERSION   1.4
  FILENAME  G4NEUTRONXS
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONXSDATA
  MD5SUM    665a12771267e3b31a08c622ba1238a7
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
  VERSION   1.0
  FILENAME  RealSurface
  EXTENSION tar.gz
  ENVVAR    G4REALSURFACEDATA
  MD5SUM    0dde95e00fcd3bcd745804f870bb6884
  )

# - SAID
geant4_add_dataset(
  NAME      G4SAIDDATA
  VERSION   1.1
  FILENAME  G4SAIDDATA
  EXTENSION tar.gz
  ENVVAR    G4SAIDXSDATA
  MD5SUM    d88a31218fdf28455e5c5a3609f7216f
  )

# - ABLA
geant4_add_dataset(
  NAME      G4ABLA
  VERSION   3.0
  FILENAME  G4ABLA
  EXTENSION tar.gz
  ENVVAR    G4ABLADATA
  MD5SUM    d7049166ef74a592cb97df0ed4b757bd
  )

# - ENSDFSTATE
geant4_add_dataset(
  NAME      G4ENSDFSTATE
  VERSION   1.2
  FILENAME  G4ENSDFSTATE
  EXTENSION tar.gz
  ENVVAR    G4ENSDFSTATEDATA
  MD5SUM    a2e88f2c626141e4be4587c838832707
  )

