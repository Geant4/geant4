# - Define datasets known and used by Geant4
# We keep this separate from the Geant4InstallData module for conveniance
# when updating and patching because datasets may change more rapidly.
# It allows us to decouple the dataset definitions from how they are
# checked/installed/configured
#

# - NDL
geant4_add_dataset(
  NAME      G4NDL
  VERSION   4.2
  FILENAME  G4NDL
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONHPDATA
  MD5SUM    2302d7ffeef1a450955302f027ee1cf4
  )

# - Low energy electromagnetics
geant4_add_dataset(
  NAME      G4EMLOW
  VERSION   6.32
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    9d3302072ba694b1d4505c330ed89d89
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
  VERSION   2.3
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    08848ebdd536280a0629d802040b70be
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
  VERSION   3.6
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    9566ff1614bea7455343328cb66abd6f
  )

# - Neutron XS
geant4_add_dataset(
  NAME      G4NEUTRONXS
  VERSION   1.2
  FILENAME  G4NEUTRONXS
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONXSDATA
  MD5SUM    092634b9258c7bc387cb83557ff1df81
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

