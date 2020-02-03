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
<<<<<<< HEAD
  VERSION   6.48
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    844064faa16a063a6a08406dc7895b68
=======
  VERSION   7.9
  FILENAME  G4EMLOW
  EXTENSION tar.gz
  ENVVAR    G4LEDATA
  MD5SUM    d28a09f0c93243522512cf2a3a733348
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  )

# - Photon evaporation
geant4_add_dataset(
  NAME      PhotonEvaporation
<<<<<<< HEAD
  VERSION   3.2
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    01d5ba17f615d3def01f7c0c6b19bd69
=======
  VERSION   5.5
  FILENAME  G4PhotonEvaporation
  EXTENSION tar.gz
  ENVVAR    G4LEVELGAMMADATA
  MD5SUM    707514c864414089af9671db0f656e35
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  )

# - Radioisotopes
geant4_add_dataset(
  NAME      RadioactiveDecay
<<<<<<< HEAD
  VERSION   4.3.2
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    ed171641682cf8c10fc3f0266c8d482e
=======
  VERSION   5.4
  FILENAME  G4RadioactiveDecay
  EXTENSION tar.gz
  ENVVAR    G4RADIOACTIVEDATA
  MD5SUM    08abe2bcc0bcd1ac4bbe09f5ae69cdbe
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  )

# - Particle XS - replaces Neutron XS
geant4_add_dataset(
<<<<<<< HEAD
  NAME      G4NEUTRONXS
  VERSION   1.4
  FILENAME  G4NEUTRONXS
  EXTENSION tar.gz
  ENVVAR    G4NEUTRONXSDATA
  MD5SUM    665a12771267e3b31a08c622ba1238a7
=======
  NAME      G4PARTICLEXS
  VERSION   2.1
  FILENAME  G4PARTICLEXS
  EXTENSION tar.gz
  ENVVAR    G4PARTICLEXSDATA
  MD5SUM    24a68bb627a95629e2edcd098131d6b3
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
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
  VERSION   2.0
  FILENAME  G4SAIDDATA
  EXTENSION tar.gz
  ENVVAR    G4SAIDXSDATA
  MD5SUM    d5d4e9541120c274aeed038c621d39da
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
  VERSION   1.2.3
  FILENAME  G4ENSDFSTATE
  EXTENSION tar.gz
  ENVVAR    G4ENSDFSTATEDATA
  MD5SUM    98fef898ea35df4010920ad7ad88f20b
  )

