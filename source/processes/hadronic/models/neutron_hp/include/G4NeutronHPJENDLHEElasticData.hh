#ifndef G4NeutronHPJENDLHEElasticData_h
#define G4NeutronHPJENDLHEElasticData_h 1

// Class Description
// Cross-section data set for a high precision (based on JENDL_HE evaluated data
// libraries) description of elastic scattering 20 MeV ~ 3 GeV;
// Class Description - End

// 15-Nov-06 First Implementation is done by T. Koi (SLAC/SCCS)

#include "G4NeutronHPJENDLHEData.hh"

class G4NeutronHPJENDLHEElasticData : public G4NeutronHPJENDLHEData
{
   public:
      G4NeutronHPJENDLHEElasticData();
};

#endif
