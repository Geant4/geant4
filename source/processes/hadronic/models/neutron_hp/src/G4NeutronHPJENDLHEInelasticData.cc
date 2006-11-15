// Class Description
// Cross-section data set for a high precision (based on JENDL_HE evaluated data
// libraries) description of elastic scattering 20 MeV ~ 3 GeV;
// Class Description - End

// 15-Nov-06 First Implementation is done by T. Koi (SLAC/SCCS)

#include "G4NeutronHPJENDLHEInelasticData.hh"
#include "G4Neutron.hh"

G4NeutronHPJENDLHEInelasticData::G4NeutronHPJENDLHEInelasticData()
:G4NeutronHPJENDLHEData( "Inelastic" , G4Neutron::Neutron() ) 
{
   ;
}
