
using namespace std;
#include "G4LENDInelasticCrossSection.hh"

G4double G4LENDInelasticCrossSection::getLENDCrossSection( GIDI4GEANT_target* target , G4double ke , G4double temperature ) 
{
   if ( target == NULL ) return 0.0;
// 090407
//   return target->getElasticCrossSection( ke/MeV , temperature )*barn;
   return target->getOthersCrossSectionAtE( ke/MeV , temperature )*barn;
}

