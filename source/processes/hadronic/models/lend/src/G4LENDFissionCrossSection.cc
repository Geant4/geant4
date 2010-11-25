using namespace std;
#include "G4LENDFissionCrossSection.hh"

G4double G4LENDFissionCrossSection::getLENDCrossSection( GIDI4GEANT_target* target , G4double ke , G4double temperature ) 
{
   if ( target == NULL ) return 0.0;
// 090407
//   return target->getElasticCrossSection( ke/MeV , temperature )*barn;
   return target->getFissionCrossSectionAtE( ke/MeV , temperature )*barn;
}

