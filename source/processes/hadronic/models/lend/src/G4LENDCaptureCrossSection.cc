#include "G4LENDCaptureCrossSection.hh"

G4double G4LENDCaptureCrossSection::getLENDCrossSection( G4GIDI_target* target , G4double ke , G4double temperature ) 
{
   if ( target == NULL ) return 0.0;
// 090407
//   return target->getElasticCrossSection( ke/MeV , temperature )*barn;
   return target->getCaptureCrossSectionAtE( ke/MeV , temperature )*barn;
}

