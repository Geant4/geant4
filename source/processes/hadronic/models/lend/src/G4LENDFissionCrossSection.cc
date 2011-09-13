#include "G4LENDFissionCrossSection.hh"

G4double G4LENDFissionCrossSection::getLENDCrossSection( G4GIDI_target* target , G4double ke , G4double temperature ) 
{
   if ( target == NULL ) return 0.0;
// 090407
//   return target->getElasticCrossSection( ke/MeV , temperature )*barn;
   //return target->getFissionCrossSectionAtE( ke/MeV , temperature )*barn;
   G4double result = target->getFissionCrossSectionAtE( ke/MeV , temperature )*barn;
   if ( result == 0.0 && ke/eV < 1.0e-4) 
   {
      G4double el = 1.0e-4*eV;
      G4double eh = 2.0e-4*eV;
      G4double xs_el = target->getFissionCrossSectionAtE( el/MeV , temperature )*barn;
      G4double xs_eh = target->getFissionCrossSectionAtE( eh/MeV , temperature )*barn;
      result = GetUltraLowEnergyExtrapolatedXS( el , eh , xs_el , xs_eh , ke );
   }
   return result;
}

