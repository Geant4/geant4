#include "G4LENDCaptureCrossSection.hh"

G4double G4LENDCaptureCrossSection::getLENDCrossSection( G4GIDI_target* target , G4double ke , G4double temperature ) 
{
   if ( target == NULL ) return 0.0;
// 090407
   //return target->getCaptureCrossSection( ke/MeV , temperature )*barn;
   //return target->getCaptureCrossSectionAtE( ke/MeV , temperature )*barn;
   G4double result = target->getCaptureCrossSectionAtE( ke/MeV , temperature )*barn;
   if ( result == 0.0 && ke/eV < 1.0e-4) 
   {
      G4double el = 1.0e-4*eV;
      G4double eh = 2.0e-4*eV;
      G4double xs_el = target->getCaptureCrossSectionAtE( el/MeV , temperature )*barn;
      G4double xs_eh = target->getCaptureCrossSectionAtE( eh/MeV , temperature )*barn;
      result = GetUltraLowEnergyExtrapolatedXS( el , eh , xs_el , xs_eh , ke );
   }
   return result;
}
