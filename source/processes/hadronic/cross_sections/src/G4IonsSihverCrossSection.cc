// 18-Sep-2003 First version is written by T. Koi

#include "G4IonsSihverCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4double G4IonsSihverCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, 
                const G4Element* anElement, G4double )
{
   G4double xsection = 0.0;

   G4int At = int ( anElement->GetN() + 0.5 );
   //Zt = int ( anElement->GetZ() + 0.5 );  // not used 

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
   //Zp = aParticle->GetDefinition()->GetPDGCharge();   // not used  
 
   G4double one_third = 1.0 / 3.0;

   G4double cubicrAt = pow ( At , one_third ); 
   G4double cubicrAp = pow ( Ap , one_third );  

   G4double b0 = 1.581 - 0.876 * ( 1.0 / cubicrAp + 1.0 / cubicrAt );

   xsection = pi * square_r0 
            * pow ( cubicrAp + cubicrAt - b0 * (  1.0 / cubicrAp + 1.0 / cubicrAt ) , 2 );
  
   return xsection; 
}
