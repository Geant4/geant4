//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

   G4double cubicrAt = std::pow ( G4double(At) , G4double(one_third) ); 
   G4double cubicrAp = std::pow ( G4double(Ap) , G4double(one_third) );  

   G4double b0 = 1.581 - 0.876 * ( 1.0 / cubicrAp + 1.0 / cubicrAt );

   xsection = pi * square_r0 
            * std::pow ( G4double(cubicrAp + cubicrAt - b0 * (  1.0 / cubicrAp + 1.0 / cubicrAt ) ), G4double(2) );
  
   return xsection; 
}
