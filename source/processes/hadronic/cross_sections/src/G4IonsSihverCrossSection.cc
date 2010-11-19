//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 18-Sep-2003 First version is written by T. Koi
// 23-Dec-2006 Isotope dependence added by D. Wright
//

#include "G4IonsSihverCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"

G4double G4IonsSihverCrossSection::
GetZandACrossSection(const G4DynamicParticle* aParticle,
                     G4int /*ZZ*/, G4int AA, G4double /*aTemperature*/)
{
   G4double xsection = 0.0;
   G4int At = AA;

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
 
   G4double one_third = 1.0 / 3.0;

   G4double cubicrAt = std::pow ( G4double(At) , G4double(one_third) ); 
   G4double cubicrAp = std::pow ( G4double(Ap) , G4double(one_third) );  

   G4double b0 = 1.581 - 0.876 * (1.0/cubicrAp + 1.0/cubicrAt);

   xsection = pi * square_r0 
            * std::pow(G4double(cubicrAp + cubicrAt - b0 * (1.0/cubicrAp + 1.0/cubicrAt)), G4double(2));
  
   return xsection; 
}


G4double G4IonsSihverCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, 
                const G4Element* anElement, G4double temperature)
{
  G4int nIso = anElement->GetNumberOfIsotopes();
  G4double xsection = 0;
    
  if (nIso) {
    G4double sig;
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    G4double* abundVector = anElement->GetRelativeAbundanceVector();
    G4int ZZ;
    G4int AA;
    
    for (G4int i = 0; i < nIso; i++) {
      ZZ = (*isoVector)[i]->GetZ();
      AA = (*isoVector)[i]->GetN();
      sig = GetZandACrossSection(aParticle, ZZ, AA, temperature);
      xsection += sig*abundVector[i];
    }
  
  } else {
    G4int ZZ = G4lrint(anElement->GetZ());
    G4int AA = G4lrint(anElement->GetN());
    xsection = GetZandACrossSection(aParticle, ZZ, AA, temperature);
  }
   
  return xsection;
}
