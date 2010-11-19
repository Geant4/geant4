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
// 12-Nov-2003 Add energy check at lower side T. Koi
// 15-Nov-2006 Above 10GeV/n Cross Section become constant T. Koi (SLAC/SCCS)
// 23-Dec-2006 Isotope dependence adde by D. Wright
//

#include "G4IonsShenCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"


G4double G4IonsShenCrossSection::
GetZandACrossSection(const G4DynamicParticle* aParticle, G4int ZZ, 
                     G4int AA, G4double /*temperature*/)
{
   G4double xsection = 0.0;

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
   G4int Zp = G4int(aParticle->GetDefinition()->GetPDGCharge()/eplus + 0.5 ); 
   G4double ke_per_N = aParticle->GetKineticEnergy() / Ap; 
   if ( ke_per_N > 10*GeV ) ke_per_N = 10*GeV;

// Apply energy check, if less than lower limit then 0 value is returned
   // if (  ke_per_N < lowerLimit ) return xsection;

   G4int At = AA;
   G4int Zt = ZZ;
 
   G4double one_third = 1.0 / 3.0;

   G4double cubicrAt = std::pow ( G4double(At) , G4double(one_third) );  
   G4double cubicrAp = std::pow ( G4double(Ap) , G4double(one_third) );  

   G4double Rt = 1.12 * cubicrAt - 0.94 * ( 1.0 / cubicrAt );
   G4double Rp = 1.12 * cubicrAp - 0.94 * ( 1.0 / cubicrAp );

   G4double r = Rt + Rp + 3.2;   // in fm
   G4double b = 1.0;   // in MeV/fm
   G4double targ_mass =
     G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Zt, At);
   G4double proj_mass = aParticle->GetMass();
   G4double proj_momentum = aParticle->GetMomentum().mag();

   G4double Ecm = calEcmValue (proj_mass, targ_mass, proj_momentum); 

   G4double B = 1.44 * Zt * Zp / r - b * Rt * Rp / ( Rt + Rp ); 
   if(Ecm <= B) return xsection;
   //G4double ke_per_N = aParticle->GetKineticEnergy() / Ap; 

   G4double c = calCeValue ( ke_per_N / MeV  );  

   G4double R1 = r0 * (cubicrAt + cubicrAp + 1.85*cubicrAt*cubicrAp/(cubicrAt + cubicrAp) - c); 

   G4double R2 = 1.0 * ( At - 2 * Zt ) * Zp / ( Ap * At );


   G4double R3 = 0.176 / std::pow(G4double(Ecm), G4double(one_third)) * cubicrAt * cubicrAp /(cubicrAt + cubicrAp);

   G4double R = R1 + R2 + R3;

   xsection = 10 * pi * R * R * ( 1 - B / Ecm );   
   xsection = xsection * millibarn;   // mulitply xsection by millibarn

   return xsection; 
}


G4double G4IonsShenCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, const G4Element* anElement,
                G4double temperature)
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


G4double
G4IonsShenCrossSection::calEcmValue(const G4double mp, const G4double mt,
                                    const G4double Plab)
{
   G4double Elab = std::sqrt ( mp * mp + Plab * Plab );
   G4double Ecm = std::sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
   G4double Pcm = Plab * mt / Ecm;
   G4double KEcm = std::sqrt ( Pcm * Pcm + mp * mp ) - mp;
   return KEcm;
}


G4double G4IonsShenCrossSection::calCeValue(const G4double ke)
{
  // Calculate c value 
  // This value is indepenent from projectile and target particle 
  // ke is projectile kinetic energy per nucleon in the Lab system 
  // with MeV unit 
  // fitting function is made by T. Koi 
  // There are no data below 30 MeV/n in Kox et al., 

   G4double Ce; 
   G4double log10_ke = std::log10 ( ke );   
   if (log10_ke > 1.5) 
   {
     Ce = -10.0/std::pow(G4double(log10_ke), G4double(5)) + 2.0;
   }
   else
   {
     Ce = (-10.0/std::pow(G4double(1.5), G4double(5) ) + 2.0) /
         std::pow(G4double(1.5) , G4double(3)) * std::pow(G4double(log10_ke), G4double(3));
   }
   return Ce;
}
