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
//

#include "globals.hh"
#include "G4ios.hh"
#include "G4KineticTrack.hh"
#include "G4XResonance.hh"
#include "Randomize.hh"
#include "G4Proton.hh"

G4XResonance::G4XResonance(const G4ParticleDefinition* in1, 
			   const G4ParticleDefinition* in2,
			   G4int iIsospinOut1, G4double spinOut1, G4double massOut1,
			   G4int iIsospinOut2, G4double spinOut2, G4double massOut2,
			   G4String subType1, G4String subType2,
			   const G4VXResonanceTable& sigmaTable) : 
  isoOut1(iIsospinOut1), iSpinOut1(spinOut1), mOut1(massOut1),
  isoOut2(iIsospinOut2), iSpinOut2(spinOut2), mOut2(massOut2)
  
{
  table = sigmaTable.CrossSectionTable();
  // Check if there is a valid cross section table for this channel
  if (table == 0)
    throw G4HadronicException(__FILE__, __LINE__, "G4XResonance::G4XResonance - no cross section table available");
  
  name = in1->GetParticleName() + in2->GetParticleName() + " -> " + subType1 + subType2;
}


G4XResonance::~G4XResonance() 
{
  delete table;
  table = 0;
}


G4bool G4XResonance::operator==(const G4XResonance &right) const
{
  return (this == (G4XResonance *) &right);
}


G4bool G4XResonance::operator!=(const G4XResonance &right) const
{
  return (this != (G4XResonance *) &right);
}


G4String G4XResonance::Name() const
{
  return name;
}

G4double G4XResonance::CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
  G4bool dummy = false;
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();

  // pp -> trk1 + trk2 cross section
  G4double sigma = table->GetValue(sqrtS,dummy);
  
  // Isospin correction
  // from G4VXResonance; 
//   G4cout << "SigmaCheck "<<sigma;
  sigma *= IsospinCorrection(trk1,trk2,isoOut1,isoOut2,iSpinOut1,iSpinOut2); // from G4VXResonance
//   G4cout << " "<<sigma<<G4endl;

  // Detailed balance
  if (trk1.GetDefinition()->IsShortLived() || trk2.GetDefinition()->IsShortLived())
  {
   sigma *= DetailedBalance(trk1,trk2, isoOut1,isoOut2, iSpinOut1,iSpinOut2, mOut1,mOut2); 
  }

  return sigma;

}

