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
    G4Exception("G4XResonance::G4XResonance - no cross section table available");
  
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

