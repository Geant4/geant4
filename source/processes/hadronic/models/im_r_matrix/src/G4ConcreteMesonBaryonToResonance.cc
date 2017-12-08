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

#include "globals.hh"
#include "G4XAnnihilationChannel.hh"
#include "G4ConcreteMesonBaryonToResonance.hh"

G4BaryonWidth*           G4ConcreteMesonBaryonToResonance::baryonWidth = nullptr;
G4BaryonPartialWidth*    G4ConcreteMesonBaryonToResonance::baryonPartialWidth = nullptr;
G4ParticleTypeConverter* G4ConcreteMesonBaryonToResonance::particleTypeConverter = nullptr;

#ifdef G4MULTITHREADED
G4Mutex G4ConcreteMesonBaryonToResonance::concreteMesonBaryonToResonanceMutex = G4MUTEX_INITIALIZER;
#endif

G4BaryonWidth & G4ConcreteMesonBaryonToResonance::theBaryonWidth()
{
  if(!baryonWidth) { G4ConcreteMesonBaryonToResonance::InitialisePointers(); }
  return *baryonWidth;
}

G4BaryonPartialWidth & G4ConcreteMesonBaryonToResonance::theBaryonPartialWidth()
{
  if(!baryonPartialWidth) { G4ConcreteMesonBaryonToResonance::InitialisePointers(); }
  return *baryonPartialWidth;
}

G4ParticleTypeConverter & G4ConcreteMesonBaryonToResonance::myConv() 
{
  if(!particleTypeConverter) { G4ConcreteMesonBaryonToResonance::InitialisePointers(); }
  return *particleTypeConverter;
}

G4ConcreteMesonBaryonToResonance::G4ConcreteMesonBaryonToResonance(const G4ParticleDefinition* aPrimary,
								   const G4ParticleDefinition* bPrimary,
								   const G4ParticleDefinition* aSecondary,
								   const G4String& partWidthLabel)
  : thePrimary1(aPrimary), thePrimary2(bPrimary), theSecondary(aSecondary)
{
  crossSectionSource = new G4XAnnihilationChannel(aSecondary, 
						  theBaryonWidth(),
						  theBaryonPartialWidth(),
						  partWidthLabel);
  InitialisePointers();
}

G4ConcreteMesonBaryonToResonance::~G4ConcreteMesonBaryonToResonance()
{ 
  delete crossSectionSource;
}

void G4ConcreteMesonBaryonToResonance::InitialisePointers()
{
  if(!baryonWidth) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&concreteMesonBaryonToResonanceMutex);
    if (!baryonWidth)  { 
#endif
      baryonWidth = new G4BaryonWidth();
      baryonPartialWidth = new G4BaryonPartialWidth();
      particleTypeConverter = new G4ParticleTypeConverter(); 
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&concreteMesonBaryonToResonanceMutex);
#endif
  }
}

G4bool G4ConcreteMesonBaryonToResonance::IsInCharge(const G4KineticTrack& trk1, 
						    const G4KineticTrack& trk2) const
{
  if (myConv().GetGenericType(trk1)==myConv().GetGenericType(thePrimary1) && 
      myConv().GetGenericType(trk2)==myConv().GetGenericType(thePrimary2)) return true;
  if (myConv().GetGenericType(trk1)==myConv().GetGenericType(thePrimary2) && 
      myConv().GetGenericType(trk2)==myConv().GetGenericType(thePrimary1)) return true;
  return false;
}

const G4ParticleDefinition* G4ConcreteMesonBaryonToResonance::GetOutgoingParticle(const G4KineticTrack& trk1, 
										  const G4KineticTrack& trk2) const
{
  G4int secondaryIso3 = trk1.GetDefinition()->GetPDGiIsospin3() + trk2.GetDefinition()->GetPDGiIsospin3();
  const G4ParticleDefinition* state;
  if ( (state = myConv().FindIso3State(myConv().GetGenericType(theSecondary), secondaryIso3)) == NULL) 
  {
    G4cerr << "for "<<static_cast<G4int>(myConv().GetGenericType(theSecondary))<<" "<<secondaryIso3<<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, "G4ConcreteMesonBaryonToResonance: Can't find iso3 state!");
  }
  return state;
}
