// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VShortLivedParticle.cc,v 1.3 1999-12-15 14:51:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
//      History: first implementation, based on object model of
//      28 June 1998 H.Kurashige
// --------------------------------------------------------------

#include "G4VShortLivedParticle.hh"

G4VShortLivedParticle::G4VShortLivedParticle(const G4String&  aName,  
               G4double         mass,     
               G4double         width,
               G4double         charge,   
               G4int            iSpin,
               G4int            iParity,
               G4int            iConjugation,
               G4int            iIsospin,   
               G4int            iIsospinZ, 
               G4int            gParity,
               const G4String&  pType,
               G4int            lepton,
               G4int            baryon,
               G4int            encoding,
               G4bool           stable,
               G4double         lifetime,
	       G4DecayTable     *decaytable)
  :G4ParticleDefinition( aName,mass,width,charge,iSpin,iParity,
           iConjugation,iIsospin,iIsospinZ,gParity,pType,
           lepton,baryon,encoding,stable,lifetime,decaytable, true)
{
   
}

void            G4VShortLivedParticle::ResetCuts()
{
  G4cout << "G4VShortLivedParticle::ResetCuts() causes no effect!!" << G4endl;
}
void            G4VShortLivedParticle::SetCuts(G4double )
{
  G4cout << "G4VShortLivedParticle::SetCuts() causes no effect!!" << G4endl;
}
void            G4VShortLivedParticle::ReCalcCuts()
{
  G4cout << "G4VShortLivedParticle::ReCalcCuts() causes no effect!!" << G4endl;
}

G4double      	G4VShortLivedParticle::GetLengthCuts() const
{
  G4cout << "G4VShortLivedParticle::GetLengthCuts() causes no effect!!" << G4endl;
  return -1.0;
}

G4double*	G4VShortLivedParticle::GetEnergyCuts() const
{
  G4cout << "G4VShortLivedParticle::GetLengthCuts() causes no effect!!" << G4endl;
  return 0;
}

G4double      	G4VShortLivedParticle::GetEnergyThreshold(const G4Material* ) const
{
  G4cout << "G4VShortLivedParticle::GetEnergyThreshold() causes no effect!!" << G4endl;
  return -1.0;
}

const G4VShortLivedParticle & G4VShortLivedParticle::operator=(const G4VShortLivedParticle& right)
{
  if (this != &right) {
  } return right;
}

G4int G4VShortLivedParticle::operator==(const G4VShortLivedParticle &right) const
{
  return (this->GetParticleName() == right.GetParticleName());
}

G4int G4VShortLivedParticle::operator!=(const G4VShortLivedParticle &right) const
{
  return (this->GetParticleName() != right.GetParticleName());
}
