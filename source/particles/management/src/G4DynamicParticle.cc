// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DynamicParticle.cc,v 1.1 1999-01-07 16:10:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//      ---------------- G4DynamicParticle  ----------------
//      first implementation by Makoto Asai, 29 January 1996
//      revised by G.Cosmo, 29 February 1996
//      revised by H.Kurashige 06 May 1996
//      revised by Hisaya Kurashige, 27 July 1996
//         modify thePreAssignedDecayProducts
//         add   void SetMomentum(G4ThreeVector &momentum)
//         add   void Set4Momentum(G4LorentzVector &momentum)
//         add   G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
//                                 G4LorentzVector &p4vector)
//      revised by Hisaya Kurashige, 19 Oct 1996
//         add    theKillProcess
//         add    ProperTime
//      revised by Hisaya Kurashige, 26 Mar 1997
//         modify destructor
//      revised by Hisaya Kurashige, 05 June 1997
//         modify DumpInfo()
//      revised by Hisaya Kurashige, 5  June 1998
//         remove    theKillProcess
//--------------------------------------------------------------
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"

G4Allocator<G4DynamicParticle> aDynamicParticleAllocator;

static const G4double EnergyMomentumRelationAllowance = keV;

G4DynamicParticle::G4DynamicParticle():
		   theParticleDefinition(NULL),
		   theMomentumDirection(),
		   theKineticEnergy(0.0),
 		   theProperTime(0.0),
                   thePreAssignedDecayProducts(NULL),
		   verboseLevel(1)
{  
   theDynamicalMass = 0.0; 
}

// -- constructors ----
G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
				     const G4ParticleMomentum& aMomentumDirection,
				     G4double aKineticEnergy):
		   theParticleDefinition(aParticleDefinition),
		   theMomentumDirection(aMomentumDirection),
		   theKineticEnergy(aKineticEnergy),
 		   theProperTime(0.0),
                   thePreAssignedDecayProducts(NULL),
		   verboseLevel(1)
{  
 theDynamicalMass = aParticleDefinition->GetPDGMass();
}

G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                                     const G4ThreeVector& aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
       		   theProperTime(0.0),
                   thePreAssignedDecayProducts(NULL),
		   verboseLevel(1)
{
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2>0.0) {
    G4double mass = theDynamicalMass;
    SetKineticEnergy(sqrt(pModule2+mass*mass)-mass);
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(aParticleMomentum.x()/pModule,
                         aParticleMomentum.y()/pModule,
                         aParticleMomentum.z()/pModule);
  } else {  
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
				     const G4LorentzVector   &aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
 		   theProperTime(0.0),
                   thePreAssignedDecayProducts(NULL),
		   verboseLevel(1)
{
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  G4double pModule2 = aParticleMomentum.x()*aParticleMomentum.x()
                       + aParticleMomentum.y()*aParticleMomentum.y()
                        + aParticleMomentum.z()*aParticleMomentum.z();
  if (pModule2>0.0) {
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(aParticleMomentum.x()/pModule,
                         aParticleMomentum.y()/pModule,
                         aParticleMomentum.z()/pModule);
    G4double totalenergy = aParticleMomentum.t();
    if (totalenergy > pModule) {
      G4double mass = sqrt(totalenergy*totalenergy - pModule2);
      theDynamicalMass = mass;
      SetKineticEnergy(totalenergy-mass);
    } else {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalenergy);
    }
  } else {  
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                                     G4double totalEnergy,  
				     const G4ThreeVector &aParticleMomentum):
                   theParticleDefinition(aParticleDefinition),
                   thePreAssignedDecayProducts(NULL),
                   theProperTime(0.0),
		   verboseLevel(1)
{
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2>0.0) {
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(aParticleMomentum.x()/pModule,
                         aParticleMomentum.y()/pModule,
                         aParticleMomentum.z()/pModule);
    if (totalEnergy > pModule2) {
      G4double mass = sqrt(totalEnergy*totalEnergy - pModule2);
      theDynamicalMass = mass;
      SetKineticEnergy(totalEnergy-mass);
    } else {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalEnergy);
    }
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

G4DynamicParticle::G4DynamicParticle(const G4DynamicParticle &right)
{
  SetMass(right.theDynamicalMass);
  SetDefinition(right.GetDefinition());
  theMomentumDirection = right.GetMomentumDirection();
  theKineticEnergy = right.GetKineticEnergy();
  thePolarization = right.GetPolarization();
  verboseLevel = right.GetVerboseLevel();

  // proper time is set to zero
  theProperTime = 0.0;

  // thePreAssignedDecayProducts must not be copied.
  thePreAssignedDecayProducts = NULL;

}

// -- destructor ----
G4DynamicParticle::~G4DynamicParticle() {

  //  delete thePreAssignedDecayProducts
  if (thePreAssignedDecayProducts != NULL) delete thePreAssignedDecayProducts;
  thePreAssignedDecayProducts = NULL;
}


G4DynamicParticle & G4DynamicParticle::operator=(const G4DynamicParticle &right)
{
  if (this != &right) {
    SetDefinition(right.GetDefinition());
    SetMass(right.theDynamicalMass);
    theMomentumDirection = right.GetMomentumDirection();
    theKineticEnergy = right.GetKineticEnergy();
    thePolarization = right.GetPolarization();
    theProperTime = right.GetProperTime();
    verboseLevel = right.GetVerboseLevel();
    
    // thePreAssignedDecayProducts must not be copied.
    thePreAssignedDecayProducts = NULL;
  }
  return *this;
}

G4int G4DynamicParticle::operator==(const G4DynamicParticle &right) const
{
  return (this == (G4DynamicParticle *) &right);
}

G4int G4DynamicParticle::operator!=(const G4DynamicParticle &right) const
{
  return (this != (G4DynamicParticle *) &right);
}


void G4DynamicParticle::SetMomentum(const G4ThreeVector &momentum)
{
  G4double pModule2 = momentum.mag2();
  if (pModule2>0.0) {
    G4double mass = theDynamicalMass;
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(momentum.x()/pModule,
                         momentum.y()/pModule,
                         momentum.z()/pModule);
    SetKineticEnergy(sqrt(pModule2 + mass*mass)-mass);
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

void G4DynamicParticle::Set4Momentum(const G4LorentzVector &momentum )
{
  G4double pModule2 = momentum.x()*momentum.x()
                       + momentum.y()*momentum.y()
                        + momentum.z()*momentum.z();
  if (pModule2>0.0) {
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(momentum.x()/pModule,
                         momentum.y()/pModule,
                         momentum.z()/pModule);
    G4double totalenergy = momentum.t();
    if (totalenergy > pModule) {
      G4double mass = sqrt(totalenergy*totalenergy - pModule2);
      theDynamicalMass = mass;
      SetKineticEnergy(totalenergy-mass);
    } else {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalenergy);
    }
  } else {  
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

void G4DynamicParticle::DumpInfo() const
{
  if (theParticleDefinition == NULL) {
    G4cout << " G4DynamicParticle::DumpInfo():: !!!Particle type not defined !!!! " << endl; 
  } else {
    G4cout << " Particle type - " << theParticleDefinition->GetParticleName() << endl
         << "   mass:        " << GetMass()/GeV <<  "[GeV]" <<endl
         << "   Direction x: " << GetMomentumDirection().x() << ", y: "
	 << GetMomentumDirection().y() << ", z: "
                             << GetMomentumDirection().z() << endl
         << "   Total Momentum = " << GetTotalMomentum() /GeV << "[GeV]" << endl
         << "   Momentum: "    << GetMomentum().x() /GeV << "[GeV]" << ", y: "
                               << GetMomentum().y() /GeV << "[GeV]" << ", z: "
                               << GetMomentum().z() /GeV << "[GeV]" << endl
         << "   Total Energy   = " << GetTotalEnergy()/GeV << "[GeV]"  << endl
         << "   Kinetic Energy = " << GetKineticEnergy() /GeV << "[GeV]" << endl
         << "   ProperTime     = " << GetProperTime() /ns <<  "[ns]" << endl;
  }
}






