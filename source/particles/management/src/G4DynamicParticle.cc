// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DynamicParticle.cc,v 1.9 2001-03-05 08:32:38 kurasige Exp $
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
//      revised by Hisaya Kurashige, 5  Mar 2001
//         fixed  SetDefinition()
//--------------------------------------------------------------
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4Allocator<G4DynamicParticle> aDynamicParticleAllocator;

static const G4double EnergyMomentumRelationAllowance = keV;

////////////////////
G4DynamicParticle::G4DynamicParticle():
		   theParticleDefinition(0),
		   theMomentumDirection(),
		   theKineticEnergy(0.0),
 		   theProperTime(0.0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1)
{  
   theDynamicalMass = 0.0; 
   theDynamicalCharge= 0.0;
   theElectronOccupancy = 0; 
}

////////////////////
// -- constructors ----
////////////////////
G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
				     const G4ThreeVector& aMomentumDirection,
				     G4double aKineticEnergy):
		   theParticleDefinition(aParticleDefinition),
		   theMomentumDirection(aMomentumDirection),
		   theKineticEnergy(aKineticEnergy),
 		   theProperTime(0.0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1)
{  
  // set dynamic charge/mass
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  theDynamicalCharge = aParticleDefinition->GetPDGCharge();
  AllocateElectronOccupancy();
}   

////////////////////
G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
                                     const G4ThreeVector& aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
       		   theProperTime(0.0),
                   thePreAssignedDecayProducts(0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1)
{
  // set dynamic charge/mass
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  theDynamicalCharge = aParticleDefinition->GetPDGCharge();
  AllocateElectronOccupancy();

  // 3-dim momentum is given
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

////////////////////
G4DynamicParticle::G4DynamicParticle(G4ParticleDefinition * aParticleDefinition,
				     const G4LorentzVector   &aParticleMomentum):
		   theParticleDefinition(aParticleDefinition),
 		   theProperTime(0.0),
                   thePreAssignedDecayTime(-1.0),
                   thePreAssignedDecayProducts(0),
		   verboseLevel(1)
{
   // set dynamic charge/mass
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  theDynamicalCharge = aParticleDefinition->GetPDGCharge();
  AllocateElectronOccupancy();
  
  // 4-momentum vector (Lorentz vecotr) is given
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
                   thePreAssignedDecayProducts(0),
                   theProperTime(0.0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1)
{
   // set dynamic charge/mass
  theDynamicalMass = aParticleDefinition->GetPDGMass();
  theDynamicalCharge = aParticleDefinition->GetPDGCharge();
  AllocateElectronOccupancy();
  
  // total energy and momentum direction are given
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2>0.0) {
    G4double pModule = sqrt(pModule2);
    SetMomentumDirection(aParticleMomentum.x()/pModule,
                         aParticleMomentum.y()/pModule,
                         aParticleMomentum.z()/pModule);
    if (totalEnergy > pModule) {
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

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4DynamicParticle &right)
{
  theDynamicalMass = right.theDynamicalMass;
  theDynamicalCharge = right.theDynamicalCharge;
  if (right.theElectronOccupancy != 0){
      theElectronOccupancy = 
	new G4ElectronOccupancy(*right.theElectronOccupancy);
  } else {
     theElectronOccupancy = 0;
  }

  theParticleDefinition = right.theParticleDefinition;
  theMomentumDirection = right.theMomentumDirection;
  theKineticEnergy = right.theKineticEnergy;
  thePolarization = right.thePolarization;
  verboseLevel = right.verboseLevel;

  // proper time is set to zero
  theProperTime = 0.0;

  // thePreAssignedDecayProducts/Time must not be copied.
  thePreAssignedDecayProducts = 0;
  thePreAssignedDecayTime = -1.0;

}

////////////////////
// -- destructor ----
////////////////////
G4DynamicParticle::~G4DynamicParticle() {

  //  delete thePreAssignedDecayProducts
  if (thePreAssignedDecayProducts != 0) delete thePreAssignedDecayProducts;
  thePreAssignedDecayProducts = 0;

  if (theElectronOccupancy != 0) delete theElectronOccupancy;
  theElectronOccupancy =0;
}


////////////////////
// -- operators ----
////////////////////
G4DynamicParticle & G4DynamicParticle::operator=(const G4DynamicParticle &right)
{
  if (this != &right) {
    theDynamicalMass = right.theDynamicalMass;
    theDynamicalCharge = right.theDynamicalCharge;
    if (right.theElectronOccupancy != 0){
      theElectronOccupancy = 
             new G4ElectronOccupancy(*right.theElectronOccupancy);
    } else {
      theElectronOccupancy = 0;
    }

    theParticleDefinition = right.theParticleDefinition;
    theMomentumDirection = right.theMomentumDirection;
    theKineticEnergy = right.theKineticEnergy;
    thePolarization = right.thePolarization;
    theProperTime = right.theProperTime;
    verboseLevel = right.verboseLevel;
    
    // thePreAssignedDecayProducts must not be copied.
    thePreAssignedDecayProducts = 0;
    thePreAssignedDecayTime = -1.0;

  }
  return *this;
}

////////////////////
void G4DynamicParticle::SetDefinition(G4ParticleDefinition * aParticleDefinition)
{
  // remove preassigned decay
  if (thePreAssignedDecayProducts != 0) {
    G4cout << " G4DynamicParticle::SetDefinition()::";
    G4cout << "!!! Pre-assigned decay products is attached !!!! " << G4endl; 
    DumpInfo(0); 
    G4cout << "!!! New Definition is " << aParticleDefinition->GetParticleName() << " !!! " << G4endl; 
    G4cout << "!!! Pre-assigned decay products will be deleted !!!! " << G4endl; 
    delete thePreAssignedDecayProducts;
  }
  thePreAssignedDecayProducts = 0;

  theParticleDefinition = aParticleDefinition;
  // set Dynamic mass/chrge
  theDynamicalMass = theParticleDefinition->GetPDGMass();
  theDynamicalCharge = theParticleDefinition->GetPDGCharge();

  // Set electron orbits
  if (theElectronOccupancy != 0) delete theElectronOccupancy;
  theElectronOccupancy =0;
  AllocateElectronOccupancy();

}

////////////////////
G4int G4DynamicParticle::operator==(const G4DynamicParticle &right) const
{
  return (this == (G4DynamicParticle *) &right);
}

////////////////////
G4int G4DynamicParticle::operator!=(const G4DynamicParticle &right) const
{
  return (this != (G4DynamicParticle *) &right);
}



////////////////////
// -- AllocateElectronOccupancy --
////////////////////
void  G4DynamicParticle::AllocateElectronOccupancy()
{
  G4ParticleDefinition* particle = GetDefinition();

  if (G4IonTable::IsIon(particle)) {
    // Only ions can have ElectronOccupancy 
    theElectronOccupancy = new G4ElectronOccupancy();

  } else {
    theElectronOccupancy = 0;

  }
}

////////////////////
// -- methods for setting Energy/Momentum  --
////////////////////
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

////////////////////
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


////////////////////
//  --- Dump Information --
////////////////////
void G4DynamicParticle::DumpInfo(G4int mode) const
{
  if (theParticleDefinition == 0) {
    G4cout << " G4DynamicParticle::DumpInfo():: !!!Particle type not defined !!!! " << G4endl; 
  } else {
    G4cout << " Particle type - " << theParticleDefinition->GetParticleName() << G4endl
         << "   mass:        " << GetMass()/GeV <<  "[GeV]" <<G4endl
         << "   charge:      " << GetCharge()/eplus <<  "[e]" <<G4endl
         << "   Direction x: " << GetMomentumDirection().x() << ", y: "
	 << GetMomentumDirection().y() << ", z: "
                             << GetMomentumDirection().z() << G4endl
         << "   Total Momentum = " << GetTotalMomentum() /GeV << "[GeV]" << G4endl
         << "   Momentum: "    << GetMomentum().x() /GeV << "[GeV]" << ", y: "
                               << GetMomentum().y() /GeV << "[GeV]" << ", z: "
                               << GetMomentum().z() /GeV << "[GeV]" << G4endl
         << "   Total Energy   = " << GetTotalEnergy()/GeV << "[GeV]"  << G4endl
         << "   Kinetic Energy = " << GetKineticEnergy() /GeV << "[GeV]" << G4endl
         << "   ProperTime     = " << GetProperTime() /ns <<  "[ns]" << G4endl;
    if (mode>0) {
      if( theElectronOccupancy != 0) {
	theElectronOccupancy->DumpInfo();
      }
    }
  }
}

////////////////////////
G4double  G4DynamicParticle::GetElectronMass() const
{
  static G4double electronMass = 0.0;

  // check if electron exits and get the mass
  if (electronMass<=0.0) {
    G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    if (electron == 0) {
      G4Exception("G4DynamicParticle: G4Electron is not defined !!"); 
    }
    electronMass = electron->GetPDGMass();
  }

  return electronMass;
}




