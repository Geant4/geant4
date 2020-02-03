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
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
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
//      revised by V.Ivanchenko,    12 June 2003
//         fixed problem of massless particles
//      revised by V.Ivanchenko,    18 June 2003
//         take into account the case of virtual photons
//      revised by M.Kelsey         12 May 2010
//	   ensure that all constructors initialize all data members
//--------------------------------------------------------------

#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4PrimaryParticle.hh"

G4Allocator<G4DynamicParticle>*& pDynamicParticleAllocator()
{
    G4ThreadLocalStatic G4Allocator<G4DynamicParticle>* _instance = nullptr;
    return _instance;
}

static const G4double EnergyMomentumRelationAllowance = 1.0e-2*CLHEP::keV;
static const G4double EnergyMRA2 = 
  EnergyMomentumRelationAllowance*EnergyMomentumRelationAllowance;

////////////////////
G4DynamicParticle::G4DynamicParticle():
		   theMomentumDirection(0.0,0.0,1.0),
		   thePolarization(0.0,0.0,0.0),
		   theParticleDefinition(nullptr),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(0.0),
                   theLogKineticEnergy(DBL_MAX),
 		   theProperTime(0.0),
		   theDynamicalMass(0.0),
		   theDynamicalCharge(0.0),
		   theDynamicalSpin(0.0),
		   theDynamicalMagneticMoment(0.0),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
}

////////////////////
// -- constructors ----
////////////////////
G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
				     const G4ThreeVector& aMomentumDirection,
				     G4double aKineticEnergy):
		   theMomentumDirection(aMomentumDirection),
		   thePolarization(0.0,0.0,0.0),
		   theParticleDefinition(aParticleDefinition),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(aKineticEnergy),
		   theLogKineticEnergy(DBL_MAX),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
}

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
				     const G4ThreeVector& aMomentumDirection,
				     G4double aKineticEnergy,
                                     const G4double dynamicalMass):
		   theMomentumDirection(aMomentumDirection),
		   thePolarization(0.0,0.0,0.0),
		   theParticleDefinition(aParticleDefinition),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(aKineticEnergy),
                   theLogKineticEnergy(DBL_MAX),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
  if (std::abs(theDynamicalMass-dynamicalMass)> EnergyMomentumRelationAllowance) {
    if (dynamicalMass>EnergyMomentumRelationAllowance) theDynamicalMass= dynamicalMass;
    else  theDynamicalMass= 0.0;
  } 
}

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
                                     const G4ThreeVector& aParticleMomentum):
		   thePolarization(0.0,0.0,0.0),
		   theParticleDefinition(aParticleDefinition),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(0.0),
                   theLogKineticEnergy(DBL_MAX),
       		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
  SetMomentum(aParticleMomentum);  // 3-dim momentum is given
}

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
				     const G4LorentzVector   &aParticleMomentum):
		   thePolarization(0.0,0.0,0.0),
		   theParticleDefinition(aParticleDefinition),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(0.0),
                   theLogKineticEnergy(DBL_MAX),
 		   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
  Set4Momentum(aParticleMomentum);  // 4-momentum vector (Lorentz vector) is given
}

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition * aParticleDefinition,
                                     G4double totalEnergy,
				     const G4ThreeVector &aParticleMomentum):
		   thePolarization(0.0,0.0,0.0),
                   theParticleDefinition(aParticleDefinition),
		   theElectronOccupancy(nullptr),
                   thePreAssignedDecayProducts(nullptr),
		   primaryParticle(nullptr),
		   theKineticEnergy(0.0),
                   theLogKineticEnergy(DBL_MAX),
                   theProperTime(0.0),
		   theDynamicalMass(aParticleDefinition->GetPDGMass()),
		   theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
		   theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
		   theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment()),
                   thePreAssignedDecayTime(-1.0),
		   verboseLevel(1),
                   thePDGcode(0)
{
  // total energy and 3-dim momentum are given
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2 > 0.0) {
    G4double mass2 = totalEnergy*totalEnergy - pModule2;
    G4double PDGmass2 = (aParticleDefinition->GetPDGMass())*(aParticleDefinition->GetPDGMass());
    SetMomentumDirection(aParticleMomentum.unit());
    if (mass2 < EnergyMRA2) {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalEnergy);
    } else {
      if (std::abs(PDGmass2-mass2)>EnergyMRA2){
        theDynamicalMass = std::sqrt(mass2);
        SetKineticEnergy(totalEnergy-theDynamicalMass);
      } else {
	SetKineticEnergy(totalEnergy-theDynamicalMass);
      }
    }
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

////////////////////
G4DynamicParticle::G4DynamicParticle(const G4DynamicParticle &right):
  theMomentumDirection(right.theMomentumDirection),
  thePolarization(right.thePolarization),
  theParticleDefinition(right.theParticleDefinition),
  theElectronOccupancy(nullptr),
  thePreAssignedDecayProducts(nullptr),	// Do not copy preassignedDecayProducts
  primaryParticle(right.primaryParticle),
  theKineticEnergy(right.theKineticEnergy),
  theLogKineticEnergy(right.theLogKineticEnergy),
  theProperTime(right.theProperTime),
  theDynamicalMass(right.theDynamicalMass),
  theDynamicalCharge(right.theDynamicalCharge),
  theDynamicalSpin(right.theDynamicalSpin),
  theDynamicalMagneticMoment(right.theDynamicalMagneticMoment),
  thePreAssignedDecayTime(-1.0),
  verboseLevel(right.verboseLevel),
  thePDGcode(right.thePDGcode)
{
  if (right.theElectronOccupancy != nullptr) {
      theElectronOccupancy =
	new G4ElectronOccupancy(*right.theElectronOccupancy);
  }
}

////////////////////
// -- destructor ----
////////////////////
G4DynamicParticle::~G4DynamicParticle()
{
  //  delete thePreAssignedDecayProducts
  if (thePreAssignedDecayProducts != nullptr) delete thePreAssignedDecayProducts;
  thePreAssignedDecayProducts = nullptr;

  if (theElectronOccupancy != nullptr) delete theElectronOccupancy;
  theElectronOccupancy =nullptr;
}

////////////////////
// -- operators ----
////////////////////
G4DynamicParticle & G4DynamicParticle::operator=(const G4DynamicParticle &right)
{
  if (this != &right) {
    theMomentumDirection = right.theMomentumDirection;
    theParticleDefinition = right.theParticleDefinition;
    thePolarization = right.thePolarization;
    theKineticEnergy = right.theKineticEnergy;
    theProperTime = right.theProperTime;

    theDynamicalMass = right.theDynamicalMass;
    theDynamicalCharge = right.theDynamicalCharge;
    theDynamicalSpin = right.theDynamicalSpin;
    theDynamicalMagneticMoment = right.theDynamicalMagneticMoment;

    if (theElectronOccupancy != nullptr) delete theElectronOccupancy;
    if (right.theElectronOccupancy != nullptr){
      theElectronOccupancy =
             new G4ElectronOccupancy(*right.theElectronOccupancy);
    } else {
      theElectronOccupancy = nullptr;
    }

    // thePreAssignedDecayProducts must not be copied.
    thePreAssignedDecayProducts = nullptr;
    thePreAssignedDecayTime = -1.0;

    verboseLevel = right.verboseLevel;

    // Primary particle information must be preserved
    //*** primaryParticle = right.primaryParticle;

    thePDGcode = right.thePDGcode;
  }
  return *this;
}

////////////////////
void G4DynamicParticle::SetDefinition(const G4ParticleDefinition * aParticleDefinition)
{
  // remove preassigned decay
  if (thePreAssignedDecayProducts != nullptr) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << " G4DynamicParticle::SetDefinition()::"
             << "!!! Pre-assigned decay products is attached !!!! " << G4endl;
      G4cout << "!!! New Definition is " << aParticleDefinition->GetParticleName() 
	     << " !!! " << G4endl;
      G4cout << "!!! Pre-assigned decay products will be deleted !!!! " << G4endl;
    }
#endif
    delete thePreAssignedDecayProducts;
  }
  thePreAssignedDecayProducts = nullptr;

  theParticleDefinition = aParticleDefinition;

  // set Dynamic mass/charge
  theDynamicalMass = theParticleDefinition->GetPDGMass();
  theDynamicalCharge = theParticleDefinition->GetPDGCharge();
  theDynamicalSpin = theParticleDefinition->GetPDGSpin();
  theDynamicalMagneticMoment = theParticleDefinition->GetPDGMagneticMoment();

  // Set electron orbits
  if (theElectronOccupancy != nullptr) {
    delete theElectronOccupancy;
    theElectronOccupancy = nullptr;
  }
}

////////////////////
G4bool G4DynamicParticle::operator==(const G4DynamicParticle &right) const
{
  return (this == (G4DynamicParticle *) &right);
}

////////////////////
G4bool G4DynamicParticle::operator!=(const G4DynamicParticle &right) const
{
  return (this != (G4DynamicParticle *) &right);
}

////////////////////
// -- AllocateElectronOccupancy --
////////////////////
void  G4DynamicParticle::AllocateElectronOccupancy()
{
  if (G4IonTable::IsIon(theParticleDefinition)) {
    // Only ions can have ElectronOccupancy
    theElectronOccupancy = new G4ElectronOccupancy();

  } else {
    theElectronOccupancy = nullptr;

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
    SetMomentumDirection(momentum.unit());
    SetKineticEnergy(pModule2/(std::sqrt(pModule2 + mass*mass)+mass));
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

////////////////////
void G4DynamicParticle::Set4Momentum(const G4LorentzVector &momentum )
{
  G4double pModule2 = momentum.vect().mag2();
  if (pModule2>0.0) {
    SetMomentumDirection(momentum.vect().unit());
    G4double totalenergy = momentum.t();
    G4double mass2 = totalenergy*totalenergy - pModule2;
    G4double PDGmass2 = (theParticleDefinition->GetPDGMass())*(theParticleDefinition->GetPDGMass());
    if (mass2 < EnergyMRA2) {
      theDynamicalMass = 0.;
    } else if (std::abs(PDGmass2-mass2)>EnergyMRA2){
      theDynamicalMass = std::sqrt(mass2);
    }
    SetKineticEnergy(totalenergy-theDynamicalMass);
  } else {
    SetMomentumDirection(1.0,0.0,0.0);
    SetKineticEnergy(0.0);
  }
}

////////////////////
//  --- Dump Information --
////////////////////
#ifdef G4VERBOSE
void G4DynamicParticle::DumpInfo(G4int mode) const
{
  if (theParticleDefinition == 0) {
    G4cout << " G4DynamicParticle::DumpInfo():: !!!Particle type not defined !!!! " << G4endl;
  } else {
    G4cout << " Particle type - " << theParticleDefinition->GetParticleName() << G4endl
         << "   mass:        " << GetMass()/CLHEP::GeV <<  "[GeV]" <<G4endl
	   << "   charge:      " << GetCharge()/CLHEP::eplus <<  "[e]" <<G4endl
         << "   Direction x: " << GetMomentumDirection().x() << ", y: "
	 << GetMomentumDirection().y() << ", z: "
                             << GetMomentumDirection().z() << G4endl
         << "   Total Momentum = " << GetTotalMomentum()/CLHEP::GeV << "[GeV]" << G4endl
         << "   Momentum: "    << GetMomentum().x() /CLHEP::GeV << "[GeV]" << ", y: "
                               << GetMomentum().y() /CLHEP::GeV << "[GeV]" << ", z: "
                               << GetMomentum().z() /CLHEP::GeV << "[GeV]" << G4endl
         << "   Total Energy   = " << GetTotalEnergy()/CLHEP::GeV << "[GeV]"  << G4endl
         << "   Kinetic Energy = " << GetKineticEnergy() /CLHEP::GeV << "[GeV]" << G4endl
	 << " MagneticMoment  [MeV/T]: " << GetMagneticMoment()/CLHEP::MeV*CLHEP::tesla << G4endl
         << "   ProperTime     = " << GetProperTime() /CLHEP::ns <<  "[ns]" << G4endl;

    if (mode>0) {
      if( theElectronOccupancy != nullptr) {
	theElectronOccupancy->DumpInfo();
      }
    }
  }
}
#else
void G4DynamicParticle::DumpInfo(G4int) const
{
  return;
}
#endif

////////////////////////
G4double G4DynamicParticle::GetElectronMass() const
{
  return CLHEP::electron_mass_c2;
}
