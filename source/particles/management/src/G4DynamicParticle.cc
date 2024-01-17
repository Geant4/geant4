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
// G4DynamicParticle class implementation
//
// History:
// - 2 December 1995, G.Cosmo - first design, based on object model.
// - 29 January 1996, M.Asai - first implementation.
// - 1996 - 2007,     H.Kurashige - revisions.
// - 15 March 2019,   M.Novak - log-kinetic energy value is computed only
//                    on demand if its stored value is not up-to-date.
//---------------------------------------------------------------------

#include "G4DynamicParticle.hh"

#include "G4DecayProducts.hh"
#include "G4IonTable.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"

G4Allocator<G4DynamicParticle>*& pDynamicParticleAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4DynamicParticle>* _instance = nullptr;
  return _instance;
}

static const G4double EnergyMomentumRelationAllowance = 1.0e-2 * CLHEP::keV;
static const G4double EnergyMRA2 =
  EnergyMomentumRelationAllowance * EnergyMomentumRelationAllowance;

G4DynamicParticle::G4DynamicParticle()
  : theMomentumDirection(0.0, 0.0, 1.0), thePolarization(0.0, 0.0, 0.0)
{}

G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                                     const G4ThreeVector& aMomentumDirection,
                                     G4double aKineticEnergy)
  : theMomentumDirection(aMomentumDirection),
    thePolarization(0.0, 0.0, 0.0),
    theParticleDefinition(aParticleDefinition),
    theKineticEnergy(aKineticEnergy),
    theDynamicalMass(aParticleDefinition->GetPDGMass()),
    theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
    theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
    theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment())
{}

G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                                     const G4ThreeVector& aMomentumDirection,
                                     G4double aKineticEnergy, const G4double dynamicalMass)
  : theMomentumDirection(aMomentumDirection),
    thePolarization(0.0, 0.0, 0.0),
    theParticleDefinition(aParticleDefinition),
    theKineticEnergy(aKineticEnergy),
    theDynamicalMass(aParticleDefinition->GetPDGMass()),
    theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
    theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
    theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment())
{
  if (std::abs(theDynamicalMass - dynamicalMass) > EnergyMomentumRelationAllowance) {
    if (dynamicalMass > EnergyMomentumRelationAllowance)
      theDynamicalMass = dynamicalMass;
    else
      theDynamicalMass = 0.0;
  }
}

G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                                     const G4ThreeVector& aParticleMomentum)
  : thePolarization(0.0, 0.0, 0.0),
    theParticleDefinition(aParticleDefinition),
    theDynamicalMass(aParticleDefinition->GetPDGMass()),
    theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
    theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
    theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment())
{
  SetMomentum(aParticleMomentum);  // 3-dim momentum is given
}

G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                                     const G4LorentzVector& aParticleMomentum)
  : thePolarization(0.0, 0.0, 0.0),
    theParticleDefinition(aParticleDefinition),
    theDynamicalMass(aParticleDefinition->GetPDGMass()),
    theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
    theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
    theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment())
{
  Set4Momentum(aParticleMomentum);  // 4-momentum vector (Lorentz vector)
}

G4DynamicParticle::G4DynamicParticle(const G4ParticleDefinition* aParticleDefinition,
                                     G4double totalEnergy, const G4ThreeVector& aParticleMomentum)
  : thePolarization(0.0, 0.0, 0.0),
    theParticleDefinition(aParticleDefinition),
    theDynamicalMass(aParticleDefinition->GetPDGMass()),
    theDynamicalCharge(aParticleDefinition->GetPDGCharge()),
    theDynamicalSpin(aParticleDefinition->GetPDGSpin()),
    theDynamicalMagneticMoment(aParticleDefinition->GetPDGMagneticMoment())
{
  // total energy and 3-dim momentum are given
  G4double pModule2 = aParticleMomentum.mag2();
  if (pModule2 > 0.0) {
    G4double mass2 = totalEnergy * totalEnergy - pModule2;
    G4double PDGmass2 = (aParticleDefinition->GetPDGMass()) * (aParticleDefinition->GetPDGMass());
    SetMomentumDirection(aParticleMomentum.unit());
    if (mass2 < EnergyMRA2) {
      theDynamicalMass = 0.;
      SetKineticEnergy(totalEnergy);
    }
    else {
      if (std::abs(PDGmass2 - mass2) > EnergyMRA2) {
        theDynamicalMass = std::sqrt(mass2);
        SetKineticEnergy(totalEnergy - theDynamicalMass);
      }
      else {
        SetKineticEnergy(totalEnergy - theDynamicalMass);
      }
    }
  }
  else {
    SetMomentumDirection(1.0, 0.0, 0.0);
    SetKineticEnergy(0.0);
  }
}

G4DynamicParticle::G4DynamicParticle(const G4DynamicParticle& right)
  : theMomentumDirection(right.theMomentumDirection),
    thePolarization(right.thePolarization),
    theParticleDefinition(right.theParticleDefinition),
    // Don't copy preassignedDecayProducts
    primaryParticle(right.primaryParticle),
    theKineticEnergy(right.theKineticEnergy),
    theLogKineticEnergy(right.theLogKineticEnergy),
    theBeta(right.theBeta),
    theProperTime(right.theProperTime),
    theDynamicalMass(right.theDynamicalMass),
    theDynamicalCharge(right.theDynamicalCharge),
    theDynamicalSpin(right.theDynamicalSpin),
    theDynamicalMagneticMoment(right.theDynamicalMagneticMoment),

    verboseLevel(right.verboseLevel),
    thePDGcode(right.thePDGcode)
{
  if (right.theElectronOccupancy != nullptr) {
    theElectronOccupancy = new G4ElectronOccupancy(*right.theElectronOccupancy);
  }
}

G4DynamicParticle::G4DynamicParticle(G4DynamicParticle&& from)
  : theMomentumDirection(from.theMomentumDirection),
    thePolarization(from.thePolarization),
    theParticleDefinition(from.theParticleDefinition),
    theElectronOccupancy(from.theElectronOccupancy),
    // Don't move preassignedDecayProducts
    primaryParticle(from.primaryParticle),
    theKineticEnergy(from.theKineticEnergy),
    theLogKineticEnergy(from.theLogKineticEnergy),
    theBeta(from.theBeta),
    theProperTime(from.theProperTime),
    theDynamicalMass(from.theDynamicalMass),
    theDynamicalCharge(from.theDynamicalCharge),
    theDynamicalSpin(from.theDynamicalSpin),
    theDynamicalMagneticMoment(from.theDynamicalMagneticMoment),

    verboseLevel(from.verboseLevel),
    thePDGcode(from.thePDGcode)
{
  // Release the data from the source object
  from.theParticleDefinition = nullptr;
  from.theElectronOccupancy = nullptr;
  from.thePreAssignedDecayProducts = nullptr;
  from.primaryParticle = nullptr;
}

G4DynamicParticle::~G4DynamicParticle()
{
  delete thePreAssignedDecayProducts;
  thePreAssignedDecayProducts = nullptr;

  delete theElectronOccupancy;
  theElectronOccupancy = nullptr;
}

G4DynamicParticle& G4DynamicParticle::operator=(const G4DynamicParticle& right)
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

    delete theElectronOccupancy;
    if (right.theElectronOccupancy != nullptr) {
      theElectronOccupancy = new G4ElectronOccupancy(*right.theElectronOccupancy);
    }
    else {
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

G4DynamicParticle& G4DynamicParticle::operator=(G4DynamicParticle&& from)
{
  if (this != &from) {
    theMomentumDirection = from.theMomentumDirection;
    thePolarization = from.thePolarization;
    theKineticEnergy = from.theKineticEnergy;
    theProperTime = from.theProperTime;

    theDynamicalMass = from.theDynamicalMass;
    theDynamicalCharge = from.theDynamicalCharge;
    theDynamicalSpin = from.theDynamicalSpin;
    theDynamicalMagneticMoment = from.theDynamicalMagneticMoment;

    delete theElectronOccupancy;
    theElectronOccupancy = from.theElectronOccupancy;
    from.theElectronOccupancy = nullptr;

    // thePreAssignedDecayProducts must not be moved.
    thePreAssignedDecayProducts = nullptr;
    from.thePreAssignedDecayProducts = nullptr;
    thePreAssignedDecayTime = -1.0;

    theParticleDefinition = from.theParticleDefinition;
    from.theParticleDefinition = nullptr;

    verboseLevel = from.verboseLevel;

    primaryParticle = from.primaryParticle;
    from.primaryParticle = nullptr;

    thePDGcode = from.thePDGcode;
  }
  return *this;
}

void G4DynamicParticle::SetDefinition(const G4ParticleDefinition* aParticleDefinition)
{
  // remove preassigned decay
  if (thePreAssignedDecayProducts != nullptr) {
#ifdef G4VERBOSE
    if (verboseLevel > 0) {
      G4cout << " G4DynamicParticle::SetDefinition()::"
             << "!!! Pre-assigned decay products is attached !!!! " << G4endl;
      G4cout << "!!! New Definition is " << aParticleDefinition->GetParticleName() << " !!! "
             << G4endl;
      G4cout << "!!! Pre-assigned decay products will be deleted !!!! " << G4endl;
    }
#endif
    delete thePreAssignedDecayProducts;
  }
  thePreAssignedDecayProducts = nullptr;

  theParticleDefinition = aParticleDefinition;

  // set Dynamic mass/charge
  SetMass(theParticleDefinition->GetPDGMass());
  theDynamicalCharge = theParticleDefinition->GetPDGCharge();
  theDynamicalSpin = theParticleDefinition->GetPDGSpin();
  theDynamicalMagneticMoment = theParticleDefinition->GetPDGMagneticMoment();

  // Set electron orbits
  if (theElectronOccupancy != nullptr) {
    delete theElectronOccupancy;
    theElectronOccupancy = nullptr;
  }
}

G4bool G4DynamicParticle::operator==(const G4DynamicParticle& right) const
{
  return (this == (G4DynamicParticle*)&right);
}

G4bool G4DynamicParticle::operator!=(const G4DynamicParticle& right) const
{
  return (this != (G4DynamicParticle*)&right);
}

void G4DynamicParticle::AllocateElectronOccupancy()
{
  if (G4IonTable::IsIon(theParticleDefinition)) {
    // Only ions can have ElectronOccupancy
    theElectronOccupancy = new G4ElectronOccupancy();
  }
  else {
    theElectronOccupancy = nullptr;
  }
}

void G4DynamicParticle::SetMomentum(const G4ThreeVector& momentum)
{
  G4double pModule2 = momentum.mag2();
  if (pModule2 > 0.0) {
    const G4double mass = theDynamicalMass;
    SetMomentumDirection(momentum.unit());
    SetKineticEnergy(pModule2 / (std::sqrt(pModule2 + mass * mass) + mass));
  }
  else {
    SetMomentumDirection(1.0, 0.0, 0.0);
    SetKineticEnergy(0.0);
  }
}

void G4DynamicParticle::Set4Momentum(const G4LorentzVector& momentum)
{
  G4double pModule2 = momentum.vect().mag2();
  if (pModule2 > 0.0) {
    SetMomentumDirection(momentum.vect().unit());
    const G4double totalenergy = momentum.t();
    const G4double mass2 = totalenergy * totalenergy - pModule2;
    const G4double PDGmass2 =
      (theParticleDefinition->GetPDGMass()) * (theParticleDefinition->GetPDGMass());
    if (mass2 < EnergyMRA2) {
      theDynamicalMass = 0.;
    }
    else if (std::abs(PDGmass2 - mass2) > EnergyMRA2) {
      theDynamicalMass = std::sqrt(mass2);
    }
    SetKineticEnergy(totalenergy - theDynamicalMass);
  }
  else {
    SetMomentumDirection(1.0, 0.0, 0.0);
    SetKineticEnergy(0.0);
  }
}

#ifdef G4VERBOSE
void G4DynamicParticle::DumpInfo(G4int mode) const
{
  if (theParticleDefinition == nullptr) {
    G4cout << " G4DynamicParticle::DumpInfo() - Particle type not defined !!! " << G4endl;
  }
  else {
    G4cout << " Particle type - " << theParticleDefinition->GetParticleName() << G4endl
           << "   mass:        " << GetMass() / CLHEP::GeV << "[GeV]" << G4endl
           << "   charge:      " << GetCharge() / CLHEP::eplus << "[e]" << G4endl
           << "   Direction x: " << GetMomentumDirection().x()
           << ", y: " << GetMomentumDirection().y() << ", z: " << GetMomentumDirection().z()
           << G4endl << "   Total Momentum = " << GetTotalMomentum() / CLHEP::GeV << "[GeV]"
           << G4endl << "   Momentum: " << GetMomentum().x() / CLHEP::GeV << "[GeV]"
           << ", y: " << GetMomentum().y() / CLHEP::GeV << "[GeV]"
           << ", z: " << GetMomentum().z() / CLHEP::GeV << "[GeV]" << G4endl
           << "   Total Energy   = " << GetTotalEnergy() / CLHEP::GeV << "[GeV]" << G4endl
           << "   Kinetic Energy = " << GetKineticEnergy() / CLHEP::GeV << "[GeV]" << G4endl
           << " MagneticMoment  [MeV/T]: " << GetMagneticMoment() / CLHEP::MeV * CLHEP::tesla
           << G4endl << "   ProperTime     = " << GetProperTime() / CLHEP::ns << "[ns]" << G4endl;

    if (mode > 0) {
      if (theElectronOccupancy != nullptr) {
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

G4double G4DynamicParticle::GetElectronMass() const
{
  return CLHEP::electron_mass_c2;
}
