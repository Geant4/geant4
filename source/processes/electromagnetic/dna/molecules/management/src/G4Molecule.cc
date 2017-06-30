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
// $Id: G4Molecule.cc 103042 2017-03-10 11:50:07Z gcosmo $
//
// ---------------------------------------------------------------------
//	GEANT 4 class header file
//
//	History: first implementation, based on G4DynamicParticle
//           New dependency : G4VUserTrackInformation
//
//      ---------------- G4Molecule  ----------------
//      first design&implementation by Alfonso Mantero, 7 Apr 2009
//      New developments Alfonso Mantero & Mathieu Karamitros
//      Oct/Nov 2009 Class Name changed to G4Molecule
//                   Removed dependency from G4DynamicParticle
//                   New constructors :
//                    copy constructor
//                    direct ionized/excited molecule
//                   New methods :
//                    Get : name,atoms' number,nb electrons,decayChannel
//                    PrintState //To get the electronic level and the
//                                 corresponding name of the excitation
//                    Kinematic :
//                    BuildTrack,GetKineticEnergy,GetDiffusionVelocity
//                    Change the way dynCharge and eNb is calculated
// ---------------------------------------------------------------------

#include "G4Molecule.hh"
#include "G4MolecularConfiguration.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VMoleculeCounter.hh"

using namespace std;

///*G4ThreadLocal*/G4double G4Molecule::fgTemperature = 310; // 310*kelvin;
// 37°C, used to shoot an energy

ITImp(G4Molecule)

G4ThreadLocal G4Allocator<G4Molecule> *aMoleculeAllocator = 0;

//______________________________________________________________________________

template<>
G4KDNode<G4Molecule>::~G4KDNode(){
  fPoint->SetNode(nullptr);
}

//______________________________________________________________________________

G4Molecule* GetMolecule(const G4Track& track)
{
  return (G4Molecule*) (GetIT(track));
}

//______________________________________________________________________________

G4Molecule* GetMolecule(const G4Track* track)
{
  return (G4Molecule*) (GetIT(track));
}

//______________________________________________________________________________

G4Molecule* G4Molecule::GetMolecule(const G4Track* track)
{
  return (G4Molecule*) (GetIT(track));
}

//______________________________________________________________________________

void G4Molecule::Print() const
{
  G4cout << "The user track information is a molecule" << G4endl;
}

//______________________________________________________________________________

G4Molecule::G4Molecule(const G4Molecule& right) :
    G4VUserTrackInformation("G4Molecule"), G4IT(right)
{
  fpMolecularConfiguration = right.fpMolecularConfiguration;
}

//______________________________________________________________________________

G4Molecule& G4Molecule::operator=(const G4Molecule& right)
{
  if (&right == this) return *this;
  fpMolecularConfiguration = right.fpMolecularConfiguration;
  return *this;
}

//______________________________________________________________________________

G4bool G4Molecule::operator==(const G4Molecule& right) const
{
  if (fpMolecularConfiguration == right.fpMolecularConfiguration)
  {
    return true;
  }
  return false;
}

//______________________________________________________________________________

G4bool G4Molecule::operator!=(const G4Molecule& right) const
{
  return !(*this == right);
}

//______________________________________________________________________________
/**	The two methods below are the most called of the simulation :
 *	compare molecules in the MoleculeStackManager or in
 *  the InteractionTable
 */

G4bool G4Molecule::operator<(const G4Molecule& right) const
{
  return fpMolecularConfiguration < right.fpMolecularConfiguration;
}

//______________________________________________________________________________
/** Default molecule builder
 */
G4Molecule::G4Molecule() :
    G4VUserTrackInformation("G4Molecule"), G4IT()

{
  fpMolecularConfiguration = 0;
}

//______________________________________________________________________________

G4Molecule::~G4Molecule()
{
  if (fpTrack != nullptr)
  {
    if (G4VMoleculeCounter::Instance()->InUse())
    {
      G4VMoleculeCounter::Instance()->
        RemoveAMoleculeAtTime(fpMolecularConfiguration,
                              fpTrack->GetGlobalTime(),
                              &(fpTrack->GetPosition()));
    }
    fpTrack = nullptr;
  }
  fpMolecularConfiguration = nullptr;
}

//______________________________________________________________________________
/** Build a molecule at ground state according to a given
 *  G4MoleculeDefinition that can be obtained from G4GenericMoleculeManager
 */
//////////////////////////
G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition) :
    G4VUserTrackInformation("G4Molecule"), G4IT()
//////////////////////////
{
  fpMolecularConfiguration =
      G4MolecularConfiguration::
        GetOrCreateMolecularConfiguration(moleculeDefinition);
}

//______________________________________________________________________________

G4Molecule::G4Molecule(G4MoleculeDefinition* moleculeDefinition, int charge)
{
  fpMolecularConfiguration =
      G4MolecularConfiguration::
        GetOrCreateMolecularConfiguration(moleculeDefinition,
                                          charge);
}

//______________________________________________________________________________
/** Build a molecule at a specific excitation/ionisation state according
 *  to a ground state that can be obtained from G4GenericMoleculeManager.
 *  Put 0 in the second option if this is a ionisation.
 */
//////////////////////////
G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition,
                       G4int OrbitalToFree,
                       G4int OrbitalToFill) :
    G4VUserTrackInformation("G4Molecule"), G4IT()
//////////////////////////
{
  if (moleculeDefinition->GetGroundStateElectronOccupancy())
  {
    G4ElectronOccupancy dynElectronOccupancy(
        *moleculeDefinition->GetGroundStateElectronOccupancy());

    if (OrbitalToFill != 0)
    {
      dynElectronOccupancy.RemoveElectron(OrbitalToFree - 1, 1);
      dynElectronOccupancy.AddElectron(OrbitalToFill - 1, 1);
      // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    if (OrbitalToFill == 0)
    {
      dynElectronOccupancy.RemoveElectron(OrbitalToFree - 1, 1);
      // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    fpMolecularConfiguration =
        G4MolecularConfiguration::GetOrCreateMolecularConfiguration(
            moleculeDefinition, dynElectronOccupancy);
  }
  else
  {
    fpMolecularConfiguration = 0;
    G4Exception(
        "G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition, "
        "G4int OrbitalToFree, G4int OrbitalToFill)",
        "G4Molecule_wrong_usage_of_constructor",
        FatalErrorInArgument,
        "If you want to use this constructor, the molecule definition has to be "
        "first defined with electron occupancies");
  }
}

//______________________________________________________________________________
/** Specific builder for water molecules to be used in Geant4-DNA,
 * the last option Excitation is true if the molecule is excited, is
 * false is the molecule is ionized.
 */

G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition,
                       G4int Level,
                       G4bool Excitation) :
    G4VUserTrackInformation("G4Molecule"), G4IT()
{
  if (moleculeDefinition->GetGroundStateElectronOccupancy())
  {
    G4ElectronOccupancy dynElectronOccupancy(
        *moleculeDefinition->GetGroundStateElectronOccupancy());

    if (Excitation == true)
    {
      dynElectronOccupancy.RemoveElectron(Level, 1);
      dynElectronOccupancy.AddElectron(5, 1);
      // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    if (Excitation == false)
    {
      dynElectronOccupancy.RemoveElectron(Level, 1);
      // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    fpMolecularConfiguration =
        G4MolecularConfiguration::GetOrCreateMolecularConfiguration(
            moleculeDefinition, dynElectronOccupancy);
  }
  else
  {
    fpMolecularConfiguration = 0;
    G4Exception(
        "G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition, "
        "G4int OrbitalToFree, G4int OrbitalToFill)",
        "G4Molecule_wrong_usage_of_constructor",
        FatalErrorInArgument,
        "If you want to use this constructor, the molecule definition has to be "
        "first defined with electron occupancies");

  }
}

//______________________________________________________________________________

G4Molecule::G4Molecule(G4MolecularConfiguration* molConf)
{
  fpMolecularConfiguration = molConf;
}

//______________________________________________________________________________

void G4Molecule::SetElectronOccupancy(const G4ElectronOccupancy* occ)
{
  fpMolecularConfiguration =
      G4MolecularConfiguration::GetOrCreateMolecularConfiguration(
          fpMolecularConfiguration->GetDefinition(), *occ);
}

//______________________________________________________________________________

/** Method used in Geant4-DNA to excite water molecules
 */
void G4Molecule::ExciteMolecule(G4int ExcitedLevel)
{
  fpMolecularConfiguration = fpMolecularConfiguration->ExciteMolecule(
      ExcitedLevel);
}

//______________________________________________________________________________

/** Method used in Geant4-DNA to ionize water molecules
 */
void G4Molecule::IonizeMolecule(G4int IonizedLevel)
{
  fpMolecularConfiguration = fpMolecularConfiguration->IonizeMolecule(
      IonizedLevel);
}

//______________________________________________________________________________

void G4Molecule::AddElectron(G4int orbit, G4int number)
{
  fpMolecularConfiguration = fpMolecularConfiguration->AddElectron(orbit,
                                                                   number);
}

//______________________________________________________________________________

void G4Molecule::RemoveElectron(G4int orbit, G4int number)
{
  fpMolecularConfiguration = fpMolecularConfiguration->RemoveElectron(orbit,
                                                                      number);
}

//______________________________________________________________________________

void G4Molecule::MoveOneElectron(G4int orbitToFree, G4int orbitToFill)
{
  fpMolecularConfiguration = fpMolecularConfiguration->MoveOneElectron(
      orbitToFree, orbitToFill);
}

//______________________________________________________________________________

const G4String& G4Molecule::GetName() const
{
  return fpMolecularConfiguration->GetName();
}

//______________________________________________________________________________

const G4String& G4Molecule::GetFormatedName() const
{
  return fpMolecularConfiguration->GetFormatedName();
}

//______________________________________________________________________________

G4int G4Molecule::GetAtomsNumber() const
{
  return fpMolecularConfiguration->GetAtomsNumber();
}

//______________________________________________________________________________

G4double G4Molecule::GetNbElectrons() const
{
  return fpMolecularConfiguration->GetNbElectrons();
}

//______________________________________________________________________________

void G4Molecule::PrintState() const
{
  fpMolecularConfiguration->PrintState();
}

//______________________________________________________________________________

G4Track * G4Molecule::BuildTrack(G4double globalTime,
                                 const G4ThreeVector& position)
{
  if (fpTrack != 0){
    G4Exception("G4Molecule::BuildTrack", "Molecule001", FatalErrorInArgument,
                "A track was already assigned to this molecule");
  }

  // Kinetic Values
  // Set a random direction to the molecule
  G4double costheta = (2 * G4UniformRand()-1);
  G4double theta = acos(costheta);
  G4double phi = 2 * pi * G4UniformRand();

  G4double xMomentum = cos(phi) * sin(theta);
  G4double yMomentum = sin(theta) * sin(phi);
  G4double zMomentum = costheta;

  G4ThreeVector MomentumDirection(xMomentum, yMomentum, zMomentum);
  G4double KineticEnergy = GetKineticEnergy();

  G4DynamicParticle* dynamicParticle = new G4DynamicParticle(
      fpMolecularConfiguration->GetDefinition(), MomentumDirection,
      KineticEnergy);

  if (G4VMoleculeCounter::InUse()){
    G4VMoleculeCounter::Instance()->
      AddAMoleculeAtTime(fpMolecularConfiguration,
                         globalTime,
                         &(fpTrack->GetPosition()));
  }

  //Set the Track
  fpTrack = new G4Track(dynamicParticle, globalTime, position);
  fpTrack->SetUserInformation(this);

  return fpTrack;
}

//______________________________________________________________________________

G4double G4Molecule::GetKineticEnergy() const
{
  ////
  // Ideal Gaz case
  double v = GetDiffusionVelocity();
  double E = (fpMolecularConfiguration->GetMass() / (c_squared)) * (v * v) / 2.;
  ////
  return E;
}

//______________________________________________________________________________

G4double G4Molecule::GetDiffusionVelocity() const
{
  double moleculeMass = fpMolecularConfiguration->GetMass() / (c_squared);

  ////
  // Different possibilities
  ////
  // Ideal Gaz case : Maxwell Boltzmann Distribution
  //    double sigma = k_Boltzmann * fgTemperature / mass;
  //    return G4RandGauss::shoot( 0, sigma );
  ////
  // Ideal Gaz case : mean velocity from equipartition theorem
  return sqrt(3 * k_Boltzmann *
              G4MolecularConfiguration::GetGlobalTemperature()/ moleculeMass);
  ////
  // Using this approximation for liquid is wrong
  // However the brownian process avoid taking
  // care of energy consideration and plays only
  // with positions
}

//______________________________________________________________________________

// added - to be transformed in a "Decay method"
const vector<const G4MolecularDissociationChannel*>*
G4Molecule::GetDecayChannel() const
{
  return fpMolecularConfiguration->GetDecayChannel();
}

//______________________________________________________________________________

G4int G4Molecule::GetFakeParticleID() const
{
  return fpMolecularConfiguration->GetFakeParticleID();
}

//______________________________________________________________________________

G4int G4Molecule::GetMoleculeID() const
{
  return fpMolecularConfiguration->GetMoleculeID();
}

//______________________________________________________________________________

void G4Molecule::SetDecayTime(G4double dynDecayTime)
{
  fpMolecularConfiguration->SetDecayTime(dynDecayTime);
}

//______________________________________________________________________________

G4double G4Molecule::GetDecayTime() const
{
  return fpMolecularConfiguration->GetDecayTime();
}

//______________________________________________________________________________

void G4Molecule::SetVanDerVaalsRadius(G4double dynVanDerVaalsRadius)
{
  fpMolecularConfiguration->SetVanDerVaalsRadius(dynVanDerVaalsRadius);
}

//______________________________________________________________________________

G4double G4Molecule::GetVanDerVaalsRadius() const
{
  return fpMolecularConfiguration->GetVanDerVaalsRadius();
}

//______________________________________________________________________________

G4int G4Molecule::GetCharge() const
{
  return fpMolecularConfiguration->GetCharge();
}

//______________________________________________________________________________

void G4Molecule::SetMass(G4double aMass)
{
  fpMolecularConfiguration->SetMass(aMass);
}

//______________________________________________________________________________

G4double G4Molecule::GetMass() const
{
  return fpMolecularConfiguration->GetMass();
}

//______________________________________________________________________________

const G4ElectronOccupancy* G4Molecule::GetElectronOccupancy() const
{
  return fpMolecularConfiguration->GetElectronOccupancy();
}

//______________________________________________________________________________

const G4MoleculeDefinition* G4Molecule::GetDefinition() const
{
  return fpMolecularConfiguration->GetDefinition();
}

//______________________________________________________________________________

void G4Molecule::SetDiffusionCoefficient(G4double dynDiffusionCoefficient)
{
  fpMolecularConfiguration->SetDiffusionCoefficient(dynDiffusionCoefficient);
}

//______________________________________________________________________________

G4double G4Molecule::GetDiffusionCoefficient() const
{
  return fpMolecularConfiguration->GetDiffusionCoefficient();
}

//______________________________________________________________________________

G4double G4Molecule::GetDiffusionCoefficient(const G4Material* mat,
                                             double temperature) const
{
  return fpMolecularConfiguration->GetDiffusionCoefficient(mat,
                                                           temperature);
}

//______________________________________________________________________________

G4MolecularConfiguration* G4Molecule::GetMolecularConfiguration() const
{
  return fpMolecularConfiguration;
}

//______________________________________________________________________________

//void G4Molecule::SetGlobalTemperature(G4double temperature)
//{
//  fgTemperature = temperature;
//}
//
////______________________________________________________________________________
//
//G4double G4Molecule::GetGlobalTemperature()
//{
//  return fgTemperature;
//}

//______________________________________________________________________________

const G4String& G4Molecule::GetLabel() const
{
  return fpMolecularConfiguration->GetLabel();
}

//______________________________________________________________________________

void G4Molecule::SetLabel(const G4String& label)
{
  fpMolecularConfiguration->SetLabel(label);
}

//______________________________________________________________________________

void G4Molecule::ChangeConfigurationToLabel(const G4String& label)
{
  // TODO check fpMolecularConfiguration already exists
  // and new one as well
  // TODO notify for stack change
  fpMolecularConfiguration =
      G4MolecularConfiguration::GetMolecularConfiguration(
          fpMolecularConfiguration->GetDefinition(), label);

  assert(fpMolecularConfiguration!=0);
}
