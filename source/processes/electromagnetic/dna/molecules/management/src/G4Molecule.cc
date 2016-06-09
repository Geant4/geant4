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
#include "G4Track.hh"

using namespace std;

double G4Molecule::fgTemperature = 310*kelvin;
// 37°C, used to shoot an energy

ITImp(G4Molecule)

G4Allocator<G4Molecule> aMoleculeAllocator;

G4Molecule* GetMolecule(const G4Track& track)
{
    return (G4Molecule*)(GetIT(track));
}

G4Molecule* GetMolecule(const G4Track* track)
{
    return (G4Molecule*)(GetIT(track));
}

void G4Molecule::Print() const
{
    G4cout<<"The user track information is a molecule"<<G4endl;
}

G4Molecule::G4Molecule(const G4Molecule& right) :
    G4VUserTrackInformation("G4Molecule"), G4IT(right)
{
    Init();
    fMolecularConfiguration = right . fMolecularConfiguration;
}

G4Molecule& G4Molecule::operator=(const G4Molecule& right)
{
    if (&right==this) return *this;
    Init();
    fMolecularConfiguration = right . fMolecularConfiguration;
    return *this;
}

G4bool G4Molecule::operator==(const G4Molecule& right) const
{
    if(fMolecularConfiguration==right.fMolecularConfiguration)
    {
        return true;
    }
    return false;
}

G4bool G4Molecule::operator!=(const G4Molecule& right) const
{
    return !(*this == right);
}

////////////////////////////////////////////////////////////////////////
///	The two methods below are the most called of the simulation :
///	compare molecules in the MoleculeStackManager or in
///	the InteractionTable

G4bool G4Molecule::operator<(const G4Molecule& right) const
{
    return fMolecularConfiguration < right.fMolecularConfiguration ;
}
////////////////////////////////////////////////////////////////////////
void G4Molecule::Init()
{
    fMolecularConfiguration = 0 ;
    fDynamicParticle = 0;
}

////////////////////////////////////////////////////////////////////////
/** Default molecule builder
 */
//////////////////////////
G4Molecule::G4Molecule() : G4VUserTrackInformation("G4Molecule"), G4IT()
//////////////////////////
{
    Init();
}

//////////////////////////
G4Molecule::~G4Molecule()
//////////////////////////
{
    if(fTrack!=NULL)
    {
        fTrack = 0;
    }
    fMolecularConfiguration = 0;
    fDynamicParticle = 0;
    // DEBUG
    // G4cout<<"Molecule killed"<<G4endl;
}

/** Build a molecule at ground state according to a given
 *  G4MoleculeDefinition that can be obtained from G4GenericMoleculeManager
 */
//////////////////////////
G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition) :
    G4VUserTrackInformation("G4Molecule"), G4IT()
//////////////////////////
{
    Init();
    fMolecularConfiguration = G4MolecularConfiguration::GetMolecularConfiguration(moleculeDefinition);
}

/** Build a molecule at a specific excitation/ionisation state according
 *  to a ground state that can be obtained from G4GenericMoleculeManager.
 *  Put 0 in the second option if this is a ionisation.
 */
//////////////////////////
G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition, G4int OrbitalToFree, G4int OrbitalToFill):
    G4VUserTrackInformation("G4Molecule"), G4IT()
//////////////////////////
{
    Init();

    G4ElectronOccupancy dynElectronOccupancy (*moleculeDefinition->GetGroundStateElectronOccupancy());

    if (OrbitalToFill != 0)
    {
        dynElectronOccupancy.RemoveElectron(OrbitalToFree-1,1);
        dynElectronOccupancy.AddElectron(OrbitalToFill-1,1);
        // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    if (OrbitalToFill == 0)
    {
        dynElectronOccupancy.RemoveElectron(OrbitalToFree-1,1);
        // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    fMolecularConfiguration = G4MolecularConfiguration::GetMolecularConfiguration(moleculeDefinition, dynElectronOccupancy);
}

/** Specific builder for water molecules to be used in Geant4-DNA,
 * the last option Excitation is true if the molecule is excited, is
 * false is the molecule is ionized.
 */

G4Molecule::G4Molecule(G4MoleculeDefinition * moleculeDefinition, G4int Level, G4bool Excitation):
    G4VUserTrackInformation("G4Molecule"), G4IT()
{
    Init();

    G4ElectronOccupancy dynElectronOccupancy (*moleculeDefinition->GetGroundStateElectronOccupancy());

    if (Excitation == true)
    {
        dynElectronOccupancy.RemoveElectron(Level,1);
        dynElectronOccupancy.AddElectron(5,1);
        // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    if (Excitation == false)
    {
        dynElectronOccupancy.RemoveElectron(Level,1);
        // dynElectronOccupancy.DumpInfo(); // DEBUG
    }

    fMolecularConfiguration = G4MolecularConfiguration::GetMolecularConfiguration(moleculeDefinition, dynElectronOccupancy);
}

void G4Molecule::SetElectronOccupancy(const G4ElectronOccupancy* occ)
{
    fMolecularConfiguration = G4MolecularConfiguration::GetMolecularConfiguration(fMolecularConfiguration->GetDefinition(), *occ);
}

/** Method used in Geant4-DNA to excite water molecules
 */
void G4Molecule::ExciteMolecule(G4int ExcitedLevel)
{
    fMolecularConfiguration = fMolecularConfiguration->ExciteMolecule(ExcitedLevel);
}

/** Method used in Geant4-DNA to ionize water molecules
 */
void G4Molecule::IonizeMolecule(G4int IonizedLevel)
{
    fMolecularConfiguration = fMolecularConfiguration->IonizeMolecule(IonizedLevel);
}

void G4Molecule::AddElectron(G4int orbit, G4int number)
{
    fMolecularConfiguration = fMolecularConfiguration->AddElectron(orbit,number);
}

void G4Molecule::RemoveElectron(G4int orbit,G4int number)
{
    fMolecularConfiguration = fMolecularConfiguration->RemoveElectron(orbit,number);
}

void G4Molecule::MoveOneElectron(G4int orbitToFree,G4int orbitToFill)
{
    fMolecularConfiguration = fMolecularConfiguration->MoveOneElectron(orbitToFree,orbitToFill);
}

const G4String& G4Molecule::GetName() const
{
    return fMolecularConfiguration->GetName();
}

G4int G4Molecule::GetAtomsNumber() const
{
    return fMolecularConfiguration->GetAtomsNumber();
}

G4double G4Molecule::GetNbElectrons() const
{
    return fMolecularConfiguration->GetNbElectrons();
}

void G4Molecule::PrintState() const
{
    fMolecularConfiguration->PrintState();
}

G4Track * G4Molecule::BuildTrack(G4double globalTime, const G4ThreeVector& Position)
{
    if(fTrack != 0)
    {
        G4Exception("G4Molecule::BuildTrack","Molecule001",
                    FatalErrorInArgument,"A track was already assigned to this molecule");
    }

    // Kinetic Values
    // Set a random direction to the molecule
    G4double costheta = (2*G4UniformRand()-1);
    G4double theta = acos (costheta);
    G4double phi = 2*pi*G4UniformRand();

    G4double xMomentum = cos(phi)* sin(theta);
    G4double yMomentum = sin(theta)*sin(phi);
    G4double zMomentum = costheta;

    G4ThreeVector MomentumDirection(xMomentum, yMomentum, zMomentum);
    G4double KineticEnergy = GetKineticEnergy();
    // G4cout << " **** KineticEnergy : " << KineticEnergy << G4endl;
    fDynamicParticle = new G4DynamicParticle(fMolecularConfiguration->GetDefinition(),
            MomentumDirection,
            KineticEnergy);

    //Set the Track
    fTrack = new G4Track(fDynamicParticle, globalTime, Position);
    fTrack -> SetUserInformation (this);

    return fTrack;
}

G4double G4Molecule::GetKineticEnergy() const
{
    ////
    // Ideal Gaz case
    double v = GetDiffusionVelocity();
    double E = (fMolecularConfiguration->GetMass()/(c_squared))*(v*v)/2.;
    ////
    return E;
}

G4double G4Molecule::GetDiffusionVelocity() const
{
    double m = fMolecularConfiguration->GetMass()/(c_squared);

    ////
    // Different possibilities
    ////
    // Ideal Gaz case : Maxwell Boltzmann Distribution
    //    double sigma = k_Boltzmann * fgTemperature / m;
    //    return G4RandGauss::shoot( 0, sigma );
    ////
    // Ideal Gaz case : mean velocity from equipartition theorem
    return sqrt(3*k_Boltzmann*fgTemperature/m);
    ////
    // Using this approximation for liquid is wrong
    // However the brownian process avoid taking
    // care of energy consideration and plays only
    // with positions
}

// added - to be transformed in a "Decay method"
const vector <const G4MolecularDecayChannel*>* G4Molecule::GetDecayChannel() const
{
    return fMolecularConfiguration->GetDecayChannel();
}

G4int G4Molecule::GetMoleculeID() const
{
    return fMolecularConfiguration->GetMoleculeID();
}

void G4Molecule::SetDecayTime(G4double dynDecayTime)
{
    fMolecularConfiguration->SetDecayTime(dynDecayTime);
}

G4double G4Molecule::GetDecayTime() const
{
    return fMolecularConfiguration->GetDecayTime();
}

void G4Molecule::SetVanDerVaalsRadius(G4double dynVanDerVaalsRadius)
{
    fMolecularConfiguration->SetVanDerVaalsRadius(dynVanDerVaalsRadius);
}

G4double G4Molecule::GetVanDerVaalsRadius() const
{
    return fMolecularConfiguration->GetVanDerVaalsRadius();
}

G4int G4Molecule::GetCharge() const
{
    return fMolecularConfiguration->GetCharge() ;
}

void G4Molecule::SetMass(G4double aMass)
{
    fMolecularConfiguration->SetMass(aMass);
}

G4double G4Molecule::GetMass() const
{
    return fMolecularConfiguration->GetMass();
}

G4ElectronOccupancy G4Molecule::GetElectronOccupancy() const
{
    return *(fMolecularConfiguration->GetElectronOccupancy());
}

const G4MoleculeDefinition* G4Molecule::GetDefinition() const
{
    return fMolecularConfiguration->GetDefinition();
}

void G4Molecule::SetDiffusionCoefficient(G4double dynDiffusionCoefficient)
{
    fMolecularConfiguration->SetDiffusionCoefficient(dynDiffusionCoefficient);
}

G4double G4Molecule::GetDiffusionCoefficient() const
{
    return fMolecularConfiguration->GetDiffusionCoefficient();
}
