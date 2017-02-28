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
// Contact: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 
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

#ifndef G4Molecule_h
#define G4Molecule_h 1

#include "G4IT.hh"
#include "G4Allocator.hh"
#include "G4MoleculeDefinition.hh"

class G4Molecule;
template<>
G4KDNode<G4Molecule>::~G4KDNode();

class G4Molecule;
class G4MolecularConfiguration;
class G4MoleculeDefinition;
class G4MolecularDissociationChannel;
class G4DynamicParticle;
class G4Material;

G4Molecule* GetMolecule(const G4Track& track);
G4Molecule* GetMolecule(const G4Track* track);

/** Class Description
 *  The dynamic molecule holds all the data that change for a molecule
 *  It has a pointer to G4MoleculeDefinition object, which holds
 *   all the "ground level" information.
 */

class G4Molecule : public G4IT
{

public:
  // With Description

ITDef(G4Molecule)

  //From G4VUserTrackInformation
  void Print() const;

  //  new/delete operators are overloded to use G4Allocator
  inline void *operator new(size_t);
#ifdef __IBMCPP__
  inline void *operator new(size_t sz, void* p)
  { return p;}
#endif
  inline void operator delete(void*);

  G4Molecule(const G4Molecule&);
  G4Molecule & operator=(const G4Molecule &right);
  G4bool operator==(const G4Molecule &right) const;
  G4bool operator!=(const G4Molecule &right) const;
  G4bool operator<(const G4Molecule &right) const;

  operator int() const
  {
    return GetMoleculeID();
  }

  virtual G4ITType GetITSubType() const
  {
    return GetMoleculeID();
  }

public:
  //------ Constructors --------------------------
  /** To build a molecule at ground state according to a given
   * G4MoleculeDefinition that can be obtained from G4GenericMoleculeManager
   */
  G4Molecule(G4MoleculeDefinition * molecule);

  G4Molecule(G4MoleculeDefinition* molDef, int charge);

  /** To build a molecule at a specific excitation/ionisation state according
   * to a ground state that can be obtained from G4GenericMoleculeManager
   */
  G4Molecule(G4MoleculeDefinition * molecule, G4int, G4int);

  /** Specific builder for water molecules to be used in Geant4-DNA,
   * the last option Excitation is true if the molecule is excited, is
   * false is the molecule is ionized.
   */
  G4Molecule(G4MoleculeDefinition * molecule, G4int, G4bool);

  G4Molecule(G4MolecularConfiguration*);

  virtual ~G4Molecule();

  //-------- Methods -------------------------------
  //Get from static definition
  /** Returns the name of the molecule
   */
  const G4String& GetName() const;

  /** Returns the formated name of the molecule
   */
  const G4String& GetFormatedName() const;

  /** Returns the nomber of atoms compouning the molecule
   */
  G4int GetAtomsNumber() const;

  /** Will set up the correct molecularConfiguration given
   * an electron configuration
   */
  void SetElectronOccupancy(const G4ElectronOccupancy*);

  /** Method used in Geant4-DNA to excite water molecules
   */
  void ExciteMolecule(G4int);

  /** Method used in Geant4-DNA to ionize water molecules
   */
  void IonizeMolecule(G4int);

  /** Add n electrons to a given orbit.
   * Note : You can add as many electrons to a given orbit, the result
   * may be unrealist.
   */
  void AddElectron(G4int orbit, G4int n = 1);

  /** Remove n electrons to a given orbit.
   */
  void RemoveElectron(G4int, G4int number = 1);

  /** Move one electron from an orbit to another.
   */
  void MoveOneElectron(G4int /*orbit*/, G4int /*orbit*/);

  /** Returns the number of electron.
   */
  G4double GetNbElectrons() const; //This method can be used to check if the electron s number is physical

  /** Show the electronic state of the molecule.
   */
  void PrintState() const;

  G4Track * BuildTrack(G4double globalTime, const G4ThreeVector& Position);

  G4double GetKineticEnergy() const;

  G4double GetDiffusionVelocity() const;

  const std::vector<const G4MolecularDissociationChannel*>* GetDecayChannel() const;

  G4int GetFakeParticleID() const;
  G4int GetMoleculeID() const;

  //-------------Inline functions ---------------------
  /**  Get molecule definition. This G4MoleculeDefinition has the ground
   * electronic state of the molecule.
   */
  const G4MoleculeDefinition* GetDefinition() const;

  //methods to set/get changing parameters

/////////////////////////////////////////////////////////////////////////////
  /** Sets the diffusion coefficient D of the molecule used in diffusion
   * processes to calculate the mean square jump distance between two
   * changes of direction. In three dimension : <x^2> = 6 D t where t is
   * the mean jump time between two changes of direction.
   */
  void SetDiffusionCoefficient(G4double);

  /** Returns the diffusion coefficient D.
   */
  G4double GetDiffusionCoefficient() const;

  /** Returns the diffusion coefficient D.
   */
  G4double GetDiffusionCoefficient(const G4Material*,
                                   double temperature) const;

  /** Set the decay time of the molecule.
   */
  void SetDecayTime(G4double);

  /** Returns the decay time of the molecule.
   */
  G4double GetDecayTime() const;

  /** The Van Der Valls Radius of the molecule
   */
  void SetVanDerVaalsRadius(G4double);
  G4double GetVanDerVaalsRadius() const;

  /** Returns the object ElectronOccupancy describing the electronic
   * configuration of the molecule.
   */
  const G4ElectronOccupancy* GetElectronOccupancy() const;

  /** Returns the charge of molecule.
   */
  G4int GetCharge() const;

  /** Set the total mass of the molecule.
   */
  void SetMass(G4double);

  /** Returns the total mass of the molecule.
   */
  G4double GetMass() const;

  /** Returns the label of the molecule configuration
   */
  const G4String& GetLabel() const;

  void SetLabel(const G4String& label);

  void ChangeConfigurationToLabel(const G4String& label);

////////////////////////////////////////////////////////////////////////

  G4MolecularConfiguration* GetMolecularConfiguration() const;

//  static void SetGlobalTemperature(G4double);
//  static G4double GetGlobalTemperature();

  static G4Molecule* GetMolecule(const G4Track*);

private:
  /** Default molecule builder
   */
  G4Molecule();

  G4MolecularConfiguration* fpMolecularConfiguration;

//  double fCachedDiffusionCoefficient;

//  static /*G4ThreadLocal*/double fgTemperature;
};

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4ThreadLocal G4Allocator<G4Molecule> *aMoleculeAllocator;
#else
extern G4DLLIMPORT G4ThreadLocal G4Allocator<G4Molecule> *aMoleculeAllocator;
#endif

//////////////////////////
inline void * G4Molecule::operator new(size_t)
//////////////////////////
{
  if (!aMoleculeAllocator) aMoleculeAllocator = new G4Allocator<G4Molecule>;
  return (void *) aMoleculeAllocator->MallocSingle();
}

//////////////////////////
inline void G4Molecule::operator delete(void * aMolecule)
//////////////////////////
{
  aMoleculeAllocator->FreeSingle((G4Molecule *) aMolecule);
}

#endif
