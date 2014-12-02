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
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

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

#ifndef G4MolecularConfiguration_
#define G4MolecularConfiguration_ 1

#include <vector>
#include <map>
#include <CLHEP/Utility/memory.h>
#include "G4Threading.hh"

class G4MolecularDissociationChannel;
class G4MoleculeDefinition;

#include "G4ElectronOccupancy.hh"

struct comparator
{
  bool operator()(const G4ElectronOccupancy& occ1,
                  const G4ElectronOccupancy& occ2) const
  {
    // Since this method is called a lot of time,
    // we retrieve only once the totOcc
    G4int totalOcc1 = occ1.GetTotalOccupancy();
    G4int totalOcc2 = occ2.GetTotalOccupancy();
    if (totalOcc1 != totalOcc2)
    {
      return totalOcc1 < totalOcc2;
    }
    else
    {
      G4int occupancy1 = -1;
      G4int occupancy2 = -1;
      const G4int sizeOrbit = occ1.GetSizeOfOrbit();
      for (G4int i = 0; i < occ1.GetSizeOfOrbit();)
      {
        // Since this method is called a lot of time,
        // we retrieve only once the Occ

        occupancy1 = occ1.GetOccupancy(i);
        occupancy2 = occ2.GetOccupancy(i);

        if (occupancy1 != occupancy2)
        {
          return occupancy1 < occupancy2;
        }
        else
        {
          i++;
          if (i >= sizeOrbit) return false;
        }
      }
    }
    return false;
  }
};

/** The pointer G4MolecularConfiguration will be shared by all the
 * molecules having the same molecule definition and the same
 * electron occupancy
 * BE CAREFUlL !!! : If you change the mass for instance of a OH^-,
 * this will affect all the OH^- molecule diffusing around
 */
class G4MolecularConfiguration
{
public:

  /////////////////
  // Static methods

  // Get ground state electronic configuration
  static G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition*);

  // Get for a given moleculeDefinition and a given electronic configuration, the mol conf
  static G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition*,
                                                             const G4ElectronOccupancy& electronOccupancy);

  // Get for a given moleculeDefinition and a given electronic configuration, the mol conf
  static G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition*,
                                                             int charge);

  // Release memory of the mol conf manager
  static void DeleteManager();
  ///////////////

  // Methods
  const G4MoleculeDefinition* GetDefinition() const;

  /** Returns the name of the molecule
   */
  const G4String& GetName() const;

  /** Returns the formated name of the molecule
   */
  const G4String& GetFormatedName() const;

  /** Returns the nomber of atoms compouning the molecule
   */
  G4int GetAtomsNumber() const;

  /** Method used in Geant4-DNA to excite water molecules
   */
  G4MolecularConfiguration* ExciteMolecule(G4int);

  /** Method used in Geant4-DNA to ionize water molecules
   */
  G4MolecularConfiguration* IonizeMolecule(G4int);

  /** Add n electrons to a given orbit.
   * Note : You can add as many electrons to a given orbit, the result
   * may be unrealist.
   */
  G4MolecularConfiguration* AddElectron(G4int orbit, G4int n = 1);

  /** Remove n electrons to a given orbit.
   */
  G4MolecularConfiguration* RemoveElectron(G4int, G4int number = 1);

  /** Move one electron from an orbit to another.
   */
  G4MolecularConfiguration* MoveOneElectron(G4int /*orbit*/, G4int /*orbit*/);

  /** Returns the number of electron.
   */
  G4double GetNbElectrons() const;

  /** Display the electronic state of the molecule.
   */
  void PrintState() const;

  const std::vector<const G4MolecularDissociationChannel*>* GetDecayChannel() const;

  G4int GetFakeParticleID() const;

  inline G4int GetMoleculeID() const;

  /** Sets the diffusion coefficient D of the molecule used in diffusion
   * processes to calculate the mean square jump distance between two
   * changes of direction. In three dimension : <x^2> = 6 D t where t is
   * the mean jump time between two changes of direction.
   *
   * Note : Diffusion Coefficient in one medium only
   * For the time being, we will consider only one diffusion
   * coefficient for the all simulation => diffusion in one medium only
   * If the user needs to use the diffusion in different materials,
   * she/he should contact the developpers/mainteners of this package
   */
  inline void SetDiffusionCoefficient(G4double);

  /** Returns the diffusion coefficient D.
   */
  inline G4double GetDiffusionCoefficient() const;

  /** Set the decay time of the molecule.
   */
  inline void SetDecayTime(G4double);

  /** Returns the decay time of the molecule.
   */
  inline G4double GetDecayTime() const;

  /** The Van Der Valls Radius of the molecule
   */
  inline void SetVanDerVaalsRadius(G4double);
  inline G4double GetVanDerVaalsRadius() const;

  /** Returns the object ElectronOccupancy describing the electronic
   * configuration of the molecule.
   */
  inline const G4ElectronOccupancy* GetElectronOccupancy() const;

  /** Returns the charge of molecule.
   */
  inline G4int GetCharge() const;

  /** Set the total mass of the molecule.
   */
  inline void SetMass(G4double);

  /** Returns the total mass of the molecule.
   */
  inline G4double GetMass() const;

protected:
  G4MolecularConfiguration(const G4MoleculeDefinition*,
                           const G4ElectronOccupancy&);
  G4MolecularConfiguration(const G4MoleculeDefinition*, int);
  G4MolecularConfiguration(const G4MolecularConfiguration&);
  G4MolecularConfiguration & operator=(G4MolecularConfiguration &right);
  ~G4MolecularConfiguration();
  G4MolecularConfiguration* ChangeConfiguration(const G4ElectronOccupancy& newElectronOccupancy);
  G4MolecularConfiguration* ChangeConfiguration(int charge);

  const G4MoleculeDefinition* fMoleculeDefinition;
  const G4ElectronOccupancy* fElectronOccupancy;

  void CheckElectronOccupancy(const char* line) const;

public:
  class G4MolecularConfigurationManager
  {
  public:
    G4MolecularConfigurationManager() :
        fMoleculeCreationMutex()
    {
      fLastMoleculeID = -1;
    }
    ~G4MolecularConfigurationManager();

    G4int SetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                    const G4ElectronOccupancy& eOcc,
                                    G4MolecularConfiguration* molConf);

    G4int SetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                    int charge,
                                    G4MolecularConfiguration* molConf);

    const G4ElectronOccupancy* FindCommonElectronOccupancy(const G4MoleculeDefinition* molDef,
                                                           const G4ElectronOccupancy& eOcc);

    G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                                        const G4ElectronOccupancy& eOcc);

    G4MolecularConfiguration* GetMolecularConfiguration(const G4MoleculeDefinition* molDef,
                                                        int charge);

    static G4Mutex fManagerCreationMutex;

    void RemoveMolecularConfigurationFromTable(G4MolecularConfiguration*);

  private:
    typedef std::map<const G4MoleculeDefinition*,
        std::map<G4ElectronOccupancy, G4MolecularConfiguration*, comparator> > MolecularConfigurationTable;
    MolecularConfigurationTable fTable;
    typedef std::map<const G4MoleculeDefinition*,
        std::map<int, G4MolecularConfiguration*, comparator> > MolChargeConfigurationTable;
    MolChargeConfigurationTable fChargeTable;
    G4int fLastMoleculeID;

    G4Mutex fMoleculeCreationMutex;
  };

protected:
  static G4MolecularConfigurationManager* fgManager;
  static G4MolecularConfigurationManager* GetManager();

  G4double fDynDiffusionCoefficient;
  G4double fDynVanDerVaalsRadius;
  G4double fDynDecayTime;
  G4double fDynMass;
  G4int fDynCharge;
  G4int fMoleculeID;
  /*mutable*/ G4String fFormatedName;
  /*mutable*/ G4String fName;
};

inline const G4MoleculeDefinition* G4MolecularConfiguration::GetDefinition() const
{
  return fMoleculeDefinition;
}

inline const G4ElectronOccupancy* G4MolecularConfiguration::GetElectronOccupancy() const
{
  return fElectronOccupancy;
}

inline void G4MolecularConfiguration::SetDiffusionCoefficient(G4double dynDiffusionCoefficient)
{
  fDynDiffusionCoefficient = dynDiffusionCoefficient;
}

inline G4double G4MolecularConfiguration::GetDiffusionCoefficient() const
{
  return fDynDiffusionCoefficient;
}

inline void G4MolecularConfiguration::SetDecayTime(G4double dynDecayTime)
{
  fDynDecayTime = dynDecayTime;
}

inline G4double G4MolecularConfiguration::GetDecayTime() const
{
  return fDynDecayTime;
}

inline void G4MolecularConfiguration::SetVanDerVaalsRadius(G4double dynVanDerVaalsRadius)
{
  fDynVanDerVaalsRadius = dynVanDerVaalsRadius;
}

inline G4double G4MolecularConfiguration::GetVanDerVaalsRadius() const
{
  return fDynVanDerVaalsRadius;
}

inline G4int G4MolecularConfiguration::GetCharge() const
{
  return fDynCharge;
}

inline void G4MolecularConfiguration::SetMass(G4double aMass)
{
  fDynMass = aMass;
}

inline G4double G4MolecularConfiguration::GetMass() const
{
  return fDynMass;
}

inline G4int G4MolecularConfiguration::GetMoleculeID() const
{
  return fMoleculeID;
}
#endif
