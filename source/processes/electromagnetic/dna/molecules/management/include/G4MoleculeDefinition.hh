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
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      21 Oct 2009 first implementation by A. Mantero and M.Karamitros
//                  Based on prototype of A.Mantero
// **********************************************************************
//
// Author: Mathieu Karamitros

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

#ifndef G4MoleculeDefinition_h
#define G4MoleculeDefinition_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "G4ElectronOccupancy.hh"
#include "G4MolecularDissociationTable.hh"
#include "G4MolecularDissociationChannel.hh"
#include "G4FakeParticleID.hh"

class G4MolecularDissociationChannel;
class G4MolecularDissociationTable;
class G4MolecularConfiguration;

// -----------------------------------------------------------------------------
// ###                          MoleculeDefinition                           ###
// -----------------------------------------------------------------------------

class G4MoleculeDefinition: public G4ParticleDefinition
{
public:
  G4MoleculeDefinition(const G4String& name,
                       G4double mass,
                       G4double diffCoeff,
                       G4int charge = 0,
                       G4int electronicLevels = 0,
                       G4double radius = -1,
                       G4int atomsNumber = -1,
                       G4double lifetime = -1,
                       G4String aType = "",
                       G4FakeParticleID ID = G4FakeParticleID::Create());

  virtual ~G4MoleculeDefinition();

  // Set the electronic configuration at ground level
  void SetLevelOccupation(G4int,
                          G4int eNb = 2);
  // set the occupation(0(def), 1 or 2) of the level specified
  //(levels numbering starts from 0)

  //methods to set/get diffusion properties
  inline void SetDiffusionCoefficient(G4double);
  inline G4double GetDiffusionCoefficient() const;

  inline void SetAtomsNumber(G4int);
  inline G4int GetAtomsNumber() const;

  inline void SetVanDerVaalsRadius(G4double);
  inline G4double GetVanDerVaalsRadius() const;

  //____________________________________________________________________________
  // Create more molecular configurations for this molecule definition
  // Other ways : through G4MolecularTable

  // Note: the userID of the created molecule configuration will be:
  //       MoleculeDefinationName_excitedStateLabel
  G4MolecularConfiguration*
  NewConfiguration(const G4String& excitedStateLabel);

  G4MolecularConfiguration*
  NewConfigurationWithElectronOccupancy(const G4String& excitedStateLabel,
                                        const G4ElectronOccupancy&,
                                        double decayTime = 0.);

  G4MolecularConfiguration*
  GetConfigurationWithLabel(const G4String& molecularConfLabel);

  //____________________________________________________________________________
  // Build the decay table
  // Version 1

  void AddDecayChannel(const G4MolecularConfiguration* molConf,
                       const G4MolecularDissociationChannel* channel);

  // Version 2

  void AddDecayChannel(const G4String& molecularConfLabel,
                       const G4MolecularDissociationChannel* channel);

  //____________________________________________________________________________
  // "Get" methods related to decay

  const std::vector<const G4MolecularDissociationChannel*>*
    GetDecayChannels(const G4MolecularConfiguration*) const;
  const std::vector<const G4MolecularDissociationChannel*>*
    GetDecayChannels(const G4String&) const;

  inline const G4MolecularDissociationTable* GetDecayTable() const;
  inline G4MolecularDissociationTable* GetDecayTable();
  inline G4double GetDecayTime() const;

  //____________________________________________________________________________
  // General "Get" methods
  inline const G4ElectronOccupancy* GetGroundStateElectronOccupancy() const;
  inline G4int GetCharge() const;
  inline const G4String& GetName() const;
  inline G4double GetMass() const;
  inline const G4String& GetType() const;
  inline G4int GetNbElectrons() const;
  inline G4int GetNbMolecularShells() const;

  inline const G4String& GetFormatedName() const
  {
    return fFormatedName;
  }

  //____________________________________________________________________________
  inline void SetFormatedName(const G4String& name)
  {
    fFormatedName = name;
  }

  void Finalize();

  static G4MoleculeDefinition* Load(std::istream&);
  void Serialize(std::ostream&);

protected:
  G4MoleculeDefinition();
  G4MoleculeDefinition(const G4MoleculeDefinition&);

private:
  const G4MoleculeDefinition & operator=(const G4MoleculeDefinition &right);

private:
  G4int fCharge;

  // Diffusion Coefficient in one medium only
  // Note : For the time being, we will consider only one diffusion
  // coefficient for the all simulation => diffusion in one medium only
  // If the user needs to use the diffusion in different materials,
  // she/he should contact the developers/maintainers of this package
  G4double fDiffusionCoefficient;

  G4int fAtomsNb;
  G4double fVanDerVaalsRadius;

  G4String fFormatedName;

  G4ElectronOccupancy* fElectronOccupancy;
  G4MolecularDissociationTable* fDecayTable;
};

inline void G4MoleculeDefinition::SetDiffusionCoefficient(G4double value)
{
  fDiffusionCoefficient = value;
}

inline G4double G4MoleculeDefinition::GetDiffusionCoefficient() const
{
  return fDiffusionCoefficient;
}

inline G4int G4MoleculeDefinition::GetCharge() const
{
  return fCharge;
}

inline G4double G4MoleculeDefinition::GetDecayTime() const
{
  return GetPDGLifeTime();
}

inline void G4MoleculeDefinition::SetAtomsNumber(G4int val)
{
  fAtomsNb = val;
}

inline G4int G4MoleculeDefinition::GetAtomsNumber() const
{
  return fAtomsNb;
}

inline void G4MoleculeDefinition::SetVanDerVaalsRadius(G4double val)
{
  fVanDerVaalsRadius = val;
}

inline G4double G4MoleculeDefinition::GetVanDerVaalsRadius() const
{
  return fVanDerVaalsRadius;
}

inline const G4ElectronOccupancy* G4MoleculeDefinition::GetGroundStateElectronOccupancy() const
{
  return fElectronOccupancy;
}

inline const G4String& G4MoleculeDefinition::GetName() const
{

  return GetParticleName();
}

inline G4double G4MoleculeDefinition::GetMass() const
{
  return GetPDGMass();
}

inline const G4String& G4MoleculeDefinition::GetType() const
{
  return GetParticleSubType();
}

inline G4int G4MoleculeDefinition::GetNbElectrons() const
{
  if (fElectronOccupancy)
  {
    return fElectronOccupancy->GetTotalOccupancy();
  }

  return 0;
  //    return fNbOfElectrons;
}

inline G4int G4MoleculeDefinition::GetNbMolecularShells() const
{
  if (fElectronOccupancy)
  {
    return fElectronOccupancy->GetSizeOfOrbit();
  }

  return 0;
}

inline const G4MolecularDissociationTable* G4MoleculeDefinition::GetDecayTable() const
{
  return fDecayTable;
}

inline G4MolecularDissociationTable* G4MoleculeDefinition::GetDecayTable()
{
  return fDecayTable;
}
#endif

