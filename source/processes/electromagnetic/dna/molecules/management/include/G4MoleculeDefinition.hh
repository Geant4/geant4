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


#ifndef G4MoleculeDefinition_h
#define G4MoleculeDefinition_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "G4ElectronOccupancy.hh"
#include "G4MolecularDecayTable.hh"
#include "G4MolecularDecayChannel.hh"
#include "G4MoleculeID.hh"

// ######################################################################
// ###                          Molecule                              ###
// ######################################################################

class G4MolecularDecayChannel;
class G4MolecularDecayTable;
class G4MolecularConfiguration;

class G4MoleculeDefinition : public G4ParticleDefinition
{
// Class Description
//  This is the base class for all molecules
//  All molecules created are objects of this class.

public: //With Description
    G4MoleculeDefinition(const G4String& name,
                         G4double mass,
                         G4int    electronsNumber,
                         G4int    electronicLevels,
                         G4double diffCoeff,
                         G4int atomsNumber = -1,
                         G4double radius = -1,
                         G4double lifetime = -1,
                         G4String aType = "",
                         G4MoleculeID ID = G4MoleculeID::Create()
                         );

    virtual ~G4MoleculeDefinition();

    // Set the electronic configuration at ground level
    void SetLevelOccupation(G4int, G4int eNb=2);
    // set the occupation(0(def), 1 or 2) of the level specified
    //(levels numbering starts from 0)

    //methods to set/get diffusion properties
    inline void SetDiffusionCoefficient(G4double);
    inline G4double GetDiffusionCoefficient() const;

    inline void SetAtomsNumber(G4int);
    inline G4int GetAtomsNumber() const;

    inline void SetVanDerVaalsRadius(G4double);
    inline G4double GetVanDerVaalsRadius() const;

    //°°°°°°°°°°°°°°°°°°°°°°°°
    // Build the decay table
    void AddExcitedState(const G4String&) ;
    const G4String& GetExcitedState(const G4ElectronOccupancy*) const;
    void AddDecayChannel(const G4String&, const G4MolecularDecayChannel*) ;
    void AddeConfToExcitedState(const G4String&,const G4ElectronOccupancy&, double decayTime = 0.);

    //°°°°°°°°°°°°°°°°°°°°°°°°
    // "Get" methods related to decay
    const std::vector<const G4MolecularDecayChannel*>* GetDecayChannels(const G4ElectronOccupancy*) const;
    const std::vector<const G4MolecularDecayChannel*>* GetDecayChannels(const G4String&) const;

    //°°°°°°°°°°°°°°°°°°°°°°°°
    // General "Get" methods
    inline const G4ElectronOccupancy*   GetGroundStateElectronOccupancy() const;
    inline const G4String&              GetName() const;
    inline G4double                     GetMass() const;
    inline const G4String&              GetType() const;
    inline G4int                        GetNbElectrons() const;
    inline G4int                        GetNbMolecularShells() const;
    inline const G4MolecularDecayTable* GetDecayTable() const ;
    inline G4MolecularDecayTable*       GetDecayTable() ;
    inline G4double                     GetDecayTime() const;

protected :
    G4MoleculeDefinition();
    G4MoleculeDefinition(const G4MoleculeDefinition&);

private :
    const G4MoleculeDefinition & operator=(const G4MoleculeDefinition &right);

private:
    G4double fMass;
    G4String fType;

    G4int fNbOfElectrons;
    G4int fNbOfMolecularShells;

    // Diffusion Coefficient in one medium only
    // Note : For the time being, we will consider only one diffusion
    // coefficient for the all simulation => diffusion in one medium only
    // If the user needs to use the diffusion in different medium,
    // he should contact the developpers/mainteners of this package
    G4double fDiffusionCoefficient;

    G4int fAtomsNb;
    G4double fVanDerVaalsRadius;

    G4ElectronOccupancy* fElectronOccupancy;
    G4MolecularDecayTable* fDecayTable;
};


inline void G4MoleculeDefinition::SetDiffusionCoefficient(G4double value)
{
    fDiffusionCoefficient = value;
}

inline G4double G4MoleculeDefinition::GetDiffusionCoefficient() const
{
    return fDiffusionCoefficient;
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
    return fMass;
}

inline const G4String& G4MoleculeDefinition::GetType() const
{

    return GetParticleSubType();
}

inline G4int G4MoleculeDefinition::GetNbElectrons() const
{

    return fNbOfElectrons;
}

inline G4int G4MoleculeDefinition::GetNbMolecularShells() const
{

    return fNbOfMolecularShells;
}

inline const G4MolecularDecayTable* G4MoleculeDefinition::GetDecayTable() const
{
    return fDecayTable;
}

inline G4MolecularDecayTable* G4MoleculeDefinition::GetDecayTable()
{
    return fDecayTable;
}
#endif








