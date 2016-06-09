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
// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Alfonso Mantero 6 Mar 2009
// ----------------------------------------------------------------
//


#ifndef G4MolecularDecayChannel_h
#define G4MolecularDecayChannel_h 1

#include <vector>
#include <map>
#include "G4VMolecularDecayDisplacer.hh"
#include "G4MoleculeHandleManager.hh"

// ######################################################################
// ###                  MolecularDecayChannel                         ###
// ######################################################################

class G4Molecule;

struct CompMoleculePointer;

class G4MolecularDecayChannel
{
    // Class Description
    // This is where are stored and can be retrieved data of one decay channel
    // of excited or ionized molecules: products. energy and probability.

public: //With Description

    G4MolecularDecayChannel();
    G4MolecularDecayChannel(G4String);
    ~G4MolecularDecayChannel();
    G4MolecularDecayChannel(const G4MolecularDecayChannel&);
    G4MolecularDecayChannel & operator=(const G4MolecularDecayChannel &right);

public:

    //Root Mean Square radial distance thermalisation of a product
    G4double GetRMSRadialDisplacementOfProduct(const G4Molecule*);

    // methods to construct decay channels "interactively"

    void  AddProduct(const G4Molecule*,G4double = 0);
    inline void  SetName(const G4String&);
    inline void  SetEnergy(G4double);
    inline void  SetProbability(G4double);
    inline void  SetDecayTime(G4double);
    inline void  SetRMSMotherMoleculeDisplacement(G4double);
    inline void  SetDisplacementType(DisplacementType);

    // get methods to retrieve data

    inline const G4String& GetName() const;
    G4int GetNbProducts() const;
    const G4Molecule* GetProduct(int) const;
    inline const std::vector<G4double>& GetRMSProductsDisplacement() const;
    inline G4double GetEnergy() const;
    inline G4double GetProbability() const;
    inline G4double GetDecayTime() const;
    inline G4double GetRMSMotherMoleculeDisplacement() const;
    inline DisplacementType GetDisplacementType() const;

private:

    DisplacementType fDisplacementType;
    G4String fName;
    std::vector<G4MoleculeHandle>* fProductsVector;
    G4double fReleasedEnergy;
    G4double fProbability;
    G4double fDecayTime; // To be taken into account in the next releases

    //Root Mean Square radial distance jump of a excited/ionised MotherMolecule molecule
    G4double fRMSMotherMoleculeDisplacement;
    std::vector<G4double> fRMSProductsDisplacementVector;
};

inline void G4MolecularDecayChannel::SetName(const G4String& value)
{
    fName = value;
}

inline void G4MolecularDecayChannel::SetEnergy(G4double value)
{
    fReleasedEnergy = value;
}


inline void G4MolecularDecayChannel::SetProbability(G4double value)
{
    fProbability = value;
}

inline void G4MolecularDecayChannel::SetDecayTime(G4double value)
{

    fDecayTime = value;
}

inline void G4MolecularDecayChannel::SetRMSMotherMoleculeDisplacement(G4double value)
{
    fRMSMotherMoleculeDisplacement = value;
}

inline const G4String& G4MolecularDecayChannel::GetName() const
{
    return fName;
}

inline const std::vector<G4double>& G4MolecularDecayChannel::GetRMSProductsDisplacement() const
{
    return fRMSProductsDisplacementVector;
}

inline G4double G4MolecularDecayChannel::GetEnergy() const
{

    return fReleasedEnergy;
}

inline G4double G4MolecularDecayChannel::GetProbability() const
{
    return fProbability;
}

inline G4double G4MolecularDecayChannel::GetDecayTime() const
{
    return fDecayTime;
}

inline G4double G4MolecularDecayChannel::GetRMSMotherMoleculeDisplacement() const
{
    return fRMSMotherMoleculeDisplacement;
}

inline void G4MolecularDecayChannel::SetDisplacementType(DisplacementType aDisplacementType)
{
    fDisplacementType = aDisplacementType;
}

inline DisplacementType G4MolecularDecayChannel::GetDisplacementType() const
{
    return fDisplacementType;
}
#endif








