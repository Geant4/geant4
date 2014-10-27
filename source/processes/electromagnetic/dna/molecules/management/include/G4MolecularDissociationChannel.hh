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
// we ask that you please cite the following reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 
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
#include "G4VMolecularDissociationDisplacer.hh"
#include "G4MoleculeHandleManager.hh"

// ######################################################################
// ###                  MolecularDecayChannel                         ###
// ######################################################################

class G4Molecule;

struct CompMoleculePointer;

class G4MolecularDissociationChannel
{
    // Class Description
    // This is where are stored and can be retrieved data of one decay channel
    // of excited or ionized molecules: products. energy and probability.

public: //With Description

    G4MolecularDissociationChannel();
    G4MolecularDissociationChannel(G4String);
    ~G4MolecularDissociationChannel();
    G4MolecularDissociationChannel(const G4MolecularDissociationChannel&);
    G4MolecularDissociationChannel & operator=(const G4MolecularDissociationChannel &right);

public:

    //Root Mean Square radial distance thermalisation of a product
    G4double GetRMSRadialDisplacementOfProduct(const G4Molecule*);

    // methods to construct decay channels "interactively"

    void AddProduct(const G4Molecule*,G4double = 0);
    void AddProduct(const G4String& molecule, G4double displacement = 0);

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
//    std::vector<G4MoleculeHandle>* fProductsVector;
    std::vector<const G4Molecule*>* fProductsVector;
    G4double fReleasedEnergy;
    G4double fProbability;
    G4double fDecayTime; // To be taken into account in the next releases

    //Root Mean Square radial distance jump of a excited/ionised MotherMolecule molecule
    G4double fRMSMotherMoleculeDisplacement;
    std::vector<G4double> fRMSProductsDisplacementVector;
};

inline void G4MolecularDissociationChannel::SetName(const G4String& value)
{
    fName = value;
}

inline void G4MolecularDissociationChannel::SetEnergy(G4double value)
{
    fReleasedEnergy = value;
}


inline void G4MolecularDissociationChannel::SetProbability(G4double value)
{
    fProbability = value;
}

inline void G4MolecularDissociationChannel::SetDecayTime(G4double value)
{

    fDecayTime = value;
}

inline void G4MolecularDissociationChannel::SetRMSMotherMoleculeDisplacement(G4double value)
{
    fRMSMotherMoleculeDisplacement = value;
}

inline const G4String& G4MolecularDissociationChannel::GetName() const
{
    return fName;
}

inline const std::vector<G4double>& G4MolecularDissociationChannel::GetRMSProductsDisplacement() const
{
    return fRMSProductsDisplacementVector;
}

inline G4double G4MolecularDissociationChannel::GetEnergy() const
{

    return fReleasedEnergy;
}

inline G4double G4MolecularDissociationChannel::GetProbability() const
{
    return fProbability;
}

inline G4double G4MolecularDissociationChannel::GetDecayTime() const
{
    return fDecayTime;
}

inline G4double G4MolecularDissociationChannel::GetRMSMotherMoleculeDisplacement() const
{
    return fRMSMotherMoleculeDisplacement;
}

inline void G4MolecularDissociationChannel::SetDisplacementType(DisplacementType aDisplacementType)
{
    fDisplacementType = aDisplacementType;
}

inline DisplacementType G4MolecularDissociationChannel::GetDisplacementType() const
{
    return fDisplacementType;
}
#endif








