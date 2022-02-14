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

#pragma once

#include <vector>
#include <map>
#include "G4VMolecularDissociationDisplacer.hh"

class G4Molecule;
class G4MolecularConfiguration;

class G4MolecularDissociationChannel
{
public:
    G4MolecularDissociationChannel();
    explicit G4MolecularDissociationChannel(const G4String&);
    ~G4MolecularDissociationChannel() = default;
    G4MolecularDissociationChannel(const G4MolecularDissociationChannel&) = default;

    G4MolecularDissociationChannel&
    operator=(const G4MolecularDissociationChannel& right) = default;

    using Product = const G4MolecularConfiguration;
    using ProductList = std::vector<Product*>;

public:

    void AddProduct(Product*, G4double displacement = 0.);

    inline void SetName(const G4String&);
    inline void SetEnergy(G4double);
    inline void SetProbability(G4double);
    inline void SetDecayTime(G4double);
    inline void SetRMSMotherMoleculeDisplacement(G4double);
    inline void SetDisplacementType(DisplacementType);


    inline const G4String& GetName() const;
    G4int GetNbProducts() const;
    Product* GetProduct(int) const;
    inline const std::vector<G4double>& GetRMSProductsDisplacement() const;
    inline G4double GetEnergy() const;
    inline G4double GetProbability() const;
    inline G4double GetDecayTime() const;
    inline G4double GetRMSMotherMoleculeDisplacement() const;
    inline DisplacementType GetDisplacementType() const;

    G4double GetRMSRadialDisplacementOfProduct(Product*);

private:
    DisplacementType fDisplacementType;
    G4String fName;
    ProductList fProductsVector;
    G4double fReleasedEnergy;
    G4double fProbability;
    G4double fDecayTime;
    // Root Mean Square radial distance jump of the mother molecule
    G4double fRMSMotherMoleculeDisplacement;
    std::vector<G4double> fRMSProductsDisplacementVector;
};

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::SetName(const G4String& value)
{
    fName = value;
}

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::SetEnergy(G4double value)
{
    fReleasedEnergy = value;
}

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::SetProbability(G4double value)
{
    fProbability = value;
}

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::SetDecayTime(G4double value)
{

    fDecayTime = value;
}

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::
SetRMSMotherMoleculeDisplacement(G4double value)
{
    fRMSMotherMoleculeDisplacement = value;
}

//______________________________________________________________________________

inline const G4String& G4MolecularDissociationChannel::GetName() const
{
    return fName;
}

//______________________________________________________________________________

inline const std::vector<G4double>&
G4MolecularDissociationChannel::GetRMSProductsDisplacement() const
{
    return fRMSProductsDisplacementVector;
}

//______________________________________________________________________________

inline G4double G4MolecularDissociationChannel::GetEnergy() const
{
    return fReleasedEnergy;
}

//______________________________________________________________________________

inline G4double G4MolecularDissociationChannel::GetProbability() const
{
    return fProbability;
}

//______________________________________________________________________________

inline G4double G4MolecularDissociationChannel::GetDecayTime() const
{
    return fDecayTime;
}

//______________________________________________________________________________

inline G4double G4MolecularDissociationChannel::
GetRMSMotherMoleculeDisplacement() const
{
    return fRMSMotherMoleculeDisplacement;
}

//______________________________________________________________________________

inline void G4MolecularDissociationChannel::
SetDisplacementType(DisplacementType aDisplacementType)
{
    fDisplacementType = aDisplacementType;
}

//______________________________________________________________________________

inline DisplacementType G4MolecularDissociationChannel::
GetDisplacementType() const
{
    return fDisplacementType;
}
