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
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation by Alfonso Mantero 4 Mar 2009
//
// **********************************************************************

#include "G4MolecularDissociationChannel.hh"
#include "G4VMolecularDissociationDisplacer.hh"

using namespace std;

//______________________________________________________________________________

G4MolecularDissociationChannel::G4MolecularDissociationChannel(const G4String& aName)
    : G4MolecularDissociationChannel()
{
    fName = aName;
}

//______________________________________________________________________________

G4MolecularDissociationChannel::G4MolecularDissociationChannel()
    : fDisplacementType(G4VMolecularDissociationDisplacer::NoDisplacement)
    , fReleasedEnergy(0.)
    , fProbability(0.)
    , fDecayTime(0.)
    , fRMSMotherMoleculeDisplacement(0.)
{
}

//______________________________________________________________________________

void G4MolecularDissociationChannel::AddProduct(Product* pProduct,
                                                G4double displacement)
{
    fProductsVector.push_back(pProduct);
    fRMSProductsDisplacementVector.push_back(displacement);
}

//______________________________________________________________________________

G4int G4MolecularDissociationChannel::GetNbProducts() const
{
    return (G4int)fProductsVector.size();
}

//______________________________________________________________________________

G4MolecularDissociationChannel::Product* G4MolecularDissociationChannel::GetProduct(int index) const
{
    return fProductsVector[index];
}

//______________________________________________________________________________

G4double
G4MolecularDissociationChannel::GetRMSRadialDisplacementOfProduct(Product* pProduct)
{
    if (fProductsVector.empty())
    {
        return -1.;
    }

    auto it = std::find_if(fProductsVector.begin(), fProductsVector.end(),
                           [pProduct](Product* _pProduct) {
                               return _pProduct == pProduct;
                           });

    if (it == fProductsVector.end())
    {
        return -1.;
    }
    auto index = std::distance(fProductsVector.begin(), it);
    return fRMSProductsDisplacementVector[index];
}

