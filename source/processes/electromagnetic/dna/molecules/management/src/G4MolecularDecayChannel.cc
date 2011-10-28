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
//
//
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation by Alfonso Mantero 4 Mar 2009
//
// **********************************************************************

#include "G4MolecularDecayChannel.hh"
#include "G4Molecule.hh"
#include "G4MoleculeHandleManager.hh"

using namespace std;

// ######################################################################
// ###                      MolecularDecayChannel                     ###
// ######################################################################

struct CompMoleculePointer
{
    bool operator()(const G4Molecule* mol1, const G4Molecule* mol2) const
    {
        return (*mol1) < (*mol2);
    }
};

// constructors and destructor

G4MolecularDecayChannel::G4MolecularDecayChannel(G4String aName) : fName(aName)
{
    fProductsVector = 0;
    fProductsVector = 0;
    fReleasedEnergy = 0;
}

G4MolecularDecayChannel::G4MolecularDecayChannel()
{
    fProductsVector = 0;
    fProductsVector = 0;
    fReleasedEnergy = 0;
}

G4MolecularDecayChannel::~G4MolecularDecayChannel()
{
    if(fProductsVector)
    {
        fProductsVector->clear();
        delete fProductsVector;
    }
}

G4MolecularDecayChannel::G4MolecularDecayChannel(const G4MolecularDecayChannel& right)
{
    *this = right ;
}

G4MolecularDecayChannel& G4MolecularDecayChannel::operator=
(const G4MolecularDecayChannel& right)
{
    fName = right.fName;
    fDisplacementType = right.fDisplacementType;
    if(right.fProductsVector)
    {
        fProductsVector = new vector<G4MoleculeHandle>(*(right.fProductsVector));
    }
    fRMSProductsDisplacementVector = right.fRMSProductsDisplacementVector;
    fReleasedEnergy = right.fReleasedEnergy;
    fProbability= right.fProbability;
    fRMSMotherMoleculeDisplacement = right.fRMSMotherMoleculeDisplacement;
    return *this;
}

void G4MolecularDecayChannel::AddProduct(const G4Molecule* molecule, G4double displacement)
{
    if(!fProductsVector) fProductsVector = new vector<G4MoleculeHandle> ;

    G4MoleculeHandle molHandle(G4MoleculeHandleManager::Instance()->GetMoleculeHandle(molecule));
    fProductsVector->push_back(molHandle);
    fRMSProductsDisplacementVector.push_back(displacement);
}

G4int G4MolecularDecayChannel::GetNbProducts() const
{
    if(fProductsVector)
        return fProductsVector->size();
    return 0;
}

const G4Molecule* G4MolecularDecayChannel::GetProduct(int index) const
{
    if(fProductsVector)
        return ((*fProductsVector)[index]).get();
//    return ((*fProductsVector)[index])();

    return 0;
}

G4double G4MolecularDecayChannel::GetRMSRadialDisplacementOfProduct(const G4Molecule* product)
{
    if(!fProductsVector) return -1.;

    G4int sz = fProductsVector->size();
    G4double value = DBL_MAX;
    for (G4int i=0; i<sz ; i++)
    {
        if(*product != *((*fProductsVector)[i]).get())
//            if(*product != *((*fProductsVector)[i])())
        {
            value = fRMSProductsDisplacementVector[i];
        }
    }
    return value;
}

