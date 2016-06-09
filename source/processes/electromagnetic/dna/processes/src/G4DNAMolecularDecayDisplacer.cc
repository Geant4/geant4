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
// $Id: G4DNAMolecularDecayDisplacer.cc 65022 2012-11-12 16:43:12Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4DNAMolecularDecayDisplacer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"

using namespace std;

const DisplacementType G4DNAMolecularDecayDisplacer::Ionisation_DissociationDecay = G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAMolecularDecayDisplacer::A1B1_DissociationDecay = G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAMolecularDecayDisplacer::B1A1_DissociationDecay = G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAMolecularDecayDisplacer::AutoIonisation = G4VMolecularDecayDisplacer::AddDisplacement();
const DisplacementType G4DNAMolecularDecayDisplacer::DissociativeAttachment = G4VMolecularDecayDisplacer::AddDisplacement();

G4DNAMolecularDecayDisplacer::G4DNAMolecularDecayDisplacer() :
    G4VMolecularDecayDisplacer()
{;}

G4DNAMolecularDecayDisplacer::~G4DNAMolecularDecayDisplacer()
{;}

G4ThreeVector G4DNAMolecularDecayDisplacer::GetMotherMoleculeDisplacement(const G4MolecularDecayChannel* theDecayChannel) const
{
    G4int decayType = theDecayChannel -> GetDisplacementType();

    G4double RMSMotherMoleculeDisplacement=0;

    if(decayType == Ionisation_DissociationDecay)
    {
        RMSMotherMoleculeDisplacement =  2.0 * nanometer ;
    }
    else if(decayType == A1B1_DissociationDecay)
    {
        RMSMotherMoleculeDisplacement = 0. * nanometer ;
    }
    else if(decayType == B1A1_DissociationDecay)
    {
        RMSMotherMoleculeDisplacement = 0. * nanometer ;
    }
    else if(decayType == AutoIonisation)
    {
        RMSMotherMoleculeDisplacement = 2.0 * nanometer ;
    }
    else if(decayType == DissociativeAttachment)
    {
        RMSMotherMoleculeDisplacement = 0. * nanometer ;
    }

    if(RMSMotherMoleculeDisplacement==0)
    {
        return G4ThreeVector(0,0,0);
    }
    G4ThreeVector RandDirection = radialDistributionOfProducts(RMSMotherMoleculeDisplacement);

    return RandDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

vector<G4ThreeVector> G4DNAMolecularDecayDisplacer::GetProductsDisplacement(const G4MolecularDecayChannel* theDecayChannel) const
{
    G4int nbProducts = theDecayChannel -> GetNbProducts();
    vector<G4ThreeVector> theProductDisplacementVector (nbProducts);

    typedef map<const G4MoleculeDefinition*,G4double> RMSmap ;
    RMSmap theRMSmap;

    G4int decayType = theDecayChannel -> GetDisplacementType();

    if(decayType == Ionisation_DissociationDecay)
    {
        if(fVerbose)
            G4cout<<"Ionisation_DissociationDecay"<<G4endl;
        G4double RdmValue = G4UniformRand();

        if(RdmValue< 0.5)
        {
            // H3O
            theRMSmap[G4H3O::Definition()] = 0.* nanometer;
            // OH
            theRMSmap[G4OH::Definition()] = 0.8* nanometer;
        }
        else
        {
            // H3O
            theRMSmap[G4H3O::Definition()] = 0.8* nanometer;
            // OH
            theRMSmap[G4OH::Definition()] = 0.* nanometer;
        }

        for(int i = 0; i < nbProducts ; i++)
        {
            G4double theRMSDisplacement;
            const G4Molecule* product = theDecayChannel->GetProduct(i);
            theRMSDisplacement = theRMSmap[product->GetDefinition()];

            if(theRMSDisplacement==0)
            {
                theProductDisplacementVector[i] = G4ThreeVector();
            }
            else
            {
                G4ThreeVector RandDirection = radialDistributionOfProducts(theRMSDisplacement);
                theProductDisplacementVector[i] = RandDirection;
            }
        }
    }
    else if(decayType == A1B1_DissociationDecay)
    {
        if(fVerbose)
            G4cout<<"A1B1_DissociationDecay"<<G4endl;
        G4double theRMSDisplacement = 2.4 * nanometer;
        G4ThreeVector RandDirection = radialDistributionOfProducts(theRMSDisplacement);

        for(G4int i =0 ; i < nbProducts ; i++)
        {
            const G4Molecule* product = theDecayChannel->GetProduct(i);
            if(product->GetDefinition()== G4OH::Definition())
            {
                theProductDisplacementVector[i] = -1./18.*RandDirection;
            }
            else if(product->GetDefinition() == G4Hydrogen::Definition())
            {
                theProductDisplacementVector[i] = +17./18.*RandDirection;
            }
        }
    }
    else if(decayType == B1A1_DissociationDecay)
    {
        if(fVerbose)
            G4cout<<"B1A1_DissociationDecay"<<G4endl;
        G4double theRMSDisplacement = 0.8 * nanometer;
        G4ThreeVector RandDirection = radialDistributionOfProducts(theRMSDisplacement);

        G4int NbOfOH = 0;
        for(G4int i =0 ; i < nbProducts ; i++)
        {
            const G4Molecule* product = theDecayChannel->GetProduct(i);
            if(product->GetDefinition() == G4H2::Definition())
            {
                theProductDisplacementVector[i] = -2./18.*RandDirection;
            }
            else if(product->GetDefinition() == G4OH::Definition())
            {
                G4ThreeVector OxygenDisplacement = +16./18.*RandDirection;
                G4double OHRMSDisplacement = 1.1 * nanometer;

                G4ThreeVector OHDisplacement = radialDistributionOfProducts(OHRMSDisplacement) ;

                if(NbOfOH==0)
                {
                    OHDisplacement = 1./2.*OHDisplacement;
                }
                else
                {
                    OHDisplacement = -1./2.*OHDisplacement;
                }

                theProductDisplacementVector[i]  = OHDisplacement + OxygenDisplacement;

                NbOfOH ++;
            }
        }
    }
    else if(decayType == AutoIonisation)
    {
        if(fVerbose)
            G4cout<<"AutoIonisation"<<G4endl;
        G4double RdmValue = G4UniformRand();

        if(RdmValue< 0.5)
        {
            // H3O
            theRMSmap[G4H3O::Definition()] = 0.* nanometer;
            // OH
            theRMSmap[G4OH::Definition()] = 0.8* nanometer;
        }
        else
        {
            // H3O
            theRMSmap[G4H3O::Definition()] = 0.8* nanometer;
            // OH
            theRMSmap[G4OH::Definition()] = 0.* nanometer;
        }

        for(G4int i =0 ; i < nbProducts ; i++)
        {
            G4double theRMSDisplacement;
            const G4Molecule* product = theDecayChannel->GetProduct(i);
            theRMSDisplacement = theRMSmap[product->GetDefinition()];

            if(theRMSDisplacement==0)
            {
                theProductDisplacementVector[i] = G4ThreeVector();
            }
            else
            {
                G4ThreeVector RandDirection = radialDistributionOfProducts(theRMSDisplacement);
                theProductDisplacementVector[i] = RandDirection;
            }
            if(product->GetDefinition() == G4Electron_aq::Definition())
            {
                theProductDisplacementVector[i]=radialDistributionOfElectron();
            }
        }
    }
    else if(decayType == DissociativeAttachment)
    {
        if(fVerbose)
            G4cout<<"DissociativeAttachment"<<G4endl;
        G4double theRMSDisplacement = 0.8 * nanometer;
        G4ThreeVector RandDirection = radialDistributionOfProducts(theRMSDisplacement);

        G4int NbOfOH = 0;
        for(G4int i =0 ; i < nbProducts ; i++)
        {
            const G4Molecule* product = theDecayChannel->GetProduct(i);
            if(product->GetDefinition() == G4H2::Definition())
            {
                theProductDisplacementVector[i] = -2./18.*RandDirection;
            }
            else if(product->GetDefinition() == G4OH::Definition())
            {
                G4ThreeVector OxygenDisplacement = +16./18.*RandDirection;
                G4double OHRMSDisplacement = 1.1 * nanometer;

                G4ThreeVector OHDisplacement = radialDistributionOfProducts(OHRMSDisplacement) ;

                if(NbOfOH==0)
                {
                    OHDisplacement = 1./2.*OHDisplacement;
                }
                else
                {
                    OHDisplacement = -1./2.*OHDisplacement;
                }

                theProductDisplacementVector[i]  = OHDisplacement + OxygenDisplacement;

                NbOfOH ++;
            }
        }
    }

    return theProductDisplacementVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4DNAMolecularDecayDisplacer::radialDistributionOfProducts(G4double Rrms) const
{
    G4double sigma = Rrms/sqrt(3.);
    G4double expectationValue = 2.*sqrt(2./3.14)*sigma;

    G4double XValueForfMax = sqrt(2.*sigma*sigma);
    G4double fMaxValue = sqrt(2./3.14) * 1./(sigma*sigma*sigma) *
            (XValueForfMax*XValueForfMax)*
            exp(-1./2. * (XValueForfMax*XValueForfMax)/(sigma*sigma));

    G4double R(-1.);

    do
    {
        G4double aRandomfValue = fMaxValue*G4UniformRand();

        G4double sign;
        if(G4UniformRand() > 0.5)
        {
            sign = +1.;
        }
        else
        {
            sign = -1;
        }

        R = expectationValue + sign*3.*sigma* G4UniformRand();
        G4double f = sqrt(2./3.14) * 1/pow(sigma, 3) * R*R * exp(-1./2. * R*R/(sigma*sigma));

        if(aRandomfValue < f)
        {
            break;
        }
    }
    while(1);

    G4double costheta = (2.*G4UniformRand()-1.);
    G4double theta = acos (costheta);
    G4double phi = 2.*pi*G4UniformRand();

    G4double xDirection = R*cos(phi)* sin(theta);
    G4double yDirection = R*sin(theta)*sin(phi);
    G4double zDirection = R*costheta;
    G4ThreeVector RandDirection(xDirection, yDirection, zDirection);

    return RandDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4DNAMolecularDecayDisplacer::radialDistributionOfElectron() const
{

    G4double sigma = 1./2.;
    G4double expectationValue = 1. ;

    G4double XValueForfMax = 1./2.;
    G4double fMaxValue = 4. * XValueForfMax *
            exp(-2. * XValueForfMax);

    G4double R(-1.);

    do
    {
        G4double aRandomfValue = fMaxValue*G4UniformRand();

        G4double sign;
        if(G4UniformRand() > 0.5)
        {
            sign = +1;
        }
        else
        {
            sign = -1;
        }

        R = (expectationValue * G4UniformRand() )+ sign*3*sigma* G4UniformRand();
        G4double f = 4* R * exp(- 2. * R);

        if(aRandomfValue < f)
        {
            break;
        }
    }
    while(1);

    G4double Rnano = R *10* nanometer;

    G4double costheta = (2*G4UniformRand()-1);
    G4double theta = acos (costheta);
    G4double phi = 2*pi*G4UniformRand();

    G4double xDirection = Rnano*cos(phi)* sin(theta);
    G4double yDirection = Rnano*sin(theta)*sin(phi);
    G4double zDirection = Rnano*costheta;
    G4ThreeVector RandDirection(xDirection, yDirection, zDirection);

    return RandDirection;
}
