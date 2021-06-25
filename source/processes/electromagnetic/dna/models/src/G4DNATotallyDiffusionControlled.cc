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
// 20/2/2019
// Author : HoangTRAN

#include "G4DNATotallyDiffusionControlled.hh"
#include "G4IRTUtils.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4VDNAReactionModel.hh"
#include "G4OctreeFinder.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4MolecularConfiguration.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4Electron_aq.hh"
#include "G4Hydrogen.hh"
#include "G4ErrorFunction.hh"
G4DNATotallyDiffusionControlled::G4DNATotallyDiffusionControlled() 
    : G4VReactionType()
{}
G4DNATotallyDiffusionControlled::~G4DNATotallyDiffusionControlled() = default;

G4double G4DNATotallyDiffusionControlled::GetTimeToEncounter(const G4Track& trackA,
                                                             const G4Track& trackB)
{
    auto pMolConfA = GetMolecule(trackA)->GetMolecularConfiguration();
    auto pMolConfB = GetMolecule(trackB)->GetMolecularConfiguration();

    G4double D = GetDiffusionCoefficient(pMolConfA, pMolConfB);
    auto reactionData = G4DNAMolecularReactionTable::Instance()
    ->GetReactionData(pMolConfA, pMolConfB);
    G4double kobs = reactionData->GetObservedReactionRateConstant();
    G4double distance = (trackA.GetPosition() - trackB.GetPosition()).mag();
    G4double Reff  = kobs / ( 4 * CLHEP::pi * D * Avogadro );

    if( distance < Reff )
    {
         G4ExceptionDescription exceptionDescription;
         exceptionDescription << "distance = "<< distance
                              << " is uncorrected with "
                              <<" Reff = "<<Reff
                              <<" for : "<<pMolConfA->GetName()
                              <<" and "<<pMolConfB->GetName();
         G4Exception("G4DNATotallyDiffusionControlled"
                     "::GetTimeToEncounter()", "G4DNATotallyDiffusionControlled02",
         FatalException, exceptionDescription);
    }
    
    G4double Winf = Reff / distance;
    G4double U = G4UniformRand();
    G4double irt = -1.0 * ps;

    if ( U < Winf )
    {
        G4double d = ( distance - Reff ) /
                    G4ErrorFunction::erfcInv( U / Winf );
        irt = ( 1.0 / ( 4 * D ) ) * d * d;
    }
    return irt;
}
G4bool G4DNATotallyDiffusionControlled::
GeminateRecombinationProbability(const G4MolecularConfiguration* pMolA,
                                 const G4MolecularConfiguration* pMolB)
{
    if(pMolA->GetDefinition() == G4Electron_aq::Definition() ||
       pMolA->GetDefinition() == G4Hydrogen::Definition())
    {
        G4bool spinA;
        G4bool spinB;
        spinA = G4UniformRand() < 0.5;
        if(spinA &&
          (pMolB->GetDefinition() == G4Electron_aq::Definition() ||
           pMolB->GetDefinition() ==  G4Hydrogen::Definition()))
        {
  	        spinB = G4UniformRand() < 0.5;
            if( !spinB )
            {
                return true;
            }
        }
        return false;
    }
    return true;
}

G4double
G4DNATotallyDiffusionControlled::GetDiffusionCoefficient(const G4MolecularConfiguration* mA,
                                                         const G4MolecularConfiguration* mB)
{
    G4double D;
     if(mA == mB)
     {    
         D = (mA->GetDiffusionCoefficient());
     }
     else
    {
        D = (mA->GetDiffusionCoefficient() +
             mB->GetDiffusionCoefficient());//
    }
    return D;
}




