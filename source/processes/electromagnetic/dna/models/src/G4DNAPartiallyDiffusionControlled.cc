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

#include "G4DNAPartiallyDiffusionControlled.hh"
#include "G4IRTUtils.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4VDNAReactionModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4MolecularConfiguration.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4ITReactionChange.hh"
#include "G4VReactionType.hh"
#include "G4Electron_aq.hh"
#include "G4ErrorFunction.hh"
G4DNAPartiallyDiffusionControlled::G4DNAPartiallyDiffusionControlled() 
	: G4VReactionType()
{}

G4DNAPartiallyDiffusionControlled::~G4DNAPartiallyDiffusionControlled() = default;

G4bool
G4DNAPartiallyDiffusionControlled::GeminateRecombinationProbability(const G4MolecularConfiguration* mA,
                                                                    const G4MolecularConfiguration* mB)
{
    auto reactionData = G4DNAMolecularReactionTable::Instance()
        ->GetReactionData(mA, mB);

    G4double D = GetDiffusionCoefficient(mA, mB);
    G4double R = mA->GetVanDerVaalsRadius() + mB->GetVanDerVaalsRadius();
    
    const G4double Rs = 0.3 * nm;
    G4double kobs = reactionData->GetObservedReactionRateConstant() / Avogadro;
    if(mA->GetCharge() * mB->GetCharge() == 0)
    {
        G4double kdif = 4 * CLHEP::pi * D * R * Avogadro;
        G4double kact = G4IRTUtils::GetKact(kobs, kdif);
        return G4UniformRand() < Rs / ( Rs + ( kdif / kact ) * ( R + Rs ));
    }
    else
    {
        G4double rc = 0.71 * nm * mA->GetCharge() *
                                  mB->GetCharge();
        G4double sigmaEff = G4IRTUtils::EffectiveDistance(rc, R);
        G4double kdif = 4 * CLHEP::pi * D * sigmaEff;
        G4double kact = G4IRTUtils::GetKact(kobs, kdif);
        G4double a = std::exp( -rc / R );
        G4double b = std::exp( -rc / ( R + Rs ) );
        G4double Preact = ( a - b ) / ( a - b - ( kdif / kact ) * ( 1 - a ) );
        
        return G4UniformRand() < Preact;
    }
}

G4double
G4DNAPartiallyDiffusionControlled::GetDiffusionCoefficient(const G4MolecularConfiguration* mA,
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
             mB->GetDiffusionCoefficient());
    }
    return D;
}

G4double G4DNAPartiallyDiffusionControlled::GetTimeToEncounter(const G4Track& trackA,
                                                               const G4Track& trackB)
{
    auto pMolConfA = GetMolecule(trackA)->GetMolecularConfiguration();
    auto pMolConfB = GetMolecule(trackB)->GetMolecularConfiguration();

    G4double D = GetDiffusionCoefficient(pMolConfA, pMolConfB);
    auto reactionData = G4DNAMolecularReactionTable::Instance()
    ->GetReactionData(pMolConfA, pMolConfB);
    G4double Reff;
    G4double kobs = reactionData->GetObservedReactionRateConstant();
    G4double distance = (trackA.GetPosition() - trackB.GetPosition()).mag();
    G4double SmoluchowskiRadius;
    G4double RVal = pMolConfA->GetVanDerVaalsRadius() + pMolConfB->GetVanDerVaalsRadius();
    
    if((pMolConfA->GetCharge() != 0) &&
       (pMolConfB->GetCharge() != 0))
    {
        G4double rc = 0.71 * nm * pMolConfA->GetCharge() *
        pMolConfB->GetCharge();
        distance = G4IRTUtils::EffectiveDistance( rc, distance );
        Reff = G4IRTUtils::EffectiveDistance( rc, RVal );
        SmoluchowskiRadius = Reff;
    }
    else
    {
        SmoluchowskiRadius = RVal;
    }
    
    G4double Winf = SmoluchowskiRadius / distance;
    G4double U1 = G4UniformRand();
    G4double U2 = G4UniformRand();
    G4double U = G4UniformRand();
    G4double X = 0;
    G4double irt_1 = -1.0 * ps;
    G4double irt_2;

    G4double kdif = 4 * CLHEP::pi * D * SmoluchowskiRadius * Avogadro;
    G4double kact = G4IRTUtils::GetKact(kobs, kdif);
    
    if ( U < Winf )
    {
        G4double d = ( distance - SmoluchowskiRadius ) /
                     G4ErrorFunction::erfcInv( U / Winf );
        irt_1 = ( 1.0 / ( 4 * D ) ) * d * d;
    }
    
    if( irt_1 < 0)
    {
        return irt_1;
    }
    else
    {
        G4double rateFactor = kact / ( kact + kdif );
        if( U1 > rateFactor )
        {
            return -1.0 * ps;
        }
        G4double Y = std::abs(G4RandGauss::shoot(0.0,std::sqrt(2)));

        if( Y > 0)
        {
            X = - ( G4Log( U2 ) ) / Y;
        }

        G4double f = X * SmoluchowskiRadius * kdif / ( kact + kdif );
        irt_2 = ( f * f ) / D ;
    }
    
    return irt_1 + irt_2;
}

