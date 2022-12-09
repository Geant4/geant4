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


#include "G4IRTUtils.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ErrorFunction.hh"
G4double G4IRTUtils::EffectiveDistance(const G4double& rc,
                                       const G4double& r0)
{
    return r0 == 0 ? 0 : - rc / (1 - std::exp( rc / r0 ) );
}


G4double G4IRTUtils::GetRCutOff()
{
    G4double tCutOff = 1000 * ns;
    
    G4double probabilityOfReaction = 0.01;
    G4double maximumReactionRadius = 1.45*CLHEP::nm;//??
    G4double maximumRelativeDiffusionCoefficient = 2.0*9.46e9 *CLHEP::nm*CLHEP::nm/CLHEP::s;//??
    G4double erfcInv = G4ErrorFunction::erfcInv(probabilityOfReaction);
    return maximumReactionRadius + 2.0 *
    std::sqrt(maximumRelativeDiffusionCoefficient * tCutOff) * erfcInv;
}


G4double G4IRTUtils::GetRCutOff(G4double tCutOff)
{
    G4double probabilityOfReaction = 0.01;
    G4double maximumReactionRadius = 1.45*CLHEP::nm;//??
    G4double maximumRelativeDiffusionCoefficient = 2.0*9.46e9 *CLHEP::nm*CLHEP::nm/CLHEP::s;//??
    G4double erfcInv = G4ErrorFunction::erfcInv(probabilityOfReaction);
    return maximumReactionRadius + 2.0 *
    std::sqrt(maximumRelativeDiffusionCoefficient * tCutOff) * erfcInv;
}

G4double G4IRTUtils::GetDNADistanceCutOff()
{
    G4double tCutOff = 100 * ps;
    
    G4double probabilityOfReaction = 0.01;
    G4double maximumReactionRadius = 1.45*CLHEP::nm;//??
    G4double maximumRelativeDiffusionCoefficient = 2.0*9.46e9 *CLHEP::nm*CLHEP::nm/CLHEP::s;//??
    G4double erfcInv = G4ErrorFunction::erfcInv(probabilityOfReaction);
    return maximumReactionRadius + 2.0 *
    std::sqrt(maximumRelativeDiffusionCoefficient * tCutOff) * erfcInv;
}