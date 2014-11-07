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
#include "G4ILawCommonTruncatedExp.hh"
#include "G4Track.hh"

#include "G4BOptnForceCommonTruncatedExp.hh"
#include "G4BiasingProcessInterface.hh"


G4ILawCommonTruncatedExp::G4ILawCommonTruncatedExp(G4String name)
  : G4VBiasingInteractionLaw(name),
    fExpInteractionLaw("expLawFor"+name)
{}

G4ILawCommonTruncatedExp::~G4ILawCommonTruncatedExp()
{}


G4double G4ILawCommonTruncatedExp::ComputeEffectiveCrossSectionAt(G4double distance) const
{
  return fExpInteractionLaw.ComputeEffectiveCrossSectionAt( distance ) * fSelectedProcessXSfraction;
}

G4double G4ILawCommonTruncatedExp::ComputeNonInteractionProbabilityAt(G4double distance) const
{
  G4double niProba = fExpInteractionLaw.ComputeNonInteractionProbabilityAt( distance );
  
  if ( niProba <= 0.0 )
    {
      G4ExceptionDescription ed;
      ed << " Negative probability for `" << GetName() << "' p = " << niProba << " distance = " << distance <<  " !!! " << G4endl;
      G4Exception(" G4ILawCommonTruncatedExp::ComputeNonInteractionProbabilityAt(...)",
		  "BIAS.GEN.08",
		  JustWarning,
		  ed);
    }

  return niProba;

}

G4double G4ILawCommonTruncatedExp::SampleInteractionLength()
{
  fInteractionDistance = fExpInteractionLaw.SampleInteractionLength();
  return fInteractionDistance;
}

G4double G4ILawCommonTruncatedExp::UpdateInteractionLengthForStep(G4double truePathLength)
{
  fInteractionDistance = fExpInteractionLaw.UpdateInteractionLengthForStep(truePathLength);
  return fInteractionDistance;
}
