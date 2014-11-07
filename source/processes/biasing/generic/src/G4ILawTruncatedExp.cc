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
#include "G4ILawTruncatedExp.hh"
#include "Randomize.hh"
#include "G4Track.hh"


G4ILawTruncatedExp::G4ILawTruncatedExp(G4String name)
  : G4VBiasingInteractionLaw(name),
    fMaximumDistance(0.0),
    fCrossSection(0.0),
    fCrossSectionDefined(false),
    fIsSingular(false)
{}

G4ILawTruncatedExp::~G4ILawTruncatedExp()
{}

void G4ILawTruncatedExp::SetForceCrossSection(G4double crossSection)
{
  if (crossSection < 0.0)
    {
      G4Exception("G4ILawTruncatedExp::SetForceCrossSection(..)",
		  "BIAS.GEN.09",
		  JustWarning,
		  "Cross-section value passed is negative. It is set to zero !");
      fIsSingular  = true;
      crossSection = 0.0;
    }
  fIsSingular          = false;
  fCrossSectionDefined = true;
  fCrossSection        = crossSection;
}

G4double G4ILawTruncatedExp::ComputeEffectiveCrossSectionAt(G4double distance) const
{
  if ( !fCrossSectionDefined )
    {
      G4Exception("G4ILawTruncatedExp::ComputeEffectiveCrossSection(..)",
		  "BIAS.GEN.10",
		  JustWarning,
		  "Cross-section value requested, but has not been defined yet. Assumes 0 !");
      // -- zero cross-section, returns the limit form of the effective cross-section:
      return 1.0 / ( fMaximumDistance - distance );
    }
  G4double denum = 1.0 - std::exp( -fCrossSection * ( fMaximumDistance - distance) );
  return fCrossSection / denum;
}

G4double G4ILawTruncatedExp::ComputeNonInteractionProbabilityAt(G4double distance) const
{
  if (!fCrossSectionDefined)
    {
      G4Exception("G4ILawTruncatedExp::ComputeNonInteractionProbability(..)",
		  "BIAS.GEN.11",
		  JustWarning,
		  "Non interaction probability value requested, but cross section has not been defined yet. Assumes it to be 0 !");
      // -- return limit case of null cross-section:
      return 1.0 - distance / fMaximumDistance;
    }
  G4double   num = 1.0 - std::exp( -fCrossSection*distance         );
  G4double denum = 1.0 - std::exp( -fCrossSection*fMaximumDistance );
  return 1.0 - num/denum;
}

G4double G4ILawTruncatedExp::SampleInteractionLength()
{
  if ( !fCrossSectionDefined )
    {
      G4Exception("G4ILawTruncatedExp::Sample(..)",
		  "BIAS.GEN.12",
		  JustWarning,
		  "Trying to sample while cross-section is not defined, assuming 0 !");
      fInteractionDistance = G4UniformRand() * fMaximumDistance;
      return fInteractionDistance;
    }
  fInteractionDistance = -std::log(1.0 - G4UniformRand()* (1.0 - std::exp(-fCrossSection*fMaximumDistance)))/fCrossSection;
  return fInteractionDistance;
}


G4double G4ILawTruncatedExp::UpdateInteractionLengthForStep(G4double       truePathLength)
{
  fInteractionDistance -= truePathLength;
  fMaximumDistance     -= truePathLength;
  
  if ( fInteractionDistance < 0 )
    {
      G4ExceptionDescription ed;
      ed << " Negative number of interaction length for `" << GetName() << "' " << fInteractionDistance << ", set it to zero !" << G4endl; 
      G4Exception("G4ILawTruncatedExp::UpdateInteractionLengthForStep(...)",
		  "BIAS.GEN.13",
		  JustWarning,
		  "Trying to sample while cross-section is not defined, assuming 0 !");
      fInteractionDistance = 0.;
    }
  
  return  fInteractionDistance;
}
