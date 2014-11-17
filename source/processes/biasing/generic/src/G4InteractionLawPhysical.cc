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
#include "G4InteractionLawPhysical.hh"
#include "Randomize.hh"

G4InteractionLawPhysical::G4InteractionLawPhysical(G4String name)
  : G4VBiasingInteractionLaw(name),
    fCrossSection(0.0),
    fCrossSectionDefined(false),
    fNumberOfInteractionLength(-1.0)
{}

G4InteractionLawPhysical::~G4InteractionLawPhysical()
{}

void G4InteractionLawPhysical::SetPhysicalCrossSection(G4double crossSection)
{
  if (crossSection < 0.0)
    {
      G4Exception("G4InteractionLawPhysical::SetPhysicalCrossSection(..)",
		  "BIAS.GEN.14",
		  JustWarning,
		  "Cross-section value passed is negative. It is set to zero !");
      crossSection = 0.0;
    }
  fCrossSectionDefined = true;
  fCrossSection        = crossSection;
}

G4double G4InteractionLawPhysical::ComputeEffectiveCrossSectionAt(G4double) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeEffectiveCrossSection(..)",
					 "BIAS.GEN.15",
					 JustWarning,
					 "Cross-section value requested, but has not been defined yet. Assumes 0 !");
  return fCrossSection;
}

G4double G4InteractionLawPhysical::ComputeNonInteractionProbabilityAt(G4double stepLength) const
{
  if (!fCrossSectionDefined) G4Exception("G4InteractionLawPhysical::ComputeNonInteractionProbability(..)",
					 "BIAS.GEN.16",
					 JustWarning,
					 "Non interaction probabitlity value requested, but cross section has not been defined yet. Assumes it to be 0 !");
  // -- allows zero cross-section case, by convention:
  if ( fCrossSection == 0.0 ) return 1.0;
  else return std::exp(-fCrossSection*stepLength);
}

G4double G4InteractionLawPhysical::SampleInteractionLength()
{
  if ( !fCrossSectionDefined || fCrossSection < 0.0 )  G4Exception("G4InteractionLawPhysical::Sample(..)",
								   "BIAS.GEN.17",
								   FatalException,
								   "Trying to sample while cross-section is not defined or < 0 !");
  if ( fCrossSection == 0.0 ) return DBL_MAX;

  fNumberOfInteractionLength =  -std::log( G4UniformRand() );
  return fNumberOfInteractionLength/fCrossSection;
}


G4double G4InteractionLawPhysical::UpdateInteractionLengthForStep(G4double       truePathLength)
{
  fNumberOfInteractionLength -= truePathLength*fCrossSection;
  
  if ( fNumberOfInteractionLength < 0 ) 
    {
      G4ExceptionDescription ed;
      ed << " Negative number of interaction length for `" << GetName() << "' " << fNumberOfInteractionLength << ", set it to zero !" << G4endl; 
      G4Exception("G4InteractionLawPhysical::UpdateInteractionLengthForStep(...)",
		  "BIAS.GEN.13",
		  JustWarning,
		  ed);
      fNumberOfInteractionLength = 0.;
    }
  return  fNumberOfInteractionLength/fCrossSection;
}
