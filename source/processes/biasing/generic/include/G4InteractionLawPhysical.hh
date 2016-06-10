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
// $Id: $
//
// ---------------------------------------------------------------
//
// G4InteractionLawPhysical
//
// Class Description:
//     A G4VBiasingInteractionLaw representing the usual
// physical exponential law, of constant cross-section of the step.
//
// ---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi

#ifndef G4InteractionLawPhysical_hh
#define G4InteractionLawPhysical_hh 1

#include "G4VBiasingInteractionLaw.hh"

class G4InteractionLawPhysical : public G4VBiasingInteractionLaw
{
public:
  G4InteractionLawPhysical(G4String name = "exponentialLaw");
  virtual ~G4InteractionLawPhysical();

public:
  void     SetPhysicalCrossSection(G4double crossSection);
  G4double GetPhysicalCrossSection() const {return fCrossSection;}
  
public:
  virtual G4double     ComputeEffectiveCrossSectionAt(G4double               length) const;
  virtual G4double ComputeNonInteractionProbabilityAt(G4double               length) const;
  // -- sample the distribution
  virtual  G4double           SampleInteractionLength();
  // -- move by true path length, this position becomes the new initial point
  virtual G4double     UpdateInteractionLengthForStep(G4double       truePathLength);


private:
  G4double fCrossSection;
  G4bool   fCrossSectionDefined;
  G4double fNumberOfInteractionLength;

};

#endif
