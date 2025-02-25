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
// G4ILawTruncatedExp
//
// Class Description:
//
// A G4VBiasingInteractionLaw representing a truncated exponential
// law : an exponential law acting on [0,L] segment. The law is
// driven by a cross-section.
//
// Author: Marc Verderi, November 2013.
// --------------------------------------------------------------------
#ifndef G4ILawTruncatedExp_hh
#define G4ILawTruncatedExp_hh 1

#include "G4VBiasingInteractionLaw.hh"

class G4ILawTruncatedExp : public G4VBiasingInteractionLaw
{
  public:

    G4ILawTruncatedExp(const G4String& name = "expForceInteractionLaw");
    virtual ~G4ILawTruncatedExp();
  
    virtual G4double ComputeEffectiveCrossSectionAt(G4double length) const;
    virtual G4double ComputeNonInteractionProbabilityAt(G4double length) const;
    // -- sample the distribution
    virtual G4double SampleInteractionLength();
    // -- move by true path length, this position becomes the new initial point
    virtual G4double UpdateInteractionLengthForStep(G4double truePathLength);
    virtual G4bool IsSingular() const { return fIsSingular; }

    void SetForceCrossSection(G4double xs);

    void SetMaximumDistance(G4double d) { fMaximumDistance = d; }
    G4double GetMaximumDistance() const { return fMaximumDistance; }
    G4double GetInteractionDistance() const { return fInteractionDistance; }

  private:

    G4double fMaximumDistance = 0.0;
    G4double fCrossSection = 0.0;
    G4bool fCrossSectionDefined = false;
    G4bool fIsSingular = false;
    G4double fInteractionDistance = 0.0;
};

#endif
