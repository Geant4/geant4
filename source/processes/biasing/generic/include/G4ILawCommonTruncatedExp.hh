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
// G4ILawCommonTruncatedExp
//
// Class Description:
//     A G4VBiasingInteractionLaw representing a truncated exponential
//  law : an exponential law acting on [0,L] segment. This law is
//  Common as it collects several process cross-sections, which sum
//  is used to drive the law.
//
// ---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi

#ifndef G4ILawCommonTruncatedExp_hh
#define G4ILawCommonTruncatedExp_hh 1

#include "G4VBiasingInteractionLaw.hh"
#include "G4ILawTruncatedExp.hh"
#include "G4ThreeVector.hh"

class G4BOptnForceCommonTruncatedExp;

class G4ILawCommonTruncatedExp : public G4VBiasingInteractionLaw
{
public:
  G4ILawCommonTruncatedExp(G4String name = "expSharedForceInteractionLaw");
  virtual ~G4ILawCommonTruncatedExp();
  
public:
  virtual G4bool                             IsSingular() const
  {return fExpInteractionLaw.IsSingular();}
  virtual G4bool IsEffectiveCrossSectionInfinite() const
  {return fExpInteractionLaw.IsEffectiveCrossSectionInfinite();}
private:
  // -- sample the distribution:
  // -- in this case, cheat by returning DBL_MAX, to let operation proposing the
  // -- step length in the alongGPIL
  virtual G4double              SampleInteractionLength();
  // -- move by true path length, this position becomes the new initial point
  // -- in this case, cheat by returning DBL_MAX, to let operation proposing the
  // -- step length in the alongGPIL
  virtual G4double       UpdateInteractionLengthForStep(G4double       truePathLength);
public:
  virtual G4double       ComputeEffectiveCrossSectionAt(G4double               length) const;
  virtual G4double   ComputeNonInteractionProbabilityAt(G4double               length) const;

public:
  void         SetNumberOfSharing(G4int n)          { fNumberOfSharing = n; }
  G4int        GetNumberOfSharing() const           { return fNumberOfSharing; }
  void       SetForceCrossSection( G4double xs )    { fExpInteractionLaw.SetForceCrossSection( xs ); }
  void                      reset();
  void         SetMaximumDistance(G4double d) { fExpInteractionLaw.SetMaximumDistance(d); }
  G4double     GetMaximumDistance() const {return fExpInteractionLaw.GetMaximumDistance();}
  G4double GetInteractionDistance() const {return fExpInteractionLaw.GetInteractionDistance();}
  void               SetOperation( G4BOptnForceCommonTruncatedExp* operation) {fOperation = operation;}

private:
  G4ILawTruncatedExp                     fExpInteractionLaw;
  //G4bool                                        fIsSingular;
  G4int                                    fNumberOfSharing;
  //G4bool                           fFirstInitializationCall;
  G4bool                                   fFirstUpdateCall;
  G4bool                                 fFirstSamplingCall;
  G4double                             fInteractionDistance;
  // TRICKY
  G4BOptnForceCommonTruncatedExp*        fOperation;
};

#endif
