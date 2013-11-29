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
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// An abstract class to model interaction laws, not necessarily.
// exponential ones.
// These laws, p(l), are caraterized by their:
//    - non-interaction probability over l:
//      P_NI(l) = 1 - integral(0,l) p(s) ds
//    - effective cross-section at l:
//      sigma_eff(l) = p(l)/P_NI(l) (several forms exist)
// If several laws compete, the resulting net law has:
//    - an non-interaction probability over a lenght l which
//      is the product of the individual non-interaction probabilities.
//    - an effective cross-section which is the
//      sum on the individual effective cross-sections.
//
//      ----------------G4VBiasingInteractionLaw ----------------
//
// Author: M.Verderi (LLR), November 2013
// - 05/11/13 : First Implementation
// --------------------------------------------------------------------


#ifndef G4VBiasingInteractionLaw_hh
#define G4VBiasingInteractionLaw_hh 1

#include "globals.hh"
#include <vector>

class G4BiasingProcessInterface;

class G4VBiasingInteractionLaw {
public:
  G4VBiasingInteractionLaw(G4String name) : fName(name), fSampledInteractionLength(DBL_MAX) {}
  virtual ~G4VBiasingInteractionLaw() {}

public:
  const G4String&       GetName() const { return fName; }
  
  // ----------------------------
  // -- Interface to sub-classes:
  // ----------------------------
protected:
  // -- Sample the distribution for point like interaction (PostStep ones)
  virtual G4double            SampleInteractionLength()                              = 0;
public:
  // -- Compute non-interaction probability and effective cross-section:
  // -- (probability of interaction over dl =  effective_cross-section* dl)
  virtual G4double ComputeNonInteractionProbabilityAt(G4double               length) const = 0;
  virtual G4double     ComputeEffectiveCrossSectionAt(G4double               length) const = 0;
protected:
  // -- Convenience method, used in many daughters classes :
  // -- update the distribution for a made step of truePathLength size:
  virtual G4double     UpdateInteractionLengthForStep(G4double   /* truePathLength */)  { return DBL_MAX; }
public:
  // -- Methods to deal with singularities : null cross sections or infinite ones.
  // -- In such cases, weight can not always be computed.
  // -- Tells if this interaction law has singularities:
  virtual G4bool                      IsSingular() const {return false;}
    // -- method interrogated only in case interaction law is IsSingular() == true:
  virtual G4bool IsEffectiveCrossSectionInfinite() const {return false;}
  
  
  // -----------------------------------------
  // -- public interface to protected methods:
  // -----------------------------------------
public:
  G4double Sample()
  {
    fSampledInteractionLength = SampleInteractionLength();
    return fSampledInteractionLength;
  }
  G4double UpdateForStep(G4double truePathLength)
  {
    fSampledInteractionLength = UpdateInteractionLengthForStep(truePathLength);
    return fSampledInteractionLength;
  }
  G4double GetSampledInteractionLength() const { return fSampledInteractionLength; }

  
private:
  G4String                     fName;
  G4double fSampledInteractionLength;
};

#endif
