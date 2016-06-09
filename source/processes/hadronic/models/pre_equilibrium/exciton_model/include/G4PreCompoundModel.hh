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
// $Id: G4PreCompoundModel.hh,v 1.3 2006/06/29 20:58:26 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// by V. Lara

// Class Description
// Model implementation for pre-equilibrium decay models in geant4. 
// To be used in your physics list, in case you neeed this kind of physics.
// Can be used as a stand-allone model, but also in conjunction with an intra-nuclear
// transport, or any of the string-parton models.
// Class Description - End

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4LorentzVector.hh"


#include "G4NucleiProperties.hh"
#include "G4PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"


//#define debug
//#define verbose

class G4PreCompoundModel : public G4VPreCompoundModel
{
public:
    
  G4PreCompoundModel(G4ExcitationHandler * const value) : 
    G4VPreCompoundModel(value), useHETCEmission(false), useGNASHTransition(false) {}

  ~G4PreCompoundModel() {}

private:
  G4PreCompoundModel() {}

  G4PreCompoundModel(const G4PreCompoundModel &) : G4VPreCompoundModel() {}

  const G4PreCompoundModel& operator=(const G4PreCompoundModel &right);
  G4bool operator==(const G4PreCompoundModel &right) const;
  G4bool operator!=(const G4PreCompoundModel &right) const;

public:
    G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus);

  G4ReactionProductVector* DeExcite(const G4Fragment& aFragment) const;

#ifdef PRECOMPOUND_TEST
  static G4Fragment GetInitialFragmentForTest()
  { return G4PreCompoundModel::theInitialFragmentForTest; }
  static std::vector<G4String*> * GetCreatorModels()
  { return &G4PreCompoundModel::theCreatorModels; }
#endif

  inline void UseHETCEmission() { useHETCEmission = true; }
  inline void UseDefaultEmission() { useHETCEmission = false; }
  inline void UseGNASHTransition() { useGNASHTransition = true; }
  inline void UseDefaultTransition() { useGNASHTransition = false; }

private:  

  void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				  G4ReactionProductVector * theResult) const;

#ifdef debug				  
  void CheckConservation(const G4Fragment & theInitialState,
			 const G4Fragment & aFragment,
			 G4ReactionProductVector * Result) const;
#endif

  //==============
  // Data Members 
  //==============

  G4bool           useHETCEmission;
  G4bool           useGNASHTransition;
    G4HadFinalState theResult;

#ifdef PRECOMPOUND_TEST
  static G4Fragment theInitialFragmentForTest;
  static std::vector<G4String*> theCreatorModels;
#endif

};
#endif

