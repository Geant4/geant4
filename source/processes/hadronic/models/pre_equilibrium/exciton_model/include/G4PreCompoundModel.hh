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
// $Id: G4PreCompoundModel.hh,v 1.6 2008/09/22 10:18:36 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// by V. Lara

// Class Description
// Model implementation for pre-equilibrium decay models in geant4. 
// To be used in your physics list, in case you neeed this kind of physics.
// Can be used as a stand-allone model, but also in conjunction with an intra-nuclear
// transport, or any of the string-parton models.
// Class Description - End
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option.(default OPTxs=3)
// JMQ (06 September 2008) Also external choices have been added for:
//                - superimposed Coulomb barrier (if useSICB=true, default false) 
//                - "never go back"  hipothesis (if useNGB=true, default false) 
//                - soft cutoff from preeq. to equlibrium (if useSCO=true, default false)
//                - CEM transition probabilities (if useCEMtr=true, default false)  

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4LorentzVector.hh"


#include "G4NucleiProperties.hh"
#include "G4PreCompoundParameters.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"

//#include "G4PreCompoundEmission.hh"

#include "G4DynamicParticle.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#define debug
//#define verbose

class G4PreCompoundModel : public G4VPreCompoundModel
{
 
public:
  G4PreCompoundModel(G4ExcitationHandler * const value) : 
    G4VPreCompoundModel(value), useHETCEmission(false), useGNASHTransition(false), 
    OPTxs(3), useSICB(false), useNGB(false), useSCO(false), useCEMtr(false) {}

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

 //for cross section selection
  inline void SetOPTxs(G4int opt) { OPTxs = opt; }
//for the rest of external choices
  inline void UseSICB() { useSICB = true; }
  inline void UseNGB()  { useNGB = true; }
  inline void UseSCO()  { useSCO = true; }
  inline void UseCEMtr() { useCEMtr = true; }
private:  

  void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				  G4ReactionProductVector * theResult) const;

private:

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

//for cross section options
  G4int OPTxs;
//for the rest of external choices
  G4bool useSICB;
  G4bool useNGB;
  G4bool useSCO;
  G4bool useCEMtr;


    G4HadFinalState theResult;

#ifdef PRECOMPOUND_TEST
  static G4Fragment theInitialFragmentForTest;
  static std::vector<G4String*> theCreatorModels;
#endif

};
#endif

