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
// $Id: G4PreCompoundModel.hh 106233 2017-09-22 21:34:41Z gcosmo $
//
// by V. Lara
//
// Class Description
// Model implementation for pre-equilibrium decay models in geant4. 
// To be used in your physics list, in case you neeed this kind of physics.
// Can be used as a stand-allone model, but also in conjunction with an intra-nuclear
// transport, or any of the string-parton models.
// Class Description - End
//
// Modified:
// 03.09.2008 J.M.Quesada added external choice of inverse 
//            cross section option.(default OPTxs=3)
// 06.09.2008 J.M.Quesada external choices have been added for:
//                - superimposed Coulomb barrier (if useSICB=true, default false) 
//                - "never go back"  hipothesis (if useNGB=true, default false) 
//                - soft cutoff from preeq. to equlibrium (if useSCO=true, default false)
//                - CEM transition probabilities (if useCEMtr=true)
// 30.10.2009 J.M.Quesada CEM transition probabilities are set as default
// 20.08.2010 V.Ivanchenko Cleanup of the code - changed data members and inline methods
// 03.01.2012 V.Ivanchenko Added pointer to G4ExcitationHandler to the 
//                         constructor

#ifndef G4PreCompoundModel_h
#define G4PreCompoundModel_h 1

#include "G4VPreCompoundModel.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ExcitationHandler.hh"

class G4PreCompoundEmission;
class G4VPreCompoundTransitions;
class G4ParticleDefinition;

class G4PreCompoundModel : public G4VPreCompoundModel
{ 
public:

  explicit G4PreCompoundModel(G4ExcitationHandler* ptr = nullptr); 

  virtual ~G4PreCompoundModel();

  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, 
					  G4Nucleus & theNucleus) final;

  virtual G4ReactionProductVector* DeExcite(G4Fragment& aFragment) final;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) final;

  virtual void InitialiseModel() final;
  
  virtual void ModelDescription(std::ostream& outFile) const final;
  virtual void DeExciteModelDescription(std::ostream& outFile) const final;

  //====== obsolete Set methods =======
  void UseHETCEmission();
  void UseDefaultEmission();
  void UseGNASHTransition();
  void UseDefaultTransition();

  //for cross section selection
  void SetOPTxs(G4int opt);

  //for the rest of external choices
  void UseSICB();
  void UseNGB();
  void UseSCO();
  void UseCEMtr();
  //======================================

private:  

  inline 
  void PerformEquilibriumEmission(const G4Fragment & aFragment, 
				  G4ReactionProductVector * theResult) const;

  void PrintWarning(const G4String& mname);

  G4PreCompoundModel(const G4PreCompoundModel &) = delete;
  const G4PreCompoundModel& operator=(const G4PreCompoundModel &right) = delete;
  G4bool operator==(const G4PreCompoundModel &right) const = delete;
  G4bool operator!=(const G4PreCompoundModel &right) const = delete;

  //==============
  // Data Members 
  //==============

  G4PreCompoundEmission*     theEmission;
  G4VPreCompoundTransitions* theTransition;

  const G4ParticleDefinition* proton;
  const G4ParticleDefinition* neutron;

  G4double fLevelDensity;
  G4double fLimitEnergy;

  //for the rest of external choices
  G4bool useSCO;
  G4bool isInitialised;
  G4bool isActive;

  G4int  minZ;
  G4int  minA;

  G4HadFinalState theResult;

};

inline void 
G4PreCompoundModel::PerformEquilibriumEmission(const G4Fragment & aFragment,
					       G4ReactionProductVector * Result) const 
{
  G4ReactionProductVector* theEquilibriumResult = 
    GetExcitationHandler()->BreakItUp(aFragment);
  Result->insert(Result->end(),theEquilibriumResult->begin(), theEquilibriumResult->end());
  delete theEquilibriumResult;
}

#endif

