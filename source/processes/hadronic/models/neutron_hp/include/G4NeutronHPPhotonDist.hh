// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPPhotonDist.hh,v 1.2 1999-06-29 18:44:12 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 
#ifndef G4NeutronHPPhotonDist_h
#define G4NeutronHPPhotonDist_h 1
#include "globals.hh"
#include <fstream.h>
#include "G4ios.hh"
#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "G4NeutronHPLegendreTable.hh"
#include "G4NeutronHPAngularP.hh"
#include "G4NeutronHPPartial.hh"
#include "G4NeutronHPFastLegendre.hh"
#include "G4NeutronHPInterpolator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4Gamma.hh"
#include "G4InterpolationManager.hh"

class G4NeutronHPPhotonDist
{
public:

  G4NeutronHPPhotonDist();

  ~G4NeutronHPPhotonDist();
  
  G4bool InitMean(ifstream & aDataFile);
    
  void InitAngular(ifstream & aDataFile);
  
  void InitEnergies(ifstream & aDataFile);
  
  void InitPartials(ifstream & aDataFile);
  
  G4ReactionProductVector * GetPhotons(G4double anEnergy);
  
  inline G4double GetTargetMass() {return targetMass;}
  
  inline G4bool NeedsCascade() {return repFlag==2;}
  
  inline G4double GetLevelEnergy() {return theBaseEnergy;}

private:

   G4int repFlag;  //representation as multiplicities or transition probability arrays.
   G4double targetMass;
   
   G4int nDiscrete;  //number of discrete photons 
   G4int * disType;  // discrete, or continuum photons
   G4double * energy;  // photon energies
   G4NeutronHPVector * theYield; // multiplicity as a function of neutron energy.
   G4NeutronHPVector theTotalXsec;
   G4NeutronHPVector * thePartialXsec;
   G4int * isPrimary;
  
   G4int isoFlag; // isotropic or not?
   G4int tabulationType;
   G4int nDiscrete2;
   G4int nIso;
   G4double * theShells;
   G4double * theGammas;
   G4int * nNeu;
   G4InterpolationManager theLegendreManager;
   G4NeutronHPLegendreTable ** theLegendre;
   G4NeutronHPAngularP ** theAngular;
   
   G4int * distribution; // not used for the moment.                                 
   G4int nPartials;
   G4NeutronHPVector *  probs; // probabilities for the partial distributions.
   G4NeutronHPPartial ** partials; // the partials, parallel to the above

   G4int * actualMult;
   
    // for transition prob arrays start
   G4int theInternalConversionFlag;
   G4int nGammaEnergies;
   G4double theBaseEnergy;
   G4double * theLevelEnergies;
   G4double * theTransitionProbabilities;
   G4double * thePhotonTransitionFraction;
    // for transition prob arrays end

   G4NeutronHPFastLegendre theLegend; // fast look-up for leg-integrals
   G4NeutronHPInterpolator theInt; // interpolation
};

#endif
