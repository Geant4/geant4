// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPPhotonDist.hh,v 1.4 1999-12-15 14:53:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 
#ifndef G4NeutronHPPhotonDist_h
#define G4NeutronHPPhotonDist_h 1
#include "globals.hh"
#include "g4std/fstream"
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

  G4NeutronHPPhotonDist()
  {
     disType = NULL;
     energy = NULL;
     theYield = NULL;
     thePartialXsec = NULL;
     isPrimary = NULL;
     theShells = NULL;
     theGammas = NULL;
     nNeu = NULL;
     theLegendre = NULL;
     theAngular = NULL;
     distribution = NULL;
     probs = NULL;
     partials = NULL;
     actualMult = NULL;

     theLevelEnergies = NULL;
     theTransitionProbabilities = NULL;
     thePhotonTransitionFraction = NULL;
  }

  ~G4NeutronHPPhotonDist()
  {
     if(disType != NULL) delete [] disType;
     if(energy != NULL) delete [] energy;
     if(theYield != NULL) delete [] theYield;
     if(thePartialXsec != NULL) delete [] thePartialXsec;
     if(isPrimary != NULL) delete [] isPrimary;
     if(theShells != NULL) delete [] theShells;
     if(theGammas != NULL) delete [] theGammas;
     if(nNeu != NULL) delete [] nNeu;
     if(theLegendre != NULL) delete [] theLegendre;
     if(theAngular != NULL) delete [] theAngular;
     if(distribution != NULL) delete [] distribution;
     if(probs != NULL) delete [] probs;
     if(partials != NULL) delete [] partials;
     if(actualMult != NULL) delete [] actualMult;

     if(theLevelEnergies != NULL) delete theLevelEnergies;
     if(theTransitionProbabilities != NULL) delete theTransitionProbabilities;
     if(thePhotonTransitionFraction != NULL) delete thePhotonTransitionFraction;
  }
  
  G4bool InitMean(G4std::ifstream & aDataFile);
    
  void InitAngular(G4std::ifstream & aDataFile);
  
  void InitEnergies(G4std::ifstream & aDataFile);
  
  void InitPartials(G4std::ifstream & aDataFile);
  
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
