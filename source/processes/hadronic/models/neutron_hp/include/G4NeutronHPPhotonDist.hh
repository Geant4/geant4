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
// $Id: G4NeutronHPPhotonDist.hh,v 1.11 2007-06-06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 
#ifndef G4NeutronHPPhotonDist_h
#define G4NeutronHPPhotonDist_h 1
#include "globals.hh"
#include <fstream>
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
     disType = 0;
     energy = 0;
     theYield = 0;
     thePartialXsec = 0;
     isPrimary = 0;
     theShells = 0;
     theGammas = 0;
     nNeu = 0;
     theLegendre = 0;
     theAngular = 0;
     distribution = 0;
     probs = 0;
     partials = 0;
     actualMult = 0;

     theLevelEnergies = 0;
     theTransitionProbabilities = 0;
     thePhotonTransitionFraction = 0;
     
     repFlag = 0;
     nDiscrete = 0;
     nPartials = 0;

  }

  ~G4NeutronHPPhotonDist()
  {
     if(disType != 0) delete [] disType;
     if(energy != 0) delete [] energy;
     if(theYield != 0) delete [] theYield;
     if(thePartialXsec != 0) delete [] thePartialXsec;
     if(isPrimary != 0) delete [] isPrimary;
     if(theShells != 0) delete [] theShells;
     if(theGammas != 0) delete [] theGammas;
     if(nNeu != 0) delete [] nNeu;
     if(theLegendre != 0) delete [] theLegendre;
     if(theAngular != 0) delete [] theAngular;
     if(distribution != 0) delete [] distribution;
     if(probs != 0) delete [] probs;
     if(partials != 0) delete [] partials;
     if(actualMult != 0) delete [] actualMult;

     if(theLevelEnergies != 0) delete theLevelEnergies;
     if(theTransitionProbabilities != 0) delete theTransitionProbabilities;
     if(thePhotonTransitionFraction != 0) delete thePhotonTransitionFraction;
  }
  
  G4bool InitMean(std::ifstream & aDataFile);
    
  void InitAngular(std::ifstream & aDataFile);
  
  void InitEnergies(std::ifstream & aDataFile);
  
  void InitPartials(std::ifstream & aDataFile);
  
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
