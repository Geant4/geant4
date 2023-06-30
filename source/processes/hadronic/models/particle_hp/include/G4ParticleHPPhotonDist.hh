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
// Hadronic Process: Very Low Energy Neutron X-Sections
// original by H.P. Wellisch, TRIUMF, 14-Feb-97
//
// 070606 fix for Valgrind error by T. Koi
// 070612 fix memory leaking by T. Koi
// 070615 fix memory leaking by T. Koi
// 080625 fix memory leaking by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPPhotonDist_h
#define G4ParticleHPPhotonDist_h 1

#include "G4Cache.hh"
#include "G4Gamma.hh"
#include "G4InterpolationManager.hh"
#include "G4ParticleHPAngularP.hh"
#include "G4ParticleHPFastLegendre.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4ParticleHPLegendreTable.hh"
#include "G4ParticleHPPartial.hh"
#include "G4ParticleHPVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>

class G4ParticleHPPhotonDist
{
  public:
    G4ParticleHPPhotonDist()
    {
      disType = nullptr;
      energy = nullptr;
      theYield = nullptr;
      thePartialXsec = nullptr;
      theReactionXsec = nullptr;
      isPrimary = nullptr;
      theShells = nullptr;
      theGammas = nullptr;
      nNeu = nullptr;
      theLegendre = nullptr;
      theAngular = nullptr;
      distribution = nullptr;
      probs = nullptr;
      partials = nullptr;
      actualMult.Put(nullptr);

      theLevelEnergies = nullptr;
      theTransitionProbabilities = nullptr;
      thePhotonTransitionFraction = nullptr;
    }

    ~G4ParticleHPPhotonDist()
    {
      delete[] disType;
      delete[] energy;
      delete[] theYield;
      delete[] thePartialXsec;
      //     delete [] theReactionXsec;
      //     DHW: not created in this class
      delete[] isPrimary;
      delete[] theShells;
      delete[] theGammas;
      delete[] nNeu;
      delete[] theAngular;
      delete[] distribution;
      delete[] probs;

      if (theLegendre != nullptr) {
        for (G4int i = 0; i < (nDiscrete2 - nIso); i++)
          if (theLegendre[i] != nullptr) delete[] theLegendre[i];

        delete[] theLegendre;
      }

      if (partials != nullptr) {
        for (G4int i = 0; i < nPartials; i++) {
          delete partials[i];
        }

        delete[] partials;
      }

      delete[] theLevelEnergies;
      delete[] theTransitionProbabilities;
      delete[] thePhotonTransitionFraction;
      if (actualMult.Get() != nullptr) delete actualMult.Get();
    }

    G4bool InitMean(std::istream& aDataFile);

    void InitAngular(std::istream& aDataFile);

    void InitEnergies(std::istream& aDataFile);

    void InitPartials(std::istream& aDataFile, G4ParticleHPVector* theXsec = nullptr);

    G4ReactionProductVector* GetPhotons(G4double anEnergy);

    inline G4double GetTargetMass() { return targetMass; }

    inline G4bool NeedsCascade() { return repFlag == 2; }

    inline G4double GetLevelEnergy() { return theBaseEnergy; }

  private:
    G4int repFlag{0};  // representation as multiplicities or transition probability arrays.
    G4double targetMass{0.0};

    G4int nDiscrete{0};  // number of discrete photons
    G4int* disType;  // discrete, or continuum photons
    G4double* energy;  // photon energies
    G4ParticleHPVector* theYield;  // multiplicity as a function of neutron energy.
    G4ParticleHPVector theTotalXsec;
    G4ParticleHPVector* thePartialXsec;
    G4ParticleHPVector* theReactionXsec;
    G4int* isPrimary;

    G4int isoFlag{0};  // isotropic or not?
    G4int tabulationType{0};
    G4int nDiscrete2{0};
    G4int nIso{0};
    G4double* theShells;
    G4double* theGammas;
    G4int* nNeu;
    G4InterpolationManager theLegendreManager;
    G4ParticleHPLegendreTable** theLegendre;
    G4ParticleHPAngularP** theAngular;

    G4int* distribution;  // not used for the moment.
    G4int nPartials{0};
    G4ParticleHPVector* probs;  // probabilities for the partial distributions.
    G4ParticleHPPartial** partials;  // the partials, parallel to the above

    G4Cache<std::vector<G4int>*> actualMult;

    // for transition prob arrays start
    G4int theInternalConversionFlag{0};
    G4int nGammaEnergies{0};
    G4double theBaseEnergy{0.0};
    G4double* theLevelEnergies;
    G4double* theTransitionProbabilities;
    G4double* thePhotonTransitionFraction;
    // for transition prob arrays end

    G4ParticleHPFastLegendre theLegend;  // fast look-up for leg-integrals
};

#endif
