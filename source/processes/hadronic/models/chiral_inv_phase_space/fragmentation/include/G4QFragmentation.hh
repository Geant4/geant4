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
// $Id$
//
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//                 History: 
//     Created by Mikhail Kossov, October 2006
//     CHIPS QGS fragmentation class
//     For comparison mirror member functions are taken from G4 classes:
//     G4QGSParticipants
//     G4QGSModels
//     G4ExcitedStringDecay
// -----------------------------------------------------------------------------
// Short description: CHIPS string fragmentation class
// -----------------------------------------------------------------------------
//
#ifndef G4QFragmentation_h
#define G4QFragmentation_h 1

#include "globals.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4QNucleus.hh"
#include "G4Quasmon.hh"
#include "G4QHadronVector.hh"
#include "G4QEnvironment.hh"
#include "G4QInteractionVector.hh"
#include "G4QProbability.hh"
#include "G4QPartonPairVector.hh" 
#include "G4QuasiFreeRatios.hh"
#include "G4QDiffractionRatio.hh"
#include "G4QStringVector.hh" 

class G4QFragmentation
{
 public:
  G4QFragmentation(const G4QNucleus& aNucleus, const G4QHadron& aPrimary);
  ~G4QFragmentation(); 

  G4QHadronVector* Fragment(); // Calls Breeder & fragments Quasmons inside ResidualNucleus

  // Static functions
  static void SetParameters(G4int nC, G4double strTens, G4double tubeDens, G4double SigPt);

 protected:
  G4bool ExciteDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const; //@@Once
  G4bool ExciteSingDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const;//@@Onc
  std::pair<G4int,G4int> ReducePair(G4int P1, G4int P2) const; // Reduce Q-pairs to singles
  void Breeder(); // String fragmentation algoritm, which makes Hadrons & Quasmons from Str
  // @@ At present always SingleDiffractive, so this is a fake member function
  G4bool IsSingleDiffractive() {G4bool r=false; if(G4UniformRand()<1.) r=true; return r;}
  G4int SumPartonPDG(G4int PDG1, G4int PFG2) const;
  G4double ChooseX(G4double Xmin, G4double Xmax) const;
  G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;
  G4int AnnihilationOrder(G4int LS, G4int MS, G4int uP, G4int mP, G4int sP, G4int nP);
  void SwapPartons(); // Try to swap partons of strings, if one of strings has negative M2
  void EvaporateResidual(G4QHadron* hadrNuc); //Evaporate hadrNucleus, NuclFrag->theResult
 private:
  enum {SOFT, DIFFRACTIVE};
  // static model parameters
  static G4int    nCutMax;                               // Maximum number of Soft Cuts 
  static G4double stringTension;                         // String Tension to absorb energy
  static G4double tubeDensity;                           // Nucleon density in the FluxTube
  static G4double widthOfPtSquare;                       // width^2 of pt(StringExcitation)

  // Body
  G4QNucleus      theNucleus;                            // TargetNucleus moving fromLStoCM
  G4QStringVector strings;                               // Vector of created strings
  G4QuasmonVector theQuasmons;                           // Strings converter to Quasmons
  G4QHadronVector* theResult;                            // Pointer to the OUTPUT Result
  G4double        maxEn;                                 // Energy absorbed by the nucleus
  G4double        maxNuc;                                // #0fNucleons in the Flux Tube
  G4QuasiFreeRatios* theQuasiElastic;                    // For CHIPS Quasi-Elastic
  G4QDiffractionRatio* theDiffraction;                   // For CHIPS Diffraction
  G4QCHIPSWorld*  theWorld;                              // Pointer to the CHIPS World
};

#endif
