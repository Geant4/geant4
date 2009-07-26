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
// $Id: G4QFragmentation.hh,v 1.11 2009-07-26 21:14:18 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4QInteractionVector.hh"
#include "G4QPomeron.hh"
#include "G4QPartonPairVector.hh" 
#include "G4QStringVector.hh" 

class G4QFragmentation
{
 public:
  G4QFragmentation();
  G4QFragmentation(const G4QFragmentation &right);
  const G4QFragmentation& operator=(const G4QFragmentation &right);
  virtual ~G4QFragmentation(); 

  int operator==(const G4QFragmentation &right) const;
  int operator!=(const G4QFragmentation &right) const;

  G4QHadronVector* Breeder(const G4QNucleus& aNucleus, const G4QHadron& aPrimary);

  G4bool ExciteDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const; //@@Once
  G4bool ExciteSingDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const;//@@Onc
  std::pair<G4int,G4int> ReducePair(G4int P1, G4int P2) const; // Reduce Q-pairs to singles

  // Static functions
  static void SetParameters(G4int nCM, G4double radNuc, G4double SigPt);

 protected:
  G4bool IsSingleDiffractive()
                  {G4bool result=false; if(G4UniformRand()<1.) result=true; return result;}
  G4int SumPartonPDG(G4int PDG1, G4int PFG2) const;
  G4double ChooseX(G4double Xmin, G4double Xmax) const;
  G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;
  G4int AnnihilationOrder(G4int LS, G4int MS, G4int uP, G4int mP, G4int sP, G4int nP);

 private:
  // static model parameters
  static G4int    nCutMax;                               // Maximum number of Soft Cuts 
  static G4double theNucleonRadius;                      // @@ ? Is that necessary? (M.K.)
  static G4double widthOfPtSquare;                       // width^2 of pt(StringExcitation)

  // Body
  enum {SOFT, DIFFRACTIVE};
};

#endif
