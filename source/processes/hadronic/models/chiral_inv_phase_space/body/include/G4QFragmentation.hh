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
// $Id: G4QFragmentation.hh,v 1.6 2009-07-02 07:17:09 mkossov Exp $
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
#include "G4QNucleus.hh"
#include "Randomize.hh"
#include "G4QHadronVector.hh"
#include "G4ShortLivedConstructor.hh"
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

  G4QHadronVector* Scatter(const G4QNucleus& aNucleus, const G4QHadron& aPrimary);

  G4bool ExciteDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const;
  G4bool ExciteSingDiffParticipants(G4QHadron* aPartner, G4QHadron* bPartner) const;

  G4QString* BuildString(G4QPartonPair* aPair)
   {return new G4QString(aPair->GetParton1(), aPair->GetParton2(), aPair->GetDirection());}

  // Static functions
  static void SetParameters(G4int nCM, G4double thresh, G4double QGSMth, G4double radNuc,
                            G4double SigPt, G4double extraM, G4double minM);
 protected:
  G4bool IsSingleDiffractive()
                  {G4bool result=false; if(G4UniformRand()<1.) result=true; return result;}
  G4bool EnergyAndMomentumCorrector(G4QHadronVector* Output, G4LorentzVector& TotaMom);   
  G4double ChooseX(G4double Xmin, G4double Xmax) const;
  G4ThreeVector GaussianPt(G4double widthSquare, G4double maxPtSquare) const;

 private:
  // static model parameters
  static G4int    nCutMax; 
  static G4double ThresholdParameter; 
  static G4double QGSMThershold; 
  static G4double theNucleonRadius;
  // Parameters of diffractional fragmentation
  static G4double widthOfPtSquare;                       // width^2 of pt(StringExcitation)
  static G4double minExtraMass;                          // minimum excitation mass 
  static G4double minmass;                               // mean pion pt-mass for Xmin

  // Body
  enum {SOFT, DIFFRACTIVE};
};

#endif
