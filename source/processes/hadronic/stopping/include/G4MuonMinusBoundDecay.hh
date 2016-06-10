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
// $Id: G4MuonMinusBoundDecay.hh 69573 2013-05-08 13:35:53Z gcosmo $
//
//-----------------------------------------------------------------------------
//
// GEANT4 Class header file 
//
// File name:  G4MuonMinusBoundDecay
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 24 April 2012 on base of G4MuMinusCaptureAtRest
//
// Class Description: 
//
// Sample probabilities of mu- nuclear capture of decay from K-shell orbit.
// In ApplyYourself method time of projectile is changed taking
//    into account delay of decay or capture. If decay is sampled
//    primary state become stopAndKill, if not - isAlive
//
// Based of reviews:
//   N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
//   B.B.Balashov, G.Ya.Korenman, P.A.Eramgan, Atomizdat, 1978. 
//
//-----------------------------------------------------------------------------
//
// Modifications: 
// 23/04/2013  K.Genser     Made GetMuonCaptureRate and 
//                          GetMuonDecayRate public static
// 04/30/2013  K.Genser     Added GetMuonZeff
//
//
//-----------------------------------------------------------------------------

#ifndef G4MuonMinusBoundDecay_h
#define G4MuonMinusBoundDecay_h 1

#include "globals.hh"
#include "G4Nucleus.hh"
#include "G4Track.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4HadFinalState.hh"
#include "G4HadronicInteraction.hh"
#include "G4DynamicParticle.hh"

class G4MuonMinusBoundDecay : public G4HadronicInteraction
{ 
public:
 
  G4MuonMinusBoundDecay();
 
  ~G4MuonMinusBoundDecay();

  G4HadFinalState* ApplyYourself(const G4HadProjectile &aTrack, 
				 G4Nucleus & targetNucleus );

  void ModelDescription(std::ostream& outFile) const; 

  static G4double GetMuonCaptureRate(G4int Z, G4int A);

  static G4double GetMuonDecayRate(G4int Z);

  static G4double GetMuonZeff(G4int Z);

private:

  inline void AddNewParticle(G4DynamicParticle* dp, G4double time);

  // hide assignment operator as private 
  G4MuonMinusBoundDecay& operator=(const G4MuonMinusBoundDecay &right);
  G4MuonMinusBoundDecay(const G4MuonMinusBoundDecay& );

  G4HadFinalState result;
  G4double fMuMass;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4MuonMinusBoundDecay::AddNewParticle(G4DynamicParticle* dp, G4double time)
{
  G4HadSecondary hs(dp);
  hs.SetTime(time);
  result.AddSecondary(hs);
}

#endif


