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
// $Id: G4PhotonEvaporationOLD.hh 89158 2015-03-21 19:15:25Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporationOLD
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 23 October 1998
//
//Modifications:
//
// 18 October 2002, Fan Lei (flei@space.qinetiq.com)
//   
//        Implementation of Internal Convertion process in discrete deexcitation
//        The following public methods have been added. 
//
//            void SetICM (G4bool);
//            void RDMForced(G4bool);
//            void SetMaxHalfLife(G4double) ;
//            void SetEOccupancy( G4ElectronOccupancy  eOccupancy) ;
//            G4ElectronOccupancy GetEOccupancy () ;
//            G4int GetVacantShellNumber () { return _vShellNumber;};
//
//        and the following priivate menbers
//
//            G4ElectronOccupancy _eOccupancy;
//            G4int _vShellNumber;
//
// 11 May 2010, V.Ivanchenko added EmittedFragment and BreakUpFragment
//                           methods
//
// 22 October 2015, V.Ivanchenko renamed to G4PhotonEvaporationOLD
// 
// -------------------------------------------------------------------

#ifndef G4PHOTONEVAPORATIONOLD_HH
#define G4PHOTONEVAPORATIONOLD_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"
#include "G4VEmissionProbability.hh"
#include "G4VGammaDeexcitation.hh"
#include "G4ElectronOccupancy.hh"

class G4Fragment;

class G4PhotonEvaporationOLD : public G4VEvaporationChannel {

public:

  G4PhotonEvaporationOLD(const G4String & aName = "");

  virtual ~G4PhotonEvaporationOLD();

  // returns one gamma or e-
  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

  // returns "false", emitted gamma and e- are added to the results
  virtual G4bool 
  BreakUpChain(G4FragmentVector* theResult, G4Fragment* theNucleus);

  // emitted gamma and e- are added to the results
  virtual G4FragmentVector* BreakUpFragment(G4Fragment* theNucleus);

  // emitted gamma, e-, and residual fragment are added to the results
  virtual G4FragmentVector * BreakItUp(const G4Fragment & nucleus);

  // emitted gamma, e-, and residual fragment are added to the results
  virtual G4FragmentVector * BreakUp(const G4Fragment & nucleus);

  virtual G4double GetEmissionProbability(G4Fragment* theNucleus);

  virtual void SetEmissionStrategy(G4VEmissionProbability * probAlgorithm);

  void SetVerboseLevel(G4int verbose);

  virtual void SetICM (G4bool);

  virtual void RDMForced (G4bool);
  
  void SetMaxHalfLife(G4double);

  void SetTimeLimit(G4double value);
 
  void SetEOccupancy( G4ElectronOccupancy  eOccupancy) ;

  inline G4ElectronOccupancy GetEOccupancy () { return eOccupancy;} ;
   
  inline G4int GetVacantShellNumber () { return vShellNumber;};

private:

  G4int verbose;
  G4bool myOwnProbAlgorithm;
  G4VEmissionProbability * probAlgorithm;
  G4VGammaDeexcitation * discrDeexcitation;
  G4VGammaDeexcitation * contDeexcitation;

  G4ElectronOccupancy eOccupancy;
  G4int vShellNumber;

  G4Fragment* nucleus;
  G4double gammaE;

  G4PhotonEvaporationOLD(const G4PhotonEvaporationOLD & right);
  const G4PhotonEvaporationOLD & operator = (const G4PhotonEvaporationOLD & right);

  G4bool operator == (const G4PhotonEvaporationOLD & right) const;
  G4bool operator != (const G4PhotonEvaporationOLD & right) const;

};

#endif



