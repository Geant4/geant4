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
// $Id: G4PhotonEvaporation.hh,v 1.8 2010-11-17 16:50:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporation
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
// -------------------------------------------------------------------

#ifndef G4PHOTONEVAPORATION_HH
#define G4PHOTONEVAPORATION_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"
#include "G4VEmissionProbability.hh"
#include "G4VGammaDeexcitation.hh"
#include "G4ElectronOccupancy.hh"

class G4Fragment;

class G4PhotonEvaporation : public G4VEvaporationChannel {

public:

    G4PhotonEvaporation();

    virtual ~G4PhotonEvaporation();

    virtual void Initialize(const G4Fragment & fragment);

    virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus);

    virtual G4FragmentVector* BreakUpFragment(G4Fragment* theNucleus);

    virtual G4FragmentVector * BreakItUp(const G4Fragment & nucleus);

    virtual G4FragmentVector * BreakUp(const G4Fragment & nucleus);

    virtual G4double GetEmissionProbability() const;

    virtual void SetEmissionStrategy(G4VEmissionProbability * probAlgorithm);

    void SetVerboseLevel(G4int verbose);

    void SetICM (G4bool);

    void RDMForced (G4bool);
  
    void SetMaxHalfLife(G4double) ;
 
    void SetEOccupancy( G4ElectronOccupancy  eOccupancy) ;

    G4ElectronOccupancy GetEOccupancy () { return _eOccupancy;} ;
   
    G4int GetVacantShellNumber () { return _vShellNumber;};

private:

    G4int _verbose;
    G4bool _myOwnProbAlgorithm;
    G4VEmissionProbability * _probAlgorithm;
    G4VGammaDeexcitation * _discrDeexcitation;
    G4VGammaDeexcitation * _contDeexcitation;

    G4ElectronOccupancy _eOccupancy;
    G4int _vShellNumber;

    G4Fragment* _nucleus;
    G4double _gammaE;

    G4PhotonEvaporation(const G4PhotonEvaporation & right);
    const G4PhotonEvaporation & operator = (const G4PhotonEvaporation & right);

    G4bool operator == (const G4PhotonEvaporation & right) const;
    G4bool operator != (const G4PhotonEvaporation & right) const;

  //#ifdef debug
  //  void CheckConservation(const G4Fragment & theInitialState, G4FragmentVector * Result) const;
  //#endif


};

#endif



