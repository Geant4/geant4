//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4PHOTONEVAPORATION_HH
#define G4PHOTONEVAPORATION_HH

#include "globals.hh"
#include "G4VPhotonEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4VEmissionProbability.hh"
#include "G4VGammaDeexcitation.hh"

//#define debug

class G4Fragment;

class G4PhotonEvaporation : public G4VPhotonEvaporation, public G4VEvaporationChannel
{

public:

  G4PhotonEvaporation();

  virtual ~G4PhotonEvaporation();
  
  virtual G4FragmentVector* BreakItUp(const G4Fragment& nucleus);

  virtual void Initialize(const G4Fragment& fragment);

  virtual G4FragmentVector* BreakUp(const G4Fragment& nucleus);

  virtual G4double GetEmissionProbability() const;

  virtual void SetEmissionStrategy(G4VEmissionProbability* probAlgorithm);

  void SetVerboseLevel(G4int verbose);

private:  

  G4int _verbose;
  G4bool _myOwnProbAlgorithm;
  G4VEmissionProbability* _probAlgorithm;
  G4VGammaDeexcitation* _discrDeexcitation;
  G4VGammaDeexcitation* _contDeexcitation;
  G4VGammaDeexcitation* _cdDeexcitation;
  G4Fragment _nucleus;
  G4double _gammaE;

  G4PhotonEvaporation(const G4PhotonEvaporation &right);
  
  const G4PhotonEvaporation& operator=(const G4PhotonEvaporation &right);

  // MGP - Check == and != multiple inheritance... must be a mess!
  G4bool operator==(const G4PhotonEvaporation &right) const;
  G4bool operator!=(const G4PhotonEvaporation &right) const;


#ifdef debug
  void CheckConservation(const G4Fragment & theInitialState,
			 G4FragmentVector * Result) const;
#endif
  

};

#endif
