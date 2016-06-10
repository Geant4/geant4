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
// $Id: G4DiscreteGammaDeexcitation.hh 85841 2014-11-05 15:35:06Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 23 October 1998
//
//      Modifications:
//
//      8 March 2001, Fan Lei (flei@space.qinetiq.com)
//      added the following as part if the IC implementation
//        void SetICM(G4bool hl) { _icm = hl; };
//        void SetRDM(G4bool hl) { _rdm = hl; };
//        void SetHL(G4double hl) { _max_hl = hl; };
//       private:
//        G4double _max_hl;
//        G4bool _icm;
//        G4bool _rdm;
// 
//
// -------------------------------------------------------------------
//

#ifndef G4DiscreteGammaDeexcitation_hh
#define G4DiscreteGammaDeexcitation_hh

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"

#include "G4DiscreteGammaTransition.hh"
#include "G4Fragment.hh"

class G4NuclearLevelManager;
class G4NuclearLevelStore;
class G4DiscreteGammaTransition;

class G4DiscreteGammaDeexcitation : public G4VGammaDeexcitation {
public:

  G4DiscreteGammaDeexcitation();

  virtual ~G4DiscreteGammaDeexcitation();

  virtual G4bool CanDoTransition(G4Fragment* nucleus);

  inline void SetICM(G4bool hl)  { icm = hl; };

  inline void SetRDM(G4bool hl)  { rdm = hl; };
  
  inline void SetHL(G4double hl) { maxhl = hl; };

private:

  G4int nucleusZ;
  G4int nucleusA;
  G4double maxhl;
  G4bool icm;
  G4bool rdm;
  G4NuclearLevelStore* store;
  G4NuclearLevelManager* levelManager;
  G4DiscreteGammaTransition* dtransition; 
};


#endif



