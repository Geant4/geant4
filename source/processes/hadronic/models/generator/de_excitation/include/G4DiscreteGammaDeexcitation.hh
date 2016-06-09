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
// * authors in the GEANT4 collaboration.                             *
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

class G4DiscreteGammaDeexcitation : public G4VGammaDeexcitation {
public:

  // Constructor
  G4DiscreteGammaDeexcitation();

  // Destructor
  ~G4DiscreteGammaDeexcitation();

  // Functions

public:

  virtual G4VGammaTransition * CreateTransition();

  virtual G4bool CanDoTransition() const;

  void SetICM(G4bool hl) { _icm = hl; };

  void SetRDM(G4bool hl) { _rdm = hl; };
  
  void SetHL(G4double hl) { _max_hl = hl; };

private:

  G4int _nucleusZ;
  G4int _nucleusA;
  G4double _tolerance;
  G4double _max_hl;
  G4bool _icm;
  G4bool _rdm;
  G4NuclearLevelManager * _levelManager;
};


#endif



