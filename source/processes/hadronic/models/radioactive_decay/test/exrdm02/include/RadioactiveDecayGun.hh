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
#ifndef RadioactiveDecayGun_h
#define RadioactiveDecayGun_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              RadioactiveDecayGun.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleGun.hh"
#include "G4Event.hh"
#include "globals.hh"

#include "Nucleus.hh"
#include "RadioactiveDecayGunmessenger.hh"

class RadioactiveDecayGunmessenger;
////////////////////////////////////////////////////////////////////////////////
//
class RadioactiveDecayGun : public G4ParticleGun
{
  // class description
  // The RadioactiveDecayGun is an inherited version of G4ParticleGun
  // to allow user to specify an isotope as the initial tracking particle.
  // class description - end

public:
  RadioactiveDecayGun();
  ~RadioactiveDecayGun();

public: // with description

  void  SetNucleus(Nucleus theIon1);
  // Sets the isotope.
  //
    inline Nucleus GetNucleus() {return theIon;}
  // Returns the specified isotope.
  //
private:

  RadioactiveDecayGunmessenger  *theRadioactiveDecayGunMessenger;

  Nucleus theIon;

};
#endif
////////////////////////////////////////////////////////////////////////////////




