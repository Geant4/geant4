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




