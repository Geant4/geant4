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
#include <iostream.h>
#include <vector.h>

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




