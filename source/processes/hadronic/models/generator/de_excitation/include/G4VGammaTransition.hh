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
//      File name:     G4VGammaTransition
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//      
// -------------------------------------------------------------------

#ifndef G4VGAMMATRANSITION_HH
#define G4VGAMMATRANSITION_HH

#include "globals.hh"

class G4VGammaTransition 
{
public:

  G4VGammaTransition() {};

  virtual ~G4VGammaTransition() {};
  
  virtual void SelectGamma() = 0;
  virtual G4double GetGammaEnergy() = 0;
  virtual G4double GetGammaCreationTime() = 0;

//--  virtual G4double GammaEnergy() = 0;

//--  virtual G4double GetEnergyTo() const = 0;

  virtual void SetEnergyFrom(const G4double energy) = 0;

private:  

  G4VGammaTransition(const G4VGammaTransition &right);
  
  const G4VGammaTransition& operator=(const G4VGammaTransition &right);
  G4bool operator==(const G4VGammaTransition &right) const;
  G4bool operator!=(const G4VGammaTransition &right) const;
  
protected:

  G4int _verbose;

};


#endif
