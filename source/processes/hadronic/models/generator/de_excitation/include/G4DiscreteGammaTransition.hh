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
//      File name:     G4DiscreteGammaTransition
//
//      Author:        Maria Grazia Pia   (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//      
// -------------------------------------------------------------------

#ifndef G4DiscreteGammaTransition_hh
#define G4DiscreteGammaTransition_hh

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4NuclearLevel.hh"

class G4DiscreteGammaTransition : public G4VGammaTransition
{
public:

  // Constructor
  G4DiscreteGammaTransition(const G4NuclearLevel& level);

  // Destructor
  ~G4DiscreteGammaTransition();

  // Functions

public:

//--  virtual G4double GammaEnergy();
//--  virtual G4double GetEnergyTo() const;
  virtual void SetEnergyFrom(const G4double energy);
  virtual G4double GetGammaEnergy();
  virtual G4double GetGammaCreationTime();
  virtual void SelectGamma();
  
private:
  G4double _gammaEnergy;
  G4NuclearLevel _level;     
  G4double _excitation;
  G4double _gammaCreationTime;

};


#endif
