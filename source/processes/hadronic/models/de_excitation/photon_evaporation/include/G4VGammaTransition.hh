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
// $Id: G4VGammaTransition.hh 96634 2016-04-27 09:31:49Z gcosmo $
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
// 15.04.1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//             Added creation time evaluation for products of evaporation
// 30.10.2010  V.Ivanchenko moved constructor and destructor to the source
//      
// -------------------------------------------------------------------

#ifndef G4VGAMMATRANSITION_HH
#define G4VGAMMATRANSITION_HH 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VGammaTransition 
{
public:

  explicit G4VGammaTransition() {};

  virtual ~G4VGammaTransition() {};
  
  virtual void SelectGamma() = 0;
  virtual G4double GetGammaEnergy() = 0;
  virtual G4double GetGammaCreationTime() = 0;
  virtual void SetEnergyFrom(G4double energy) = 0;

private:  

  G4VGammaTransition(const G4VGammaTransition &right) = delete;
  const G4VGammaTransition& operator=(const G4VGammaTransition &right) = delete;
  G4bool operator==(const G4VGammaTransition &right) const = delete;
  G4bool operator!=(const G4VGammaTransition &right) const = delete;
 
};


#endif





