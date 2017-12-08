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
// $Id: ShieldingLEND.hh 83151 2014-08-06 16:22:30Z klg $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 7 Nov 2017  Tatsumi Koi
//   created from Shielding
//
// Modified:
//
//----------------------------------------------------------------------------
//
#ifndef TShieldingLEND_h
#define TShieldingLEND_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TShieldingLEND: public T
{
public:
  //TShieldingLEND(G4int ver = 1);
  //TShieldingLEND(G4int ver = 1, G4bool rad=true, G4bool lend=false);
  explicit TShieldingLEND( G4int verbose = 1 , G4String low_energy_neutron_model = "LEND", 
                       G4String HadrPhysVariant = "");
  virtual ~TShieldingLEND();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "ShieldingLEND.icc"
typedef TShieldingLEND<G4VModularPhysicsList> ShieldingLEND;

#endif



