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
// $Id: Shielding.hh 83413 2014-08-21 15:21:43Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2010  Tatsumi Koi, Gunter Folger
//   created from FTFP_BERT
//
// Modified:
//
// 2014.08.05 K.L.Genser added provision for Hadronic Physics Variant
//
//----------------------------------------------------------------------------
//
#ifndef TShielding_h
#define TShielding_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TShielding: public T
{
public:
  //TShielding(G4int ver = 1);
  //TShielding(G4int ver = 1, G4bool rad=true, G4bool lend=false);
  explicit TShielding( G4int verbose = 1 , G4String low_energy_neutron_model = "HP", 
                       G4String HadrPhysVariant = "");
  virtual ~TShielding();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
};
#include "Shielding.icc"
typedef TShielding<G4VModularPhysicsList> Shielding;

#endif



