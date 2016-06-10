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
//                       Rachid Mazini    Rachid.Mazini@cern.ch  
//
//   Language:           C++
//   Tested on:          g++
//   Prerequisites:      None
//   Purpose:            Header file for the Material Expert for the Atlas FCAL
//                       See FCALMaterialConsultant for more details.
//   Developped on:      10-March-2000  M.F. R.M.
//   History:
//
//----------------------------------------------------------------------------

#ifndef FCALMaterialConsultant_h
#define FCALMaterialConsultant_h 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "globals.hh"

//#include "FCALParameters.hh"

class FCALMaterialConsultant
{
public:

  ~FCALMaterialConsultant() {;};
  static FCALMaterialConsultant * GetInstance();
  G4Material * Material(G4String);

public:
  static FCALMaterialConsultant * theFCALMaterialConsultant;
  FCALMaterialConsultant();

private:

  G4Element *elH, *elD, *elHe, *elLi, *elBe, *elC, *elN;
  G4Element *elNe, *elAl, *elFe, *elCu, *elW, *elPb, *elU, *elO;
  G4Element *elCa, *elNa, *elSi;

  G4Material *Aluminium, *Iron, *Copper, *Tungsten, *Lead;
  G4Material *Air, *RhoaCell, *Vacuum, *CO2, *ArgonGas, *ShieldingConcrete;
  G4Material *Polystyrene, *StainlessSteel, *Nickel, *LiquidArgon;
  G4Material *Kapton, *FCAL1CuArKap, *FCAL1CuAr, *FCAL2CuArKap;
  G4Material *FCAL2WFeNi, *FCAL2WFeNiCuAr, *MWPCArCO2;
 
};

#endif   /* FCALMaterialConsultant.hh */




