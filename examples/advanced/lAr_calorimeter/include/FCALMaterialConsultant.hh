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

  FCALMaterialConsultant();

public:

  ~FCALMaterialConsultant() {theFCALMaterialConsultant = NULL;};
  static FCALMaterialConsultant * construct();
  G4Material * Material(G4String);

public:

  static FCALMaterialConsultant * theFCALMaterialConsultant;

private:

  G4Element *elH, *elD, *elHe, *elLi, *elBe, *elC, *elN, *elNi;
  G4Element *elNe, *elAl, *elFe, *elCu, *elW, *elPb, *elU, *elO;
  G4Element *elCa, *elNa, *elSi;

  G4Material *Aluminium, *Iron, *Copper, *Tungsten, *Lead;
  G4Material *Air, *RhoaCell, *Vacuum, *CO2, *ArgonGas, *ShieldingConcrete;
  G4Material *Polystyrene, *StainlessSteel, *Nickel, *LiquidArgon;
  G4Material *Kapton, *FCAL1CuArKap, *FCAL1CuAr, *FCAL2CuArKap;
  G4Material *FCAL2WFeNi, *FCAL2WFeNiCuAr, *MWPCArCO2;
 
};

#endif   /* FCALMaterialConsultant.hh */




