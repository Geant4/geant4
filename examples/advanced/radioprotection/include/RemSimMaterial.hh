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
// $Id: RemSimMaterial.hh,v 1.5 2004/05/22 12:57:04 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//    **********************************
//    *                                *
//    *      BrachyMaterial.hh          *
//    *                                *
//    **********************************
//
//Code developed by: Susanna Guatelli, guatelli@ge.infn.it
//
//This class manages the elements and materials needed by the simulation
// set-up ...
//
#ifndef RemSimMaterial_H
#define RemSimMaterial_H 1
#include "globals.hh"
class G4Material;

class RemSimMaterial
{ 
public:
  RemSimMaterial();
  ~ RemSimMaterial();

public:
  void  DefineMaterials();
  G4Material* GetMaterial(G4String); //returns the material

private:
  G4Material* matPb;
  G4Material* matAir;
  G4Material* matH2O; 
  G4Material* Al;
  G4Material* nylon;
  G4Material* mylar;
  G4Material* beta;
  G4Material* nextel;
  G4Material* kevlar;
  G4Material* vacuum;
  G4Material* betaCloth; 
  G4Material* eterogeneousNextel;
  G4Material* kevlarVacuum;
  G4Material* polyethylene;
  G4Material* polyacrylate;
  G4Material* evoh;
  G4Material* nomex; 
  G4Material* nomexAir;
  G4Material* kevlarAir;
  G4Material* moon;
};
#endif
