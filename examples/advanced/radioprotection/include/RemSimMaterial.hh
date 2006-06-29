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
//
// $Id: RemSimMaterial.hh,v 1.6 2006-06-29 16:22:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
