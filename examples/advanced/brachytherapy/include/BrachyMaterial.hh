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
//
//    **********************************
//    *                                *
//    *      BrachyMaterial.hh          *
//    *                                *
//    **********************************
//
//Code developed by: Susanna Guatelli, Albert Le
//
// This class manages the elements and materials
//
#ifndef BrachyMaterial_H
#define BrachyMaterial_H 1

#include "globals.hh"
class G4Material;

class BrachyMaterial
{ 
public:
  BrachyMaterial();
  ~ BrachyMaterial();

public:
  void  DefineMaterials();
  G4Material* GetMat(G4String); //returns the material

private:
  G4Material* fW; 
  G4Material* fPlexiglass;
  G4Material* fPb;
  G4Material* fIr192;
  G4Material* fTi;
  G4Material* fAir;
  G4Material* fH2O;
  G4Material* fSoft;
  G4Material* fSteel;
  G4Material* fMat304steel;
  G4Material* fAu;
  G4Material* fI; 
  G4Material* fCeramic;
  G4Material* fVacuum; 
  G4Material* fBone;
  G4Material* fMuscle;
  G4Material* fAg;
};
#endif
