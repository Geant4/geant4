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
// $Id: RemSimMaterial.hh,v 1.2 2004-02-03 09:16:45 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//    **********************************
//    *                                *
//    *      BrachyMaterial.hh          *
//    *                                *
//    **********************************
//
//Code developed by: Susanna Guatelli
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
  G4Material* matW; 
  G4Material* matplexiglass;
  G4Material* matPb;
  G4Material* matir192;
  G4Material* Titanium;
  G4Material* matAir;
  G4Material* matH2O;
  G4Material* soft;
  G4Material* matsteel;
  G4Material* gold;
  G4Material* matI; 
  G4Material* ceramic;
  G4Material*Vacuum; 
  G4Material* bone;
  G4Material* muscle;
  G4Material* Ar;
  G4Material* Al;
  G4Material* nylon;
  G4Material* mylar;
  G4Material* beta;
  G4Material* kevlar;
};
#endif
