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
// $Id: HadrontherapyMaterial.hh; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyMaterial_H
#define HadrontherapyMaterial_H 1
#include "globals.hh"
class G4Material;

class HadrontherapyMaterial
{ 
public:
  HadrontherapyMaterial();
  ~ HadrontherapyMaterial();

public:
  void  DefineMaterials();
  G4Material* GetMat(G4String); //returns the material
  
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
G4Material* Vacuum; 
G4Material* bone;
G4Material* muscle;
G4Material* Ta;
G4Material* Brass;
G4Material* Kapton;
G4Material* matAl;
G4Material* matTa;
G4Material* matCu;
G4Material* matSi;
};
#endif
