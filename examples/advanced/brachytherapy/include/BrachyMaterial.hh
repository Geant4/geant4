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
// $Id: BrachyMaterial.hh,v 1.3 2003-05-22 17:20:41 guatelli Exp $
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
#ifndef BrachyMaterial_H
#define BrachyMaterial_H 1
#include "globals.hh"
class G4Material;

class BrachyMaterial
{ public:
  BrachyMaterial();
  ~ BrachyMaterial();

public:
  void  DefineMaterials();

public:
  G4Material* GetMat(G4String); //returns the material
};

#endif
