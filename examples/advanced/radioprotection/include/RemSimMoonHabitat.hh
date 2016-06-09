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
//    **************************************
//    *                                    *
//    *    RemSimMoonHabitat.hh            *
//    *                                    *          
//    **************************************
//
// $Id: RemSimMoonHabitat.hh,v 1.4 2005/05/27 14:21:42 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#ifndef RemSimMoonHabitat_h
#define RemSimMoonHabitat_h 1

class RemSimVGeometryComponent;
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class G4VisAttributes;
class RemSimMoonHabitat: public RemSimVGeometryComponent
{
public:
  RemSimMoonHabitat();
  ~RemSimMoonHabitat();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 

private:
  RemSimMaterial* pMaterial;
  G4VPhysicalVolume* moonPhys;
};
#endif
