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
// $Id: XrayFluoDetectorConstruction.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  20 Aug 2001  Alfonso Mantero   Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoMaterials_hh
#define XrayFluoMaterials_hh 1

#include "globals.hh"
#include "G4Material.hh"

class XrayFluoMaterials
{
public:

  ~XrayFluoMaterials();
  
  static XrayFluoMaterials* GetInstance();
  
  G4Material* GetMaterial(G4String);

  
private:
  
  XrayFluoMaterials();
  
  void CreateMaterials();

private:

  static XrayFluoMaterials* instance;

  G4Material*        dolorite;
  G4Material*        FeMaterial;
  G4Material*        Al;
  G4Material*        Si;
  G4Material*        Cu;
  G4Material*        HPGe;
  G4Material*        materialMg;
  G4Material*        materialNd;
  G4Material*        Sn;
  G4Material*        Ti;
  G4Material*        mars1;
  G4Material*        anorthosite;
  G4Material*        Pb;
  G4Material*        basalt;
  G4Material*        Air;
  G4Material*        Sci;
  G4Material*        Vacuum;

};

#endif
