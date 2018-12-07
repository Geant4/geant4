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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  20 Aug 2001  Alfonso Mantero   Created
//
// -------------------------------------------------------------------

#ifndef XrayFluoNistMaterials_hh
#define XrayFluoNistMaterials_hh 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

class XrayFluoNistMaterials
{
public:

  ~XrayFluoNistMaterials();
  
  static XrayFluoNistMaterials* GetInstance();
  
  G4Material* GetMaterial(G4String);

  
private:
  
  XrayFluoNistMaterials();
  
  void CreateMaterials();

private:

  static XrayFluoNistMaterials* instance;

  G4NistManager*     nistMan;

  G4Material*        copper;
  G4Material*        SiLi;
  G4Material*        dolorite;
  G4Material*        HPGe;
  G4Material*        mars1;
  G4Material*        hawaiianWD;
  G4Material*        hawaiianRF;
  G4Material*        anorthosite;
  G4Material*        basalt;
  G4Material*        gabbro;
  G4Material*        gabbroWD;
  G4Material*        gabbroRF;
  G4Material*        Air;
  G4Material*        Sci;
  G4Material*        Vacuum;
  G4Material*        madaBasalt;
  G4Material*        icelandicBasalt;
  G4Material*        icelandicWD;
  G4Material*        icelandicRF;
  G4Material*        GaAs;
  G4Material*        galactic;

};

#endif
