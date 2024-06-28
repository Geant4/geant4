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
//  27/03/2004 : G.Barrand : creation.

#ifndef HEPVis_SoStyleCache_h
#define HEPVis_SoStyleCache_h

#include <Inventor/nodes/SoGroup.h>

//#include <HEPVis/SbStyle.h>

class SoMaterial;
class SoDrawStyle;
class SoLightModel;
class SoResetTransform;
class SbColor;
 
#define SoStyleCache Geant4_SoStyleCache

class SoStyleCache : public SoGroup {
public:
  SoStyleCache();
protected:
  ~SoStyleCache();
public:
  SoMaterial* getMaterial(const SbColor&,float = 0); 
  SoMaterial* getMaterial(float,float,float,float = 0); 
  //SoDrawStyle* getLineStyle(SbLineStyle,float = 0); 
  SoDrawStyle* getLineStyle(unsigned short = 0xFFFF,float = 0); 
  SoLightModel* getLightModelPhong(); 
  SoLightModel* getLightModelBaseColor(); 
  SoResetTransform* getResetTransform();
private:
  SoGroup* fMaterials;
  SoGroup* fLineStyles;
  SoGroup* fLightModels;
  SoResetTransform* fResetTransform;
};

#endif
