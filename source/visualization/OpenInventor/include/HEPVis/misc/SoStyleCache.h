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
