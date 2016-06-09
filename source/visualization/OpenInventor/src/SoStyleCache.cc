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
#ifdef G4VIS_BUILD_OI_DRIVER

/*----------------------------HEPVis----------------------------------------*/
/*                                                                          */
/* Node:             SoStyleCache                                           */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/misc/SoStyleCache.h>

#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoResetTransform.h>

//////////////////////////////////////////////////////////////////////////////
SoStyleCache::SoStyleCache(
) 
:fMaterials(0)
,fLineStyles(0)
,fLightModels(0)
,fResetTransform(0)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fMaterials = new SoGroup;
  addChild(fMaterials);
  fLineStyles = new SoGroup;
  addChild(fLineStyles);
  fLightModels = new SoGroup;
  addChild(fLightModels);
  fResetTransform = new SoResetTransform;
  addChild(fResetTransform);
}
//////////////////////////////////////////////////////////////////////////////
SoStyleCache::~SoStyleCache(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
SoMaterial* SoStyleCache::getMaterial(
 const SbColor& aRGB
,float aTransparency
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  int number = fMaterials->getNumChildren();
  for(int index=0;index<number;index++) { 
    SoMaterial* material = (SoMaterial*)fMaterials->getChild(index);
    if( (material->diffuseColor[0]==aRGB) &&
        (material->transparency[0]==aTransparency) ) {
      return material;
    }
  }
  SoMaterial* material = new SoMaterial;
  material->diffuseColor.setValue(aRGB);
  material->transparency.setValue(aTransparency);
  fMaterials->addChild(material);
  return material;
}
//////////////////////////////////////////////////////////////////////////////
SoMaterial* SoStyleCache::getMaterial(
 float aRed
,float aGreen
,float aBlue
,float aTransparency
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SbColor aRGB(aRed,aGreen,aBlue);
  int number = fMaterials->getNumChildren();
  for(int index=0;index<number;index++) { 
    SoMaterial* material = (SoMaterial*)fMaterials->getChild(index);
    if( (material->diffuseColor[0]==aRGB) &&
        (material->transparency[0]==aTransparency) ) {
      return material;
    }
  }
  SoMaterial* material = new SoMaterial;
  material->diffuseColor.setValue(aRGB);
  material->transparency.setValue(aTransparency);
  fMaterials->addChild(material);
  return material;
}
/*
//////////////////////////////////////////////////////////////////////////////
SoDrawStyle* SoStyleCache::getLineStyle(
 SbLineStyle aStyle
,float aWidth
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  unsigned short pattern = 0xFFFF;
  switch(aStyle) {
  case SbLineDashed:
    pattern = 0x00FF;
    break;
  case SbLineDotted:
    pattern = 0x0101;
    break;
  case SbLineDashDotted:
    pattern = 0x1C47;
    break;
  default: //SbLineSolid:
    pattern = 0xFFFF;
    break;
  }
  int number = fLineStyles->getNumChildren();
  for(int index=0;index<number;index++) { 
    SoDrawStyle* drawStyle = (SoDrawStyle*)fLineStyles->getChild(index);
    if( (drawStyle->style.getValue()==SoDrawStyle::LINES) &&
        (drawStyle->lineWidth.getValue()==aWidth) &&
        (drawStyle->linePattern.getValue()==pattern) ) {
      return drawStyle;
    }
  }
  SoDrawStyle* drawStyle = new SoDrawStyle;
  drawStyle->style.setValue(SoDrawStyle::LINES);
  drawStyle->lineWidth.setValue(aWidth);
  drawStyle->linePattern.setValue(pattern);
  fLineStyles->addChild(drawStyle);
  return drawStyle;
}
*/
//////////////////////////////////////////////////////////////////////////////
SoDrawStyle* SoStyleCache::getLineStyle(
 unsigned short aPattern
,float aWidth
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  int number = fLineStyles->getNumChildren();
  for(int index=0;index<number;index++) { 
    SoDrawStyle* drawStyle = (SoDrawStyle*)fLineStyles->getChild(index);
    if( (drawStyle->style.getValue()==SoDrawStyle::LINES) &&
        (drawStyle->lineWidth.getValue()==aWidth) &&
        (drawStyle->linePattern.getValue()==aPattern) ) {
      return drawStyle;
    }
  }
  SoDrawStyle* drawStyle = new SoDrawStyle;
  drawStyle->style.setValue(SoDrawStyle::LINES);
  drawStyle->lineWidth.setValue(aWidth);
  drawStyle->linePattern.setValue(aPattern);
  fLineStyles->addChild(drawStyle);
  return drawStyle;
}
//////////////////////////////////////////////////////////////////////////////
SoLightModel* SoStyleCache::getLightModelPhong(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SoLightModel* lightModel = new SoLightModel;
  lightModel->model.setValue(SoLightModel::PHONG);
  fLightModels->addChild(lightModel);
  return lightModel;
}
//////////////////////////////////////////////////////////////////////////////
SoLightModel* SoStyleCache::getLightModelBaseColor(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SoLightModel* lightModel = new SoLightModel;
  lightModel->model.setValue(SoLightModel::BASE_COLOR);
  fLightModels->addChild(lightModel);
  return lightModel;
}
//////////////////////////////////////////////////////////////////////////////
SoResetTransform* SoStyleCache::getResetTransform(
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fResetTransform;
}

#endif
