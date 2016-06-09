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


// this :
#include <HEPVis/actions/SoAlternateRepAction.h>

// Inventor :
#include <Inventor/nodes/SoNode.h>
#include <Inventor/elements/SoSwitchElement.h>
#include <Inventor/elements/SoCoordinateElement.h>

SO_ACTION_SOURCE(SoAlternateRepAction)

//////////////////////////////////////////////////////////////////////////////
void SoAlternateRepAction::initClass(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_INIT_CLASS(SoAlternateRepAction,SoAction);

  SO_ENABLE(SoAlternateRepAction,SoSwitchElement); //SGI wants that.
  SO_ENABLE(SoAlternateRepAction,SoCoordinateElement); //For SoMarkerSet.

  SO_ACTION_ADD_METHOD(SoNode,nodeAction);
}
//////////////////////////////////////////////////////////////////////////////
SoAlternateRepAction::SoAlternateRepAction(
)
:fGenerate(FALSE)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  SO_ACTION_CONSTRUCTOR(SoAlternateRepAction);
}
//////////////////////////////////////////////////////////////////////////////
SoAlternateRepAction::~SoAlternateRepAction(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
}
//////////////////////////////////////////////////////////////////////////////
void SoAlternateRepAction::setGenerate(
 SbBool aGenerate
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fGenerate = aGenerate;
}
//////////////////////////////////////////////////////////////////////////////
SbBool SoAlternateRepAction::getGenerate(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  return fGenerate;
}

// Inventor :
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoCoordinate3.h>

// HEPVisUtils :
//#include <HEPVis/nodes/SoHighlightMaterial.h>
//#include <HEPVis/nodes/SoSplineCurve.h>

#define IF_CLASS(aClass) \
  if(aNode->isOfType(aClass::getClassTypeId())) {\
    /*printf("debug : SoAlternateRepAction::nodeAction : %s generate %d\n",#aClass,This->fGenerate);*/\
    aClass* node = (aClass*)aNode;\
    if(This->fGenerate==TRUE) {\
      if(node->alternateRep.getValue()==NULL) {\
        node->generateAlternateRep();\
        /* Then go down if needed : */\
        SoNode* altRep = node->alternateRep.getValue();\
        if((altRep!=NULL) && altRep->isOfType(SoGroup::getClassTypeId()))\
          altRep->doAction(This);\
      }\
    } else {\
      /* First go down if needed : */\
      SoNode* altRep = node->alternateRep.getValue();\
      if((altRep!=NULL) && altRep->isOfType(SoGroup::getClassTypeId()))\
        altRep->doAction(This);\
      /* Then clear : */\
      node->clearAlternateRep();\
    }\
  }

//////////////////////////////////////////////////////////////////////////////
void SoAlternateRepAction::nodeAction(
 SoAction* aThis
,SoNode* aNode
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  /*printf("debug : SoAlternateRepAction::nodeAction : %s %s\n",
    aNode->getName().getString(),
    aNode->getTypeId().getName().getString());*/
  //SoAlternateRepAction* This = (SoAlternateRepAction*)aThis;

  if(false) {}
  // HEPVisUtils :
  //else IF_CLASS(SoHighlightMaterial)
  //else IF_CLASS(SoSplineCurve)

  // Inventor :
  else if(aNode->isOfType(SoGroup::getClassTypeId())) {
    aNode->doAction(aThis);
  }
  else if(aNode->isOfType(SoCoordinate3::getClassTypeId())) { //For SoMakerSet.
    aNode->doAction(aThis);
  }
  // Else :
  else if(aNode->getField("alternateRep")!=NULL) {
    aNode->doAction(aThis);
  }

}

#endif
