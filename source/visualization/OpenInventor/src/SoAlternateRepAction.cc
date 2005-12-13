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
