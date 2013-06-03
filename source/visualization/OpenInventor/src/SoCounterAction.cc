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
/* Node:             SoCounterAction                                        */
/* Author:           Guy Barrand                                            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

// this :
#include <HEPVis/actions/SoCounterAction.h>

#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodekits/SoBaseKit.h>
#include <Inventor/elements/SoSwitchElement.h>

SO_ACTION_SOURCE(SoCounterAction)

void SoCounterAction::initClass(void){
  SO_ACTION_INIT_CLASS(SoCounterAction,SoAction);

  SO_ENABLE(SoCounterAction,SoSwitchElement);

  SO_ACTION_ADD_METHOD(SoNode,SoCounterAction::actionMethod);
}
SoCounterAction::SoCounterAction()
:fCount(0),fLookFor(NODE),fType(),fCheckDerived(TRUE) {
  SO_ACTION_CONSTRUCTOR(SoCounterAction);
}
SoCounterAction::~SoCounterAction(){}
void SoCounterAction::beginTraversal(SoNode* node){
  fCount = 0;
  SoAction::beginTraversal(node);
}
void SoCounterAction::actionMethod(SoAction* aThis,SoNode* aNode) {
  //printf("debug : begin : %s %s\n",aNode->getName().getString(),
    //aNode->getTypeId().getName().getString());
  SoCounterAction* This = (SoCounterAction*)aThis;
  if(This->fLookFor==NODE) {
    This->fCount++;
  } else if(This->fLookFor==TYPE) {
    if(This->fCheckDerived==TRUE) {
      if(aNode->getTypeId().isDerivedFrom(This->fType)) This->fCount++;
    } else {
      if(aNode->getTypeId()==This->fType) This->fCount++;
    }
  } else if(This->fLookFor==NAME) {
    if(This->fName==aNode->getName()) This->fCount++;
  }
  if(aNode->isOfType(SoSwitch::getClassTypeId())) {
    SoSwitch* sw = (SoSwitch*)aNode;
    SbBool flag = sw->whichChild.enableNotify(FALSE);
    int old = sw->whichChild.getValue();
    sw->whichChild.setValue(SO_SWITCH_ALL);
    aNode->doAction(aThis);
    sw->whichChild.setValue(old);
    sw->whichChild.enableNotify(flag);
  } else if(aNode->isOfType(SoGroup::getClassTypeId())) {
    aNode->doAction(aThis);
  } else if(aNode->isOfType(SoBaseKit::getClassTypeId())) {
    aNode->doAction(aThis);
  }
}
void SoCounterAction::setLookFor(LookFor aLookFor) {
  fLookFor = aLookFor;
}
void SoCounterAction::setType(const SoType aType,SbBool aCheckDerived) {
  fType = aType;
  fCheckDerived = aCheckDerived;
}
void SoCounterAction::setName(const SbName aName){
  fName = aName;
}
int SoCounterAction::getCount() const {
  return fCount;
}

#endif
