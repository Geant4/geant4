// this :
#include <HEPVis/actions/SoCounterAction.h>

#include <Inventor/nodes/SoNode.h>

SO_ACTION_SOURCE(SoCounterAction);

void SoCounterAction::initClass(void){
  SO_ACTION_INIT_CLASS(SoCounterAction,SoAction);
  SO_ACTION_ADD_METHOD(SoNode,SoCounterAction::actionMethod);
}
SoCounterAction::SoCounterAction()
:fCount(0),fLookFor(NODE),fCheckDerived(TRUE){
  SO_ACTION_CONSTRUCTOR(SoCounterAction);
}
SoCounterAction::~SoCounterAction(){}
void SoCounterAction::beginTraversal(SoNode* node){
  fCount = 0;
  SoAction::beginTraversal(node);
}
void SoCounterAction::actionMethod(SoAction* aThis,SoNode* aNode) {
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
