#ifndef HEPVis_SoAlternateRepAction_h
#define HEPVis_SoAlternateRepAction_h 

// Inheritance :
#include <Inventor/actions/SoAction.h>

#include <Inventor/actions/SoSubAction.h>

#define SoAlternateRepAction Geant4_SoAlternateRepAction

class SoAlternateRepAction : public SoAction {
  SO_ACTION_HEADER(SoAlternateRepAction);
public:
  SoAlternateRepAction();
  virtual ~SoAlternateRepAction();  
  static void initClass(void);
  void setGenerate(SbBool);
  SbBool getGenerate() const;
private:
  static void nodeAction(SoAction*,SoNode*);
private:
  SbBool fGenerate;
};

#define SO_ALTERNATEREP_DO_ACTION(aAction) \
  if(aAction->isOfType(SoAlternateRepAction::getClassTypeId())) {\
    if(((SoAlternateRepAction*)aAction)->getGenerate()==TRUE) {\
      if(alternateRep.getValue()==NULL) {\
        generateAlternateRep();\
        SoNode* altRep = alternateRep.getValue();\
        if((altRep!=NULL) && altRep->isOfType(SoGroup::getClassTypeId()))\
          altRep->doAction(aAction);\
      }\
    } else {\
      SoNode* altRep = alternateRep.getValue();\
      if((altRep!=NULL) && altRep->isOfType(SoGroup::getClassTypeId()))\
        altRep->doAction(aAction);\
      clearAlternateRep();\
    }\
    return;\
  }


#endif
