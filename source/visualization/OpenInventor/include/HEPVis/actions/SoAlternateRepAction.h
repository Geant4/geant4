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
