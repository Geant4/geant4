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
