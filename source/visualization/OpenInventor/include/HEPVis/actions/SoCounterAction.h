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
#ifndef HEPVis_SoCounterAction_h
#define HEPVis_SoCounterAction_h 

// Inheritance :
#include <Inventor/actions/SoAction.h>

#include <Inventor/actions/SoSubAction.h>

#define SoCounterAction Geant4_SoCounterAction

class SoCounterAction : public SoAction {
  SO_ACTION_HEADER(SoCounterAction);
public:
  SoCounterAction();
  virtual ~SoCounterAction();  
  static void initClass(void);
  enum LookFor { NODE = 1,TYPE = 2, NAME = 3 };
  void setLookFor(LookFor);
  void setType(const SoType,SbBool = TRUE);
  void setName(const SbName);
public:
  int getCount() const;
protected:
  virtual void beginTraversal(SoNode*);
private:
  static void actionMethod(SoAction*,SoNode*);
private:
  int fCount;
  int fLookFor;
  SbName fName;
  SoType fType;
  SbBool fCheckDerived;
};

#endif
