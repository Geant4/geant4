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
