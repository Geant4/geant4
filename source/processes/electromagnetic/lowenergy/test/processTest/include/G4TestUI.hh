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
//
// $Id: G4TestUI.hh,v 1.1 2001-10-28 18:00:34 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// File name:     G4TestUI
//
// Author:        Maria Grazia Pia
// 
// Creation date: 1 October 2001
//
// Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#ifndef G4TESTUI_HH
#define G4TESTUI_HH

class G4ProcessTest;

class G4TestUI
{
  public:

  G4TestUI();

  virtual ~G4TestUI();
  
  void selectMaterial();
  
  void selectProcess();
  
  void selectTestTopic();
  
  void selectNumberOfIterations();
  
  const G4Material* getSelectedMaterial() const;
  
  G4ProcessTest* getSelectedProcess() const;

  const G4String& getTestTopic() const;
  
  G4int getNIterations() const;
  
 private:

  const G4String& GetTestTopic();
  void operator=(const G4TestUI& right);

  void selectProcessType();
  void selectProcessCategory();
  void isPolarised();

  G4int materialId;            
  G4int type;
  G4int category;
  G4bool polarised;
  G4int nIterations;
  G4int topic;
  vector<G4String> types;
  vector<G4String> topics;
  vector<G4String> categories;

};

#endif 
