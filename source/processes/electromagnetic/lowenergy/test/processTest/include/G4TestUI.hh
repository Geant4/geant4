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
// $Id: G4TestUI.hh,v 1.3 2001-11-01 17:26:18 pia Exp $
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

#ifndef G4TESTUI_HH
#define G4TESTUI_HH

#include "globals.hh"
#include "g4std/vector"

class G4ProcessTest;
class G4Material;
class G4ParticleDefinition;

class G4TestUI
{
  public:

  G4TestUI();
  virtual ~G4TestUI();
  
  void configure();

  void selectNumberOfIterations();

  void selectMaterial(); 

  void selectProcess();

  void selectTestTopic();  

  void selectEnergyRange();

  G4int getNumberOfIterations() const;

  const G4Material* getSelectedMaterial() const;

  const G4String& getProcessType() const;

  const G4String& getProcessCategory() const;

  const G4String& getTestTopic() const ;

  G4bool getPolarisationSelection() const;

  G4ParticleDefinition* getParticleDefinition() const;

  G4double getMinEnergy() const { return eMin; }

  G4double getMaxEnergy() const { return eMax; }
  
 private:

  void operator=(const G4TestUI& right);

  void selectProcessType();
  void selectProcessCategory();
  void isPolarised();

  G4int nIterations;
  G4int materialId;            
  G4int type;
  G4int category;
  G4int topic;
  G4bool polarised;
  G4double eMin;
  G4double eMax;
  G4std::vector<G4String> types;
  G4std::vector<G4String> topics;
  G4std::vector<G4String> categories;

};

#endif 
