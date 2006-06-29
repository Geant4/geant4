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
//
// $Id: G4TestUI.hh,v 1.6 2006-06-29 19:48:36 gunter Exp $
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
#include <vector>

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
  std::vector<G4String> types;
  std::vector<G4String> topics;
  std::vector<G4String> categories;

};

#endif 
