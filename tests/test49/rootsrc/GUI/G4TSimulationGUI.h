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
// Class G4TSimulationGUI
//
// Class description:
//
// The GUI for the simulation tool
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------
#ifndef G4TSimulationGUI_H_
#define G4TSimulationGUI_H_

#include "../CommonHeaders.h"
#include "../Database/G4TCatalog.h"
#include "../Helpers/G4TSimHelper.h"
#include "../G4TSimulationTool.h"

#include "Riostream.h"


class G4TSimulationGUI : public TGMainFrame {

  private:
  TGTextButton* bRun;
  TGTextButton* bClose;
  TGComboBox*  cPublication;
  TGComboBox*  cModel;
  TGNumberEntry* nCrossSection;
  TGNumberEntry* nRuns;
  TGCheckButton*  cExisting;

  void Initialize();

  void Run_Click();
  void Close_Click();

  public:


 G4TSimulationGUI(const TGWindow *p, UInt_t w, UInt_t h);
 virtual ~G4TSimulationGUI ();

 virtual void CloseWindow();
 virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);


 ClassDef(G4TSimulationGUI, 1)  //The class for Geant4 Simulation
};



#endif




