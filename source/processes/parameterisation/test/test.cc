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
// $Id: test.cc,v 1.2 2006-06-29 21:09:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------

#include <iostream.h>

#include "G4GlobalFastSimulationManager.hh"
#include "G4FastSimulationManager.hh"
#include "DummyModel.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

int main()
{
  G4Box box("box",10,10,10);
  G4LogicalVolume vol1(&box,0,"vol1");
  G4LogicalVolume vol2(&box,0,"vol2");
  G4LogicalVolume vol3(&box,0,"vol3");

  DummyModel dummy11 ("dummy11", &vol1);
  DummyModel dummy12 ("dummy12", &vol1);
  DummyModel to_find1("to_find", &vol1);
  DummyModel dummy13 ("dummy13", &vol1);
  
  DummyModel to_find2("to_find", &vol2);
  DummyModel dummy21 ("dummy21", &vol2);
  DummyModel dummy22 ("dummy22", &vol2);
  DummyModel dummy23 ("dummy23", &vol2);
  
  DummyModel dummy31 ("dummy31", &vol3);
  DummyModel to_find3("to_find", &vol3);
  DummyModel to_find4("to_find", &vol3);
  DummyModel to_find5("to_find", &vol3);
  DummyModel dummy32 ("dummy32", &vol3);

  G4GlobalFastSimulationManager* gbl = G4GlobalFastSimulationManager::GetGlobalFastSimulationManager();

  cout << "model addresses to retrieve:" << endl;
  cout << &to_find1 << endl;
  cout << &to_find2 << endl;
  cout << &to_find3 << endl;
  cout << &to_find4 << endl;
  cout << &to_find5 << endl;

  cout << "-- individual : " << endl;
  cout << gbl->GetFastSimulationModel("to_find") << endl;
  cout << gbl->GetFastSimulationModel("to_find") << endl;
  
  cout << "-- loop: " << endl;
  G4VFastSimulationModel *model, *previousModel(0);
  while ((model = gbl->GetFastSimulationModel("to_find", previousModel)))
    {
      cout << "->" << model << endl;
      previousModel = model;
    }
  

}
