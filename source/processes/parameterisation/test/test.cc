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
// $Id: test.cc,v 1.1 2002-11-02 00:10:35 mverderi Exp $
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
