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
// $Id: VisActionWithAttDefs.cc,v 1.3 2005-03-28 10:33:27 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "VisActionWithAttDefs.hh"

#include "globals.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"

VisActionWithAttDefs::VisActionWithAttDefs ()
{
  fpBox = new G4Box("boxAtts",1*m,2*m,3*m);

  fpAttDefs = new std::map<G4String,G4AttDef>;
  (*fpAttDefs)["name"] = G4AttDef("name","Name of box","","","G4String");
  (*fpAttDefs)["dimX"] = G4AttDef("dimX","Half x-side","","","G4BestUnit");
  (*fpAttDefs)["dimY"] = G4AttDef("dimY","Half y-side","","","km");
  (*fpAttDefs)["dimZ"] = G4AttDef("dimZ","Half z-side","","RubBish","G4BestUni");
  (*fpAttDefs)["dims"] =
    G4AttDef("dims","Half sides","","G4ThreeVector","G4BestUnit");

  fpAttValues = new std::vector<G4AttValue>;
  fpAttValues->push_back(G4AttValue("name","Cyan box with attributes",""));
  fpAttValues->push_back
    (G4AttValue("dimX",G4BestUnit(fpBox->GetXHalfLength(),"Length"),""));
  fpAttValues->push_back
    (G4AttValue("dimY",
		G4UIcommand::ConvertToString(fpBox->GetYHalfLength()/km),
		""));
  fpAttValues->push_back
    (G4AttValue("dimZ",G4BestUnit(fpBox->GetZHalfLength(),"Length"),""));
  fpAttValues->push_back
    (G4AttValue
     ("dims",
      G4BestUnit
      (G4ThreeVector
       (fpBox->GetXHalfLength(),
	fpBox->GetYHalfLength(),
	fpBox->GetZHalfLength()),
       "Length"),""));
  fpAttValues->push_back
    (G4AttValue("ScoobyDoo","Rubbish",""));

  G4AttCheck(fpAttValues,fpAttDefs).Check();  // Check only.

  // Print...
  G4cout << "\nVisActionWithAttDefs: constructor: att values:\n"
	 << G4AttCheck(fpAttValues,fpAttDefs) << G4endl;

  G4AttCheck standard = G4AttCheck(fpAttValues,fpAttDefs).Standard();
  // Creates new AttValues and AttDefs on the heap...
  G4cout << "\nVisActionWithAttDefs: constructor: standardised versions:\n"
	 << G4AttCheck(standard.GetAttValues(),
		       standard.GetAttDefs())
	 << G4endl;
  // ...so don't forget to delete them...
  delete standard.GetAttValues();
  delete standard.GetAttDefs();

  fpVisAtts = new G4VisAttributes(G4Colour(1,0,1));
  fpVisAtts->SetAttDefs(fpAttDefs);
  fpVisAtts->SetAttValues(fpAttValues);

  fpTransform = new G4Translate3D(-6*m,6*m,0);
}

VisActionWithAttDefs::~VisActionWithAttDefs ()
{
  delete fpVisAtts;
  delete fpAttValues;
  delete fpAttDefs;
  delete fpBox;
}

void VisActionWithAttDefs::Draw()
{
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    pVisManager->Draw(*fpBox,*fpVisAtts,*fpTransform);
  }
}
