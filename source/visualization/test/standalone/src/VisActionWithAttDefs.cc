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
// $Id: VisActionWithAttDefs.cc,v 1.4 2005-03-28 19:15:33 allison Exp $
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
  (*fpAttDefs)["name"] =
    G4AttDef("name","Name of box","Bookkeeping","","G4String");
  (*fpAttDefs)["dimX"] =
    G4AttDef("dimX","Half x-side","Physics","G4BestUnit","G4double");
  (*fpAttDefs)["dimY"] =
    G4AttDef("dimY","Half y-side","Physics","km","G4double");
  (*fpAttDefs)["dimZ"] =
    G4AttDef("dimZ","Half z-side","Fizziks","RubBish","G4BestUnit");
  (*fpAttDefs)["diag"] =
    G4AttDef("diag","Diagonal","Physics","mmm","G4double");
  (*fpAttDefs)["dims"] =
    G4AttDef("dims","Half sides","Physics","G4BestUnit","G4ThreeVector");

  fpAttValues = new std::vector<G4AttValue>;
  fpAttValues->push_back(G4AttValue("name","Cyan box with attributes",""));
  G4double x = fpBox->GetXHalfLength();
  G4double y = fpBox->GetYHalfLength();
  G4double z = fpBox->GetZHalfLength();
  fpAttValues->push_back
    (G4AttValue("dimX",G4BestUnit(x,"Length"),""));
  fpAttValues->push_back
    (G4AttValue("dimY",G4UIcommand::ConvertToString(y/km),""));
  fpAttValues->push_back
    (G4AttValue("dimZ",G4BestUnit(z,"Length"),""));
  fpAttValues->push_back
    (G4AttValue("diag",G4BestUnit(sqrt(x*x+y*y+z*z),"Length"),""));
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

  std::vector<G4AttValue>* pStandardAttValues = new std::vector<G4AttValue>;
  std::map<G4String,G4AttDef>* pStandardAttDefs =
    new std::map<G4String,G4AttDef>;
  G4AttCheck(fpAttValues,fpAttDefs).Standard
    (pStandardAttValues,pStandardAttDefs);
  G4cout << "\nVisActionWithAttDefs: constructor: standardised versions:\n"
	 << G4AttCheck(pStandardAttValues,pStandardAttDefs)
	 << G4endl;
  delete pStandardAttDefs;
  delete pStandardAttValues;

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
