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
// $Id: VisActionWithAttDefs.cc,v 1.1 2005-03-23 17:43:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "VisActionWithAttDefs.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4AttCheck.hh"
#include "G4UnitsTable.hh"

VisActionWithAttDefs::VisActionWithAttDefs ()
{
  fpBox = new G4Box("boxAtts",1*m,2*m,3*m);

  fpAttDefs = new std::map<G4String,G4AttDef>;
  (*fpAttDefs)["name"] = G4AttDef("name","Name of box","","","G4String");
  (*fpAttDefs)["dimX"] = G4AttDef("dimX","Half x-side","","","G4double");
  (*fpAttDefs)["dimY"] = G4AttDef("dimY","Half y-side","","","G4double");
  (*fpAttDefs)["dimZ"] = G4AttDef("dimZ","Half z-side","","","G4double");

  fpAttValues = new std::vector<G4AttValue>;
  fpAttValues->push_back(G4AttValue("name","Cyan box with attributes",""));
  fpAttValues->push_back
    (G4AttValue("dimX",G4BestUnit(fpBox->GetXHalfLength(),"Length"),""));
  fpAttValues->push_back
    (G4AttValue("dimY",G4BestUnit(fpBox->GetYHalfLength(),"Length"),""));
  fpAttValues->push_back
    (G4AttValue("dimZ",G4BestUnit(fpBox->GetZHalfLength(),"Length"),""));

  G4cout << "\nVisActionWithAttDefs: constructor: att values check:\n"
	 << G4AttCheck(fpAttValues,fpAttDefs) << G4endl;

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
