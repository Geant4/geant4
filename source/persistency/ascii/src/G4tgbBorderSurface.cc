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
// G4tgbBorderSurface
//
#include "G4tgbBorderSurface.hh"
#include "G4tgrMaterialPropertiesTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4tgbMaterialPropertiesTable.hh"

// ---------------------------------------------------------------------
G4tgbBorderSurface::G4tgbBorderSurface(G4tgrBorderSurface* tgrbrdr)
{
  theTgrBorderSurface = tgrbrdr;
  theName = tgrbrdr->GetName();

  fOptic = theTgrBorderSurface->GetOptic();

  G4tgrMaterialPropertiesTable *tgrmpt = 
      theTgrBorderSurface->GetTgrMaterialPropertiesTable();
  if (tgrmpt)
  {
    G4tgbMaterialPropertiesTable* tgbmpt = new 
                      G4tgbMaterialPropertiesTable(tgrmpt);
    fOptic->SetMaterialPropertiesTable(
                      tgbmpt->BuildG4MaterialPropertiesTable());
  }
}

// ---------------------------------------------------------------------
G4tgbBorderSurface::~G4tgbBorderSurface()
{
}

// ---------------------------------------------------------------------
void G4tgbBorderSurface::BuildG4BorderSurface()
{
  if (theTgrBorderSurface->GetV1Name() == "skin")
  {
  	BuildG4SkinBorderSurface();
  } else
  {
  	BuildG4LogicalBorderSurface();
  }
}

// ---------------------------------------------------------------------
void G4tgbBorderSurface::BuildG4SkinBorderSurface()
{

  G4LogicalVolume *v1_log = fTgbVolmgr->FindG4LogVol(
                              theTgrBorderSurface->GetV2Name());

  if (v1_log) { // add some error getting tag

    new G4LogicalSkinSurface(theTgrBorderSurface->GetName(), 
      v1_log, fOptic);
    G4cout<<"Border surface "<<theTgrBorderSurface->GetName()<<" around "
          << theTgrBorderSurface->GetV2Name() << ":" << " added." << G4endl;
  }
}

// ---------------------------------------------------------------------
void G4tgbBorderSurface::BuildG4LogicalBorderSurface()
{

  G4LogicalVolume *m1 = fTgbVolmgr->FindG4PhysVol(
                    theTgrBorderSurface->GetV1Name())->GetMotherLogical();
  G4LogicalVolume *m2 = fTgbVolmgr->FindG4PhysVol(
                    theTgrBorderSurface->GetV2Name())->GetMotherLogical();

  // search for physics volumes on the sides of the border
  G4VPhysicalVolume *v1=0, *v2=0;
  for (int i=0; i<(int)m1->GetNoDaughters(); i++) {
    v1 = m1->GetDaughter(i);
    if (v1->GetCopyNo()==theTgrBorderSurface->GetCopyNo1()) break;
  }
  for (int i=0; i<(int)m2->GetNoDaughters(); i++) {
    v2 = m2->GetDaughter(i);
    if (v2->GetCopyNo()==theTgrBorderSurface->GetCopyNo2()) break;
  }
  if (v1 && v2) {
    new G4LogicalBorderSurface(theTgrBorderSurface->GetName(),
      v1, v2, fOptic);
    G4cout<<"Border surface "<<theTgrBorderSurface->GetName()<<" in between "
        << theTgrBorderSurface->GetV1Name()  << ":"       << 
           theTgrBorderSurface->GetCopyNo1() << " and "   << 
           theTgrBorderSurface->GetV2Name()  << ":"       << 
           theTgrBorderSurface->GetCopyNo2() << " added." << G4endl;
  }
}
