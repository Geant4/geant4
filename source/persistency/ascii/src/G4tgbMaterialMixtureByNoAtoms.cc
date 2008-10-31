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
// $Id: G4tgbMaterialMixtureByNoAtoms.cc,v 1.2 2008-10-31 18:33:30 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4tgbMaterialMixtureByNoAtoms

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbMaterialMixtureByNoAtoms.hh"
#include "G4tgbMaterial.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"


// -------------------------------------------------------------------------
G4tgbMaterialMixtureByNoAtoms::G4tgbMaterialMixtureByNoAtoms()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialMixtureByNoAtoms::~G4tgbMaterialMixtureByNoAtoms()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialMixtureByNoAtoms::
G4tgbMaterialMixtureByNoAtoms( G4tgrMaterial* hg)
{
  theTgrMate = hg;
}


// -------------------------------------------------------------------------
G4Material* G4tgbMaterialMixtureByNoAtoms::BuildG4Material()
{ 
  //----- construct new G4Material with components materials (a mixture)
  G4Material* mate = new G4Material( theTgrMate->GetName(),
                                     theTgrMate->GetDensity(),
                                     theTgrMate->GetNumberOfComponents(),
                                     kStateUndefined, STP_Temperature );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbMaterialMixtureByNoAtoms::BuildG4Material() -"
           << " Constructing new G4Material: " << theTgrMate->GetName()
           << " " << theTgrMate->GetDensity() << G4endl;
  }
#endif
  //--- add components
  G4Element* compElem;
  G4Material* compMate;
  compAreElements = 0;
  compAreMaterials = 0;
  G4double totalWeight = 0.; //use for number of atoms -> weight conversion
  std::vector<G4Material*> compMateList; 
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  for( G4int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++)
  {
    //look if this component is an element
    G4tgbElement* hselem = mf->FindG4tgbElement( GetComponent(ii) );
    if( hselem != 0 )
    {
      compElem = mf->FindOrBuildG4Element( GetComponent(ii) );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMixtureByNoAtoms::BuildG4Material() -"
               << " Adding component element ..." << G4endl;
      }
#endif
      if( compAreMaterials )
      {
        G4String ErrMessage = "Material with some component elements and "
                            + G4String("some materials: ")
                            + theTgrMate->GetName();
        G4Exception("G4tgbMaterialMixtureByNoAtoms::BuildG4Material()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
      compAreElements = 1;
      //add it by number of atoms
      mate->AddElement( compElem, G4int(GetFraction(ii)) );
      //if it is not an element look if it is a material
    }
    else
    { 
      compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMixtureByNoAtoms::BuildG4Material() -"
               << " compMate: " << GetFraction(ii) << G4endl;
      }
#endif
      if( compMate != 0 )
      { 
        if( compAreElements )
        { 
          G4String ErrMessage = "Material with some component materials and "
                              + G4String("some elements: ")
                              + theTgrMate->GetName();
          G4Exception("G4tgbMaterialMixtureByNoAtoms::BuildG4Material()",
                      "InvalidSetup", FatalException, ErrMessage);
        }
        compAreMaterials = 1;
/*
        G4String ErrMessage = "Adding materials by atoms is not supported: "
                            + theTgrMate->GetName() + ", sorry ...";
        G4Exception("G4tgbMaterialMixtureByNoAtoms::buildG4Material()",
                    "NotImplemented", FatalException, ErrMessage);
*/
        // If it is a material transform No Atoms to Weight
	//        if( fr > 1.0 ) { fr = 1.; }
        G4double fr = GetFraction(ii);
	fr *= compMate->GetDensity();
	totalWeight += fr;
	//Do it after normalization  mate->AddMaterial( compMate, GetFraction( ii ) );
	compMateList.push_back(compMate); 
      }
      else
      {
        G4String ErrMessage = "Component " + GetComponent(ii)
                            + " of material " +  theTgrMate->GetName()
                            + "\n" + "is not an element nor a material !";
        G4Exception("G4tgbMaterialMixtureByWeight::buildG4Material()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
    } 
  }

  //renormalize after the conversion number of atoms -> weight
  if( compAreMaterials ) {
    for( G4int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++) {
      G4double fr = GetFraction(ii);
      fr *= compMateList[ii]->GetDensity()/totalWeight; 
      mate->AddMaterial( compMateList[ii], fr ); 
    }
  }
      
  return mate;

}


// -------------------------------------------------------------------------
void G4tgbMaterialMixtureByNoAtoms::TransformToFractionsByWeight() 
{
}
