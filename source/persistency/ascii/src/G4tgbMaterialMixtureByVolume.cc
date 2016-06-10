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
// $Id: G4tgbMaterialMixtureByVolume.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbMaterialMixtureByVolume

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbMaterialMixtureByVolume.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgbMaterial.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgbMaterialMixtureByVolume::G4tgbMaterialMixtureByVolume()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialMixtureByVolume::~G4tgbMaterialMixtureByVolume()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialMixtureByVolume::G4tgbMaterialMixtureByVolume( G4tgrMaterial* hg )
{
  theTgrMate = hg;
}


// -------------------------------------------------------------------------
G4Material* G4tgbMaterialMixtureByVolume::BuildG4Material()
{
 
  //----- construct new G4Material with components materials (a mixture)
  G4Material* mate = new G4Material( theTgrMate->GetName(),
                                     theTgrMate->GetDensity(),
                                     theTgrMate->GetNumberOfComponents(),
                                     theTgrMate->GetState(),
                                     theTgrMate->GetTemperature(),
                                     theTgrMate->GetPressure() );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbMaterialMixtureByVolume::buildG4Material() -"
           << " Constructing new G4Material:"
           << " " << theTgrMate->GetName() 
           << " " << theTgrMate->GetDensity()/g*cm3
           << " " << theTgrMate->GetNumberOfComponents()
           << " " << theTgrMate->GetState()
           << " " << theTgrMate->GetTemperature()
           << " " << theTgrMate->GetPressure() << G4endl;
  }
#endif
 
  //----- Transform fractions by volume to fractions by weight
  TransformToFractionsByWeight();

  //----- Add components
  G4Material* compMate = 0;
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  for( G4int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++)
  {
    // Look if this component is a material
    compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
    if( compMate != 0 )
    { 
      // If it is a material add it by weight fraction 
      mate->AddMaterial( compMate, theFractionsByWeight[ii] );
    }
    else
    {
      G4String ErrMessage = "Component " + GetComponent(ii)
                          + " of material " + theTgrMate->GetName()
                          + "\n" + "is not an element nor a material !";
      G4Exception("G4tgbMaterialMixtureByVolume::BuildG4Material()",
                  "InvalidSetup", FatalException, ErrMessage);
    }
  } 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Constructing new G4Material by volume: " << *mate << G4endl; 
  }
#endif      

  return mate;
}


// -------------------------------------------------------------------------
void G4tgbMaterialMixtureByVolume::TransformToFractionsByWeight() 
{ 
  G4tgbMaterialMgr* mf = G4tgbMaterialMgr::GetInstance();
  G4Material* compMate = 0;
  G4double totalfd = 0.;
  for( G4int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++ )
  {
    compMate = mf->FindOrBuildG4Material( GetComponent(ii) );
    if( compMate != 0 )
    {
      // If it is a material add it by weight fraction 
      theFractionsByWeight.push_back( GetFraction(ii)*compMate->GetDensity() );
      totalfd += theFractionsByWeight[ii];
    }
    else
    {
      G4String ErrMessage = "Component " + GetComponent(ii)
                          + " of material " + theTgrMate->GetName()
                          + "\n" + "is not a material !";
      G4Exception("G4tgbMaterialMixtureByVolume::BuildG4Material()",
                  "InvalidSetup", FatalException, ErrMessage);
    }   
  }
  for( G4int ii = 0; ii < theTgrMate->GetNumberOfComponents(); ii++ )
  {
    theFractionsByWeight[ii] /= totalfd;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMixtureByVolume::TransformToFractionsByWeight()"
             << " Component " << ii  << " : "
             << mf->FindOrBuildG4Material( GetComponent(ii) )->GetName()
             << " FractionByVolume= " << GetFraction(ii)
             << " FractionByWeight= " << theFractionsByWeight[ii] 
             << G4endl;
    }
#endif
  }
}
