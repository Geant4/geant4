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
// $Id: G4tgbMaterialMgr.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgbMaterialMgr

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbMaterialMgr.hh"
#include "G4tgbMaterialMixtureByWeight.hh"
#include "G4tgbMaterialMixtureByVolume.hh"
#include "G4tgbMaterialMixtureByNoAtoms.hh"
#include "G4tgbMaterialSimple.hh"

#include "G4tgrMaterialFactory.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMaterialMixture.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "G4NistManager.hh"


G4ThreadLocal G4tgbMaterialMgr* G4tgbMaterialMgr::theInstance = 0;


// -------------------------------------------------------------------------
G4tgbMaterialMgr::G4tgbMaterialMgr()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialMgr* G4tgbMaterialMgr::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgbMaterialMgr;
    theInstance->CopyIsotopes();
    theInstance->CopyElements();
    theInstance->CopyMaterials();
  }
  return theInstance;
}


// -------------------------------------------------------------------------
G4tgbMaterialMgr::~G4tgbMaterialMgr()
{
  G4mstgbisot::const_iterator isotcite;
  for( isotcite = theG4tgbIsotopes.begin();
       isotcite != theG4tgbIsotopes.end(); isotcite++)
  {
    delete (*isotcite).second;
  }
  theG4tgbIsotopes.clear();

  G4mstgbelem::const_iterator elemcite;
  for( elemcite = theG4tgbElements.begin();
       elemcite != theG4tgbElements.end(); elemcite++)
  {
    delete (*elemcite).second;
  }
  theG4tgbElements.clear();

  G4mstgbmate::const_iterator matcite;
  for( matcite = theG4tgbMaterials.begin();
       matcite != theG4tgbMaterials.end(); matcite++)
  {
    delete (*matcite).second;
  }
  theG4tgbMaterials.clear();

  delete theInstance;
}


// -------------------------------------------------------------------------
void G4tgbMaterialMgr::CopyIsotopes()
{
  const G4mstgrisot tgrIsots
        = G4tgrMaterialFactory::GetInstance()->GetIsotopeList();
  G4mstgrisot::const_iterator cite;
  for( cite = tgrIsots.begin(); cite != tgrIsots.end(); cite++ )
  {
    G4tgrIsotope* tgr = (*cite).second;
    G4tgbIsotope* tgb = new G4tgbIsotope( tgr );
    theG4tgbIsotopes[tgb->GetName()] = tgb;
  }
}


// -------------------------------------------------------------------------
void G4tgbMaterialMgr::CopyElements()
{
  const G4mstgrelem tgrElems
        = G4tgrMaterialFactory::GetInstance()->GetElementList();
  G4mstgrelem::const_iterator cite;
  for( cite = tgrElems.begin(); cite != tgrElems.end(); cite++ )
  {
    G4tgrElement* tgr = (*cite).second;
    G4tgbElement* tgb = new G4tgbElement( tgr );
    theG4tgbElements[tgb->GetName()] = tgb;
  }
}


// -------------------------------------------------------------------------
void G4tgbMaterialMgr::CopyMaterials()
{
  const G4mstgrmate tgrMates
        = G4tgrMaterialFactory::GetInstance()->GetMaterialList();
  G4mstgrmate::const_iterator cite;
  for( cite = tgrMates.begin(); cite != tgrMates.end(); cite++ )
  {
    G4tgrMaterial* tgr = (*cite).second;
    G4tgbMaterial* tgb = 0;
    if( tgr->GetType() == "MaterialSimple" ) {
      tgb = new G4tgbMaterialSimple( tgr );
    } else if( tgr->GetType() == "MaterialMixtureByWeight" ) {
      tgb = new G4tgbMaterialMixtureByWeight( tgr );
    } else if( tgr->GetType() == "MaterialMixtureByNoAtoms" ) {
      tgb = new G4tgbMaterialMixtureByNoAtoms( tgr );
    } else if( tgr->GetType() == "MaterialMixtureByVolume" ) {
      tgb = new G4tgbMaterialMixtureByVolume( tgr );
    } else {
      return;
    }
    theG4tgbMaterials[tgb->GetName()] = tgb;
  }
}


// -------------------------------------------------------------------------
G4Isotope* G4tgbMaterialMgr::FindOrBuildG4Isotope(const G4String & name) 
{
  G4Isotope* g4isot = FindBuiltG4Isotope( name );
  if( g4isot == 0 )
  {
    G4tgbIsotope* tgbisot = FindG4tgbIsotope( name );
    // FindG4tgbIsotope never returns 0, otherwise if not found, crashes
    g4isot = tgbisot->BuildG4Isotope();
    // Register it
    G4String isotname = g4isot->GetName();
    theG4Isotopes[isotname] = g4isot;
  }
  else
  { 
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " G4tgbMaterialMgr::FindOrBuildG4Isotope() -"
             << " G4Isotope already built: " << g4isot->GetName() << G4endl;
    }
#endif
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbMaterialMgr::FindOrBuildG4Isotope() - Isotope: "
           << name << G4endl;
  }
#endif
  return g4isot;
} 


// -------------------------------------------------------------------------
G4Isotope* G4tgbMaterialMgr::FindBuiltG4Isotope(const G4String & name) const 
{
  G4Isotope* g4isot = 0;
 
  G4msg4isot::const_iterator cite = theG4Isotopes.find( name );
  if( cite != theG4Isotopes.end() )
  {
    g4isot = (*cite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMgr::FindBuiltG4Isotope() - Isotope: "
               << name << " = " << g4isot << G4endl;
      }
#endif

  } 

  return g4isot;
}


// -------------------------------------------------------------------------
G4tgbIsotope* G4tgbMaterialMgr::FindG4tgbIsotope(const G4String & name,
                                                       G4bool bMustExist ) const 
{
  G4tgbIsotope* isot = 0;

  G4mstgbisot::const_iterator cite = theG4tgbIsotopes.find( name ); 
  if( cite != theG4tgbIsotopes.end() )
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMgr::FindG4tgbIsotope() -"
             << " G4tgbIsotope found: " << ( (*cite).second )->GetName()
             << G4endl;
    }
#endif
    isot = (*cite).second;
  }
  if( (isot == 0) && bMustExist )
  {
    G4String ErrMessage = "Isotope " + name + " not found !";
    G4Exception("G4tgbMaterialMgr::FindG4tgbIsotope()",
                "InvalidSetup", FatalException, ErrMessage);
  }

  return isot;
}


// -------------------------------------------------------------------------
G4Element* G4tgbMaterialMgr::FindOrBuildG4Element(const G4String & name,
                                                  G4bool bMustExist ) 
{
  G4Element* g4elem = FindBuiltG4Element( name );
  if( g4elem == 0 )
  {
    G4tgbElement* tgbelem = FindG4tgbElement( name, false );
    if( tgbelem == 0)
    {
      // If FindG4tgbElement returns 0, look for a G4NISTElement
      G4cout << "  G4NistManager::Instance()->FindOrBuildElement( " << G4endl;
      g4elem = G4NistManager::Instance()->FindOrBuildElement(name);
    }
    else
    {
      if( tgbelem->GetType() == "ElementSimple" )
      {
        g4elem = tgbelem->BuildG4ElementSimple();
      }
      else if( tgbelem->GetType() == "ElementFromIsotopes" )
      {
        g4elem = tgbelem->BuildG4ElementFromIsotopes();
      }
      else
      {
        G4String ErrMessage = "Element type " + tgbelem->GetType()
          + " does not exist !";
        G4Exception("G4tgbMaterialMgr::GetG4Element()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
    }
    // Register it
    if( (g4elem != 0) )
    {
      theG4Elements[g4elem->GetName()] = g4elem;
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
        {
          G4cout << " G4tgbMaterialMgr::FindOrBuildG4Element() - Element: "
                 << name << G4endl;
        }
#endif
    }
    else
    {
      if( bMustExist )
      {
        G4String ErrMessage = "Element " + name + " not found !";
        G4Exception("G4tgbMaterialMgr::FindOrBuildG4Element()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMgr::FindOrBuildG4Element() - Element: "
               << name << " not found  " << G4endl;
      }
#endif
    }
  }
  else
  { 
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " G4tgbMaterialMgr::GetG4Element() -"
             << " G4Element already built: " << g4elem->GetName() << G4endl; 
    }
#endif
  }

  return g4elem;
} 


// -------------------------------------------------------------------------
G4Element* G4tgbMaterialMgr::FindBuiltG4Element(const G4String & name) const 
{
  G4Element* g4elem = 0;
 
  G4msg4elem::const_iterator cite = theG4Elements.find( name );
  if( cite != theG4Elements.end() )
  {
    g4elem = (*cite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMgr::FindBuiltG4Element() - Element: "
             << name << " = " << g4elem << G4endl;
    }
#endif
  } 

  return g4elem;
}


// -------------------------------------------------------------------------
G4tgbElement* G4tgbMaterialMgr::FindG4tgbElement(const G4String & name,
                                                       G4bool bMustExist ) const
{
  G4tgbElement* elem = 0;

  G4mstgbelem::const_iterator cite = theG4tgbElements.find( name ); 
  if( cite != theG4tgbElements.end() )
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMgr::FindG4tgbElement() -"
             << " G4tgbElement found: " << ( (*cite).second )->GetName()
             << G4endl;
    }
#endif
    elem = (*cite).second;
  }
  if( (elem == 0) && bMustExist )
  {
    G4String ErrMessage = "Element " + name + "  not found !";
    G4Exception("G4tgbMaterialMgr::FindG4tgbElement()",
                "InvalidSetup", FatalException, ErrMessage);
  }

  return elem;
}


// -------------------------------------------------------------------------
G4Material* G4tgbMaterialMgr::FindOrBuildG4Material(const G4String & name,
                                                    G4bool bMustExist )  
{
  G4Material* g4mate = FindBuiltG4Material( name );
  if( g4mate == 0)
  {
    G4tgbMaterial* tgbmate = FindG4tgbMaterial( name, false );
   
    if( tgbmate == 0)
    {
      // if FindG4tgbMaterial() returns 0, look for a G4NISTMaterial
      g4mate = G4NistManager::Instance()->FindOrBuildMaterial(name);
    }
    else
    {
      g4mate = tgbmate->BuildG4Material();

      if( tgbmate->GetTgrMate()->GetIonisationMeanExcitationEnergy() != -1. )
      {
        g4mate->GetIonisation()->SetMeanExcitationEnergy(tgbmate->
                GetTgrMate()->GetIonisationMeanExcitationEnergy());
      }
    }

    // Register it
    if( g4mate != 0 )
    {
      theG4Materials[g4mate->GetName()] = g4mate;
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMgr::FindOrBuildG4Material() - Material: "
               << name << G4endl;
      }
#endif
    }
    else
    {
      if( bMustExist )
      {
        G4String ErrMessage = "Material " + name + "  not found !";
        G4Exception("G4tgbMaterialMgr::FindOrBuildG4Material()",
                    "InvalidSetup", FatalException, ErrMessage);
      }  
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbMaterialMgr::FindOrBuildG4Material() - Element: "
               << name << " not found  " << G4endl;
      }
#endif
    }
  }
  else
  { 
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " G4tgbMaterialMgr::FindOrBuildG4Material() -"
             << " G4Material already built: " << g4mate->GetName() << G4endl;
    }
#endif
  }

  return g4mate;
}


// -------------------------------------------------------------------------
G4Material* G4tgbMaterialMgr::FindBuiltG4Material(const G4String & name) const 
{
  G4Material* g4mate = 0;
  //---------- look for an existing G4Material
  G4msg4mate::const_iterator cite = theG4Materials.find( name );
  if( cite != theG4Materials.end() )
  {
    g4mate = (*cite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMgr::FindBuiltG4Material() - Material: "
             << name << " = " << g4mate << G4endl;
    }
#endif
  }

  return g4mate;
}


// -------------------------------------------------------------------------
G4tgbMaterial* G4tgbMaterialMgr::FindG4tgbMaterial(const G4String & name,
                                                   G4bool bMustExist ) const 
{
  G4tgbMaterial* mate = 0;
  G4mstgbmate::const_iterator cite = theG4tgbMaterials.find( name );
  if( cite != theG4tgbMaterials.end() )
  {
    mate = (*cite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbMaterialMgr::FindG4tgbMaterial() -"
             << " G4tgbMaterial found: " << ( (*cite).second )->GetName()
             << " type " << ( (*cite).second )->GetName() << G4endl;
    }
#endif
  }

  if( (mate == 0) && bMustExist )
  {
    G4String ErrMessage = "Material " + name + "  not found !";
    G4Exception("G4tgbMaterialMgr::FindG4tgbMaterial()",
                "InvalidSetup", FatalException, ErrMessage);
  }

  return mate;
}
