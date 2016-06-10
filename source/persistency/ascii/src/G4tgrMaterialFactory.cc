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
// $Id: G4tgrMaterialFactory.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrMaterialFactory

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrMaterialFactory.hh"
#include "G4tgrUtils.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMaterialMixture.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"


G4ThreadLocal G4tgrMaterialFactory* G4tgrMaterialFactory::theInstance = 0;


//-------------------------------------------------------------
G4tgrMaterialFactory::G4tgrMaterialFactory()
{
}


//-------------------------------------------------------------
G4tgrMaterialFactory* G4tgrMaterialFactory::GetInstance()
{
  if( !theInstance )
  {
    theInstance = new G4tgrMaterialFactory;
  }
  return theInstance;
}


//-------------------------------------------------------------
G4tgrMaterialFactory::~G4tgrMaterialFactory()
{
  G4mstgrisot::iterator isotcite;
  for( isotcite = theG4tgrIsotopes.begin();
       isotcite != theG4tgrIsotopes.end(); isotcite++)
  {
    delete (*isotcite).second;
  }
  theG4tgrIsotopes.clear();

  G4mstgrelem::iterator elemcite;
  for( elemcite = theG4tgrElements.begin();
       elemcite != theG4tgrElements.end(); elemcite++)
  {
    delete (*elemcite).second;
  }
  theG4tgrElements.clear();

  G4mstgrmate::iterator matcite;
  for( matcite = theG4tgrMaterials.begin();
       matcite != theG4tgrMaterials.end(); matcite++)
  {
    delete (*matcite).second;
  }
  theG4tgrMaterials.clear();
  delete theInstance;
}


//-------------------------------------------------------------
G4tgrIsotope*
G4tgrMaterialFactory::AddIsotope( const std::vector<G4String>& wl )
{
  //---------- Look if isotope exists
  if( FindIsotope( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    ErrorAlreadyExists("isotope", wl );
  }
  
  G4tgrIsotope* isot = new G4tgrIsotope( wl );
  theG4tgrIsotopes[isot->GetName()] = isot;

  return isot;
}

//-------------------------------------------------------------
G4tgrElementSimple*
G4tgrMaterialFactory::AddElementSimple( const std::vector<G4String>& wl )
{
  //---------- Look if element exists
  if( FindElement( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    ErrorAlreadyExists("element", wl );
  }
  
  G4tgrElementSimple* elem = new G4tgrElementSimple( wl );
  theG4tgrElements[elem->GetName()] = elem;

  return elem;
}


//-------------------------------------------------------------
G4tgrElementFromIsotopes*
G4tgrMaterialFactory::AddElementFromIsotopes( const std::vector<G4String>& wl )
{
  //---------- Look if element exists
  if( FindElement( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    ErrorAlreadyExists("element", wl );
  }
  
  G4tgrElementFromIsotopes* elem = new G4tgrElementFromIsotopes( wl );
  theG4tgrElements[elem->GetName()] = elem;

  return elem;
}


//-------------------------------------------------------------
G4tgrMaterialSimple*
G4tgrMaterialFactory::AddMaterialSimple( const std::vector<G4String>& wl )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrMaterialFactory::AddMaterialSimple" << wl[1] << G4endl;
  }
#endif

  //---------- Look if material exists
  if( FindMaterial( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    ErrorAlreadyExists("material simple", wl );
  }

  G4tgrMaterialSimple* mate = new G4tgrMaterialSimple("MaterialSimple", wl );

  //---------- register this material
  theG4tgrMaterials[ mate->GetName() ] = mate;
  
  return mate;
}


//-------------------------------------------------------------
G4tgrMaterialMixture*
G4tgrMaterialFactory::AddMaterialMixture( const std::vector<G4String>& wl,
                                          const G4String& mixtType )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrMaterialFactory::AddMaterialMixture " << wl[1] << G4endl;
  }
#endif

  //---------- Look if material already exists
  if( FindMaterial( G4tgrUtils::GetString(wl[1]) ) != 0 )
  {
    ErrorAlreadyExists("material mixture", wl );
  }

  G4tgrMaterialMixture* mate; 
  mate = new G4tgrMaterialMixture( mixtType, wl );
  
  //---------- register this material
  theG4tgrMaterials[ mate->GetName() ] = mate;
  
  return mate;
}


//-------------------------------------------------------------
G4tgrIsotope* G4tgrMaterialFactory::FindIsotope(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
     G4cout << " G4tgrMaterialFactory::FindIsotope() - " << name << G4endl;
  }
#endif

  G4mstgrisot::const_iterator cite;
  cite = theG4tgrIsotopes.find( name ); 
  if( cite == theG4tgrIsotopes.end() )
  {
    return 0;
  }
  else
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 3 )
    {
      G4cout << " G4tgrIsotope found: "
             << ( (*cite).second )->GetName() << G4endl;
    }
#endif
    return (*cite).second;
  }
}


//-------------------------------------------------------------
G4tgrElement* G4tgrMaterialFactory::FindElement(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgrMaterialFactory::FindElement() - " << name << G4endl;
  }
#endif
  G4mstgrelem::const_iterator cite;
  cite = theG4tgrElements.find( name ); 
  if( cite == theG4tgrElements.end() )
  {
    return 0;
  }
  else
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 3 )
    {
      DumpElementList();
      G4cout << " G4tgrElement found: "
             << ( (*cite).second )->GetName() << G4endl;
    }
#endif
    return (*cite).second;
  }
}


//-------------------------------------------------------------
G4tgrMaterial* G4tgrMaterialFactory::FindMaterial(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 3 )
  {
    G4cout << " G4tgrMaterialFactory::FindMaterial() - " << name << G4endl;
  }
#endif
  G4mstgrmate::const_iterator cite;
  cite = theG4tgrMaterials.find( name );
  if( cite == theG4tgrMaterials.end() )
  {
    return 0;
  }
  else
  {
    return (*cite).second;
  }
}


//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpIsotopeList() const
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrIsotope's List " << G4endl;
  G4mstgrisot::const_iterator cite;
  for(cite = theG4tgrIsotopes.begin(); cite != theG4tgrIsotopes.end(); cite++)
  {
    G4cout << " ISOT: " << (*cite).second->GetName() << G4endl;
  }
}


//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpElementList() const 
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrElement's List " << G4endl;
  G4mstgrelem::const_iterator cite;
  for(cite = theG4tgrElements.begin(); cite != theG4tgrElements.end(); cite++)
  {
    G4cout << " ELEM: " << (*cite).second->GetName() << G4endl;
  }
}


//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpMaterialList() const
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrMaterial's List " << G4endl;
  G4mstgrmate::const_iterator cite;
  for(cite = theG4tgrMaterials.begin(); cite != theG4tgrMaterials.end(); cite++)
  {
    G4tgrMaterial* mate = (*cite).second;
    G4cout << " MATE: " << mate->GetName() << " Type: " << mate->GetType() 
           << " NoComponents= " << mate->GetNumberOfComponents() << G4endl;
  }
}
 

//-------------------------------------------------------------
void G4tgrMaterialFactory::
ErrorAlreadyExists(const G4String& object,
                   const std::vector<G4String>& wl, const G4bool bNoRepeating )
{
  G4String msg = object + G4String(" repeated");
  if( bNoRepeating )
  {
    G4tgrUtils::DumpVS( wl, (G4String("!!!! EXITING: ") + msg).c_str() );
    G4Exception("G4tgrMaterialFactory", "FatalError",
                FatalException, "Aborting...");
  }
  else
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4tgrUtils::DumpVS( wl, (G4String("!! WARNING: ") + msg).c_str() ); 
    }
#endif
  }
}
