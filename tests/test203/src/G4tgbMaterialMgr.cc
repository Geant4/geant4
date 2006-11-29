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

//#include <stdlib.h>

G4tgbMaterialMgr* G4tgbMaterialMgr::theInstance = 0;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialMgr* G4tgbMaterialMgr::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgbMaterialMgr;
    theInstance->CopyIsotopes();
    theInstance->CopyElements();
    theInstance->CopyMaterials();
  }
  return theInstance;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialMgr::~G4tgbMaterialMgr()
{
  mstgbisot::const_iterator isotcite;
  for( isotcite = theG4tgbIsotopes.begin(); isotcite != theG4tgbIsotopes.end(); isotcite++) {
    delete (*isotcite).second;
  }
  theG4tgbIsotopes.clear();

  mstgbelem::const_iterator elemcite;
  for( elemcite = theG4tgbElements.begin(); elemcite != theG4tgbElements.end(); elemcite++) {
    delete (*elemcite).second;
  }
  theG4tgbElements.clear();

  mstgbmate::const_iterator matcite;
  for( matcite = theG4tgbMaterials.begin(); matcite != theG4tgbMaterials.end(); matcite++) {
    delete (*matcite).second;
  }
  theG4tgbMaterials.clear();

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbMaterialMgr::CopyIsotopes()
{
  const mstgrisot tgrIsots = (G4tgrMaterialFactory::GetInstance())->GetIsotopeList();
  mstgrisot::const_iterator cite;
  for( cite = tgrIsots.begin(); cite != tgrIsots.end(); cite++ ){
    G4tgrIsotope* tgr = (*cite).second;
    G4tgbIsotope* tgb = new G4tgbIsotope( tgr );
    theG4tgbIsotopes[tgb->GetName()] = tgb;
  }

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbMaterialMgr::CopyElements()
{
  const mstgrelem tgrElems = (G4tgrMaterialFactory::GetInstance())->GetElementList();
  mstgrelem::const_iterator cite;
  for( cite = tgrElems.begin(); cite != tgrElems.end(); cite++ ){
    G4tgrElement* tgr = (*cite).second;
    G4tgbElement* tgb = new G4tgbElement( tgr );
    theG4tgbElements[tgb->GetName()] = tgb;
  }

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbMaterialMgr::CopyMaterials()
{
  const mstgrmate tgrMates = (G4tgrMaterialFactory::GetInstance())->GetMaterialList();
  mstgrmate::const_iterator cite;
  for( cite = tgrMates.begin(); cite != tgrMates.end(); cite++ ){
    G4tgrMaterial* tgr = (*cite).second;
    G4tgbMaterial* tgb;
    if( tgr->GetType() == "MaterialSimple" ) {
      tgb = new G4tgbMaterialSimple( tgr );
    }else if( tgr->GetType() == "MaterialMixtureByWeight" ) {
      tgb = new G4tgbMaterialMixtureByWeight( tgr );
    }else if( tgr->GetType() == "MaterialMixtureByNoAtoms" ) {
      tgb = new G4tgbMaterialMixtureByNoAtoms( tgr );
    }else if( tgr->GetType() == "MaterialMixtureByVolume" ) {
      tgb = new G4tgbMaterialMixtureByVolume( tgr );
    } 
    theG4tgbMaterials[tgb->GetName()] = tgb;
  }

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Isotope* G4tgbMaterialMgr::FindOrBuildG4Isotope(const G4String & name) 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindOrBuildG4Isotope " << name << endl;
#endif
  G4Isotope* g4isot = FindG4Isotope( name );
  if( g4isot == 0 ) {
    G4tgbIsotope* tgbisot = FindG4tgbIsotope( name );
    //FindG4tgbIsotope never returns 0, because if it is not found, it will crash
    g4isot = tgbisot->BuildG4Isotope();
    //!register it
    G4String isotname = g4isot->GetName();
    theG4Isotopes[isotname] = g4isot;
  } else { 
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
      cout << " G4Isotope already built: " << g4isot->GetName() << endl; 
#endif
  }
  return g4isot;
} 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Isotope* G4tgbMaterialMgr::FindG4Isotope(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4Isotope " << name << endl;
#endif
  G4Isotope* g4isot = 0;
 
  msg4isot::const_iterator cite = theG4Isotopes.find( name );
  if( cite != theG4Isotopes.end() ) {
    g4isot = (*cite).second;
  } 

  return g4isot;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbIsotope* G4tgbMaterialMgr::FindG4tgbIsotope(const G4String & name, G4bool bMustExist ) const 
{
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4tgbIsotope " << name << endl;
#endif
 
  G4tgbIsotope* isot = 0;

  mstgbisot::const_iterator cite = theG4tgbIsotopes.find( name ); 
  if( cite != theG4tgbIsotopes.end() ) {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbIsotope found " << ( (*cite).second )->GetName() << endl;
#endif
    isot = (*cite).second;
  }
  if( isot == 0 && bMustExist ){
    G4Exception("G4tgbMaterialMgr::FindG4tgbIsotope. Isotope " + name + "  not found ");
  }

  return isot;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Element* G4tgbMaterialMgr::FindOrBuildG4Element(const G4String & name) 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindOrBuildG4Element " << name << endl;
#endif
  G4Element* g4elem = FindG4Element( name );
  if( g4elem == 0 ) {
    G4tgbElement* tgbelem = FindG4tgbElement( name );
    //FindG4tgbElement never returns 0, because if it is not found, it will crash
    if( tgbelem->GetType() == "ElementSimple" ){
      g4elem = tgbelem->BuildG4ElementSimple();
    } else if( tgbelem->GetType() == "ElementFromIsotopes" ){
      g4elem = tgbelem->BuildG4ElementFromIsotopes();
    } else {
      G4Exception("G4tgbMaterialMgr::GetG4Element element type does not exist " +  tgbelem->GetType());
    }

    //!register it
    G4String elemname = g4elem->GetName();
    theG4Elements[elemname] = g4elem;
  } else { 
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
    cout << " G4Element already built: " << g4elem->GetName() << endl; 
#endif
  }
  return g4elem;
} 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Element* G4tgbMaterialMgr::FindG4Element(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4Element " << name << endl;
#endif
  G4Element* g4elem = 0;
 
  msg4elem::const_iterator cite = theG4Elements.find( name );
  if( cite != theG4Elements.end() ) {
    g4elem = (*cite).second;
  } 

  return g4elem;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbElement* G4tgbMaterialMgr::FindG4tgbElement(const G4String & name, G4bool bMustExist ) const 
{
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4tgbElement " << name << endl;
#endif
 
  G4tgbElement* elem = 0;

  mstgbelem::const_iterator cite = theG4tgbElements.find( name ); 
  if( cite != theG4tgbElements.end() ) {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbElement found " << ( (*cite).second )->GetName() << endl;
#endif
    elem = (*cite).second;
  }
  if( elem == 0 && bMustExist ){
    G4Exception("!!! EXITING  G4tgbMaterialMgr::FindG4tgbElement. Element " + name + "  not found ");
  }

  return elem;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialMgr::FindOrBuildG4Material(const G4String & name) 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindOrBuildG4Material " << name << endl;
#endif
  G4Material* g4mate = FindG4Material( name );
  if( g4mate == 0) {
    G4tgbMaterial* tgbmate = FindG4tgbMaterial( name, 1 );
    //FindG4tgbMaterial never returns 0, because if it is not found, it will crash
    g4mate = tgbmate->BuildG4Material();
    //! register it
    theG4Materials[g4mate->GetName()] = g4mate;
  } else { 
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
      cout << " G4Material already built: " << g4mate->GetName() << endl;
#endif
  }
  return g4mate;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialMgr::FindG4Material(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4Material " << name << endl;
#endif
  G4Material* g4mate = 0;
  //---------- look for an existing G4Material
  msg4mate::const_iterator cite = theG4Materials.find( name );
  if( cite != theG4Materials.end() ) {
    g4mate = (*cite).second;
  }

  return g4mate;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterial* G4tgbMaterialMgr::FindG4tgbMaterial(const G4String & name, G4bool bMustExist ) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbMaterialMgr::FindG4tgbMaterial " << name << endl;
#endif
 
  G4tgbMaterial* mate = 0;
  mstgbmate::const_iterator cite = theG4tgbMaterials.find( name );
  if( cite != theG4tgbMaterials.end() ) {
    mate = (*cite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
      G4cout << " G4tgbMaterial found: " << ( (*cite).second )->GetName() << " type " << ( (*cite).second )->GetName() << G4endl;
#endif
  }

  if( mate == 0 && bMustExist ){
    G4Exception("G4tgbMaterialMgr::FindG4tgbMaterial. Material " + name + "  not found ");
  }
  
  /*-  for( cite =  theG4tgbMaterials.begin(); cite != 
	 theG4tgbMaterials.end(); cite++ ) {
    cout << "tgb material" << (*cite).first << endl;
  }
  */

  return mate;
}
