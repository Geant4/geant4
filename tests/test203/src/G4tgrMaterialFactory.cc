#include "G4tgrMaterialFactory.hh"
#include "G4tgrUtils.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMaterialMixture.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"

//#include <stdlib.h>

G4tgrMaterialFactory* G4tgrMaterialFactory::theInstance = 0;


//-------------------------------------------------------------
G4tgrMaterialFactory* G4tgrMaterialFactory::GetInstance()
{
 if( !theInstance ) {
    theInstance = new G4tgrMaterialFactory;
  }
  return theInstance;
}


//-------------------------------------------------------------
G4tgrMaterialFactory::~G4tgrMaterialFactory()
{
  mstgrisot::iterator isotcite;
  for( isotcite = theG4tgrIsotopes.begin(); isotcite != theG4tgrIsotopes.end(); isotcite++) {
    delete (*isotcite).second;
  }
  theG4tgrIsotopes.clear();

  mstgrelem::iterator elemcite;
  for( elemcite = theG4tgrElements.begin(); elemcite != theG4tgrElements.end(); elemcite++) {
    delete (*elemcite).second;
  }
  theG4tgrElements.clear();

  mstgrmate::iterator matcite;
  for( matcite = theG4tgrMaterials.begin(); matcite != theG4tgrMaterials.end(); matcite++) {
    delete (*matcite).second;
  }
  theG4tgrMaterials.clear();

}


//-------------------------------------------------------------
G4tgrIsotope* G4tgrMaterialFactory::AddIsotope( const vector<G4String>& wl )
{
  //---------- Look if isotope exists
  if( FindIsotope( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
    ErrorAlreadyExists("isotope", wl );
  }
  
  G4tgrIsotope* isot = new G4tgrIsotope( wl );
  theG4tgrIsotopes[isot->GetName()] = isot;

  return isot;
}

//-------------------------------------------------------------
G4tgrElementSimple* G4tgrMaterialFactory::AddElementSimple( const vector<G4String>& wl )
{
  //---------- Look if element exists
  if( FindElement( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
    ErrorAlreadyExists("element", wl );
  }
  
  G4tgrElementSimple* elem = new G4tgrElementSimple( wl );
  theG4tgrElements[elem->GetName()] = elem;

  return elem;
}


//-------------------------------------------------------------
G4tgrElementFromIsotopes* G4tgrMaterialFactory::AddElementFromIsotopes( const vector<G4String>& wl )
{
  //---------- Look if element exists
  if( FindElement( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
    ErrorAlreadyExists("element", wl );
  }
  
  G4tgrElementFromIsotopes* elem = new G4tgrElementFromIsotopes( wl );
  theG4tgrElements[elem->GetName()] = elem;

  return elem;
}


//-------------------------------------------------------------
G4tgrMaterialSimple* G4tgrMaterialFactory::AddMaterialSimple( const vector<G4String>& wl )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << " add material simple " << wl[1] << G4endl;
#endif

  //---------- Look if material exists
  if( FindMaterial( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
    ErrorAlreadyExists("material simple", wl );
  }

  G4tgrMaterialSimple* mate = new G4tgrMaterialSimple("MaterialSimple", wl );
  //---------- register this material
  theG4tgrMaterials[ mate->GetName() ] = mate;
  
  return mate;
}


//-------------------------------------------------------------
G4tgrMaterialMixture* G4tgrMaterialFactory::AddMaterialMixture( const vector<G4String>& wl, const G4String& mixtType )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " adding material mixture " << wl[1] << G4endl;
#endif

  //---------- Look if material already exists
  if( FindMaterial( G4tgrUtils::SubQuotes(wl[1]) ) != 0 ) {
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
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
     G4cout << " G4tgrMaterialFactory::FindIsotope " << name << G4endl;
#endif

  mstgrisot::const_iterator cite;
  cite = theG4tgrIsotopes.find( name ); 
  if( cite == theG4tgrIsotopes.end() ) {
    return 0;
  } else {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4cout << " G4tgrIsotope found " << ( (*cite).second )->GetName() << G4endl;
#endif
   return (*cite).second;
  }
}

//-------------------------------------------------------------
G4tgrElement* G4tgrMaterialFactory::FindElement(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgrMaterialFactory::FindElement " << name << G4endl;
#endif
  mstgrelem::const_iterator cite;
  cite = theG4tgrElements.find( name ); 
  if( cite == theG4tgrElements.end() ) {
    return 0;
  } else {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    DumpElementList();
    G4cout << " G4tgrElement found " << ( (*cite).second )->GetName() << G4endl;
#endif
   return (*cite).second;
  }
}


//-------------------------------------------------------------
G4tgrMaterial* G4tgrMaterialFactory::FindMaterial(const G4String & name) const 
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgrMaterialFactory::FindMaterial " << name << G4endl;
#endif
  mstgrmate::const_iterator cite;
  cite = theG4tgrMaterials.find( name );
  if( cite == theG4tgrMaterials.end() ) {
    return 0;
  } else {
    return (*cite).second;
  }

}


//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpIsotopeList() const
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrIsotope's List " << G4endl;
  mstgrisot::const_iterator cite;
  for(cite = theG4tgrIsotopes.begin(); cite != theG4tgrIsotopes.end(); cite++) {
    G4cout << " ISOT: " << (*cite).second->GetName() << G4endl;
  }
}

//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpElementList() const 
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrElement's List " << G4endl;
  mstgrelem::const_iterator cite;
  for(cite = theG4tgrElements.begin(); cite != theG4tgrElements.end(); cite++) {
    G4cout << " ELEM: " << (*cite).second->GetName() << G4endl;
  }
}


//-------------------------------------------------------------
void G4tgrMaterialFactory::DumpMaterialList() const
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrMaterial's List " << G4endl;
  mstgrmate::const_iterator cite;
  for(cite = theG4tgrMaterials.begin(); cite != theG4tgrMaterials.end(); cite++) {
    G4tgrMaterial* mate = (*cite).second;
    G4cout << " MATE: " << mate->GetName() << " Type: " << mate->GetType() 
         << " NoComponents= " << mate->GetNumberOfComponents() << G4endl;
  }
}
 

//-------------------------------------------------------------
void G4tgrMaterialFactory::ErrorAlreadyExists(const G4String& object,  const vector<G4String>& wl, const G4bool bNoRepeating )
{
  G4String msg = object + G4String(" repeated");
  if( bNoRepeating ) {
    G4tgrUtils::DumpVS( wl, (G4String("!!!! EXITING: ") + msg).c_str() );
    G4Exception("");
  } else {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4tgrUtils::DumpVS( wl, (G4String("!! WARNING: ") + msg).c_str() ); 
#endif
  }
}
  
