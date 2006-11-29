#include "G4tgbMaterialSimple.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMessenger.hh"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbMaterialSimple::G4tgbMaterialSimple( G4tgrMaterial* hgmate)
{
  theTgrMate = hgmate;
  G4tgrMaterialSimple* matesimp = static_cast<G4tgrMaterialSimple*>(hgmate);
  theZ = matesimp->GetZ();
  theA = matesimp->GetA();
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4Material* G4tgbMaterialSimple::BuildG4Material()
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgbMaterialSimple::BuildG4Material() " << G4endl;
#endif
  //----- construct new G4Material with no components (only itself)
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) {
    cout << " building new G4MaterialSimple " <<  " " << theTgrMate->GetName() << " " << " " << *this << endl;
    G4cout << " Z " << theTgrMate->GetZ() << " A " << theTgrMate->GetA() << G4endl; 
  }
#endif
  G4Material* mate = new G4Material( GetName(), GetZ(), GetA(), theTgrMate->GetDensity() );
  
  return mate;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ostream& operator<<(ostream& os, const G4tgbMaterialSimple& mate) 
{
  os << "MATERIAL SIMPLE: " << mate.GetName() << endl
     << " Z= " << mate.GetZ() 
     << " A= " << mate.GetA() 
     << " density= " << mate.GetDensity() << endl;
  return os;
}

