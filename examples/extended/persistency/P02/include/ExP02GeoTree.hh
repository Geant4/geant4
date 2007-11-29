// $Id: ExP02GeoTree.hh,v 1.1 2007-11-29 17:05:22 witoldp Exp $
#ifndef INCLUDE_EXP02GEOTREE_H 
#define INCLUDE_EXP02GEOTREE_H 1

// Include files
#include "G4VPhysicalVolume.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

class ExP02GeoTree {
public: 

  ExP02GeoTree( ); 
  ExP02GeoTree(G4VPhysicalVolume* vol, const G4ElementTable* et, const G4MaterialTable* mt);
  
  virtual ~ExP02GeoTree( );

  G4VPhysicalVolume* TopVol();

private:

  G4VPhysicalVolume* topV;
  const G4ElementTable* eltab;
  const G4MaterialTable* mattab;
  
};
#endif // INCLUDE_EXP02GEOTREE_H
