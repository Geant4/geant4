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
// $Id: ExP02GeoTree.hh,v 1.2 2007-12-10 16:29:11 gunter Exp $
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
