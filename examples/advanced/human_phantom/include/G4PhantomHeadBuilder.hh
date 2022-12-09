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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//

#ifndef G4PhantomHeadBuilder_h
#define G4PhantomHeadBuilder_h 1

#include "G4VPhysicalVolume.hh"
#include "G4BasePhantomBuilder.hh"
#include "globals.hh"

class G4BasePhantomBuilder;
class G4VPhysicalVolume;
class G4VBodyFactory;
class G4PhantomHeadBuilder: public G4BasePhantomBuilder
{
public:
  G4PhantomHeadBuilder();
  ~G4PhantomHeadBuilder();

  void BuildHead(const G4String&,G4bool,G4bool);
  void BuildSkull(const G4String&,G4bool,G4bool);
  void BuildBrain(const G4String&,G4bool,G4bool);

  void SetModel(G4String);
  void SetMotherVolume(G4VPhysicalVolume*);


G4VPhysicalVolume* GetPhantom();

protected: 
  G4VBodyFactory* fBody;
  G4String fModel;

  G4VPhysicalVolume* fMotherVolume;
  G4VPhysicalVolume* fHeadVolume;
};
#endif
