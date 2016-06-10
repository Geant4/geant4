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
// $Id: G4VFlavoredParallelWorld.hh 67965 2013-03-13 09:35:29Z gcosmo $
//
// 
// Abstract interface for GEANT4 Flavored Parallel World.
// P. Mora de Freitas & M. Verderi 14/April/1999.
//

#ifndef G4VFLAVOREDPARALLELWORLD_HH
#define G4VFLAVOREDPARALLELWORLD_HH

class G4VPhysicalVolume;

class G4VFlavoredParallelWorld {

public:

  virtual ~G4VFlavoredParallelWorld () {}

  // G4VFlavoredParallelWorld Interface for visualisation.

  virtual 
  G4VPhysicalVolume* GetThePhysicalVolumeWorld() const =0;
};
#endif
