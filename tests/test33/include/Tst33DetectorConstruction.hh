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
// $Id: Tst33DetectorConstruction.hh,v 1.3 2006-06-29 21:59:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33DetectorConstruction
//
// Class description:
// This class only takes a pointer to the world volume and gives it
// to the Geant4 kernel when asked for it. The detector construction
// is done elsewhere.

#ifndef Tst33DetectorConstruction_hh
#define Tst33DetectorConstruction_hh Tst33DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class Tst33DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  explicit Tst33DetectorConstruction();
  virtual ~Tst33DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  void SetWorldVolume(G4VPhysicalVolume *v);
  
private:
  Tst33DetectorConstruction(const Tst33DetectorConstruction &);
  Tst33DetectorConstruction &operator=(const Tst33DetectorConstruction &);
  G4VPhysicalVolume* fWorldVolume;
};

#endif
