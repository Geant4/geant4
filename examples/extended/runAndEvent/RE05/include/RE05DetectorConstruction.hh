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
<<<<<<< HEAD
// $Id: RE05DetectorConstruction.hh 69764 2013-05-14 09:59:36Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file RE05/include/RE05DetectorConstruction.hh
/// \brief Definition of the RE05DetectorConstruction class
//

#ifndef RE05DetectorConstruction_h
#define RE05DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;
class G4Element;

class RE05DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    RE05DetectorConstruction();
    virtual ~RE05DetectorConstruction();

  public:
     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();

  private:
     void DefineMaterials();

#include "RE05DetectorParameterDef.hh"

  G4Material* Air;
  G4Material* Ar;
  G4Material* Silicon;
  G4Material* Scinti;
  G4Material* Lead;

  G4Element* H;
  G4Element* C;
  G4Element* N;
  G4Element* O;

};

#endif

