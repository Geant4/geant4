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
/// \file DetConstruction.hh
/// \brief Definition of the DetConstruction class
///
/// \author J. Yarba; FNAL
///

#ifndef DetConstruction_H
#define DetConstruction_H

#include "G4VUserDetectorConstruction.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class G4VPhysicalVolume;

class DetConstruction : public G4VUserDetectorConstruction
{

   public:

      // ctor & dtor
      //
      DetConstruction() : G4VUserDetectorConstruction(), 
                          fWorld(0), fDet(0)  {}
      ~DetConstruction() {} /* we don't need to delete fWorld or fDet 
                               as it'll be done "globally" */
    
      // methods/functions
      // 
      virtual G4VPhysicalVolume* Construct();
   
   private:

     // data members
     //
     G4VPhysicalVolume* fWorld;
     G4VPhysicalVolume* fDet;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
