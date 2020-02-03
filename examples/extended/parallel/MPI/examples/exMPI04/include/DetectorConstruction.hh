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
<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4InteractionCode.hh
// $Id: G4InteractionCode.hh 67999 2013-03-13 11:14:32Z gcosmo $
//
#ifndef G4InteractionCode_h
#define G4InteractionCode_h 1
=======
/// @file DetectorConstruction.hh
/// @brief Define geometry
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parallel/MPI/examples/exMPI04/include/DetectorConstruction.hh

#ifndef DETECTOR_CONSTRUCTION_H
#define DETECTOR_CONSTRUCTION_H

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4InteractionCode.hh
class G4InteractionCode
{

  public:
    G4InteractionCode(G4String & aCode);
=======
#include "G4VUserDetectorConstruction.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parallel/MPI/examples/exMPI04/include/DetectorConstruction.hh

class G4LogicalVolume;

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4InteractionCode.hh
  private:

    G4String theCode;
};
=======
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction();
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parallel/MPI/examples/exMPI04/include/DetectorConstruction.hh

  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

private:
  G4LogicalVolume* flv_voxel;

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/management/include/G4InteractionCode.hh
=======
};

>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/parallel/MPI/examples/exMPI04/include/DetectorConstruction.hh
#endif
