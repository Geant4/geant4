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
/// \file ExN04DetectorParameterDef.hh
/// \brief Definition of the ExN04DetectorParameterDef class
//

     G4double expHall_x;
     G4double expHall_y;
     G4double expHall_z;

     G4double trkTubs_rmax;
     G4double trkTubs_rmin;
     G4double trkTubs_dz;
     G4double trkTubs_sphi;
     G4double trkTubs_dphi;

     G4int notrkLayers;
     G4double tracker_radius[5];
     G4double tracker_thick;
     G4double tracker_length[5];

     G4double caloTubs_rmax;
     G4double caloTubs_rmin;
     G4double caloTubs_dz;
     G4double caloTubs_sphi;
     G4double caloTubs_dphi;

     G4int nocaloLayers;
     G4double absorber_thick;
     G4double scinti_thick;

     G4int segmentsinZ;
     G4double caloRing_rmax;
     G4double caloRing_rmin;
     G4double caloRing_dz;
     G4double caloRing_sphi;
     G4double caloRing_dphi;

     G4int segmentsinPhi;
     G4double caloCell_rmax;
     G4double caloCell_rmin;
     G4double caloCell_dz;
     G4double caloCell_sphi;
     G4double caloCell_dphi;

     G4int nomucounter;
     G4double muBox_radius;
     G4double muBox_width;
     G4double muBox_thick;
     G4double muBox_length;


