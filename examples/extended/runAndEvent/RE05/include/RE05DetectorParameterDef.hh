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
/// \file RE05/include/RE05DetectorParameterDef.hh
/// \brief Definition of the RE05DetectorParameterDef class
//

     G4double fExpHall_x;
     G4double fExpHall_y;
     G4double fExpHall_z;

     G4double fTrkTubs_rmax;
     G4double fTrkTubs_rmin;
     G4double fTrkTubs_dz;
     G4double fTrkTubs_sphi;
     G4double fTrkTubs_dphi;

     G4int fNotrkLayers;
     G4double fTracker_radius[5];
     G4double fTracker_thick;
     G4double fTracker_length[5];

     G4double fCaloTubs_rmax;
     G4double fCaloTubs_rmin;
     G4double fCaloTubs_dz;
     G4double fCaloTubs_sphi;
     G4double fCaloTubs_dphi;

     G4int fNocaloLayers;
     G4double fAbsorber_thick;
     G4double fScinti_thick;

     G4int fSegmentsinZ;
     G4double fCaloRing_rmax;
     G4double fCaloRing_rmin;
     G4double fCaloRing_dz;
     G4double fCaloRing_sphi;
     G4double fCaloRing_dphi;

     G4int fSegmentsinPhi;
     G4double fCaloCell_rmax;
     G4double fCaloCell_rmin;
     G4double fCaloCell_dz;
     G4double fCaloCell_sphi;
     G4double fCaloCell_dphi;

     G4int fNomucounter;
     G4double fMuBox_radius;
     G4double fMuBox_width;
     G4double fMuBox_thick;
     G4double fMuBox_length;
