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
/// \file eventgenerator/HepMC/HepMCEx01/include/ExN04DetectorParameterDef.hh
/// \brief Definition of the ExN04DetectorParameterDef class
//
//

G4double fexpHall_x;
G4double fexpHall_y;
G4double fexpHall_z;

G4double ftrkTubs_rmax;
G4double ftrkTubs_rmin;
G4double ftrkTubs_dz;
G4double ftrkTubs_sphi;
G4double ftrkTubs_dphi;

G4int fnotrkLayers;
G4double ftracker_radius[5];
G4double ftracker_thick;
G4double ftracker_length[5];

G4double fcaloTubs_rmax;
G4double fcaloTubs_rmin;
G4double fcaloTubs_dz;
G4double fcaloTubs_sphi;
G4double fcaloTubs_dphi;

G4int fnocaloLayers;
G4double fabsorber_thick;
G4double fscinti_thick;

G4int fsegmentsinZ;
G4double fcaloRing_rmax;
G4double fcaloRing_rmin;
G4double fcaloRing_dz;
G4double fcaloRing_sphi;
G4double fcaloRing_dphi;

G4int fsegmentsinPhi;
G4double fcaloCell_rmax;
G4double fcaloCell_rmin;
G4double fcaloCell_dz;
G4double fcaloCell_sphi;
G4double fcaloCell_dphi;

G4int fnomucounter;
G4double fmuBox_radius;
G4double fmuBox_width;
G4double fmuBox_thick;
G4double fmuBox_length;
