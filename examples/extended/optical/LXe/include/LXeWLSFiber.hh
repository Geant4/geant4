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
/// \file optical/LXe/include/LXeWLSFiber.hh
/// \brief Definition of the LXeWLSFiber class
//
#ifndef LXeWLSFiber_H
#define LXeWLSFiber_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeWLSFiber : public G4PVPlacement
{
  public:

    LXeWLSFiber(G4RotationMatrix *pRot,
                const G4ThreeVector &tlate,
                G4LogicalVolume *pMotherLogical,
                G4bool pMany,
                G4int pCopyNo,
                LXeDetectorConstruction* c);

  private:

    void CopyValues();

    static G4LogicalVolume* fClad2_log;

    G4double fFiber_rmin;
    G4double fFiber_rmax;
    G4double fFiber_z;
    G4double fFiber_sphi;
    G4double fFiber_ephi;

    G4double fClad1_rmin;
    G4double fClad1_rmax;    
    G4double fClad1_z;
    G4double fClad1_sphi;
    G4double fClad1_ephi; 

    G4double fClad2_rmin;
    G4double fClad2_rmax;    
    G4double fClad2_z;
    G4double fClad2_sphi;
    G4double fClad2_ephi;

    LXeDetectorConstruction* fConstructor;
};

#endif
