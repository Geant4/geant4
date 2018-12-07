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
/// \file optical/LXe/include/LXeWLSSlab.hh
/// \brief Definition of the LXeWLSSlab class
//
//
#ifndef LXeWLSSlab_H
#define LXeWLSSlab_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeWLSSlab : public G4PVPlacement
{
  public:

    LXeWLSSlab(G4RotationMatrix *pRot,
               const G4ThreeVector &tlate,
               G4LogicalVolume *pMotherLogical,
               G4bool pMany,
               G4int pCopyNo,
               LXeDetectorConstruction* c);

  private:

    void CopyValues();

    LXeDetectorConstruction* fConstructor;

    static G4LogicalVolume* fScintSlab_log;

    G4int fNfibers;
    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;
    G4double fSlab_z;
};

#endif
