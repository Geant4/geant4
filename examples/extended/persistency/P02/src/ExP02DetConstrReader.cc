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
/// \file persistency/P02/src/ExP02DetConstrReader.cc
/// \brief Implementation of the ExP02DetConstrReader class
//
//
//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"

//G4
#include "G4Element.hh"
#include "G4PVPlacement.hh"
#include "ExP02GeoTree.hh"

// local
#include "ExP02DetConstrReader.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP02DetConstrReader::ExP02DetConstrReader()
 : G4VUserDetectorConstruction()
{  
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libExP02ClassesDict");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP02DetConstrReader::~ExP02DetConstrReader()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExP02DetConstrReader::Construct()
{
  ExP02GeoTree* geotree;

  TFile fo("geo.root");
  fo.GetObject("my_geo", geotree);

  return geotree->TopVol();
}
