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
// $Id: Tst33CellScorer.cc,v 1.2 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33CellScorer.cc
//
// ----------------------------------------------------------------------

#include "Tst33CellScorer.hh"
#include "G4GeometryCell.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"


Tst33CellScorer::Tst33CellScorer()
{}

Tst33CellScorer::~Tst33CellScorer()
{}

void Tst33CellScorer::ScoreAnExitingStep(const G4Step &aStep,
				       const G4GeometryCell &pre_gCell){
  fG4CellScorer.ScoreAnExitingStep(aStep, pre_gCell);
}

void Tst33CellScorer::ScoreAnEnteringStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnEnteringStep(aStep, post_gCell);
  return;
}

void Tst33CellScorer::ScoreAnInVolumeStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnInVolumeStep(aStep, post_gCell);
}

