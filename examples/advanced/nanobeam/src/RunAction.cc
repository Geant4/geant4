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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include <fstream>
#include <fstream>
#include <vector>
#include <cmath>

// #define MATRIX_BOUND_CHECK

#include "RunAction.hh"

#include "G4ios.hh"
#include "Randomize.hh"

#include "G4SteppingManager.hh"
#include "G4Run.hh"
#include "G4Material.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* pri,
HistoManager* hi)
:detector(det),primary(pri),hist(hi)
{   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
   // Vector initialization 
   row=0;
     
   xVector = CLHEP::HepVector(32);
   yVector = CLHEP::HepVector(32);
   thetaVector = CLHEP::HepVector(32);
   phiVector = CLHEP::HepVector(32);
   
  // Histograms
  hist->book();

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{

if (detector->GetCoef()==1)
{

	CLHEP::HepMatrix m;
	
	// VECTOR READING
	// VECTOR READING

	m = CLHEP::HepMatrix(32,32);
	m = primary->GetMatrix();

	G4cout << G4endl;
	G4cout << "===> NANOBEAM LINE INTRINSIC ABERRATION COEFFICIENTS (units of micrometer and mrad) :" << G4endl;
	G4cout << G4endl;

	int inv;

	m.invert(inv);
	CLHEP::HepVector tmp(32,0);
	tmp=m*xVector;
	CLHEP::HepVector b;
	b=tmp.sub(2,2);   G4cout << "<x|theta>=" << b << G4endl;
	b=tmp.sub(8,8);   G4cout << "<x|theta*delta>=" << b << G4endl;
	b=tmp.sub(10,10); G4cout << "<x|theta^3>=" << b << G4endl;
	b=tmp.sub(12,12); G4cout << "<x|theta*phi^2>=" << b << G4endl;
	m.invert(inv);

	m.invert(inv);
	tmp = m*thetaVector;
	m.invert(inv);
	b=tmp.sub(2,2); G4cout << "<x|x>=" << b << G4endl;

	m.invert(inv);
	tmp=m*yVector;
	b=tmp.sub(3,3);   G4cout << "<y|phi>=" << b << G4endl;
	b=tmp.sub(9,9);   G4cout << "<y|phi*delta>=" << b << G4endl;
	b=tmp.sub(11,11); G4cout << "<y|theta^2*phi>=" << b << G4endl;
	b=tmp.sub(13,13); G4cout << "<y|phi^3>=" << b << G4endl;
	m.invert(inv);

	m.invert(inv);
	tmp = m*phiVector;
	m.invert(inv);
	b=tmp.sub(3,3); G4cout << "<y|y>=" << b << G4endl;

}

  //save histograms      
  hist->save();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
