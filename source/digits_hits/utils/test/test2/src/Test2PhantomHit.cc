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

#include "Test2PhantomHit.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"


G4Allocator<Test2PhantomHit> Test2PhantomHitAllocator;

Test2PhantomHit::Test2PhantomHit() {
fSegment[0]=1;
fSegment[1]=1;
fSegment[2]=1;
}

Test2PhantomHit::Test2PhantomHit(G4int & x, G4int & y, G4int & z,
				 G4int segment[3])
  : fXCellID(x), fYCellID(y), fZCellID(z){
  fSegment[0]=segment[0];
  fSegment[1]=segment[1];
  fSegment[2]=segment[2];

  fEdep = 0.0;
  fDose = 0.0;
  fTrackLength = 0.0;
  flatSurfCurr=0.0;
  passageCellCurr = 0.0;
  flatSurfaceFlux = 0.0;
  cellFlux = 0.0;
  passageCellFlux = 0.0;
  noOfSecondary = 0.0;
  cellCharge = 0.0;

}

Test2PhantomHit::~Test2PhantomHit() {
  ;
}

Test2PhantomHit::Test2PhantomHit(const Test2PhantomHit & right)
  : G4VHit() {

  fXCellID = right.fXCellID;
  fYCellID = right.fYCellID;
  fZCellID = right.fZCellID;
  fSegment[0] = right.fSegment[0];
  fSegment[1] = right.fSegment[1];
  fSegment[2] = right.fSegment[2];
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
  fParticleName = right.fParticleName;
}

const Test2PhantomHit & Test2PhantomHit::operator=(const Test2PhantomHit & right) {

  fXCellID = right.fXCellID;
  fYCellID = right.fYCellID;
  fZCellID = right.fZCellID;
  fSegment[0] = right.fSegment[0];
  fSegment[1] = right.fSegment[1];
  fSegment[2] = right.fSegment[2];
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
  fParticleName = right.fParticleName;

  return *this;
}


G4int Test2PhantomHit::operator==(const Test2PhantomHit &right) const {

  return ((fXCellID == right.fXCellID) &&
	  (fYCellID == right.fYCellID) &&
	  (fZCellID == right.fZCellID));
}

void Test2PhantomHit::Print() {
  ;
}


