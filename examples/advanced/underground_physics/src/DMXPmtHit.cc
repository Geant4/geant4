//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PmtHit (sensitive detector) program
// --------------------------------------------------------------

#include "DMXPmtHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "g4std/iomanip"

G4Allocator<DMXPmtHit> DMXPmtHitsAllocator;


DMXPmtHit::DMXPmtHit() {
  
  pos=(0., 0., 0.);
  time=0.;

}


DMXPmtHit::~DMXPmtHit() {;}


DMXPmtHit::DMXPmtHit(const DMXPmtHit& right) {

  pos  = right.pos;
  time = right.time;

}



const DMXPmtHit& DMXPmtHit::operator=(const DMXPmtHit& right) {

  pos  = right.pos;
  time = right.time;

  return *this;

}


int DMXPmtHit::operator==(const DMXPmtHit& right) const {

  return 0;
}



void DMXPmtHit::Draw()  {

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    G4Circle circle(pos);
    circle.SetScreenSize(0.002);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,1.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }


}


void DMXPmtHit::Print() {

  //  G4cout << "      PMT hit: " << G4BestUnit(pos,"Length") << G4endl;

  G4cout << "      PMT hit: " << G4std::setw(5) << G4BestUnit(time,"Time") 
	 << ", at " << G4BestUnit(pos,"Length") << G4endl;


}



