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
// $Id: G4PBox.cc,v 1.3 2001/07/11 10:02:24 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// 
//
// Implementation for G4PBox class
//
// History:
// 19.06.98 A.Kimura Converted G4Box.cc


#include "G4VSolid.hh"
#include "G4PBox.hh"
#include "G4Box.hh"

// Constructor - check & set half widths
G4PBox::G4PBox(const G4Box* theBox) : G4PCSGSolid(theBox->GetName()) {
    G4double pX = theBox->GetXHalfLength();
    G4double pY = theBox->GetYHalfLength();
    G4double pZ = theBox->GetZHalfLength();

    if (pX>0&&pY>0&&pZ>0) {
	fDx=pX; fDy=pY; fDz=pZ;
    } else {
	G4Exception("Error in G4PBox::Box - negative parameters");
    }	

}

// Destructor
G4PBox::~G4PBox() {
    ;
}

// make a transient object
G4VSolid* G4PBox::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Box(GetName(),
				       fDx, fDy, fDz);
    return transientObject;
}

