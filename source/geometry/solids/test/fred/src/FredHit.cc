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
// FredHit.cc
//
// Implementation of FredHit: intersection of track with a test volume
//
// For more comments, you should really take a look at FredHit.hh
//

#include "FredHit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

//
// Implement FredHit's allocator, via template
G4Allocator<FredHit> FredHitAllocator;


//
// Constructor (position and enters specified)
//
FredHit::FredHit( const G4ThreeVector aPos, G4bool aEnters, const G4Track *aTrack )
{
	pos = aPos;
	enters = aEnters;
	track = aTrack;
}


//
// Copy (initialization)
//
FredHit::FredHit( const FredHit &right )
  : G4VHit(right)
{
	// Assign position
	pos = right.pos;
	enters = right.enters;
	track = right.track;
}

//
// Copy (assignment operator)
//
const FredHit& FredHit::operator=( const FredHit &right )
{
	pos = right.pos;
	enters = right.enters;
	track = right.track;
	return *this;
}

//
// Draw
// Our standard drawing
//
void FredHit::Draw()
{
	DrawPrim(0);
}

//
// DrawShadowHit
//
// Special draw mode for drawing tracks that miss the test volume
//
void FredHit::DrawShadowHit()
{
	DrawPrim(1);
}

//
// DrawShadowHit
//
// Special draw mode for drawing tracks that hit the test volume
//
void FredHit::DrawMaskHit()
{
	DrawPrim(3);
}

//
// DrawErrorHit
//
// Special draw mode for drawing tracks with tracking errors (bugs)
//
void FredHit::DrawErrorHit()
{
	DrawPrim(2);
}

//
// DrawPrim
// With all this wonderful c++ technology, why is drawing still
// such a pain?? 
//

void FredHit::DrawPrim(G4int imode)
{
	//
	// Make sure there is something to draw first
	//
	G4VVisManager *visManager = G4VVisManager::GetConcreteInstance();
	
	if (visManager) {
		//
		// Heh. Sometimes something bad happens, and we get
		// a point at infinity. Don't draw such a beasty, as it
		// messes up VRML viewing.
		//
		if (pos.mag() > 200.0*m) return;
	
	
		// How about a little circle?
		//
		// I took this from example N04. Note the temporary
		// nature of the "G4VisAttribute". We can do this
		// because circle.SetVisAttributes copies it.
		
		G4Circle circle( pos );
		if (imode==0) {
			//circle.SetScreenSize( 3.0 );
			circle.SetWorldSize( 5*cm );
			G4Color color( enters ? 0:1, enters ? 0:1, 1 );
			G4VisAttributes attribs( color );
			circle.SetVisAttributes( attribs );
		}
		else if (imode==1) {
			circle.SetScreenSize( 2.0 );
			G4Color color( 0.25, 0.25, 1.0 );
			G4VisAttributes attribs( color );
			circle.SetVisAttributes( attribs );
		}
		else if (imode==2) {
			circle.SetScreenSize( 3.0 );
			G4Color color( 1.0, 1.0, 1.0 );
			G4VisAttributes attribs( color );
			circle.SetVisAttributes( attribs );
		}
		else if (imode==3) {
			circle.SetScreenSize( 3.0 );
			G4Color color( 1.0, 0.25, 0.25 );
			G4VisAttributes attribs( color );
			circle.SetVisAttributes( attribs );
		}
		
		circle.SetFillStyle( G4Circle::filled );
		visManager->Draw( circle );
	}
}

//
// Print
//
void FredHit::Print()
{

	// Print the position
	;
}
