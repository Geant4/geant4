// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08ChamberParameterisation.cc,v 1.1 1999-01-08 16:35:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "T08ChamberParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

T08ChamberParameterisation::T08ChamberParameterisation(  
        G4int    NoChambers, 
        G4double startZ,          //  Z of center of first 
        G4double spacingZ,        //  Z spacing of centers
        G4double widthChamber, 
        G4double lengthInitial, 
        G4double lengthFinal )
{
   fNoChambers=  NoChambers; 
   fStartZ=      startZ; 
   fHalfWidth=   widthChamber*0.5;
   fSpacing=     spacingZ;
   fHalfLengthFirst= 0.5 * lengthInitial; 
   // fHalfLengthLast=  lengthFinal;
   if( NoChambers > 1 ){
      fHalfLengthIncr=  0.5 * (lengthFinal-lengthInitial)/(NoChambers-1);

      if( spacingZ < widthChamber ) {
         G4Exception( "T08ChamberParameterisation construction: Width > Spacing" );
      }
   }
   
}

T08ChamberParameterisation::~T08ChamberParameterisation()
{}

void T08ChamberParameterisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4double      Zposition= fStartZ + copyNo * fSpacing;
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

void T08ChamberParameterisation::ComputeDimensions
(G4Box & trackerChamber, const G4int copyNo,
 const G4VPhysicalVolume * physVol) const
{
  G4double  halfLength= fHalfLengthFirst + (copyNo-1) * fHalfLengthIncr;
  trackerChamber.SetXHalfLength(halfLength);
  trackerChamber.SetYHalfLength(halfLength);
  trackerChamber.SetZHalfLength(fHalfWidth);
}
