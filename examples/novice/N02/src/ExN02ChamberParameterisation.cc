// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02ChamberParameterisation.cc,v 1.3 2000-12-04 16:24:07 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "ExN02ChamberParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02ChamberParameterisation::ExN02ChamberParameterisation(  
        G4int    NoChambers, 
        G4double startZ,          //  Z of center of first 
        G4double spacingZ,        //  Z spacing of centers
        G4double widthChamber, 
        G4double lengthInitial, 
        G4double lengthFinal )
{
   fNoChambers =  NoChambers; 
   fStartZ     =  startZ; 
   fHalfWidth  =  widthChamber*0.5;
   fSpacing    =  spacingZ;
   fHalfLengthFirst = 0.5 * lengthInitial; 
   // fHalfLengthLast = lengthFinal;
   if( NoChambers > 1 ){
      fHalfLengthIncr =  0.5 * (lengthFinal-lengthInitial)/(NoChambers-1);

      if( spacingZ < widthChamber ) {
         G4Exception( "ExN02ChamberParameterisation construction: Width > Spacing" );
      }
   }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02ChamberParameterisation::~ExN02ChamberParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExN02ChamberParameterisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4double      Zposition= fStartZ + copyNo * fSpacing;
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExN02ChamberParameterisation::ComputeDimensions
(G4Box & trackerChamber, const G4int copyNo,
 const G4VPhysicalVolume * physVol) const
{
  G4double  halfLength= fHalfLengthFirst + (copyNo-1) * fHalfLengthIncr;
  trackerChamber.SetXHalfLength(halfLength);
  trackerChamber.SetYHalfLength(halfLength);
  trackerChamber.SetZHalfLength(fHalfWidth);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
