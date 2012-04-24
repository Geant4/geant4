
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "BeamTestCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamTestCellParameterisation::BeamTestCellParameterisation(  
        G4int    NoChambers, 
        G4double startZ,          //  Z of center of first 
        G4double spacingZ,        //  Z spacing of centers
        G4double widthChamber )
{
	fNoChambers =  NoChambers; 
	fStartZ     =  startZ; 
	fHalfWidth  =  widthChamber*0.5;
	fSpacing    =  spacingZ;
	; 
	// fHalfLengthLast = lengthFinal;
	/*if( NoChambers > 0 )
	  {
	  fHalfLengthIncr =  0.5 * (lengthFinal-lengthInitial)/NoChambers;
	  if (spacingZ < widthChamber) 
	  {
	  G4Exception("ExN02ChamberParameterisation construction: Width>Spacing");
	  }
	  }*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamTestCellParameterisation::~BeamTestCellParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamTestCellParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
	/*if(copyNo==0)
	{
		G4double Zposition= fStartZ + fHalfWidth*2 ;
		G4ThreeVector origin(0,0,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(0);
	}*/
	//else
//	{
    G4double Zposition= fStartZ + (copyNo) * (fSpacing + fHalfWidth*2);
		G4ThreeVector origin(0,0,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(0);

//	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void BeamTestCellParameterisation::ComputeDimensions
  (G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const
  {
  G4double  halfLength= fHalfLengthFirst + copyNo * fHalfLengthIncr;
  trackerChamber.SetXHalfLength(halfLength);
  trackerChamber.SetYHalfLength(halfLength);
  trackerChamber.SetZHalfLength(fHalfWidth);
  }
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
