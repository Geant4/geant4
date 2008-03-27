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
// FredMessenger.hh
//
// Declaration of fred's options
//

#ifndef FredMessenger_HH
#define FredMessenger_HH

#include "G4UImessenger.hh"

enum VolumeType {
	NATALIA,
	VOXEL,
	BOX,
	CONE,
	CONE2,
        ORB,
	PARA,
	SPHERE,
	TORUS1,
	TORUS2,
	TRAP,
	TRD,
	TUBS,
        ELLIPS,
        ELCONE,
        ELTUBE,
        EXTRU1,
        EXTRU2,
        EXTRU3,
        EXTRU4,
        HYPE,
	PCON,
	PCON2,
	PCON3,
        PCON4,
	PGON,
	BOOL1,
	PGON2,
	PGON3,
        TESSEL1,
        TESSEL2,
        TESSEL3,
        TESSEL4,
        TET,
        TWBOX,
        TWTRAP1,
        TWTRAP2,
        TWTRD,
        TWTUBS
};
#define FRED_VOLUMETYPE_NUM 39

enum GunType {
	SPRAY,
	GRID,
	G4
};
#define FRED_GUNTYPE_NUM  3

enum DrawType {
	NORMAL,
	SHADOW
};
#define FRED_DRAWTYPE_NUM  2

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

class G4VSolid;
class FredTest3Messenger;
class FredVoxelTestMessenger;

class FredMessenger : public G4UImessenger
{
	public:
	FredMessenger();
	~FredMessenger();
	
        void SetNewValue( G4UIcommand *command, G4String newValues );
        G4String GetCurrentValue( G4UIcommand *command );

	VolumeType	SelectedVolume();
	GunType		SelectedGun();
	DrawType	SelectedDrawing();
	
	G4double	StartPhi() { return startPhi; }
	G4double	DeltaPhi() { return deltaPhi; }
	G4int		NumSide()  { return numSide; }
	
	void SetTestVolume( const G4VSolid *theTestVolume );
	inline const G4VSolid *GetTestVolume() const { return testVolume; }
	
	private:
	void PauseInput();

	private:
	VolumeType	testVolumeType;
	GunType		gunType;
	DrawType	drawType;
	G4double	startPhi, 
			deltaPhi;
	G4int		numSide;
			
	G4String	volumeNames[FRED_VOLUMETYPE_NUM];
	
	G4String	gunNames[FRED_GUNTYPE_NUM];
	
	G4String	drawNames[FRED_DRAWTYPE_NUM];

	G4UIdirectory		*fredDirectory;
	G4UIcmdWithAString	*volumeTypeNameCmd;
	G4UIcmdWithAString	*gunTypeNameCmd;
	G4UIcmdWithAString	*drawTypeNameCmd;
	G4UIcmdWithoutParameter	*pauseCmd;
	G4UIcmdWithADouble	*startPhiCmd;
	G4UIcmdWithADouble	*deltaPhiCmd;
	G4UIcmdWithAnInteger	*numSideCmd;
	
	FredTest3Messenger	*test3Messenger;
	FredVoxelTestMessenger  *voxelTestMessenger;
 	const G4VSolid		*testVolume;
};

#endif
