//
// FredMessenger.hh
//
// Declaration of fred's options
//

#ifndef FredMessenger_HH
#define FredMessenger_HH

#include "G4UImessenger.hh"

enum VolumeType {
	TUBS,
	SPHERE,
	BOX,
	PCON,
	PCON2,
	PCON3,
	PCON4,
	PGON,
	CONE,
	CONE2,
	VOXEL,
	NATALIA,
	BOOL1,
	PGON2,
	PGON3,
	PGON4,
	TRAP,
	PARA,
	TORUS1,
	TORUS2,
	TRD
};
#define FRED_VOLUMETYPE_NUM 21

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
	
	const G4double	StartPhi() { return startPhi; }
	const G4double	DeltaPhi() { return deltaPhi; }
	const G4int	NumSide()  { return numSide; }
	
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
