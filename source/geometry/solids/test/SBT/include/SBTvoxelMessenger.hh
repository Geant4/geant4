//
// SBTvoxelMessenger.hh
//
// Definition of the UI messenger for SBTvoxel
//

#ifndef SBTvoxelMessenger_hh
#define SBTvoxelMessenger_hh

#include "G4UImessenger.hh"
#include "SBTvoxel.hh"
#include "G4RotationMatrix.hh"

class ifstream;
class G4VSolid;
class G4SolidQuery;

class G4AffineTransform;
class G4VoxelLimits;

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;

class G4UIcmdWithPargs;
class G4UIcmdParg;

class SBTvoxelMessenger : public G4UImessenger {
	public:
	SBTvoxelMessenger( const G4String prefix, const G4SolidQuery *solidQuery, SBTVisManager *visManager );
	~SBTvoxelMessenger();
	
	void SetNewValue( G4UIcommand *command, G4String newValues );
	G4String GetCurrentValue( G4UIcommand *command );
	
	inline const G4SolidQuery *GetSolidQuery() const { return solidQuery; }

        void InvokeTest();
	void Draw();
	void Debug();

	private:
	SBTvoxel		tester;
	const G4SolidQuery 	*solidQuery;
	SBTVisManager	*visManager;
	
	G4VoxelLimits		*voxel;
	G4ThreeVector		*translate;
	G4RotationMatrix	*rotate;
	G4ThreeVector		*point;
	G4double		limits[2];
	EAxis			axis;
	
	G4String		errorFile;
	
	G4UIdirectory		*voxelDirectory;
	G4UIcmdWith3VectorAndUnit	*targetCmd;
	G4UIcmdWith3VectorAndUnit	*widthsCmd;
	G4UIcmdWithAnInteger		*maxVoxelsCmd;
	G4UIcmdWithAnInteger		*maxErrorsCmd;
	G4UIcmdWithAString		*errorFileCmd;
	G4UIcmdWithoutParameter		*runCmd;
	
	G4UIdirectory			*pictDirectory;
	G4UIcmdParg			*pictVoxelArgs[6];
	G4UIcmdWithPargs		*pictVoxelCmd;
	G4UIcmdWith3Vector		*pictTranCmd;
	G4UIcmdParg			*pictRotArgs[4];
	G4UIcmdWithPargs		*pictRotCmd;
	G4UIcmdWith3Vector	        *pictPointCmd;
	G4UIcmdParg			*pictLimitArgs[3];
	G4UIcmdWithPargs		*pictLimitCmd;
	G4UIcmdWithoutParameter		*pictDrawCmd;
	G4UIcmdWithoutParameter		*pictDebugCmd;
};

#endif
