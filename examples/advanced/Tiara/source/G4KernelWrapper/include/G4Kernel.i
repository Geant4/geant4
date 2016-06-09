# $Id: G4Kernel.i,v 1.5 2003/06/20 12:41:06 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-05-02 $
# -------------------------------------------------------------------

%module G4Kernel
%{
#include <memory>
#include <string>
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4FieldManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4UserLimits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VScorer.hh"
#include "G4VCellScorer.hh"
#include "G4CellScorer.hh"
#include "G4CellStoreScorer.hh"
#include "G4VCellScorerStore.hh"
#include "G4ScoreTable.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "G4GeometryCell.hh"
#include "G4CellScoreValues.hh"
#include "G4CellScorerStore.hh"
#include "G4Scorer.hh"
#include "G4VSampler.hh"
#include "G4MassGeometrySampler.hh"
#include "G4ParallelGeometrySampler.hh"
#include "G4VImportanceAlgorithm.hh"
#include "G4CellScoreComposer.hh"
#include "G4VIStore.hh"
#include "G4IStore.hh"
#include "G4UserSteppingAction.hh"
#include <strstream>
#include "G4VisManager.hh"
#include "G4VUserDetectorConstruction.hh"
%}

%include std_string.i
%include G4String.i

%import CLHEP.i
%inline %{
  typedef bool G4bool;
  typedef double G4double;
  typedef int G4int;
%}


%include G4VIStore.hh

%include G4GeometryCell.hh

%include G4IStore.hh

%include G4CellScoreValues.hh

%include G4CellScoreComposer.hh

%include G4VScorer.hh

%include G4VCellScorer.hh
%include G4CellScorer.hh

%include G4VCellScorerStore.hh
%include G4CellScorerStore.hh

%include G4CellStoreScorer.hh
%include G4Scorer.hh

%include G4VImportanceAlgorithm.hh

%include G4VSampler.hh

%include G4MassGeometrySampler.hh

%include G4ParallelGeometrySampler.hh

%include G4ScoreTable.hh
%extend G4ScoreTable {
  const char *Write(const G4MapGeometryCellCellScorer &cs){
    std::ostrstream tmpout;
    self->Print(cs, &tmpout);
    std::string *value = new std::string(tmpout.str());
    return value->c_str();
  };
}


class G4VisManager {
public:
  void Initialize();
};



class G4VPhysicalVolume{
};

class G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct();
};


class G4VUserPrimaryGeneratorAction {
};

class G4VUserPhysicsList {
public:
  void Construct();  
};

class G4VModularPhysicsList: public G4VUserPhysicsList {
public:
  void Construct();  
};

class G4UserSteppingAction{
};

class G4UserEventAction{
};

class G4FieldManager {
};

class G4VSensitiveDetector {
};

class G4UserLimits {
};

class G4Material {
};


class G4VSolid {
};

class G4CSGSolid : public G4VSolid{
};

class G4Box : public G4CSGSolid
{
public:  
  G4Box(const G4String &pName, G4double pX, G4double pY, G4double pZ);
};


class G4Tubs : public G4CSGSolid
{
public:  
  G4Tubs(const G4String &pName, 
	 G4double pRMin,
	 G4double pRMax,
	 G4double pDz,
	 G4double pSPhi,
	 G4double pDPhi );
};



class G4LogicalVolume
{
public:  
  G4LogicalVolume(G4VSolid *pSolid,
		  G4Material *pMaterial,
		  const G4String &pName,
		  G4FieldManager* pFieldMgr=0,
		  G4VSensitiveDetector* pSDetector=0,
		  G4UserLimits* pULimits=0,
		  G4bool optimise=true);
  // Constructor. The solid and material pointer must be non null.
  // The parameters for field, detector and user limits are optional.
  // The volume also enters itself into the logical volume Store.
  // Optimisation of the geometry (voxelisation) for the volume
  // hierarchy is applied by default. For parameterised volumes in
  // the hierarchy, optimisation is -always- applied.

};




class G4RotationMatrix{
public:
  G4RotationMatrix();
};

class G4PVPlacement : public G4VPhysicalVolume
{
public:  
  G4PVPlacement(G4RotationMatrix *pRot,
		const Hep3Vector &tlate,
		G4LogicalVolume *pCurrentLogical,
		const G4String &pName,
		G4LogicalVolume *pMotherLogical = 0,
		G4bool pMany = false,
		G4int pCopyNo = 0);
};



