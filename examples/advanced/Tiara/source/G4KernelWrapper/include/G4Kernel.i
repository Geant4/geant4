# $Id: G4Kernel.i,v 1.2 2003-06-16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4String.hh"
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
#include "g4std/strstream"
#include "G4Dimensions.hh"
#include "G4VisManager.hh"
#include "G4VUserDetectorConstruction.hh"
%}

%include std_string.i

//%include G4ThreeVector.hh

%import CLHEP.i
%inline %{
  typedef bool G4bool;
  typedef double G4double;
  typedef int G4int;
%}
typedef std::string G4String;

namespace G4Dimensions {
  const G4double cm;
  const G4double MeV;
  const G4double eV;
}


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
    G4std::ostrstream tmpout;
    self->Print(cs, &tmpout);
    G4std::string *value = new G4std::string(tmpout.str());
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
public:  // with description
  %extend {
    G4Box(const char* pName, G4double pX, G4double pY, G4double pZ){
      G4String name(pName);
      return new G4Box(name, pX, pY, pZ);
    }
  }
      // Construct a box with name, and half lengths pX,pY,pZ
};


class G4Tubs : public G4CSGSolid
{
public:  // with description
  %extend {
    G4Tubs(const char* pName, 
	   G4double pRMin,
	   G4double pRMax,
	   G4double pDz,
	   G4double pSPhi,
	   G4double pDPhi ){
      G4String name(pName);
      return new G4Tubs(name,
			pRMin,
			pRMax,
			pDz,
			pSPhi,
			pDPhi );
    }
  }
};



class G4LogicalVolume
{
public:  // with description
  %extend {
    G4LogicalVolume(G4VSolid *pSolid,
                    G4Material *pMaterial,
		    const char *pName,
                    G4FieldManager* pFieldMgr=0,
                    G4VSensitiveDetector* pSDetector=0,
                    G4UserLimits* pULimits=0,
                    G4bool optimise=true){
      G4String name(pName);
      return new G4LogicalVolume(pSolid,
				 pMaterial,
				 name,
				 pFieldMgr,
				 pSDetector,
				 pULimits,
				 optimise);
      // Constructor. The solid and material pointer must be non null.
      // The parameters for field, detector and user limits are optional.
      // The volume also enters itself into the logical volume Store.
      // Optimisation of the geometry (voxelisation) for the volume
      // hierarchy is applied by default. For parameterised volumes in
      // the hierarchy, optimisation is -always- applied.
    }  
  }  
  

};




class G4RotationMatrix{
public:
  G4RotationMatrix();
};

class G4PVPlacement : public G4VPhysicalVolume
{
public:  
  %extend {
    G4PVPlacement(G4RotationMatrix *pRot,
		  const Hep3Vector &tlate,
		  G4LogicalVolume *pCurrentLogical,
		  const char *pName,
		  G4LogicalVolume *pMotherLogical = 0,
		  G4bool pMany = false,
		  G4int pCopyNo = 0){
      return new G4PVPlacement(pRot,
			       tlate,
			       pCurrentLogical,
			       pName,
			       pMotherLogical,
			       pMany,
			       pCopyNo);
    }
  }
};



