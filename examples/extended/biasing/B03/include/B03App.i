%module B03App
%{
#include "globals.hh"
#include "B03App.hh"
#include "G4RunManager.hh"
#include "G4IStore.hh"
#include "G4VCellScorerStore.hh"
#include "G4VScorer.hh"
#include "G4ParallelGeometrySampler.hh"
#include "G4CellScorerStore.hh"
#include "G4CellStoreScorer.hh"
#include "G4GeometryCell.hh"
#include "G4CellScoreComposer.hh"
#include "G4CellScoreValues.hh"
#include "G4CellScorer.hh"
#include "G4ScoreTable.hh"
#include "G4VIStore.hh"
#include <string>
#include "g4std/strstream"
#include <memory>
#include "G4VPhysicalVolume.hh"
#include "B03ImportanceDetectorConstruction.hh"
#include "G4VImportanceAlgorithm.hh"
%}

class G4RunManager {
public:
  void BeamOn(int e);
};

class G4VIStore{
};

class G4VPhysicalVolume {
};

class G4GeometryCell{
public:
  G4GeometryCell(const G4VPhysicalVolume &aVolume, int RepNum);
  ~G4GeometryCell();
  const G4VPhysicalVolume &GetPhysicalVolume() const;
  int GetReplicaNumber() const;
};

class G4IStore : public G4VIStore{
public:
  G4IStore(const G4VPhysicalVolume &worldvolume);
  ~G4IStore();
  void AddImportanceGeometryCell(double importance,
			   const G4GeometryCell &gCell);
  void ChangeImportance(double importance,
			const G4GeometryCell &gCell);
  double GetImportance(const G4GeometryCell &gCell) const;
  bool IsKnown(const G4GeometryCell &gCell) const;
  const G4VPhysicalVolume &GetWorldVolume() const;
};


struct G4CellScoreValues {
  G4CellScoreValues();
  double fSumSL;
  double fSumSLW;
  double fSumSLW_v;
  double fSumSLWE;
  double fSumSLWE_v;
  int fSumTracksEntering;
  int fSumPopulation;
  int fSumCollisions;
  double fSumCollisionsWeight;
  double fNumberWeightedEnergy;
  double fFluxWeightedEnergy;
  double fAverageTrackWeight;
  double fImportance;
};


class G4CellScoreComposer {
public:
  G4CellScoreComposer();
  ~G4CellScoreComposer();
  G4CellScoreValues GetStandardCellScoreValues() const;
};


class G4CellScorer{
public:
  G4CellScoreComposer GetCellScoreComposer() const;
  G4CellScoreValues GetCellScoreValues() const;
};

typedef G4std::map<G4GeometryCell, G4CellScorer *, G4PTkComp> G4MapGeometryCellCellScorer;

class G4VCellScorerStore{
};

class G4CellScorerStore : public G4VCellScorerStore{
public:	
  G4CellScorerStore();
  ~G4CellScorerStore();
  void SetAutoScorerCreate();
  G4CellScorer *AddCellScorer(const G4GeometryCell &gCell);
  const G4MapGeometryCellCellScorer &GetMapGeometryCellCellScorer() const;    
};

class G4VScorer{
};

class G4CellStoreScorer : public G4VScorer{
public:
  G4CellStoreScorer(G4VCellScorerStore &csc);
};

class G4VImportanceAlgorithm{
};

class G4ParallelGeometrySampler{
public:
  G4ParallelGeometrySampler(G4VPhysicalVolume &worldvolume,
			    const char *);
  ~G4ParallelGeometrySampler();
  void PrepareScoring(G4VScorer *Scorer);
  void PrepareImportanceSampling(G4VIStore *istore,
				 const G4VImportanceAlgorithm *ialg = 0);
  void PrepareWeightRoulett(G4double wsurvive, 
			    G4double wlimit,
			    G4double isource);
  
  void Configure();

  void ClearSampling();
  G4bool IsConfigured() const;
};

class G4ScoreTable
{
public:
  G4ScoreTable(const G4VIStore *aIStore = 0);
  ~G4ScoreTable(){}
  void Print(const G4MapGeometryCellCellScorer &cs, 
	     G4std::ostream *out = 0);
  %addmethods {
    const char *Write(const G4MapGeometryCellCellScorer &cs){
      G4std::ostrstream tmpout;
      self->Print(cs, &tmpout);
      string *value = new string(tmpout.str());
      return value->c_str();
    };
  }
};

class B03ImportanceDetectorConstruction{
public:
  B03ImportanceDetectorConstruction();
  ~B03ImportanceDetectorConstruction();

  const G4VPhysicalVolume &GetPhysicalVolumeByName(const char *name) const;
  G4VPhysicalVolume &GetWorldVolume() const;
  %addmethods {
    const char *ListPhysNames(){
      G4String *value = new G4String(self->ListPhysNamesAsG4String());
      return value->c_str();
    }
  }
  G4GeometryCell GetGeometryCell(int i);
};

class B03AppBase{
public:
  ~B03AppBase();
  static B03AppBase &GetB03AppBase();
  G4RunManager &GetRunManager();
private:
  B03AppBase();

};



