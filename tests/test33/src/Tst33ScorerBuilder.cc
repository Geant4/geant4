#include "Tst33ScorerBuilder.hh"
#include "G4CellScorerStore.hh"
#include "G4CellScorer.hh"
#include "Tst33VGeometry.hh"

Tst33ScorerBuilder::Tst33ScorerBuilder()
{}

Tst33ScorerBuilder::~Tst33ScorerBuilder()
{}


G4CellScorerStore *Tst33ScorerBuilder::
CreateScorer(Tst33VGeometry *samplegeo, 
	     const G4CellScorer **specialCellScorer){
  G4GeometryCell gWorldCell(samplegeo->GetWorldVolume(), -1);
  
  G4CellScorerStore *cs_store = new G4CellScorerStore();
  if (!cs_store) {
    G4std::G4Exception("Tst33ScorerBuilder::CreateScorer: new failed to create G4CellScorerStore!");
  }
  cs_store->AddCellScorer(gWorldCell);
  
  
  G4int i = 1;
  for (i=1; i <= 19; i++) {
    G4String volname = samplegeo->GetCellName(i);
    G4double imp = pow(2,i-1);
    if (i==19) {
      imp = pow(2,17);
    }
    
    const G4VPhysicalVolume *pvol = samplegeo->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      const G4CellScorer *s = cs_store->AddCellScorer(gCell);
      if (i==18) {
	*specialCellScorer = s;
      }
    }
  }
  return cs_store;
}
