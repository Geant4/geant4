#ifndef Tst33ScorerBuilder_hh
#define Tst33ScorerBuilder_hh Tst33ScorerBuilder_hh

class G4CellScorerStore;
class Tst33VGeometry;
class G4CellScorer;

class Tst33ScorerBuilder {
public:
  Tst33ScorerBuilder();
  ~Tst33ScorerBuilder();
  G4CellScorerStore *CreateScorer(Tst33VGeometry *samplegeo,
				  const G4CellScorer **specialCellScorer);  
};


#endif
