#ifndef B01MassScoring_hh
#define B01MassScoring_hh B01MassScoring_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class B01Run;
class G4MassScoreSampler;
class G4Scorer;


class B01MassScoring : public B01VSimulation {
public:
  B01MassScoring();
  ~B01MassScoring();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01MassScoring(const B01MassScoring &);
  B01MassScoring &operator=(const B01MassScoring &);

  G4String fName;
  B01VGeometry *fGeometry;
  B01Run *fRun;
  G4Scorer *fScorer;
  G4MassScoreSampler *fSampler;
};

#endif
