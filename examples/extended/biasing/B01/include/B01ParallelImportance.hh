#ifndef B01ParallelImportance_hh
#define B01ParallelImportance_hh B01ParallelImportance_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class B01Run;
class G4IStore;
class G4ParallelImportanceSampler;
class B01ParallelGeometry;


class B01ParallelImportance : public B01VSimulation {
public:
  B01ParallelImportance();
  ~B01ParallelImportance();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01ParallelImportance(const B01ParallelImportance &);
  B01ParallelImportance &operator=
  (const B01ParallelImportance &);

  G4String fName;
  B01VGeometry *fMassGeometry;
  B01Run *fRun;
  B01ParallelGeometry *fParallelGeometry;
  G4IStore *fIStore;
  G4ParallelImportanceSampler *fSampler;
};

#endif
