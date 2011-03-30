#ifndef analysis_h
#define analysis_h 1

#ifdef G4ANALYSIS_USE

class G4Run;
class G4Event;
class G4Step;

class TFile;
class TNtuple;
class TH1F;

class analysis {
public:
  analysis();
  virtual ~analysis();
public:
  virtual void BeginOfRun(const G4Run*); 
  virtual void EndOfRun(const G4Run*); 
  virtual void BeginOfEvent(const G4Event*); 
  virtual void EndOfEvent(const G4Event*); 
  virtual void Step(const G4Step*);

private:
  int runId;

  TFile* hfile;
  TNtuple* ntuple;
  TH1F* histoTrajectories;
  TH1F* histoE;
};

#else

class analysis;

#endif

#endif
