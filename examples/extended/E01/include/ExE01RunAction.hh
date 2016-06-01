// $Id: ExE01RunAction.hh,v 1.1 1998/10/14 15:19:39 allison Exp $

#ifndef ExE01RunAction_h
#define ExE01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "HepODBMS/clustering/HepDbApplication.h"
#include "histograms.h"
class G4Run;

class ExE01RunAction : public G4UserRunAction
{
  public:
    ExE01RunAction() { runIDcounter = 0; }
    ~ExE01RunAction() {}

    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

    static HepRef(Histo1D) get_1d() { return myHisto1D; }
    static HepRef(Histo2D) Get2d() { return myHisto2D; }

  private:
    G4int runIDcounter;

    static HepRef(Histo1D) myHisto1D;
    static HepRef(Histo2D) myHisto2D;

    HepDbApplication dbApp;
};

#endif

