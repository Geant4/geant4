// $Id: ExE01RunAction.hh,v 1.3 1999/04/22 21:49:17 asaim Exp $

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
    ExE01RunAction() { }
    virtual ~ExE01RunAction() {}

    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

    static HepRef(Histo1D) get_1d() { return myHisto1D; }
    static HepRef(Histo2D) Get2d() { return myHisto2D; }

  private:
    static HepRef(Histo1D) myHisto1D;
    static HepRef(Histo2D) myHisto2D;

    HepDbApplication dbApp;
};

#endif

