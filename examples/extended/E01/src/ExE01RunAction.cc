// $Id: ExE01RunAction.cc,v 1.1 1998/10/14 15:20:12 allison Exp $

#include "ExE01RunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

HepRef(Histo1D) ExE01RunAction::myHisto1D;
HepRef(Histo2D) ExE01RunAction::myHisto2D;

void ExE01RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/event/Verbose 0");
  UI->ApplyCommand("/tracking/Verbose 0");

  dbApp.Init();
  dbApp.startUpdate();

  HepDatabaseRef dbH = dbApp.db("ExE01DB");

  // Number of histogram bins  
  int nBins1d = 100;
  int nBins2d = 40;

  // Limits for values
  double  low = 0,  high = 1;
  double xlow = 0, xhigh = 1, ylow = 0, yhigh = 1;

  myHisto1D = new(dbH) 
    Histo1D("HepJamesRandom1", nBins1d, low, high);

  myHisto2D = new(dbH) 
    Histo2D("HepJamesVsDRand48", nBins2d, xlow, xhigh, nBins2d, ylow, yhigh);
}

void ExE01RunAction::EndOfRunAction(G4Run* aRun)
{
  // Print histograms
  HistPrintout p(G4cout);

  p.print(*myHisto1D);
  p.print(*myHisto2D);

  // commit changes to Objectivty database
  dbApp.commit();
}


