#include <sstream>

#include "RootIO.hh"
//
#include "Cintex/Cintex.h"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
//

static RootIO* instance = 0;

RootIO::RootIO():Nevents(0)
{
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libClassesDict");

  ROOT::Cintex::Cintex::SetDebug(2);
  ROOT::Cintex::Cintex::Enable();
  gDebug = 1;

  fo = new TFile("hits.root","RECREATE");
}

RootIO::~RootIO()
{}

RootIO* RootIO::GetInstance()
{
  if (instance == 0 )
  {
    instance = new RootIO();
  }
  return instance;
}

void RootIO::Write(std::vector<ExP01TrackerHit*>* hcont)
{
  Nevents++;

  std::ostringstream os;
  os << Nevents;
  std::string stevt = "Event_" + os.str(); 
  const char* chevt = stevt.c_str();

  std::cout << "writing " << stevt << std::endl;


  fo->WriteObject(hcont, chevt);

}
