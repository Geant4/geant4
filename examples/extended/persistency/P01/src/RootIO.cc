//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
