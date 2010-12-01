//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

void RootIO::Close()
{
  fo->Close();
}
