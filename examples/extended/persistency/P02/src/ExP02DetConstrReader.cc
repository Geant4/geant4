//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "Cintex/Cintex.h"

//G4
#include "G4Element.hh"
#include "G4PVPlacement.hh"
#include "ExP02GeoTree.hh"

// local
#include "ExP02DetConstrReader.hh"

ExP02DetConstrReader::ExP02DetConstrReader()
{  
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libClassesDict");

  //  ROOT::Cintex::Cintex::SetDebug(2);
  ROOT::Cintex::Cintex::Enable();
  //  gDebug = 1;

}

ExP02DetConstrReader::~ExP02DetConstrReader()
{}

G4VPhysicalVolume* ExP02DetConstrReader::Construct()
{
  ExP02GeoTree* geotree;

  TFile fo("geo.root");
  fo.GetObject("my_geo", geotree);

  return geotree->TopVol();
}
