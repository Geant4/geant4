// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
//*****************************************************************************
// To run this macro in cint do (after replacing the location_of_your_libraries below):
// .include  $G4INCLUDE
// gSystem->Load("libCintex.so");
// gSystem->Load("<location_of_your_libraries>/libClassesDict.so");
// ROOT::Cintex::Cintex::Enable();
// .L hits.C++
// hits();
//*****************************************************************************
#include "include/ExP01TrackerHit.hh"

void hits()
{
  TFile fo("hits.root");
   
  std::vector<ExP01TrackerHit*>* hits;
  fo.GetListOfKeys()->Print();
 
  TIter next(fo.GetListOfKeys());
  TKey *key;
  double tot_en;
  while ((key=(TKey*)next()))
  {
    fo.GetObject(key->GetName(), hits);
 
    tot_en = 0;
    cout << "Collection: " << key->GetName() << endl;
    cout << "Number of hits: " << hits->size() << endl;
    for (int i=0;i!=hits->size();i++)
    {
      (*hits)[i]->Print();
    }         
  }
}