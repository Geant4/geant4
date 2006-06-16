// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
//
#include "Cintex/Cintex.h"
//
#include "include/ExP01TrackerHit.hh"


int main(int argc,char** argv) 
{
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCintex");
  gSystem->Load("libClassesDict");
  //  ROOT::Cintex::Cintex::SetDebug(2);
  ROOT::Cintex::Cintex::Enable();
  if(argc<2) G4cout << "Missing name of the file to read!" << G4endl;
 
  TFile fo(argv[1]);
   
  std::vector<ExP01TrackerHit*>* hits;
  fo.GetListOfKeys()->Print();
 
  TIter next(fo.GetListOfKeys());
  TKey *key;
  double tot_en;
  while ((key=(TKey*)next()))
  {
    fo.GetObject(key->GetName(), hits);
 
    tot_en = 0;
    G4cout << "Collection: " << key->GetName() << G4endl;
    G4cout << "Number of hits: " << hits->size() << G4endl;
    for (int i=0;i!=hits->size();i++)
    {
      (*hits)[i]->Print();
    }         
  }
}


