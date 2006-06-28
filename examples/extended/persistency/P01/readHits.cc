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


