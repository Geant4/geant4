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
//
/// \file persistency/P01/readHits.cc
/// \brief Main program of the persistency/P01 example
//
// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
//
#include "include/ExP01TrackerHit.hh"

int main(int argc,char** argv) 
{
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libExP01ClassesDict");
  if(argc<2) G4cout << "Missing name of the file to read!" << G4endl;
 
  TFile fo(argv[1]);
   
  std::vector<ExP01TrackerHit*>* hits;
  fo.GetListOfKeys()->Print();
 
  TIter next(fo.GetListOfKeys());
  TKey *key;
  //double tot_en;
  while ((key=(TKey*)next()))
  {
    fo.GetObject(key->GetName(), hits);
 
    //tot_en = 0;
    G4cout << "Collection: " << key->GetName() << G4endl;
    G4cout << "Number of hits: " << hits->size() << G4endl;
    for (size_t i=0;i!=hits->size();i++)
    {
      (*hits)[i]->Print();
    }         
  }
}

