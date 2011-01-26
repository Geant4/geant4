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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TCatalog.h"

G4TCatalog *gCatalog = new G4TCatalog();

ClassImp(G4TCatalog)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
vector<TString> G4TCatalog::GetPublications()
{
  vector<TString> files = this->Load();
  vector<TString> result;
  for(UInt_t i=0; i < files.size(); ++i)
  {
    TString str = files[i];
    if(str.Contains("data")) result.push_back(str);
  }
  return result;
}

//______________________________________________________________________________
vector<TString> G4TCatalog::Load()
{
  vector<TString> catalog;
  ifstream asciiFile;
  string   asciiLine;
  // open the file
  asciiFile.open(fFileName.Data());
  if (!asciiFile)
  {
    cout<< "Warning>>G4TCatalog::Load: Error in opening file="<< fFileName.Data()<< endl;
    return catalog;
  }
  // process
  while(getline(asciiFile, asciiLine)) catalog.push_back(TString(asciiLine));

  asciiFile.close();
  return catalog;
}

//______________________________________________________________________________
Bool_t G4TCatalog::ContainsLine(TString const& line)
{
  ifstream asciiFile;
  string   asciiLine;
  // open the file
  asciiFile.open(fFileName.Data());
  if (!asciiFile)
  {
    cout<<"Warning>G4TCatalog::ContainsLine: ErrorInOpening file="<<fFileName.Data()<<endl;
    return false;
  }
  // process
  while(getline(asciiFile, asciiLine))
  {
    if(line == TString(asciiLine))
    {
      asciiFile.close();
      return true; // already there
    }
  }
  asciiFile.close();
  return false;
}

//______________________________________________________________________________
void G4TCatalog::Insert(TString const& name)
{
  ofstream  asciiFile;
  asciiFile.open(fFileName.Data(), std::ios::out | std::ios::app );
  if (!asciiFile || asciiFile.fail())
  {
    cout<<"!Warning->G4TCatalog::Insert: Error in opening file="<<fFileName.Data()
        <<", file="<<asciiFile<<", fail="<<asciiFile.fail()<<endl;
    return;
  }
  if(ContainsLine(name))
  {
    cout<<"Warning>G4TCatalog::Insert:File="<<name<<" exists -> rewrite"<<endl;
    return; // This name is is already in the Catalog, so just rewriting the file itself
  }
  else cout<<"G4TCatalog::Insert: File="<<name<<" is new"<<endl;
  asciiFile << name << endl;
  asciiFile.close();
  return;
}

