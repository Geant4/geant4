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
//
// $Id: G4StandardScoreTable.cc,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4StandardScoreTable.cc
//
// ----------------------------------------------------------------------

#include "G4StandardScoreTable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VIStore.hh"
#include "g4std/set"
#include "G4Pstring.hh"


G4StandardScoreTable::G4StandardScoreTable(const G4VIStore *aIStore) :
  fIStore(aIStore)
{
  FieldName = 25;
  FieldValue = 12;
}

void G4StandardScoreTable::
Print(const G4MapPtkStandardCellScorer &mptkscorer,
      G4std::ostream *out){
  PrintHeader(out);
  PrintTable(mptkscorer, out);
}

void G4StandardScoreTable::PrintHeader(G4std::ostream *out)
{
  G4std::vector<G4String> vecScoreName;
  vecScoreName.push_back("Tr.Entering");
  vecScoreName.push_back("Population");
  vecScoreName.push_back("Collisions");
  vecScoreName.push_back("Coll*WGT");
  vecScoreName.push_back("NumWGTedE");
  vecScoreName.push_back("FluxWGTedE");
  vecScoreName.push_back("Av.Tr.WGT");

  // head line
  G4std::string vname = FillString("Volume name", ' ', FieldName+1);
  *out << vname << '|';
  for (G4std::vector<G4String>::iterator it = vecScoreName.begin();
       it != vecScoreName.end(); it++) {
    vname = FillString((*it),
		       ' ', 
		       FieldValue+1, 
		       false);
    *out << vname << '|';
  }
  *out << G4endl;  
}

G4String G4StandardScoreTable::CreateName(G4PTouchableKey ptk) {
  G4String name(ptk.fVPhysiclaVolume->GetName());
  name += "_rep:" + str(ptk.fRepNum);
  return name;
}
   

void G4StandardScoreTable::PrintTable(const G4MapPtkStandardCellScorer 
				      &mptkscorer,
				      G4std::ostream *out) {
  
  // this lines sort the ScoreValues according to the volume 
  // name they belong to
  
  G4std::map<G4String , G4StandardCellScorer> MapStringSCScorer;
  for (G4MapPtkStandardCellScorer::const_iterator mit = 
	 mptkscorer.begin();
       mit != mptkscorer.end(); mit++) {
    G4PTouchableKey ptk = (*mit).first; // get a key identifying a volume
    G4String name(CreateName(ptk)); 

    G4double importance = 1;
    if (fIStore) {
      importance = fIStore->GetImportance(ptk);
    }
    MapStringSCScorer[name] = (*mit).second;
    MapStringSCScorer[name].SetImportnace(importance);
  }
  // now do the printing
  for ( G4std::map<G4String , G4StandardCellScorer>::iterator
	  it = MapStringSCScorer.begin();
	it != MapStringSCScorer.end(); it++) {
    G4String name((*it).first);
    PrintLine(name, (*it).second.GetStandardCellScoreValues(), out); 
  }
}

void G4StandardScoreTable::PrintLine(G4String &name,
				     G4StandardCellScoreValues sc_scores,
				     G4std::ostream *out) 
{
  G4std::string fname = FillString(name, '.', FieldName);
  *out << fname << " |";

  *out << G4std::setw(FieldValue) << sc_scores.fSumTracksEntering 
       << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fSumPopulation << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fSumCollisions << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fSumCollisionsWeight 
       << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fNumberWeightedEnergy 
       << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fFluxWeightedEnergy 
       << " |"; 
  *out << G4std::setw(FieldValue) << sc_scores.fAverageTrackWeight*
    sc_scores.fImportance
       << " |"; 
  
  *out << G4endl;
}


G4std::string G4StandardScoreTable::FillString(const G4std::string &name, 
                                          char c, G4int n, G4bool back)
{
  G4std::string fname;
  G4int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += G4std::string(k,c);
    }
    else {
      fname = G4std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}









