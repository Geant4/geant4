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
// $Id: G4ScoreTable.cc,v 1.2 2002-11-04 10:52:39 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ScoreTable.cc
//
// ----------------------------------------------------------------------

#include "G4ScoreTable.hh"
#include "g4std/strstream"

#include "G4VPhysicalVolume.hh"
#include "G4VIStore.hh"
#include "g4std/set"
#include "G4CellScorer.hh"


G4ScoreTable::G4ScoreTable(const G4VIStore *aIStore) :
  fIStore(aIStore),
  FieldName(25),
  FieldValue(12)
{}

G4ScoreTable::~G4ScoreTable()
{}
void G4ScoreTable::
Print(const G4MapGeometryCellCellScorer &cs,
      G4std::ostream *out){
  if (!out) {
    out = &G4cout;
  }
  PrintHeader(out);
  PrintTable(cs, out);
}

void G4ScoreTable::PrintHeader(G4std::ostream *out)
{
  G4std::vector<G4String> vecScoreName;
  vecScoreName.push_back("Importance");
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

G4String G4ScoreTable::CreateName(const G4GeometryCell &gCell) {
  
  char st[200];
  G4std::ostrstream os(st,200);
  os << gCell.GetPhysicalVolume().GetName()
     << "_rep:" << gCell.GetReplicaNumber()
     << '\0';
  G4String name(st);

  return name;
}
   

void G4ScoreTable::PrintTable(const G4MapGeometryCellCellScorer &mcs,
				      G4std::ostream *out) {
  
  // this lines sort the ScoreValues according to the volume 
  // name they belong to
  G4std::map<G4String , G4CellScoreComposer> MapStringSCScorer;
  for (G4MapGeometryCellCellScorer::const_iterator mit = 
	 mcs.begin();
       mit != mcs.end(); ++mit) {
    G4GeometryCell gCell = (*mit).first; // get a key identifying a volume
    G4String name(CreateName(gCell)); 

    G4double importance = 1;
    if (fIStore) {
      if (fIStore->IsKnown(gCell)) {
	importance = fIStore->GetImportance(gCell);
      }
    }
    MapStringSCScorer[name] = (*mit).second->GetCellScoreComposer();
    MapStringSCScorer[name].SetImportnace(importance);
  }
  // now do the printing
  for ( G4std::map<G4String , G4CellScoreComposer>::iterator
	  it = MapStringSCScorer.begin();
	it != MapStringSCScorer.end(); ++it) {
    G4String name((*it).first);
    PrintLine(name, (*it).second.GetStandardCellScoreValues(), out); 
  }
  *out << '\0';
  out->flush();
}

void G4ScoreTable::PrintLine(const G4String &name,
			     const G4CellScoreValues &sc_scores,
			     G4std::ostream *out) 
{
  G4std::string fname = FillString(name, '.', FieldName);
  *out << fname << " |";
  *out << G4std::setw(FieldValue) << sc_scores.fImportance 
       << " |";
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


G4std::string G4ScoreTable::FillString(const G4std::string &name, 
				       char c, G4int n, G4bool back)
{
  G4std::string fname("");
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









