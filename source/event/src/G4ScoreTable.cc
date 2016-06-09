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
// $Id: G4ScoreTable.cc,v 1.3 2003/06/16 16:50:36 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ScoreTable.cc
//
// ----------------------------------------------------------------------

#include "G4ScoreTable.hh"
#include <strstream>

#include "G4VPhysicalVolume.hh"
#include "G4VIStore.hh"
#include <set>
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
      std::ostream *out){
  if (!out) {
    out = &G4cout;
  }
  PrintHeader(out);
  PrintTable(cs, out);
}

void G4ScoreTable::PrintHeader(std::ostream *out)
{
  std::vector<G4String> vecScoreName;
  vecScoreName.push_back("Importance");
  vecScoreName.push_back("Tr.Entering");
  vecScoreName.push_back("Population");
  vecScoreName.push_back("Collisions");
  vecScoreName.push_back("Coll*WGT");
  vecScoreName.push_back("NumWGTedE");
  vecScoreName.push_back("FluxWGTedE");
  vecScoreName.push_back("Av.Tr.WGT");

  // head line
  std::string vname = FillString("Volume name", ' ', FieldName+1);
  *out << vname << '|';
  for (std::vector<G4String>::iterator it = vecScoreName.begin();
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
  std::ostrstream os(st,200);
  os << gCell.GetPhysicalVolume().GetName()
     << "_rep:" << gCell.GetReplicaNumber()
     << '\0';
  G4String name(st);

  return name;
}
   

void G4ScoreTable::PrintTable(const G4MapGeometryCellCellScorer &mcs,
				      std::ostream *out) {
  
  // this lines sort the ScoreValues according to the volume 
  // name they belong to
  std::map<G4String , G4CellScoreComposer> MapStringSCScorer;
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
  for ( std::map<G4String , G4CellScoreComposer>::iterator
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
			     std::ostream *out) 
{
  std::string fname = FillString(name, '.', FieldName);
  *out << fname << " |";
  *out << std::setw(FieldValue) << sc_scores.fImportance 
       << " |";
  *out << std::setw(FieldValue) << sc_scores.fSumTracksEntering 
       << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fSumPopulation << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fSumCollisions << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fSumCollisionsWeight 
       << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fNumberWeightedEnergy 
       << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fFluxWeightedEnergy 
       << " |"; 
  *out << std::setw(FieldValue) << sc_scores.fAverageTrackWeight*
    sc_scores.fImportance
       << " |"; 
  
  *out << G4endl;
}


std::string G4ScoreTable::FillString(const std::string &name, 
				       char c, G4int n, G4bool back)
{
  std::string fname("");
  G4int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += std::string(k,c);
    }
    else {
      fname = std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}









