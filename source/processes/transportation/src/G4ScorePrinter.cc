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
// $Id: G4ScorePrinter.cc,v 1.1 2002-06-13 07:12:33 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScorePrinter.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Sigma.hh"
#include "G4VIStore.hh"
#include "g4std/set"
#include "G4Pstring.hh"
#include "G4StandardTally.hh"
#include "G4ImpStandardTally.hh"
#include "G4ImportanceTally.hh"

G4ScorePrinter::G4ScorePrinter(const G4VIStore *aIStore) :
  fIStore(aIStore)
{
  FieldName = 25;
  FieldValue = 12;
  fVecCombiTally.push_back(new G4StandardTally("Histo. Ent.",
				     "HistorysEnteringWeighted",
					       "Xsum"));
  
  fVecCombiTally.push_back(new G4StandardTally("sum hi.ent*w",
				      "HistorysEnteringWeighted",
					       "SumOfWeights"));
  fVecCombiTally.push_back(new G4StandardTally("Collisions",
					       "CollisionsWeighted",
					       "Xsum"));
  fVecCombiTally.push_back(new G4StandardTally("sum col*w",
					       "CollisionsWeighted",
					       "SumOfWeights"));
  fVecCombiTally.push_back(new G4StandardTally(" AV E/ent.Tr",
				       "EnergyEnteringHistoryWeighted",
					       "Mean"));
  fVecCombiTally.push_back(new G4StandardTally(" AV E/col",
					       "CollisionEnergyWeighted",
					       "Mean"));
  fVecCombiTally.push_back(new G4StandardTally("W*SL",
					       "StepLentgh",
					       "WeightedXsum"));
  fVecCombiTally.push_back(new G4StandardTally("W*SL*E",
					       "SteplengthTimesEnergy",
					       "WeightedXsum"));

  if (fIStore) {
    fVecImpCombiTally.push_back(new G4ImportanceTally("importance"));
    fVecImpCombiTally.push_back(new G4ImpStandardTally("av w/tr rel.",
					 "WeighteOfHistorysEntering",
						 "Mean"));
    fVecImpCombiTally.push_back(new G4ImpStandardTally("av w/co rel.",
						 "WeighteOfCollisions",
						 "Mean"));
  }
}

void G4ScorePrinter::PrintHeader(G4std::ostream *out)
{
  // head line
  G4std::string vname = FillString("Volume name", ' ', FieldName+1);
  *out << vname << '|';
  for (G4VecCombiTally::iterator it = fVecCombiTally.begin();
       it != fVecCombiTally.end(); it++) {
    vname = FillString((*it)->GetName(),
		       ' ', 
		       FieldValue+1, 
		       false);
    *out << vname << '|';
  }
  for (G4VecImpCombiTally::iterator itI = fVecImpCombiTally.begin();
       itI != fVecImpCombiTally.end(); itI++) {
    vname = FillString((*itI)->GetName(),
		       ' ', 
		       FieldValue+1, 
		       false);
    *out << vname << '|';
  }
  *out << G4endl;  
  
}

G4String G4ScorePrinter::CreateName(G4PTouchableKey ptk) {
  G4String name(ptk.fVPhysiclaVolume->GetName());
  name += "_rep:" + str(ptk.fRepNum);
  return name;
}
   

void G4ScorePrinter::PrintTable(const G4PMapPtkTallys &aMapPtkTallys, 
				G4std::ostream *out) {
  // this lines sort the tallies according to the volume 
  // name they belong to
  
  G4std::map<G4String , G4TallyImp> MapStringTallies;
  for (G4PMapPtkTallys::const_iterator mit = aMapPtkTallys.begin();
       mit != aMapPtkTallys.end(); mit++) {
    G4PTouchableKey ptk = (*mit).first; // get a key identifying a volume
    G4String name(CreateName(ptk)); 
    G4double importance = 1;
    if (fIStore) {
      importance = fIStore->GetImportance(ptk);
    }
    G4PMapNameTally mtallies = (*mit).second;
    G4TallyImp tim;
    tim.tally = mtallies;
    tim.imp = importance;
    MapStringTallies[name] = tim;
  }
  // now do the printing
  for ( G4std::map<G4String , G4TallyImp>::iterator
	  it = MapStringTallies.begin();
	it != MapStringTallies.end(); it++) {
    G4String name((*it).first);
    CalculateTalliesForAVolume((*it).second.tally, 
			       (*it).second.imp); // calculate
    PrintLine(name, out); // and print line
  }
}

void G4ScorePrinter::
CalculateTalliesForAVolume(G4PMapNameTally &mtallies, 
			   G4double importance){
  for (G4PMapNameTally::iterator mt = mtallies.begin();
       mt != mtallies.end(); mt++) {
    G4String rawtallyname((*mt).first);
    for (G4VecCombiTally::iterator it = fVecCombiTally.begin();
	 it != fVecCombiTally.end(); it++) {
      (*it)->Tally(rawtallyname, (*mt).second);
    }
    for (G4VecImpCombiTally::iterator itI = fVecImpCombiTally.begin();
	 itI != fVecImpCombiTally.end(); itI++) {
      (*itI)->Tally(rawtallyname, (*mt).second, importance);
    }
  } 
}

void G4ScorePrinter::PrintLine(G4String &name,
			       G4std::ostream *out) 
{
  G4std::string fname = FillString(name, '.', FieldName);
  *out << fname << " |";
  for (G4VecCombiTally::iterator it = fVecCombiTally.begin();
       it != fVecCombiTally.end(); it++) {
    *out << G4std::setw(FieldValue) << (*it)->GetValue() << " |"; 
    (*it)->Reset();
  }
  for (G4VecImpCombiTally::iterator itI = fVecImpCombiTally.begin();
       itI != fVecImpCombiTally.end(); itI++) {
    *out << G4std::setw(FieldValue) << (*itI)->GetValue() << " |"; 
    (*itI)->Reset();
  }
   

  *out << G4endl;
  if (name == "imp_world_phys") G4cout << G4endl;

}


G4std::string G4ScorePrinter::FillString(const G4std::string &name, 
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
