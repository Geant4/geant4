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
// $Id: B06ScorePrinter.cc,v 1.5 2002/04/19 10:54:33 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#include "B06ScorePrinter.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Sigma.hh"

B06ScorePrinter::B06ScorePrinter()
{
  FieldName = 25;
  FieldValue = 12;
}

void B06ScorePrinter::PrintHeader(G4std::ostream *out)
{
  // head line
  G4std::string vname = FillString("Volume name", ' ', FieldName+1);
  *out << vname << '|';
  vname = FillString("AV w/Ent.Hi.", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("    sigma ", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString(" AV w/Coll ", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("    sigma ", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("Coll_Ent.Tr", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString(" AV E/Track ", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString(" sigma", ' ', FieldValue+1);
  *out << vname << '|';
  *out << G4endl;
}

void B06ScorePrinter::PrintTable(const G4PMapPtkTallys &aMapPtkTallys, 
                                 G4std::ostream *out)
{
  for (G4PMapPtkTallys::const_iterator mit = aMapPtkTallys.begin();
       mit != aMapPtkTallys.end(); mit++) {
    G4PTouchableKey ptk = (*mit).first; // get a key identifying a volume
    G4PMapNameTally mtallies = (*mit).second; // get tallies of the volume
    G4String name(ptk.fVPhysiclaVolume->GetName()); // print volume name
    
    G4double avW_EntHist = -1;
    G4double avW_EntHist_sigma = -1;
    G4double avW_Coll = -1;
    G4double avW_Coll_sigma = -1;

    G4double sumTrWeight = -1;
    G4double sumCollWeight = -1;
    G4double Coll_Ent_Tr = -1;

    G4double meanTrackEnergy = -1, sigmaTrackEnergy = -1;

    for (G4PMapNameTally::iterator mt = mtallies.begin();
	 mt != mtallies.end(); mt++) {
      G4String tmp((*mt).first);
      if (tmp == "WeighteOfHistorysEntering") {
	 avW_EntHist = (*mt).second.GetMean();
	 avW_EntHist_sigma = (*mt).second.GetSigma();
      }
      if (tmp == "WeighteOfCollisions") {
	 avW_Coll = (*mt).second.GetMean();
	 avW_Coll_sigma = (*mt).second.GetSigma();
      }
      if (tmp == "HistorysEnteringWeighted") {
	sumTrWeight = (*mt).second.GetSumOfWeights();
      }
      if (tmp == "CollisionsWeighted") {
	sumCollWeight = (*mt).second.GetSumOfWeights();
      }
      if (tmp == "EnergyEnteringHistoryWeighted") {
	meanTrackEnergy =  (*mt).second.GetMean();
	sigmaTrackEnergy = (*mt).second.GetSigma();
      }
    }

    if (sumTrWeight!=0.&&sumCollWeight!=-1) {
      Coll_Ent_Tr = sumCollWeight / sumTrWeight;
    }
    // print values

    G4std::string fname = FillString(name, '.', FieldName);
    *out << fname << " |";
    *out << G4std::setw(FieldValue) << avW_EntHist << " |"; 
    *out << G4std::setw(FieldValue) << avW_EntHist_sigma << " |";
    *out << G4std::setw(FieldValue) << avW_Coll << " |";
    *out << G4std::setw(FieldValue) << avW_Coll_sigma << " |";
    *out << G4std::setw(FieldValue) << Coll_Ent_Tr << " |";
    *out << G4std::setw(FieldValue) << meanTrackEnergy << " |";
    *out << G4std::setw(FieldValue) << sigmaTrackEnergy << " |";
    *out << G4endl;
  }
}

G4std::string B06ScorePrinter::FillString(const G4std::string &name, 
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
