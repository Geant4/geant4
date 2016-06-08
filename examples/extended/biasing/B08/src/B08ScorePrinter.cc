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
// $Id: B08ScorePrinter.cc,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "B08ScorePrinter.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Sigma.hh"
#include "G4VIStore.hh"

B08ScorePrinter::B08ScorePrinter(const G4VIStore *aIStore) :
  fIStore(aIStore)
{
  FieldName = 25;
  FieldValue = 12;
}

void B08ScorePrinter::PrintHeader(G4std::ostream *out)
{
  // head line
  G4std::string vname = FillString("Volume name", ' ', FieldName+1);
  *out << vname << '|';
  vname = FillString("Histo. Ent.", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("sum hi.ent*w", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("Collisions", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("sum col*w", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString(" AV E/ent.Tr", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString(" AV E/col", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("W*SL", ' ', FieldValue+1, false);
  *out << vname << '|';
  vname = FillString("W*SL*E", ' ', FieldValue+1, false);
  *out << vname << '|';
  if (fIStore) {
    vname = FillString(" importance ", ' ', FieldValue+1, false);
    *out << vname << '|';
    vname = FillString("av w/tr rel.", ' ', FieldValue+1, false);
    *out << vname << '|';
    vname = FillString("av w/co rel.", ' ', FieldValue+1, false);
    *out << vname << '|';
  }
  *out << G4endl;
}

void B08ScorePrinter::PrintTable(const G4PMapPtkTallys &aMapPtkTallys, 
                                 G4std::ostream *out)
{
  for (G4PMapPtkTallys::const_iterator mit = aMapPtkTallys.begin();
       mit != aMapPtkTallys.end(); mit++) {
    G4PTouchableKey ptk = (*mit).first; // get a key identifying a volume
    G4PMapNameTally mtallies = (*mit).second; // get tallies of the volume
    G4String name(ptk.fVPhysiclaVolume->GetName()); // print volume name
    G4double importance = 1;
    if (fIStore) {
      importance = fIStore->GetImportance(ptk);
    }
    G4double histoEntering = 0;
    G4double collSum = 0;
    G4double collwsum = 0;
    G4double hinsum = 0;
    G4double collEnergyWeightedSum = 0;
    G4double WeightEnergy = 0;
    G4double entTrackEnergyWeightedSum = 0;
    G4double wsl = 0.;
    G4double wsle = 0.;
    G4double wsl_tr = 0.;
    G4double wsle_tr = 0.;
    G4double avWE_Coll = 0.;
    

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
	hinsum = (*mt).second.GetEntries();
	avW_EntHist = (*mt).second.GetMean();
	avW_EntHist_sigma = (*mt).second.GetSigma();
      }
      if (tmp == "WeighteOfCollisions") {
	avW_Coll = (*mt).second.GetMean();
	avW_Coll_sigma = (*mt).second.GetSigma();
      }
      if (tmp == "HistorysEnteringWeighted") {
	sumTrWeight = (*mt).second.GetSumOfWeights();
	histoEntering = (*mt).second.GetXsum();
      }
      if (tmp == "CollisionsWeighted") {
	sumCollWeight = (*mt).second.GetSumOfWeights();
	collSum = (*mt).second.GetXsum();
        collwsum = (*mt).second.GetWeightedXsum();
      }
      if (tmp == "EnergyEnteringHistoryWeighted") {
	meanTrackEnergy =  (*mt).second.GetMean();
	sigmaTrackEnergy = (*mt).second.GetSigma();
	entTrackEnergyWeightedSum = (*mt).second.GetWeightedXsum();
      }

      if (tmp == "CollisionEnergyWeighted") {
	collEnergyWeightedSum = (*mt).second.GetWeightedXsum();
	avWE_Coll = (*mt).second.GetMean();
       }
      
      if (tmp == "StepLentgh") {
	wsl = (*mt).second.GetWeightedXsum();
      }
      if (tmp == "SteplengthTimesEnergy") {
	wsle = (*mt).second.GetWeightedXsum();
      }

    }
    
    
    G4double sumtrcoll = sumTrWeight + sumCollWeight;
    if (sumtrcoll!=0) {
      WeightEnergy = (collEnergyWeightedSum + entTrackEnergyWeightedSum) /
	(sumtrcoll);
      wsl_tr = wsl / sumtrcoll;
      wsle_tr = wsle / sumtrcoll;
    }
    
    if (sumTrWeight!=0.&&sumCollWeight!=-1) {
      Coll_Ent_Tr = sumCollWeight / sumTrWeight;
    }
    // print values

    G4std::string fname = FillString(name, '.', FieldName);
    *out << fname << " |";
    *out << G4std::setw(FieldValue) << histoEntering << " |"; 
    *out << G4std::setw(FieldValue) << sumTrWeight << " |"; 
    *out << G4std::setw(FieldValue) << collSum << " |";
    if (hinsum==0) hinsum = 1;
    *out << G4std::setw(FieldValue) << collwsum << " |"; 
    *out << G4std::setw(FieldValue) << meanTrackEnergy << " |";
    *out << G4std::setw(FieldValue) << avWE_Coll << " |";
    *out << G4std::setw(FieldValue) << wsl << " |";
    *out << G4std::setw(FieldValue) << wsle << " |";
    if (fIStore) {
      *out << G4std::setw(FieldValue) << importance << " |"; 
      *out << G4std::setw(FieldValue) << avW_EntHist*importance << " |"; 
      *out << G4std::setw(FieldValue) << avW_Coll*importance << " |";
    }
    *out << G4endl;
    if (name == "imp_world_phys") G4cout << G4endl;
  }
}

G4std::string B08ScorePrinter::FillString(const G4std::string &name, 
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
