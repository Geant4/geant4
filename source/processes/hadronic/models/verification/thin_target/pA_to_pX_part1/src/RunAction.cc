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
// $Id: RunAction.cc,v 1.3 2003-06-19 14:44:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "RunAction.hh"
#include "TargetConstruction.hh"
#include "EnergyAngleCrossSection.hh"

#include <iomanip>
#include "G4Run.hh"
#include "G4Material.hh"


RunAction::RunAction(TargetConstruction* targ)
  : loBinEdge(0.), hiBinEdge(800.), binSize(10.), 
    ehists(loBinEdge,hiBinEdge,binSize), 
    target(targ), theData(0), material(0), twopi(6.283185)
{
  // Create histograms
  ehists.CreateHists();
}


RunAction::~RunAction()
{
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  material =  target->GetTargetMaterial();
  const G4String& matName =  material->GetName();  

  // Load the data

  if (matName == "Carbon") {
    theData = new EnergyAngleCrossSection("data/PPprime800_C.dat");

  } else if (matName == "Calcium") {
    theData = new EnergyAngleCrossSection("data/PPprime800_Ca.dat");

  } else if (matName == "Lead") {
    theData = new EnergyAngleCrossSection("data/PPprime800_Pb.dat");

  } else {
    G4Exception("No such material");
  }

  G4cout << matName << " data loaded " << G4endl;

  // Set up array of 4 degree wide angle bins

  G4double deg_to_rad = twopi/360.;
  std::vector<G4double>& angles = theData->GetAngles();
  binEdges.clear();
  std::pair<G4double, G4double> bin_edge_pair;

  for (G4int i = 0; i < (G4int)angles.size(); i++) {
    bin_edge_pair.first = cos(deg_to_rad*(angles[i] - 2.0));
    bin_edge_pair.second = cos(deg_to_rad*(angles[i] + 2.0));
    binEdges.push_back(bin_edge_pair);
  }
}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double thick = target->GetTargetThickness()/cm;
  G4double Natoms = material->GetTotNbOfAtomsPerVolume()*cm3;
  G4double Nevents = aRun->GetNumberOfEvent();
  G4double scale_factor = Nevents*Natoms*thick*binSize*1e-27;
  G4cout << " Natoms = " << Natoms << ", Nevents = " << Nevents << ", scale factor = " << scale_factor << G4endl;
  for (G4int i = 0; i < (G4int)binEdges.size(); i++) {
    G4double dOmega = twopi*abs(binEdges[i].first - binEdges[i].second);
    G4double Weight = 1./(scale_factor*dOmega);
    ehists.ScaleHists(Weight,i);
  }
}


void RunAction::FillHists(G4double charge, G4double ke, G4double cost)
{
  std::vector<std::pair<G4double,G4double> > bin_edges; 

  for (G4int i = 0; i < (G4int)binEdges.size(); i++) {
    if(cost < binEdges[i].first && cost > binEdges[i].second ) {
      ehists.FillHists(charge, ke, 1.0, i);
      /*
      if (charge > 0 ) {
        G4cout <<"**** Spectrum " << i << G4endl;
        G4cout <<"solid angle = " << setprecision(6) << dOmega << G4endl;
        G4cout << " cost: " << binEdges[i].first << " > " 
	     << cost << " > " << binEdges[i].second << G4endl;
        G4cout << " KE = " << ke << G4endl;
      }
      */
    }
  }
}




