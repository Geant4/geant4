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
// Guy Barrand 25th September 2021

#include "G4Plotter.hh"

G4Plotter::G4Plotter() = default;

G4Plotter::G4Plotter (const G4Plotter& a_from) = default;

G4Plotter& G4Plotter::operator=(const G4Plotter& a_from) = default;

void G4Plotter::SetLayout(unsigned int a_cols,unsigned int a_rows) {
  fColumns = a_cols;
  fRows = a_rows;
}
void G4Plotter::AddStyle(const G4String& a_style) {
  fStyles.push_back(a_style);
}
void G4Plotter::AddRegionStyle(unsigned int a_region,const G4String& a_style) {
  fRegionStyles.emplace_back(a_region,a_style);
}
void G4Plotter::AddRegionParameter(unsigned int a_region,const G4String& a_parameter,const G4String& a_value) {
  fRegionParameters.emplace_back(a_region,Parameter(a_parameter,a_value));
}

void G4Plotter::AddRegionHistogram(unsigned int a_region,tools::histo::h1d* a_h) {
  fRegion_h1ds.emplace_back(a_region,a_h);
}
void G4Plotter::AddRegionHistogram(unsigned int a_region,tools::histo::h2d* a_h) {
  fRegion_h2ds.emplace_back(a_region,a_h);
}
void G4Plotter::AddRegionH1(unsigned int a_region,int a_id) {
  fRegion_h1s.emplace_back(a_region,a_id);
}
void G4Plotter::AddRegionH2(unsigned int a_region,int a_id) {
  fRegion_h2s.emplace_back(a_region,a_id);
}

void G4Plotter::Reset() {
  fColumns = 1;
  fRows = 1;
  fStyles.clear();
  fRegionStyles.clear();
  fRegionParameters.clear();
  fRegion_h1ds.clear();
  fRegion_h2ds.clear();
  fRegion_h1s.clear();
  fRegion_h2s.clear();
}
void G4Plotter::Clear() {
  fRegion_h1ds.clear();
  fRegion_h2ds.clear();
  fRegion_h1s.clear();
  fRegion_h2s.clear();
}
void G4Plotter::ClearRegion(unsigned int a_region) {
 {std::vector<Region_h1d>::iterator it;
  for(it=fRegion_h1ds.begin();it!=fRegion_h1ds.end();) {
    if((*it).first==a_region) {
      it = fRegion_h1ds.erase(it);
    } else {
      ++it;
    }
  }}
 {std::vector<Region_h2d>::iterator it;
  for(it=fRegion_h2ds.begin();it!=fRegion_h2ds.end();) {
    if((*it).first==a_region) {
      it = fRegion_h2ds.erase(it);
    } else {
      ++it;
    }
  }}

 {std::vector<Region_h1>::iterator it;
  for(it=fRegion_h1s.begin();it!=fRegion_h1s.end();) {
    if((*it).first==a_region) {
      it = fRegion_h1s.erase(it);
    } else {
      ++it;
    }
  }}
 {std::vector<Region_h2>::iterator it;
  for(it=fRegion_h2s.begin();it!=fRegion_h2s.end();) {
    if((*it).first==a_region) {
      it = fRegion_h2s.erase(it);
    } else {
      ++it;
    }
  }}
}
