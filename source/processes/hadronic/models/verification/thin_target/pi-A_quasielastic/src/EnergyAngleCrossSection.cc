//
////////////////////////////////////////////////////////////////////
//
//  EnergyAngleCrossSection
//  Author: D.H. Wright (SLAC)
//  Date: 5 March 2003
//
///////////////////////////////////////////////////////////////////
//
// $Id: EnergyAngleCrossSection.cc,v 1.1 2003-07-31 01:21:00 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "EnergyAngleCrossSection.hh"

EnergyAngleCrossSection::EnergyAngleCrossSection()
{}


EnergyAngleCrossSection::EnergyAngleCrossSection(std::ifstream datafile)
{
  if (!datafile.good()) G4Exception(" Cannot find file ");

  // Skip header lines

  G4int tag = -1;
  char fline[128];
  while ((tag == -1) && !datafile.eof()) {
    datafile.getline(fline,128);
    std::string sline(fline);
    tag = sline.find("startdata");
  }
  if (tag == -1) G4Exception(" Failed to find startdata tag ");

  // Read the data 

  G4int i = 0;
  while (!datafile.eof()) {
    datafile >> ddcs[i][0] >> ddcs[i][1] >> ddcs[i][2] >> ddcs[i][3] ;
    i++;
  }

  nDataPts = i-1;
  G4cout << nDataPts << " data points read in " << G4endl;
}


EnergyAngleCrossSection::~EnergyAngleCrossSection()
{}


std::vector<G4double>&
EnergyAngleCrossSection::GetEnergies()
{
  Energies.clear();
  G4int i;
  for (i=0; i < nDataPts; i++){
    Energies.push_back(ddcs[i][0]);
  }
  return Energies;
}

 
std::vector<G4double>&
EnergyAngleCrossSection::GetAngles()
{
  Angles.clear();
  G4double ang;
  G4bool newangle;

  for (G4int i = 0; i < nDataPts; i++){
    ang = ddcs[i][1];
    newangle = true;
    for (G4int j = 0; j < (G4int)Angles.size(); j++) {
      if (ang == Angles[j]) {
        newangle = false;
        break;
      }
    }
    if (newangle) Angles.push_back(ang);
  }
  return Angles;
}
 
// GetEnergySpectrum without argument

std::vector<std::pair<G4double, G4double> >& 
EnergyAngleCrossSection::GetEnergySpectrum()
{
  EnergySpectrum.clear();
  G4int i;
  for (i=0; i < nDataPts; i++){
    std::pair<G4double, G4double> epair;
    epair.first = ddcs[i][0];
    epair.second = ddcs[i][2];
    EnergySpectrum.push_back(epair);
  }
  return EnergySpectrum;
}
 
// GetEnergySpectrum with angle argument

std::vector<std::pair<G4double, G4double> >& 
EnergyAngleCrossSection::GetEnergySpectrum(const G4double& angle)
{
  EnergySpectrum.clear();
  G4int i;
  for (i=0; i < nDataPts; i++){
    if (ddcs[i][1] == angle) {
      std::pair<G4double, G4double> epair;
      epair.first = ddcs[i][0];
      epair.second = ddcs[i][2];
      EnergySpectrum.push_back(epair);
    }
  }
  return EnergySpectrum;
}


std::vector<std::pair<G4double, G4double> >&
EnergyAngleCrossSection::GetAngularDistribution(const G4double& energy)
{
  AngleSpectrum.clear();
  G4int i;
  for (i=0; i < nDataPts; i++){
    if (ddcs[i][0] == energy) {
      std::pair<G4double, G4double> apair;
      apair.first = ddcs[i][1];
      apair.second = ddcs[i][2];
      AngleSpectrum.push_back(apair);
    }
  }
  return AngleSpectrum;
}

