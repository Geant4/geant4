#ifndef EnergyAngleCrossSection_h
#define EnergyAngleCrossSection_h

/////////////////////////////////////////////////////////////////////
//
// Author: D.H. Wright (SLAC)
// Date: 5 March 2003
//
/////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include <vector>
#include <fstream>

class EnergyAngleCrossSection
{
  public:
    EnergyAngleCrossSection();
    EnergyAngleCrossSection(std::ifstream datafile);
    ~EnergyAngleCrossSection();

    std::vector<G4double>& GetEnergies();
    std::vector<G4double>& GetAngles();

    std::vector<std::pair<G4double, G4double> >& 
          GetEnergySpectrum();
 
    std::vector<std::pair<G4double, G4double> >& 
          GetEnergySpectrum(const G4double& angle);
 
    std::vector<std::pair<G4double, G4double> >& 
          GetAngularDistribution(const G4double& energy);

  private:

    std::vector<G4double> Angles; 
    std::vector<G4double> Energies; 
    std::vector<std::pair<G4double, G4double> > EnergySpectrum; 
    std::vector<std::pair<G4double, G4double> > AngleSpectrum;

    G4int nDataPts; 
    G4double ddcs[300][4];
};

#endif



