#ifndef EnergyAngleCrossSection_h
#define EnergyAngleCrossSection_h

/////////////////////////////////////////////////////////////////////
//
// Author: D.H. Wright (SLAC)
// Date: 5 March 2003
//
/////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include "g4std/vector"
#include "g4std/fstream"

class EnergyAngleCrossSection
{
  public:
    EnergyAngleCrossSection();
    EnergyAngleCrossSection(G4std::ifstream datafile);
    ~EnergyAngleCrossSection();

    G4std::vector<G4double>& GetEnergies();
    G4std::vector<G4double>& GetAngles();

    G4std::vector<G4std::pair<G4double, G4double> >& 
          GetEnergySpectrum();
 
    G4std::vector<G4std::pair<G4double, G4double> >& 
          GetEnergySpectrum(const G4double& angle);
 
    G4std::vector<G4std::pair<G4double, G4double> >& 
          GetAngularDistribution(const G4double& energy);

  private:

    G4std::vector<G4double> Angles; 
    G4std::vector<G4double> Energies; 
    G4std::vector<G4std::pair<G4double, G4double> > EnergySpectrum; 
    G4std::vector<G4std::pair<G4double, G4double> > AngleSpectrum;

    G4int nDataPts; 
    G4double ddcs[200][4];
};

#endif



