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
//
// $Id: G4FissionLevelDensityParameterINCLXX.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitation
// by D. Mancusi (6th October 2014)
//



#ifndef G4FissionLevelDensityParameterINCLXX_h
#define G4FissionLevelDensityParameterINCLXX_h 1


#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"


/** \class G4FissionLevelDensityParameterINCLXX
 *  \brief Revised level-density parameter for fission after INCL++
 *  \author Davide Mancusi
 *  \date 6th October 2014
 *
 * This class contains a revised level-density parameter that works better than
 * the standard one with the Li\`ege Intranuclear Cascade model (INCL++). The
 * fit parameter is the ratio of the level-density parameters in the fission
 * channel and in the neutron-evaporation channel. This is commonly known as
 * af/an and is usually very close to 1.0. Variations of a few percent are
 * likely to induce large factors in the fission probability, because of the
 * exponential growth of the Fermi level density. The best values of af/an were
 * empirically found to be about 1.02 for 1-GeV p+208Pb and 1.04 for 1-GeV
 * p+U238. A linear interpolation was adopted between the two extreme values.
 */

class G4FissionLevelDensityParameterINCLXX : public G4VLevelDensityParameter
{
public:
  G4FissionLevelDensityParameterINCLXX();
  virtual ~G4FissionLevelDensityParameterINCLXX();

private:  
  G4FissionLevelDensityParameterINCLXX(const G4FissionLevelDensityParameterINCLXX &right);

  const G4FissionLevelDensityParameterINCLXX & operator=(const G4FissionLevelDensityParameterINCLXX &right);
  G4bool operator==(const G4FissionLevelDensityParameterINCLXX &right) const;
  G4bool operator!=(const G4FissionLevelDensityParameterINCLXX &right) const;
  
public:
  G4double LevelDensityParameter(G4int A, G4int Z, G4double U) const;

  void setAfanLow(const double a) { afanLow = a; UpdateAfanSlope(); }
  void setAfanHigh(const double a) { afanHigh = a; UpdateAfanSlope(); }
  void setZLow(const int z) { ZLow = z; UpdateAfanSlope(); }
  void setZHigh(const int z) { ZHigh = z; UpdateAfanSlope(); }
  double getAfanLow() const { return afanLow; }
  double getAfanHigh() const { return afanHigh; }
  int getZLow() const { return ZLow; }
  int getZHigh() const { return ZHigh; }

private:
  
  void UpdateAfanSlope();

  G4EvaporationLevelDensityParameter theEvaporationLevelDensityParameter;

  double afanLow, afanHigh;
  int ZLow, ZHigh;
  double afanSlope;

};


#endif
