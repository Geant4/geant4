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
// Author: Mathieu Karamitros
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include <algorithm>
#include "G4DNAOneStepThermalizationModel.hh"
#include "globals.hh"
#include "G4Exp.hh"
#include "G4RandomDirection.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"

//------------------------------------------------------------------------------

namespace DNA {
namespace Penetration {
  
const double
Meesungnoen2002::gCoeff[13] =
{ -4.06217193e-08,  3.06848412e-06,  -9.93217814e-05,
  1.80172797e-03,  -2.01135480e-02,   1.42939448e-01,
  -6.48348714e-01,  1.85227848e+00,  -3.36450378e+00,
  4.37785068e+00,  -4.20557339e+00,   3.81679083e+00,
  -2.34069784e-01 };
// fit from Meesungnoen, 2002

const double
Terrisol1990::gEnergies_T1990[11] =
{ 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7,
  // The two last are not in the dataset
  8, 9}; // eV

const double
Terrisol1990::gStdDev_T1990[11] =
{ 17.68*CLHEP::angstrom,
  22.3*CLHEP::angstrom,
  28.49*CLHEP::angstrom,
  45.35*CLHEP::angstrom,
  70.03*CLHEP::angstrom,
  98.05*CLHEP::angstrom,
  120.56*CLHEP::angstrom,
  132.73*CLHEP::angstrom,
  142.60*CLHEP::angstrom,
  // the above value as given in the paper's table does not match
  // b=27.22 nm nor the mean value. 129.62*CLHEP::angstrom could be
  // a better fit.
  //
  // The two last are made up
  137.9*CLHEP::angstrom,
  120.7*CLHEP::angstrom
}; // angstrom

//----------------------------------------------------------------------------

double Meesungnoen2002::GetRmean(double k){
  G4double k_eV = k/eV;

  if(k_eV>0.1){ // data until 0.2 eV
    G4double r_mean = 0;
    for(int8_t i=12; i!=-1 ; --i){
      r_mean+=gCoeff[12-i]*std::pow(k_eV,i);
    }
    r_mean*=CLHEP::nanometer;
    return r_mean;
  }
  return 0;
}

void GetGaussianPenetrationFromRmean3D(G4double r_mean,
                                       G4ThreeVector& displacement)
{
  if(r_mean == 0)
  {
    // rare events:
    // prevent H2O and secondary electron from being placed at the same position
    displacement = G4RandomDirection() * (1e-3*CLHEP::nanometer);
    return;
  }
  
  static constexpr double convertRmean3DToSigma1D = 0.62665706865775006;
  // = sqrt(CLHEP::pi)/pow(2,3./2.)
    
  // Use r_mean to build a 3D gaussian
  const double sigma1D = r_mean * convertRmean3DToSigma1D;
  displacement = G4ThreeVector(G4RandGauss::shoot(0, sigma1D),
                               G4RandGauss::shoot(0, sigma1D),
                               G4RandGauss::shoot(0, sigma1D));
}
  
void Meesungnoen2002::GetPenetration(G4double k,
                                     G4ThreeVector& displacement)
{
  GetGaussianPenetrationFromRmean3D(GetRmean(k), displacement);
}

//----------------------------------------------------------------------------

void Ritchie1994::GetPenetration(G4double k,
                                 G4ThreeVector& displacement)
{
  GetGaussianPenetrationFromRmean3D(k/eV * 1.8 * nm, // r_mean
                                    displacement);
}

//----------------------------------------------------------------------------

double Terrisol1990::Get3DStdDeviation(double energy){
  G4double k_eV = energy/eV;
  if(k_eV < 0.2){
    // rare events:
    //  prevent H2O and secondary electron to be at the spot
    return 1e-3*CLHEP::nanometer;
  }
  else if(k_eV == 9.){
    return gStdDev_T1990[10];
  }
  else if(k_eV > 9.){
    G4ExceptionDescription description;
    description << "Terrisol1990 is not tabulated for energies greater than 9eV";
    G4Exception("Terrisol1990::Get3DStdDeviation",
                "INVALID_ARGUMENT",
                FatalErrorInArgument,
                description);
  }

  size_t lowBin, upBin;

  if(k_eV >= 1.){
    lowBin=std::floor(k_eV)+1;
    upBin=std::min(lowBin+1, size_t(10));
  }
  else{
    auto it=std::lower_bound(&gEnergies_T1990[0],
                             &gEnergies_T1990[2],
                             k_eV);
    lowBin = it-&gEnergies_T1990[0];
    upBin = lowBin+1;
  }

  double lowE = gEnergies_T1990[lowBin];
  double upE = gEnergies_T1990[upBin];

  double lowS = gStdDev_T1990[lowBin];
  double upS = gStdDev_T1990[upBin];

  double tanA = (lowS-upS)/(lowE-upE);
  double sigma3D = lowS + (k_eV-lowE)*tanA;
  return sigma3D;
}

double Terrisol1990::GetRmean(double energy){
  double sigma3D=Get3DStdDeviation(energy);

  static constexpr double s2r=1.595769121605731; // = pow(2,3./2.)/sqrt(CLHEP::pi)

  double r_mean=sigma3D*s2r;
  return r_mean;
}

void Terrisol1990::GetPenetration(G4double energy,
                                  G4ThreeVector& displacement){
  double sigma3D = Get3DStdDeviation(energy);

  static constexpr double factor = 2.20496999539; // = 1./(3. - 8./CLHEP::pi);

  double sigma1D = std::sqrt(std::pow(sigma3D, 2.)*factor);

  displacement = G4ThreeVector(G4RandGauss::shoot(0, sigma1D),
                               G4RandGauss::shoot(0, sigma1D),
                               G4RandGauss::shoot(0, sigma1D));
}

} // Penetration
} // DNA

//------------------------------------------------------------------------------

G4VEmModel* G4DNASolvationModelFactory::Create(const G4String& penetrationModel)
{
  G4String modelNamePrefix("DNAOneStepThermalizationModel_");
  
  if(penetrationModel == "Terrisol1990")
  {
    return new G4TDNAOneStepThermalizationModel<DNA::Penetration::Terrisol1990>(G4Electron::Definition(), modelNamePrefix + penetrationModel);
  }
  else if(penetrationModel == "Meesungnoen2002")
  {
    return new G4TDNAOneStepThermalizationModel<DNA::Penetration::Meesungnoen2002>(G4Electron::Definition(), modelNamePrefix + penetrationModel);
  }
  else if(penetrationModel == "Ritchie1994")
  {
    return new G4TDNAOneStepThermalizationModel<DNA::Penetration::Ritchie1994>(G4Electron::Definition(), modelNamePrefix + penetrationModel);
  }
  else
  {
    G4ExceptionDescription description;
    description << penetrationModel + " is not a valid model name.";
    G4Exception("G4DNASolvationModelFactory::Create",
                "INVALID_ARGUMENT",
                FatalErrorInArgument,
                description,
                "Options are: Terrisol1990, Meesungnoen2002, Ritchie1994.");
  }
  return nullptr;
}

//------------------------------------------------------------------------------
G4VEmModel* G4DNASolvationModelFactory::GetMacroDefinedModel()
{
  auto dnaSubType = G4EmParameters::Instance()->DNAeSolvationSubType();
  
  switch(dnaSubType)
  {
  case fRitchie1994eSolvation:
    return Create("Ritchie1994");
  case fTerrisol1990eSolvation:
    return Create("Terrisol1990");
  case fMeesungnoen2002eSolvation:
  case fDNAUnknownModel:
    return Create("Meesungnoen2002");
  default:
    G4Exception("G4DNASolvationModelFactory::GetMacroDefinedModel",
                "DnaSubType",
                FatalErrorInArgument,
                "The solvation parameter stored in G4EmParameters is unknown. Supported types are: fRitchie1994eSolvation, fTerrisol1990eSolvation, fMeesungnoen2002eSolvation.");
  }

  return nullptr;
}
