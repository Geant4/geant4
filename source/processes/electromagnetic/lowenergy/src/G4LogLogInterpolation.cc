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
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//         Sebastian Incerti (incerti@cenbg.in2p3.fr)
//         Nicolas A. Karakatsanis (knicolas@mail.ntua.gr)
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 27 Jun 2008   SI         Add check to avoid FPE errors
// 08 Dec 2008   NAK        Log-Log interpolation math formula streamlined, 
//                            self-test function
// 14 Jun 2008   NAK        New implementation for log-log interpolation 
//                           after directly loading
//                           logarithmic values from G4EMLOW dataset
// -------------------------------------------------------------------

#include "G4LogLogInterpolation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor
G4LogLogInterpolation::G4LogLogInterpolation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Destructor
G4LogLogInterpolation::~G4LogLogInterpolation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VDataSetAlgorithm* G4LogLogInterpolation::Clone() const 
{ return new G4LogLogInterpolation; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LogLogInterpolation::Calculate(G4double x, G4int bin, 
					  const G4DataVector& points, 
					  const G4DataVector& data) const
{
  G4int nBins = G4int(data.size() - 1);
  G4double value = 0.;
  if (x < points[0])
    {
      value = 0.;
    }
  else if (bin < nBins)
    {
      G4double e1 = points[bin];
      G4double e2 = points[bin+1];
      G4double d1 = data[bin];
      G4double d2 = data[bin+1];
      // Check of e1, e2, d1 and d2 values to avoid floating-point errors when estimating the 
      // interpolated value below 
      if ((d1 > 0.) && (d2 > 0.) && (e1 > 0.) && (e2 > 0.))
        {

	  value = std::log10(d1)+(std::log10(d2/d1)/std::log10(e2/e1)*std::log10(x/e1));
	  value = std::pow(10.,value);
        }
      else value = 0.;
    }
  else
    {
      value = data[nBins];
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Nicolas A. Karakatsanis: New implementation of log-log interpolation after directly loading 
//                          logarithmic values from G4EMLOW dataset
G4double G4LogLogInterpolation::Calculate(G4double x, G4int bin, 
					  const G4DataVector& points,
                                          const G4DataVector& data,
                                          const G4DataVector& log_points, 
					  const G4DataVector& log_data) const
{ 
  G4int nBins = G4int(data.size() - 1);
  G4double value = 0.;
  G4double log_x = std::log10(x);
  if (x < points[0])
    {
      value = 0.;
    }
  else if (bin < nBins)
    {
      G4double log_e1 = log_points[bin];
      G4double log_e2 = log_points[bin+1];
      G4double log_d1 = log_data[bin];
      G4double log_d2 = log_data[bin+1];
      
      // Values e1, e2, d1 and d2 are the log values of the corresponding
      // original energy and data values. Simple linear interpolation performed
      // on loagarithmic data should be equivalent to log-log interpolation
      value = log_d1 + (log_d2 - log_d1)*(log_x - log_e1)/(log_e2 - log_e1);

      // Delogarithmize to obtain interpolated value
      value = std::pow(10.,value);
   }
  else
    {
      value = data[nBins];
    }
  return value;
}
