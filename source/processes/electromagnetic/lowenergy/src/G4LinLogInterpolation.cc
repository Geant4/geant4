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
// Author:  Vladimir Ivanchenko (Vladimir.Ivantchenko@cern.ch)
//
// History:
// -----------
// 27 May 2002   VI        Created
//
// -------------------------------------------------------------------

#include "G4LinLogInterpolation.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Constructor
G4LinLogInterpolation::G4LinLogInterpolation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// Destructor
G4LinLogInterpolation::~G4LinLogInterpolation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4LinLogInterpolation::Calculate(G4double x, G4int bin,
					  const G4DataVector& points,
					  const G4DataVector& data) const
{
  //G4cout << "G4LinLogInterpolation is performed on dataset (2 arguments) " << G4endl;
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
      G4double d1 = std::log(data[bin]);
      G4double d2 = std::log(data[bin+1]);
      value = G4Exp(d1 + (d2 - d1)*(x - e1)/ (e2 - e1));
    }
  else
    {
      value = data[nBins];
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LinLogInterpolation::Calculate(G4double x, G4int bin,
					  const G4DataVector& points,
                                          const G4DataVector& data,
                                          const G4DataVector& /*log_points*/,
					  const G4DataVector& log_data) const
{
  //G4cout << "G4LinLogInterpolation is performed on dataset (4 arguments) " << G4endl;
  G4int nBins = G4int(data.size() - 1);
  G4double value = 0.;
  //G4double log_x = std::log10(x);
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
      //G4double log_e1 = log_points[bin];
      //G4double log_e2 = log_points[bin+1];
      G4double log_d1 = log_data[bin];
      G4double log_d2 = log_data[bin+1];
      if (d1 > 0.0 && d2 > 0.0)
        {
	  // Values e1, e2, d1 and d2 are the log values of the corresponding
	  // original energy and data values. 
	  // Simple linear interpolation performed
	  // on logarithmic data should be equivalent to log-log interpolation
          //value = log_d1 + (log_d2 - log_d1)*(log_x - log_e1)/(log_e2 - log_e1);
          value = G4Exp(log_d1 + (log_d2 - log_d1)*(x - e1)/(e2 - e1));
        }
      else
        {
          if (d1 == 0.0) log_d1 = -300;
          if (d2 == 0.0) log_d2 = -300;
          value = G4Exp(log_d1 + (log_d2 - log_d1)*(x - e1)/(e2 - e1));
        }
   }
 else
   {
     value = data[nBins];
   }
  return value;
}
