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
// $Id: G4SemiLogInterpolation.cc 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4SemiLogInterpolation.hh"

// Constructor

G4SemiLogInterpolation::G4SemiLogInterpolation()
{ }


// Destructor

G4SemiLogInterpolation::~G4SemiLogInterpolation()
{ }

G4VDataSetAlgorithm* G4SemiLogInterpolation::Clone() const 
{ return new G4SemiLogInterpolation; }


G4double G4SemiLogInterpolation::Calculate(G4double x, G4int bin, 
					  const G4DataVector& points, 
					  const G4DataVector& data) const
{
  //G4cout << "G4SemiLogInterpolation is performed(2 arguments) " << G4endl;
  G4int nBins = data.size() - 1;
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
      value = (d1*std::log10(e2/x) + d2*std::log10(x/e1)) / std::log10(e2/e1);
    }
  else
    {
      value = data[nBins];
    }
  return value;
}

G4double G4SemiLogInterpolation::Calculate(G4double x, G4int bin, 
				           const G4DataVector& points,
                                           const G4DataVector& data,
                                           const G4DataVector& log_points, 
				           const G4DataVector& /*log_data*/) const
{
//A combination of logarithmic interpolation on energy set and 
//linear Interpolation on data set
  //G4cout << "G4SemiLogInterpolation is performed (4 arguments)" << G4endl;
  G4int nBins = data.size() - 1;
  G4double value = 0.;
  G4double log_x = std::log10(x);
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
      G4double log_e1 = log_points[bin];
      G4double log_e2 = log_points[bin+1];
      //G4double log_d1 = log_data[bin];
      //G4double log_d2 = log_data[bin+1];
      //G4cout << "x = " << x << " , log_x = " << log_x << " , bin = " << bin << G4endl;
      //G4cout << "e1 = " << e1 << " , d1 = " << d1 << G4endl;
      //G4cout << "e2 = " << e2 << " , d2 = " << d2 << G4endl;
// Values log_e1 and log_e2 are the log values of the corresponding
// original energy actual values. Original d1 and d2 values are used.
// Simple linear interpolation performed on loagarithmic data 
// should be equivalent to semi log-log interpolation
      if (e1 == 0.0) log_e1 = -300;
      if (e2 == 0.0) log_e2 = -300;
      value = d1 + (d2 - d1)*(log_x - log_e1)/(log_e2 - log_e1);
      //G4cout << "G4SemiLogInterpolation - Final Interpolated Value: " << value << G4endl << G4endl;
   }
 else
   {
     value = data[nBins];
   }
  return value;
}

