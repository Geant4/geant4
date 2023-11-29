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
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 14 JUn 2009   NAK        New Calculation method implemented to which logarithmic values
//                          from the G4EMLOW dataset are loaded directly to enhance performance
//
// -------------------------------------------------------------------

#include "G4LinInterpolation.hh"

// Constructor
// -------------------------------------------------------------------
G4LinInterpolation::G4LinInterpolation()
{ }

// -------------------------------------------------------------------
// Destructor
G4LinInterpolation::~G4LinInterpolation()
{ }

// -------------------------------------------------------------------
G4VDataSetAlgorithm* G4LinInterpolation::Clone() const 
{ return new G4LinInterpolation; }

// -------------------------------------------------------------------
G4double G4LinInterpolation::Calculate(G4double x, G4int bin, 
				       const G4DataVector& points, 
				       const G4DataVector& data) const
{
  //G4cout << "G4LinInterpolation is performed (2 arguments)" << G4endl;
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
      value = d1 + (d2 - d1)*(x - e1)/(e2 - e1);
    }
  else
    {
      value = data[nBins];
    }
  return value;
}

// -------------------------------------------------------------------
//Nicolas A. Karakatsanis: New Calculation method implemented to which logarithmic values
//                         from the G4EMLOW dataset are loaded directly to enhance performance

G4double G4LinInterpolation::Calculate(G4double x, G4int bin, 
				       const G4DataVector& points,
                                       const G4DataVector& data,
                                       const G4DataVector& log_points, 
				       const G4DataVector& log_data) const
{
//Linear Interpolation is performed on loagarithmic data set
//Equivalent to log-log interpolation on non-loagarithmic data set
  //G4cout << "G4LinInterpolation is performed - (4 arguments)" << G4endl;
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
