//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

// Constructor

G4LinLogInterpolation::G4LinLogInterpolation()
{ }


// Destructor

G4LinLogInterpolation::~G4LinLogInterpolation()
{ }


G4double G4LinLogInterpolation::Calculate(G4double x, G4int bin, 
					  const G4DataVector& points, 
					  const G4DataVector& data) const
{
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
      G4double d1 = log(data[bin]);
      G4double d2 = log(data[bin+1]);
      value = exp(d1 + (d2 - d1)*(x - e1)/ (e2 - e1));
    }
  else
    {
      value = data[nBins];
    }
  return value;
}
