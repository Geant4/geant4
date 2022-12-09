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
// 31 Jul 2008   MGP        Revised and renamed to G4LinInterpolator
//
// -------------------------------------------------------------------

#include "G4LinInterpolator.hh"

// Constructor

G4LinInterpolator::G4LinInterpolator()
{ }

// Destructor

G4LinInterpolator::~G4LinInterpolator()
{ }


G4IInterpolator* G4LinInterpolator::Clone() const 
{ return new G4LinInterpolator; }

G4double G4LinInterpolator::Calculate(G4double x, G4int bin, 
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
      value = d1 + (d2 - d1)*(x - e1)/(e2 - e1);
    }
  else
    {
      value = data[nBins];
    }
  return value;
}
