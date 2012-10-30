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
//      File name:     G4DistributionGenerator
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 8 May 1998
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DistributionGenerator.hh"
#include "G4HadronicDeprecate.hh"
#include "G4ios.hh"
#include <assert.h>

// Constructor

G4DistributionGenerator::G4DistributionGenerator(std::vector<G4double>& x,
						 std::vector<G4double>& values)
  
{
  G4HadronicDeprecate("G4DistributionGenerator");
  _x = x;

  // Check boundaries: must be size(x) = size(values) + 1
  if (x.size() != (values.size() + 1))
    { G4cout << " Inconsistent parameters in G4DistributionGenerator "
	   << G4endl;
    }
  assert (x.size() == (values.size() + 1));

  G4double tot = 0.;
  unsigned int i;
  for (i=0; i<values.size(); i++) { tot += values[i]; }
  assert (tot > 0.);
  
  _cumProb.push_back(0.);    
  //  _cumProb.push_back(values[0] / tot);
  G4double sum = 0.;
  for (i=0; i<values.size(); i++) 
    { 
      sum += values[i];
      _cumProb.push_back(sum / tot); }

  // Debugging
  /*
  for (i=0; i<values.size(); i++)
    { G4cout << values[i] << "  " ; }
  G4cout << "  Integral = " << tot << G4endl;
  for (i=0; i<_cumProb.size(); i++)
    { 
      G4cout << "Variable " << _x[i]  
	   << " --- cumProb = " << _cumProb[i] << G4endl;
    }
  */
  // End of debugging

}

// Destructor

G4DistributionGenerator::~G4DistributionGenerator()
{
}


G4double G4DistributionGenerator::Generate(G4double ranflat)
{
  G4double xRandom = _x[0];

  G4int bin = _cumProb.size() - 1;
  unsigned int i;
  for (i=1; i<_cumProb.size(); i++)
    {
      if (ranflat >= _cumProb[i-1] && ranflat < _cumProb[i])
	{
	  bin = i - 1;
	}
    }

  if (bin >= 0 && bin < static_cast<G4int>(_cumProb.size()-1) && bin < static_cast<G4int>(_x.size()-1)) 
    {
      G4double coeff = (ranflat - _cumProb[bin]) *  (_x[bin+1] - _x[bin]) / 
                   (_cumProb[bin+1] - _cumProb[bin]);
      xRandom = _x[bin] + coeff;

      // Deugging
      /*
      G4cout << "Random = " << ranflat << " Generated " << xRandom << G4endl;
      */
      // Endo of Debugging

    }
  else
    { 	
      // Debugging
      /*
      G4cout << "Bin " << bin << " "
	   << _cumProb.size() << " " 
	   << _x.size()
	   << G4endl;
      */
      // End of debugging
    }

  return xRandom;

}


