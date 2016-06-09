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
// $Id: G4E1Probability.hh,v 1.6 2010-11-17 17:59:04 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------
//
// Geant4 header G4E1Probability
//
// by V. Lara (May 2003)
//
// Modifications:
// 18.05.2010 V.Ivanchenko trying to speedup the most slow method
//            by usage of G4Pow, integer A and introduction of const members
//
//

#ifndef G4E1Probability_hh
#define G4E1Probability_hh

#include "globals.hh"
#include "G4VEmissionProbability.hh"
#include "G4Fragment.hh"

class G4Pow;

class G4E1Probability : public G4VEmissionProbability
{

public:

  G4E1Probability();

  virtual ~G4E1Probability();

  G4double EmissionProbability(const G4Fragment& frag, G4double excite);
  G4double EmissionProbDensity(const G4Fragment& frag, G4double ePhoton);

private:

  G4E1Probability(const G4E1Probability& right);

  const G4E1Probability& operator=(const G4E1Probability& right);
  G4bool operator==(const G4E1Probability& right) const;
  G4bool operator!=(const G4E1Probability& right) const;

  // Integrator (simple Gaussian quadrature)
  G4double EmissionIntegration(const G4Fragment& frag,
                               G4double lowLim, G4double upLim);

  // members
  G4Pow*   fG4pow;
  G4double theLevelDensityParameter;
  G4double normC;

};

#endif
