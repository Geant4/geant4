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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4VAngularDistribution
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 2000
//
//      Modifications: 
//      
// Abstract class for angular distribution strategy pattern
//
// Id: G4VAngularDistribution.hh,v 1.16 2000/05/11 19:07:29 pia Exp $ //
//
// -------------------------------------------------------------------

#ifndef G4VANGULARDISTRIBUTION_HH
#define G4VANGULARDISTRIBUTION_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4HadronicException.hh"
#include "Randomize.hh"

class G4VAngularDistribution 
{

public:

  // Constructors
  G4VAngularDistribution() { }

  virtual ~G4VAngularDistribution() { }

  virtual G4double CosTheta(G4double s, G4double m1, G4double m2) const = 0;

  virtual G4double Phi() const { return 2.*CLHEP::pi*G4UniformRand(); }

protected:

private:  

};

#endif
