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
// $Id: TiaraIsotropicDirections.cc,v 1.5 2006/06/29 15:45:12 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

#include "TiaraIsotropicDirections.hh"
#include "TiaraDimensions.hh"
#include "Randomize.hh"
#include <cmath>

TiaraIsotropicDirections::
TiaraIsotropicDirections(G4double colWidth,
			 const TiaraDimensions &tiaraDimensions) :
  fMinCos(MinimumCosine(colWidth,
			tiaraDimensions))
{}

TiaraIsotropicDirections::~TiaraIsotropicDirections()
{}

TiaraIsotropicDirections::
TiaraIsotropicDirections(const TiaraIsotropicDirections& rhs)
  : TiaraVDirectionGenerator()
{
  *this = rhs;
}

G4double TiaraIsotropicDirections::
MinimumCosine(G4double colWidth,
	      const TiaraDimensions &tiaraDimensions) {
  G4double l(tiaraDimensions.distTargetExperiment + colWidth);
  G4double r(tiaraDimensions.pipeRadius);
  G4double R(std::sqrt(std::pow(l, 2) + std::pow(r, 2)));
  return l/R; 
}

TiaraIsotropicDirections& TiaraIsotropicDirections::
operator=(const TiaraIsotropicDirections& rhs){
  if (this != &rhs) {
    fMinCos = rhs.fMinCos;
  }
  return *this;
}

TiaraVDirectionGenerator *TiaraIsotropicDirections::Clone() const {
  return new TiaraIsotropicDirections(*this);
}


G4ThreeVector TiaraIsotropicDirections::GetDirection(){
  G4double cosTh(fMinCos + G4UniformRand() * (1 - fMinCos));
  G4double sinTh(std::sqrt(1-std::pow(cosTh,2)));
  G4double phi(G4UniformRand() * 2*pi);
  G4double x(sinTh * std::cos(phi));
  G4double y(sinTh * std::sin(phi));
  return G4ThreeVector(x, y, cosTh);
}

