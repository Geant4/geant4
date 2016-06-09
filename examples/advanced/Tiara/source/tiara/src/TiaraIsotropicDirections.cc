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
// $Id: TiaraIsotropicDirections.cc,v 1.4 2003/12/04 11:08:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

