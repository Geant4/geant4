// $Id: TiaraIsotropicDirections.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
TiaraIsotropicDirections(const TiaraIsotropicDirections& rhs){
  *this = rhs;
};

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

