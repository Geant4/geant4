#ifndef TiaraIsotropicDirections_hh
#define TiaraIsotropicDirections_hh TiaraIsotropicDirections_hh

#include "TiaraVDirectionGenerator.hh"
#include "G4ThreeVector.hh"

class TiaraDimensions;

class TiaraIsotropicDirections : public TiaraVDirectionGenerator {
public:
  TiaraIsotropicDirections(G4double colWidth,
			   const TiaraDimensions &tiaraDimensions);
  ~TiaraIsotropicDirections();
  TiaraIsotropicDirections(const TiaraIsotropicDirections& rhs);

  virtual G4ThreeVector GetDirection();
  virtual TiaraVDirectionGenerator *Clone() const;

  G4double MinimumCosine(G4double colWidth,
			 const TiaraDimensions &tiaraDimensions);
  
  TiaraIsotropicDirections& operator=(const TiaraIsotropicDirections& rhs);
private:
  G4double fMinCos;
};

#endif
