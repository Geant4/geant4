#ifndef included_G4SurfaceOfLinearExtrusion
#define included_G4SurfaceOfLinearExtrusion

// surface of linear extrusion

#include "G4Surface.hh"

class G4SurfaceOfLinearExtrusion: public G4Surface
{

public:

  G4SurfaceOfLinearExtrusion();
  virtual ~G4SurfaceOfLinearExtrusion();

private:

  G4SurfaceOfLinearExtrusion(const G4SurfaceOfLinearExtrusion &);
  G4SurfaceOfLinearExtrusion& operator=(const G4SurfaceOfLinearExtrusion &);

};

#endif



