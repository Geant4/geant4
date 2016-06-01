#ifndef included_G4SurfaceOfRevolution
#define included_G4SurfaceOfRevolution

// surface of linear extrusion

#include "G4Surface.hh"


class G4SurfaceOfRevolution: public G4Surface
{

public:

  G4SurfaceOfRevolution();
  virtual ~G4SurfaceOfRevolution();


private:

  G4SurfaceOfRevolution(const G4SurfaceOfRevolution &);
  G4SurfaceOfRevolution& operator=(const G4SurfaceOfRevolution &);

};

#endif



