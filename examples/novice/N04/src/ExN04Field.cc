
#include "ExN04Field.hh"

ExN04Field::ExN04Field()
{
  Bz = 3.0*tesla;
  rmax_sq = sqr(50.*cm);
  zmax = 100.*cm;
}

ExN04Field::~ExN04Field()
{;}

void ExN04Field::GetFieldValue(const double Point[3],double *Bfield) const
{
  Bfield[0] = 0.;
  Bfield[1] = 0.;
  if(abs(Point[2])<zmax && (sqr(Point[0])+sqr(Point[1]))<rmax_sq)
  { Bfield[2] = Bz; }
  else
  { Bfield[2] = 0.; }
}

