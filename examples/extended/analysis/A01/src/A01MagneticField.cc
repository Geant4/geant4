// $Id: A01MagneticField.cc,v 1.1 2002-11-13 07:23:47 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01MagneticField.hh"
#include "A01MagneticFieldMessenger.hh"

A01MagneticField::A01MagneticField()
{
  messenger = new A01MagneticFieldMessenger(this);
  By = 1.0*tesla;
  rmax_sq = sqr(1.*m);
  ymax = 100.*cm;
}

A01MagneticField::~A01MagneticField()
{ delete messenger; }

void A01MagneticField::GetFieldValue(const double Point[3],double *Bfield) const
{
  Bfield[0] = 0.;
  Bfield[2] = 0.;
  if(abs(Point[2])<ymax && (sqr(Point[0])+sqr(Point[2]))<rmax_sq)
  { Bfield[1] = By; }
  else
  { Bfield[1] = 0.; }
}

