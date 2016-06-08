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

