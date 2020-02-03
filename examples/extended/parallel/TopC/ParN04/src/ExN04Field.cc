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
/// \file ExN04Field.cc
/// \brief Implementation of the ExN04Field class
//

#include "ExN04Field.hh"
#include "G4SystemOfUnits.hh"

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
  if(std::abs(Point[2])<zmax && (sqr(Point[0])+sqr(Point[1]))<rmax_sq)
  { Bfield[2] = Bz; }
  else
  { Bfield[2] = 0.; }
}

