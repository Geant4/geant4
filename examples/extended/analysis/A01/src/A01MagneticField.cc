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
// $Id: A01MagneticField.cc,v 1.6 2006/06/29 16:32:57 gunter Exp $
// --------------------------------------------------------------
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
  if(std::abs(Point[1])<ymax && (sqr(Point[0])+sqr(Point[2]))<rmax_sq)
  { Bfield[1] = By; }
  else
  { Bfield[1] = 0.; }
}

