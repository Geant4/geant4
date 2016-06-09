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
// $Id: G4OsloMatrix.cc,v 1.7 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4OsloMatrix.cc
//
// ----------------------------------------------------------------------

#include "G4OsloMatrix.hh"

G4OsloMatrix::G4OsloMatrix()
  : next(0), offset(0), osize(0), o_vec(0)
{
}

G4OsloMatrix::G4OsloMatrix(G4int vec_size, G4int offsetparam, G4int osizeparam)
{
  next   = 0;
  o_vec  = new G4KnotVector(vec_size);
  offset = offsetparam;
  osize  = osizeparam;
}

G4OsloMatrix::~G4OsloMatrix()
{
  delete o_vec;
}

G4PointMatrix::G4PointMatrix()
{
  nr=nc=0;
  data=0;
}


G4PointMatrix::G4PointMatrix(G4int rows, G4int columns)
{
  nr=rows; nc=columns; 
  data = new G4double[nr*nc];
  for(G4int a =0; a<nr*nc;a++) 
    { data[a]=0; }
}


G4PointMatrix::G4PointMatrix(G4double vec[])
{
  nr = nc = 4; 
  data = new G4double[nr*nc];
  for(G4int a=0;a<nr*nc;a++)
    { data[a]=vec[a]; }
}

G4PointMatrix::~G4PointMatrix()
{
  delete [] data;
}
