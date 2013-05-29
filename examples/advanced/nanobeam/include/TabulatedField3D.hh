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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#include "G4MagneticField.hh"
#include "G4ios.hh"

#include "globals.hh"
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

class TabulatedField3D

#ifndef STANDALONE
 : public G4MagneticField
#endif

{
    
public:
  TabulatedField3D(G4float gr1, G4float gr2, G4float gr3, G4float gr4, G4int quadModel);
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;

private:
  vector< vector< vector< double > > > xField;
  vector< vector< vector< double > > > yField;
  vector< vector< vector< double > > > zField;
  G4int nx,ny,nz; 
  G4double minx, maxx, miny, maxy, minz, maxz;
  G4double dx, dy, dz;
  
  G4float gradient1, gradient2, gradient3, gradient4;
  G4int model;

};

