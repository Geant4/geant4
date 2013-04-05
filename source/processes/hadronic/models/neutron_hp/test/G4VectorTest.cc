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
// $Id$
//
#include "g4templates.hh"
#include "G4NeutronHPVector.hh"
#include "../src/G4NeutronHPVector.cc"

main()
{

  G4NeutronHPVector aVector;
  
  aVector.SetData(0, .000000E+00, .000000E+00);
  aVector.SetData(1, .500000E+07, .000000E+00);
  aVector.SetData(2, .550000E+07, .465857E-07);
  aVector.SetData(3, .550000E+07, .465857E-07);
  aVector.SetData(4, .600000E+07, .727875E-07);
  aVector.SetData(5, .700000E+07, .107729E-06);
  aVector.SetData(6, .900000E+07, .116469E-06);
  aVector.SetData(7, .140000E+08, .110639E-06);
  aVector.SetData(8, .149000E+08, .582307E-07);
  aVector.SetData(9, .149100E+08, .000000E+00);

  for(G4int i=0; i<100000000000; i++)
  {
    G4double it = aVector.Sample();
    G4cout << it<<G4endl;
  }
}
