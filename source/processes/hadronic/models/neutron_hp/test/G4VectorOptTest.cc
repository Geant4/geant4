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
//
// $Id: G4VectorOptTest.cc,v 1.2 2001-07-11 10:07:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
  int i;;
  for( i=9; i<10000; i++)
  {
    aVector.SetData(i,.149000E+08 + double(i)/10000.*1.E+08, (.582307+double(i)/1000.)*1.E-7);
  }
  for( i=0; i<100000; i++)
  {
    G4double it = aVector.Sample();
    G4cout << it<<G4endl;
  }
}
