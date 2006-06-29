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
//////////////////////////////////////////////////////////////////////////
// $Id: G4BREPSolidConeTest.cc,v 1.11 2006-06-29 18:43:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
//
// BREP solid test, create by L. Broglia, 20/10/98
// modification of old G4Gerep test
//

#include "G4Timer.hh"
#include <cmath>
#include <fstream>
#include "G4ios.hh" 
#include "G4BREPSolid.hh"
#include "G4BREPSolidCone.hh"


int main()
{

  G4ThreeVector tStart(19000,0,10.1);
  G4ThreeVector tDir(-1,0,0);
  G4ThreeVector tStart2(0,0,10);
  G4ThreeVector pt(0,0,50);
  G4ThreeVector pt2(100000,0,50);
  double Dist;


  G4cout << "\n ============   Cone test   ================";

  G4BREPSolidCone *MyConBox = new G4BREPSolidCone ("MyConBox"          ,
						   G4ThreeVector(0,0,0),
						   G4ThreeVector(0,0,1),
						   G4ThreeVector(1,0,0),     
						   1000.0              ,
						   0.0                 ,
						   101.0                );

  G4cout << "\n\nCone created ! ";
  G4cout << "\nDir =  -1,0,0";
  G4cout << "\nStart 19000,1,0";
  Dist = MyConBox->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist;
  MyConBox->Reset();

  Dist = MyConBox->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyConBox->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyConBox->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist << G4endl;


  return EXIT_SUCCESS;
}

