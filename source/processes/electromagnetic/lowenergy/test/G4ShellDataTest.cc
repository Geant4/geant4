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
// $Id: G4ShellDataTest.cc,v 1.1 2001-08-20 17:28:51 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4DataSetManagerTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 6 August 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4ShellData.hh"

int main()
{
  G4cout.setf( ios::scientific, ios::floatfield );

  G4String file = "fluor/binding";

  G4ShellData* dataSet = new G4ShellData(file);
  
  G4cout << "Enter Z" << G4endl;
  G4int Z;
  G4cin >> Z;

  G4int n = dataSet->NumberOfShells(Z);
  G4cout << "Z = " << Z << " has " << n << " shells" << G4endl;
  
  G4cout << "Enter shell index " << G4endl;
  G4int i;
  G4cin >> i;

  G4int id = dataSet->ShellId(Z,i);
  G4double e = dataSet->BindingEnergy(Z,i) / keV;

  G4cout << "Shell id = " << id 
	 << ", Binding energy = " << e << " keV" <<G4endl;

  G4DataVector idVector = dataSet->ShellIdVector(Z);
  for (G4int ind=0; ind<n; ind++)
    {
      G4int idx = (G4int) idVector[ind];
      G4cout << "Id vector(" << ind << ") = " << idx << G4endl; 
    }

  G4cout << "Dump data (1) or stop (2)" << G4endl;
  G4int k;
  G4cin >> k;

  if (k == 1) dataSet->PrintData();  

  delete dataSet;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}








