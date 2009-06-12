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
// --------------------------------------------------------------------
//
// GEANT4     Test file
//
// File name: testExp.cc
//
// Author:    Simona Saliceti (simona.saliceti@ge.infn.it)
// 
// History:
// --------
// 
// 13/04/2004 Simona Saliceti 1st implementation
// --------------------------------------------------------------------
//
// Test Description: 
// ------------------
// Test of second implementation of the Empiric Model for shell cross sections in proton ionisation
// --------------------------------------------------------------------
// $Id: PIXEtest.cc,v 1.1 2009-06-12 10:43:45 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "globals.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "G4VhShellCrossSection.hh"
#include "G4hShellCrossSection.hh"
#include "G4hShellCrossSectionDoubleExp.hh"
#include "G4teoCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Proton.hh"
//#include "AIDA/AIDA.h"

int main()
{ 
   G4int Z; 
   G4double incidentEnergy;
   G4double mass;
   G4double deltaEnergy;
   size_t shellNumber;
   G4String fileName;

   //    G4VhShellCrossSection* shellExp = new G4hShellCrossSectionDoubleExp();
   //  here you deciude the implementation: G4teoCrossSection is ECPSSR for K + Orlic for L shells
   //  G4hShellCrossSectionDoubleExp is previous work with ownmade fitting functions to Pauli

   G4VhShellCrossSection* shellCS = new G4teoCrossSection("ecpssr");

   //proton mass in Kg
   //mass = 1.67262158e-27 * kg;   
   G4Proton* aProtone = G4Proton::Proton();
   mass = aProtone->GetPDGMass();

   std::vector<G4double> energies;
   //Energy from 0.005 MeV to 500 MeV
   for (G4double i=-5.5; i<6.5; i=i+0.125)
     {
       energies.push_back(std::exp(i)*MeV);
     } 

   G4cout << "Enter shell number: " << G4endl;
   G4cin >> shellNumber;

   G4cout << "Enter filename: " << G4endl;
   G4cin >> fileName;

   if (fileName == "" ) {fileName = "PIXECrossSection";}           

   G4String fileNameTxt = fileName;
   char buffer[3];
   std::ofstream myfile;

   //Z is the atomic number
   for (Z = 6; Z<=92; Z++)
     { 
       G4cout << "Z = " << Z << G4endl;
       snprintf(buffer, 3, "%d", Z);

       
       fileNameTxt = fileName + buffer + ".dat";
       
       myfile.open (fileNameTxt);
       
       
       //Cross section for each incident energy
       for (size_t k=0; k<energies.size();k++)
	 {
	   incidentEnergy = energies[k]*MeV;
	   deltaEnergy = 0.0;
	   	   
	   std::vector<G4double> CS = shellCS->GetCrossSection(Z,incidentEnergy,mass,deltaEnergy,false);
	   
	   
	   myfile << incidentEnergy << "\t" << CS[shellNumber] << G4endl;  //barn  //error in ecpssr/orlic class: correct units management!
	   
	   
	 }
       
       myfile.close();
       fileNameTxt = fileName;     
       
     } 
   delete shellCS;
   
   G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}
