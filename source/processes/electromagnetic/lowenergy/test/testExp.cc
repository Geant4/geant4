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
// Author:     Simona Saliceti (simona.saliceti@ge.infn.it)
// 
// History:
// --------
// Structor from AtomicDexecitation's test
// 13/04/2004 Simona Saliceti 1st implementation
// --------------------------------------------------------------------
//
// Test Description: 
// ------------------
// Test of Empiric Model for shell cross sections in proton ionisation
// --------------------------------------------------------------------
// $Id: testExp.cc,v 1.4 2010-11-10 00:29:46 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "globals.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "G4VhShellCrossSection.hh"
#include "G4hShellCrossSection.hh"
#include "G4hShellCrossSectionExp.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4Proton.hh"
#include "AIDA/AIDA.h"

int main()
{ 
   G4int Z; 
   G4double incidentEnergy;
   G4double mass;
   G4double deltaEnergy;
   size_t shellNumber;

    G4VhShellCrossSection* shellExp = new G4hShellCrossSectionExp();
 
   //proton mass in Kg
   //mass = 1.67262158e-27 * kg;   
   G4Proton* aProtone = G4Proton::Proton();
   mass = aProtone->GetPDGMass();

   std::vector<G4double> energies;
   //Energy from 0.005 MeV to 500 MeV
   for (G4double i=-5.5; i<6.5; i=i+0.25)
     {
       energies.push_back(std::exp(i));
     } 

   G4cout << "Enter shell number: " << G4endl;
   G4cin >> shellNumber;

   // Creating the analysis factory and  the tree factory
   AIDA::IAnalysisFactory* analysisFactory = AIDA_createAnalysisFactory();
   AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();

   G4String fileName = "crossSectionElement";           
   char buffer[3];

   //Z is the atomic number
   for (Z = 6; Z<=92; Z++)
     { 
       G4cout << "Z = " << Z << G4endl;
       snprintf(buffer, 3, "%d", Z);
       fileName = fileName + buffer + ".xml";

       // Creating a tree in uncompress XML
       AIDA::ITree* tree = treeFactory->create(fileName,"xml",0,1,""); // output file
       // Creating a data point set  factory, which will be handled by the tree
       AIDA::IDataPointSetFactory* dataPointSetFactory = analysisFactory->createDataPointSetFactory(*tree);
       // Creating a 2 D data point set
       AIDA::IDataPointSet* dataPointSetElement = dataPointSetFactory->create("CS","Cross section", 2); 
       fileName = "crossSectionElement";

       //Cross section for each incident energy
       for (size_t k=0; k<energies.size();k++)
	 {
	   incidentEnergy = energies[k]*MeV;
	   deltaEnergy = 0.0;

	   // *************************************************************** //
	   // From Empiric Model                                              //
	   // *************************************************************** //
 	
	   std::vector<G4double> CS = shellExp->GetCrossSection(Z,incidentEnergy,mass,deltaEnergy,true);

	   // Write in a file.xml
	   // Fill the two dimensional IDataPointSet
	   dataPointSetElement->addPoint();
	   AIDA::IDataPoint & point = *(dataPointSetElement->point(k));
	   AIDA::IMeasurement& xPoint = *(point.coordinate(0));
	   xPoint.setValue(incidentEnergy);
	   AIDA::IMeasurement& yPoint = *(point.coordinate(1));
	   yPoint.setValue(CS[0]/barn);	     
	 }

       // Flushing the histograms into the file
       tree->commit();
       // Explicitly closing the tree
       tree->close();
     } 
   delete shellExp;

   G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}
