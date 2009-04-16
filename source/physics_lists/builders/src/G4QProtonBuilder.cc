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
// $Id: G4QProtonBuilder.cc,v 1.1 2009-04-16 09:24:05 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4QProtonBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
 #include "G4QProtonBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 
 void G4QProtonBuilder::Build()
 {
   wasActivated = true;
   std::vector<G4VProtonBuilder *>::iterator i;
   for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
   {
     (*i)->Build(theProtonInelastic);
   }
   G4ProcessManager * theProcMan = G4Proton::Proton()->GetProcessManager();
   theProcMan->AddDiscreteProcess(theProtonInelastic);
 }

 G4QProtonBuilder::
 G4QProtonBuilder(): wasActivated(false)  
 {
   theProtonInelastic=new G4ProtonInelasticProcess;
 }

 G4QProtonBuilder::
 ~G4QProtonBuilder() 
 {
   delete theProtonInelastic;
 }

 // 2009 by M. Kosov
