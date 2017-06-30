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
// $Id: G4DeuteronBuilder.cc,v 1.2 2013/02/26 10:34:03 arce Exp $
// GEANT4 tag $Name:  $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4DeuteronBuilder
//
// Author: 2013 P. Arce
//
// Modified
// 12.04.2017 A.Dotti move to new design with base class
//----------------------------------------------------------------------------
//
 #include "G4DeuteronBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 
 void G4DeuteronBuilder::Build()
 {
   wasActivated = true;
   std::vector<G4VDeuteronBuilder *>::iterator i;
   for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
   {
     (*i)->Build(theDeuteronInelastic);
   }
   G4ProcessManager * theProcMan = G4Deuteron::Deuteron()->GetProcessManager();
   theProcMan->AddDiscreteProcess(theDeuteronInelastic);
 }

 G4DeuteronBuilder::
 G4DeuteronBuilder(): wasActivated(false)  
 {
   theDeuteronInelastic=new G4DeuteronInelasticProcess;
 }

 void G4DeuteronBuilder::RegisterMe(G4PhysicsBuilderInterface* aB) {
   auto bld = dynamic_cast<G4VDeuteronBuilder*>(aB);
   if ( bld != nullptr ) {
       theModelCollections.push_back(bld);
   } else {
       G4PhysicsBuilderInterface::RegisterMe(aB);
   }
 }
