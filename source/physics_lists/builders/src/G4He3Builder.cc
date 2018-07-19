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
// $Id: G4He3Builder.cc,v 1.2 2013/02/26 10:34:08 arce Exp $
// GEANT4 tag $Name:  $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4He3Builder
//
// Author: 2013 P. Arce
//
// Modified
// 12.04.2017 A.Dotti move to new design with base class
//
//----------------------------------------------------------------------------
//
 #include "G4He3Builder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 
 void G4He3Builder::Build()
 {
   wasActivated = true;
   std::vector<G4VHe3Builder *>::iterator i;
   for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
   {
     (*i)->Build(theHe3Inelastic);
   }
   G4ProcessManager * theProcMan = G4He3::He3()->GetProcessManager();
   theProcMan->AddDiscreteProcess(theHe3Inelastic);
 }

 G4He3Builder::
 G4He3Builder(): wasActivated(false)  
 {
   theHe3Inelastic=new G4He3InelasticProcess;
 }

 void G4He3Builder::RegisterMe(G4PhysicsBuilderInterface* aB) {
   auto bld = dynamic_cast<G4VHe3Builder*>(aB);
   if ( bld != nullptr ) {
       theModelCollections.push_back(bld);
   } else {
       G4PhysicsBuilderInterface::RegisterMe(aB);
   }
 }
