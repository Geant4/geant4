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
// $Id: G4ProtonBuilder.cc,v 1.4 2006-06-15 14:15:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
 #include "G4ProtonBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 
 void G4ProtonBuilder::Build()
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

 G4ProtonBuilder::
 G4ProtonBuilder(): wasActivated(false)  
 {
   theProtonInelastic=new G4ProtonInelasticProcess;
 }

 G4ProtonBuilder::
 ~G4ProtonBuilder() 
 {
   delete theProtonInelastic;
 }

 // 2002 by J.P. Wellisch
