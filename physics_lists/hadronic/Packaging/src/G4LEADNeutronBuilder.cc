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
 #include "G4LEADNeutronBuilder.hh"
 #include "G4ParticleDefinition.hh"
 #include "G4ParticleTable.hh"
 #include "G4ProcessManager.hh"

 G4LEADNeutronBuilder::
 G4LEADNeutronBuilder() 
 {
   theMin = 0;
   theModel = new G4Mars5GeV;
 }

 G4LEADNeutronBuilder::
 ~G4LEADNeutronBuilder() {}

 void G4LEADNeutronBuilder::
 Build(G4NeutronInelasticProcess & aP)
 {
   aP.AddDataSet(&theXSec);  
   theModel->SetMinEnergy(theMin);
   aP.RegisterMe(theModel);
 }

 void G4LEADNeutronBuilder::
 Build(G4HadronElasticProcess & ){}

 void G4LEADNeutronBuilder::
 Build(G4HadronFissionProcess & ){}

 void G4LEADNeutronBuilder::
 Build(G4HadronCaptureProcess & ){}

 // 2002 by J.P. Wellisch
