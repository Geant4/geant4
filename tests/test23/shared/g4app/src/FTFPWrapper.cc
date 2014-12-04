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

#include "FTFPWrapper.hh"

#include "G4TheoFSGenerator.hh"

#include "G4FTFModel.hh"

#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"

#include "G4SystemOfUnits.hh"

FTFPWrapper::FTFPWrapper( const G4String& name, G4ProcessType type )
   : ProcessWrapper(name,type)
{

   // override
   fUseLundStrFragm = true;

}

void FTFPWrapper::Compose()
{

   G4TheoFSGenerator* gen = new G4TheoFSGenerator();

   // fInteractionModel = new G4TheoFSGenerator();

   fStringModel      = new G4FTFModel();
   
   // Note-1: one can explicitly instantiate it with G4PreCompoundModel,
   //         but if not, the ctor will either fish it out of a registry,
   //         or will create a new one on the spot
   //
   // Note-2: need the "trick" because G4VIntraNuclTransport stuff 
   //         does NOT have the SetCaptureThreshold method
   //         
   G4GeneratorPrecompoundInterface* preco = new G4GeneratorPrecompoundInterface();
   preco->SetCaptureThreshold(10.*MeV);
   fCascade = preco;
   
   // there're these 2 options for modeling string fragmentation,
   // the Lund one is "traditionally" used with FTF, and 
   // the QGSM one (older code) is typically used with QGS
   //
   if ( fUseLundStrFragm )
   {
      fStringDecay      = new G4ExcitedStringDecay(new G4LundStringFragmentation());
   }
   else
   {
      fStringDecay      = new G4ExcitedStringDecay( new G4QGSMFragmentation() );
   }
   fStringModel->SetFragmentationModel( fStringDecay );
   
   gen->SetTransport( fCascade );
   gen->SetHighEnergyGenerator( fStringModel );
   
   gen->SetMinEnergy(GeV);
   gen->SetMaxEnergy(500.*GeV);
      
   RegisterMe( gen );
   
   return;  

}
