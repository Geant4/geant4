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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LENDBertiniGammaElectroNuclearBuilder
//
// Author: 2017 Oct. T. Koi
//
// Modified: 
//----------------------------------------------------------------------------
//

#include "G4LENDBertiniGammaElectroNuclearBuilder.hh"
#include "G4LENDorBERTModel.hh"
#include "G4LENDCombinedCrossSection.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4LENDBertiniGammaElectroNuclearBuilder::G4LENDBertiniGammaElectroNuclearBuilder(G4bool eNucl) : 
G4BertiniElectroNuclearBuilder( eNucl )
{
}

G4LENDBertiniGammaElectroNuclearBuilder::~G4LENDBertiniGammaElectroNuclearBuilder() 
{
/*
  DHW 13 Jan 2020 - fix double deletion error; these deletes are already done in the base class dtor
                    (Coverity bugs 101609 and 101727) 
  if ( wasActivated ) {
     delete theFragmentation;
     delete theStringDecay;
  }
*/   
}

void G4LENDBertiniGammaElectroNuclearBuilder::Build()
{
   //G4cout << "G4LENDBertiniGammaElectroNuclearBuilder::Build()" << G4endl;

   base::Build();

   if ( !G4FindDataDir("G4LENDDATA") ) {
      G4String message = "\n Skipping activation of Low Energy Nuclear Data (LEND) model for gamma nuclear interactions.\n The LEND model needs data files and they are available from ftp://gdo-nuclear.ucllnl.org/GND_after2013/GND_v1.3.tar.gz.\n Please set the environment variable G4LENDDATA to point to the directory named v1.3 extracted from the archive file.\n"; 
      G4Exception( "G4LENDBertiniGammaElectroNuclearBuilder::Build()"
                 , "G4LENDBertiniGammaElectroNuclearBuilder001"
                 , JustWarning , message);
      return;
   }
   
   theGammaReaction->SetMinEnergy(20*MeV);
   G4LENDorBERTModel* theGammaReactionLowE = new G4LENDorBERTModel( G4Gamma::Gamma() );
   theGammaReactionLowE->DumpLENDTargetInfo(true);
   G4LENDCombinedCrossSection* theGammaCrossSectionLowE = new G4LENDCombinedCrossSection( G4Gamma::Gamma() );
   theGammaReactionLowE->SetMaxEnergy(20*MeV);
   thePhotoNuclearProcess->RegisterMe(theGammaReactionLowE);
   thePhotoNuclearProcess->AddDataSet(theGammaCrossSectionLowE);
   
}
