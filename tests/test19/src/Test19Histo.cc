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

#include "TSystem.h"

#include "Test19Histo.hh"
#include "TstReader.hh"

#include "TestNA49Histo.hh"
#include "TestNA61Histo.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


Test19Histo::Test19Histo( const TstReader* pset )
   : TstHisto(pset)
{

   if ( pset->GetExpDataSet() == "NA61" )
   {
      fHistoSet = new TestNA61Histo( fHistoTitle );
      fHistoDirName  = "na61-histo";
   }
   else if ( pset->GetExpDataSet() == "NA49" )
   {
      fHistoSet = new TestNA49Histo( fHistoTitle );
      fHistoDirName  = "na49-histo";
   }
   
   // fHistoSet->SetDoResDecay(pset->ForseResDecay());
   fHistoSet->SetDoResDecay(fDoResDecay);   

}

TFile* Test19Histo::OpenHistoFile()
{

   // AccessPathName actually returns o (false) is directory ** exists **
   // and 1 (true) if it does NOT
   //
   if ( gSystem->AccessPathName( fHistoDirName.c_str() ) )
   {   
      gSystem->mkdir( fHistoDirName.c_str() );
   }
   
   G4String fname = fHistoDirName + "/" + fBeam + fTarget + fBeamMomentum + fModel;
   
   if ( fJobID > -1 )
   {
      char buf[5];
      sprintf( buf, "%i", fJobID );
      fname += "-";
      fname.append( buf ); 
   }  
   if ( fDoResDecay )
   {
      fname += "-with-decays";
   }
   fname += ".root";

   return new TFile( fname.c_str(), "recreate" );

}
