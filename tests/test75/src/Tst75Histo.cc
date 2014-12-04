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

//#include "TSystem.h"

#include "Tst75Histo.hh"
//#include "TstReader.hh"

#include "Tst75HistoSet.hh"

//#include "G4SystemOfUnits.hh"

//#include <iostream>
//#include <fstream>
//#include <iomanip>


Tst75Histo::Tst75Histo( const TstReader* pset )
   : TstHisto(pset)
{
  
   // NOTE: This methid can be successfully called only AFTER G4ParticleTable is filled up.
   //       In principle, the way things are implemented isn't safe, because nothing guarantee 
   //       that the particle table is alset... 
   //       Will need to go over this aspect some more, and make it safer.
   //
   pset->SyncKinematics();
   
   std::ostringstream ossmom;
   ossmom << pset->GetBeamMomentum();
   
   fBeamMomentum = "";
   fBeamMomentum += ossmom.str();
   fBeamMomentum += "MeV";
     
   fHistoTitle = fBeam + "(" + ossmom.str() + "MeV) + " + fTarget + " -> X "; 
     
   fHistoSet = new Tst75HistoSet( fHistoTitle, pset->GetBeamKineticEnergy() );
   
}

TFile* Tst75Histo::OpenHistoFile()
{
   
   G4String fname = fBeam + fBeamMomentum + fTarget + fModel;
   
   if ( fJobID > -1 )
   {
      char buf[5];
      sprintf( buf, "%i", fJobID );
      fname += "-";
      fname.append( buf ); 
   }  
   fname += ".root";

   return new TFile( fname.c_str(), "recreate" );

}

