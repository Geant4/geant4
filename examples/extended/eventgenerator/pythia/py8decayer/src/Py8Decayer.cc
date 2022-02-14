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
/// \file eventgenerator/pythia/pythia8decayer/src/Py8Decayer.cc
/// \brief Implementation of the Py8Decayer class
///
/// \author J. Yarba; FNAL

#include "Py8Decayer.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

using namespace Pythia8;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Py8Decayer::Py8Decayer()
   : G4VExtDecayer("Py8Decayer"),
   fDecayer(nullptr)
{

   // use default path to the xml/configs but do NOT print banner
   //
   fDecayer = new Pythia( "../share/Pythia8/xmldoc", false );

   // this is the trick to make Pythia8 work as "decayer"
   //
   fDecayer->readString("ProcessLevel:all = off"); 

   fDecayer->readString("ProcessLevel:resonanceDecays=on");
    
   // shut off Pythia8 (default) verbosty
   //
   fDecayer->readString("Init:showAllSettings=false");
   fDecayer->readString("Init:showChangedSettings=false");
   fDecayer->readString("Init:showAllParticleData=false");
   fDecayer->readString("Init:showChangedParticleData=false");
   //
   // specify how many Py8 events to print out, at either level
   // in this particular case print out a maximum of 10 events
   //
   fDecayer->readString("Next:numberShowProcess = 0" );
   fDecayer->readString("Next:numberShowEvent = 10");
           
   fDecayer->init();
   
   // shut off decays of pi0's as we want Geant4 to handle them
   // if other immediate decay products should be handled by Geant4,
   // their respective decay modes should be shut off as well
   //
   fDecayer->readString("111:onMode = off");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Py8Decayer::~Py8Decayer()
{

   delete fDecayer;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DecayProducts* Py8Decayer::ImportDecayProducts(const G4Track& track)
{

   fDecayer->event.reset();
   
   G4DecayProducts* dproducts = nullptr;   
   
   G4ParticleDefinition* pd = track.GetDefinition();
   int    pdgid   = pd->GetPDGEncoding();
   
   // check if pdgid is consistent with Pythia8 particle table
   //   
   if ( !fDecayer->particleData.findParticle( pdgid ) )
   {
      G4cout << " can NOT find pdgid = " << pdgid 
             << " in Pythia8::ParticleData" << G4endl;
      return dproducts;
   }
   
   if ( !fDecayer->particleData.canDecay(pdgid) )
   {
      G4cout << " Particle of pdgid = " << pdgid 
             << " can NOT be decayed by Pythia8" << G4endl;
      return dproducts;
   }
   
   // NOTE: Energy should be in GeV 

   fDecayer->event.append( pdgid, 1, 0, 0, 
                           track.GetMomentum().x() / CLHEP::GeV, 
                           track.GetMomentum().y() / CLHEP::GeV,  
                           track.GetMomentum().z() / CLHEP::GeV,
                           track.GetDynamicParticle()->GetTotalEnergy() / CLHEP::GeV,
                           pd->GetPDGMass() / CLHEP::GeV );

   // specify polarization, if any
   
   // NOTE: while in Py8 polarization is a double variable , 
   //       in reality it's expected to be -1, 0., or 1 in case of "external" tau's, 
   //       similar to LHA SPINUP; see Particle Decays, Hadron and Tau Decays in docs at
   //       https://pythia.org/manuals/pythia8305/Welcome.html
   //       so it's not able to handle anything like 0.99, thus we're rounding off    
   fDecayer->event.back().pol( 
      round( std::cos( track.GetPolarization().angle( track.GetMomentumDirection() ) )
      ) 
   );

   int npart_before_decay = fDecayer->event.size();
   
   fDecayer->next();
   
   int npart_after_decay = fDecayer->event.size();
   
   // create & fill up decay products
   //
   dproducts = new G4DecayProducts(*(track.GetDynamicParticle()));
   
   // create G4DynamicParticle out of each fDecayer->event entry (except the 1st one)
   // and push into dproducts
   
   for ( int ip=npart_before_decay; ip<npart_after_decay; ++ip )
   {
      
      // only select final state decay products (direct or via subsequent decays);
      // skip all others
      //
      // NOTE: in general, final state decays products will have 
      //       positive status code between 91 and 99 
      //       (in case such information could be of interest in the future)
      //
      if ( fDecayer->event[ip].status() < 0 ) continue;
            
      // should we also skip neutrinos ???
      // if so, skip products with abs(fDecayer->event[ip].id()) of 12, 14, or 16
            
      G4ParticleDefinition* pddec = 
         G4ParticleTable::GetParticleTable()->FindParticle( fDecayer->event[ip].id() );
      if ( !pddec ) continue; // maybe we should print out a warning !
      G4ThreeVector momentum = G4ThreeVector( fDecayer->event[ip].px() * CLHEP::GeV,
                                              fDecayer->event[ip].py() * CLHEP::GeV,
                                              fDecayer->event[ip].pz() * CLHEP::GeV ); 
      dproducts->PushProducts( new G4DynamicParticle( pddec, momentum) ); 
   }
   
   return dproducts;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

