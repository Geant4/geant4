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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "Analysis.hh"
#include "EventActionMessenger.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():G4UserEventAction()
{
  //default parameter values
  //
  fThresEdepForSSB=8.22*eV;
  fThresDistForDSB=10;
  fTotalEnergyDeposit=0;

  //create commands
  //
  fpEventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete fpEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event*)
{
  // Initialization of parameters
  //
  fTotalEnergyDeposit=0.;
  fEdepStrand1.clear();
  fEdepStrand2.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event*)
{
  // At the end of an event, compute the number of strand breaks
  //
  G4int sb[2] = {0,0};
  ComputeStrandBreaks(sb);
  // Fill histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if ( fTotalEnergyDeposit>0. )
  {
    analysisManager->FillH1(1,fTotalEnergyDeposit);
  }
  if ( sb[0]>0 )
  {
    analysisManager->FillH1(2,sb[0]);
  }
  if ( sb[1]>0 )
  {
    analysisManager->FillH1(3,sb[1]);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ComputeStrandBreaks(G4int* sb)
{
  // sb quantities
  //
  G4int ssb1=0;
  G4int ssb2=0;
  G4int dsb=0;

  // nucleotide id and energy deposit for each strand
  G4int nucl1;
  G4int nucl2;
  G4double edep1;
  G4double edep2;

  //Read strand1
  //
  while ( !fEdepStrand1.empty() )
  {
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    fEdepStrand1.erase( fEdepStrand1.begin() );

    // SSB in strand1
    //
    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }

    // Look at strand2
    //
    if ( !fEdepStrand2.empty() )
    {
      do
      {
        nucl2 = fEdepStrand2.begin()->first;
        edep2 = fEdepStrand2.begin()->second;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2++;
        }
        fEdepStrand2.erase( fEdepStrand2.begin() );
      } while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fEdepStrand2.empty()) );

      // no dsb
      //
      if ( nucl2-nucl1 > fThresDistForDSB )
      {
        fEdepStrand2[nucl2]=edep2;
        if ( edep2 >= fThresEdepForSSB/eV )
        {
          ssb2--;
        }
      }

      // one dsb
      //
      if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
      {
        if ( ( edep2 >= fThresEdepForSSB/eV ) &&
            ( edep1 >= fThresEdepForSSB/eV ) )
        {
          ssb1--;
          ssb2--;
          dsb++;
        }
      }
    }
  }

  // End with not processed data
  //
  while ( !fEdepStrand1.empty() )
  {
    nucl1 = fEdepStrand1.begin()->first;
    edep1 = fEdepStrand1.begin()->second;
    if ( edep1 >= fThresEdepForSSB/eV )
    {
      ssb1++;
    }
    fEdepStrand1.erase( fEdepStrand1.begin() );
  }

  while ( !fEdepStrand2.empty() )
  {
    nucl2 = fEdepStrand2.begin()->first;
    edep2 = fEdepStrand2.begin()->second;
    if ( edep2 >= fThresEdepForSSB/eV )
    {
      ssb2++;
    }
    fEdepStrand2.erase( fEdepStrand2.begin() );
  }

  sb[0]=ssb1+ssb2;
  sb[1]=dsb;
}
