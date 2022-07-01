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
// G4HEPEvtInterface class implementation
//
// Author: Makoto Asai, 1997
// --------------------------------------------------------------------

#include "G4HEPEvtInterface.hh"

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4HEPEvtParticle.hh"
#include "G4Event.hh"

G4HEPEvtInterface::G4HEPEvtInterface(const char* evfile, G4int vl)
  : vLevel(vl)
{
  inputFile.open((char*)evfile);
  if (inputFile.is_open())
  {
    fileName = evfile;
    if(vl>0)
      G4cout << "G4HEPEvtInterface - " << fileName << " is open." << G4endl;
  }
  else
  {
    G4Exception("G4HEPEvtInterface::G4HEPEvtInterface","Event0201",
                FatalException, "G4HEPEvtInterface:: cannot open file.");
  }
  G4ThreeVector zero;
  particle_position = zero;
  particle_time = 0.0;
}

void G4HEPEvtInterface::GeneratePrimaryVertex(G4Event* evt)
{
  G4int NHEP = 0;  // number of entries
  if (inputFile.is_open())
  {
    inputFile >> NHEP;
  }
  else
  {
    G4Exception("G4HEPEvtInterface::G4HEPEvtInterface","Event0201",
                FatalException, "G4HEPEvtInterface:: cannot open file.");
  }
  if( inputFile.eof() ) 
  {
    G4Exception("G4HEPEvtInterface::GeneratePrimaryVertex", "Event0202",
                RunMustBeAborted,
                "End-Of-File: HEPEvt input file -- no more event to read!");
    return;
  }

  if(vLevel > 0)
  {
    G4cout << "G4HEPEvtInterface - reading " << NHEP
           << " HEPEvt particles from " << fileName << "." << G4endl;
  }
  for( G4int IHEP=0; IHEP<NHEP; ++IHEP )
  {
    G4int ISTHEP;   // status code
    G4int IDHEP;    // PDG code
    G4int JDAHEP1;  // first daughter
    G4int JDAHEP2;  // last daughter
    G4double PHEP1; // px in GeV
    G4double PHEP2; // py in GeV
    G4double PHEP3; // pz in GeV
    G4double PHEP5; // mass in GeV

    inputFile >> ISTHEP >> IDHEP >> JDAHEP1 >> JDAHEP2
              >> PHEP1 >> PHEP2 >> PHEP3 >> PHEP5;
    if( inputFile.eof() ) 
    {
      G4Exception("G4HEPEvtInterface::GeneratePrimaryVertex", "Event0203",
                  FatalException,
                  "Unexpected End-Of-File in the middle of an event");
    }
    if(vLevel > 1)
    {
      G4cout << " " << ISTHEP << " " << IDHEP << " " << JDAHEP1
             << " " << JDAHEP2 << " " << PHEP1 << " " << PHEP2
             << " " << PHEP3 << " " << PHEP5 << G4endl;
    }

    // Create G4PrimaryParticle object
    //
    auto* particle = new G4PrimaryParticle( IDHEP );
    particle->SetMass( PHEP5*GeV );
    particle->SetMomentum(PHEP1*GeV, PHEP2*GeV, PHEP3*GeV );

    // Create G4HEPEvtParticle object
    //
    auto* hepParticle
      = new G4HEPEvtParticle( particle, ISTHEP, JDAHEP1, JDAHEP2 );

    // Store
    //
    HPlist.push_back( hepParticle );
  }

  // Check if there is at least one particle
  //
  if( HPlist.empty() ) return; 

  // Make connection between daughter particles decayed from the same mother
  //
  for(auto & i : HPlist)
  {
    if( i->GetJDAHEP1() > 0 ) //  it has daughters
    {
      G4int jda1 = i->GetJDAHEP1()-1; // FORTRAN index starts from 1
      G4int jda2 = i->GetJDAHEP2()-1; // but C++ starts from 0.
      G4PrimaryParticle* mother = i->GetTheParticle();
      for( G4int j=jda1; j<=jda2; ++j )
      {
        G4PrimaryParticle* daughter = HPlist[j]->GetTheParticle();
        if(HPlist[j]->GetISTHEP()>0)
        {
          mother->SetDaughter( daughter );
          HPlist[j]->Done();
        }
      }
    }
  }

  // Create G4PrimaryVertex object
  //
  auto* vertex = new G4PrimaryVertex(particle_position,particle_time);

  // Put initial particles to the vertex
  //
  for(auto & ii : HPlist)
  {
    if( ii->GetISTHEP() > 0 ) // ISTHEP of daughters had been 
                                      // set to negative
    {
      G4PrimaryParticle* initialParticle = ii->GetTheParticle();
      vertex->SetPrimary( initialParticle );
    }
  }

  // Clear G4HEPEvtParticles
  //
  for(auto & iii : HPlist)
  {
    delete iii;
  }
  HPlist.clear();

  // Put the vertex to G4Event object
  //
  evt->AddPrimaryVertex( vertex );
}
