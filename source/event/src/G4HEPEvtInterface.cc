// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEPEvtInterface.cc,v 1.1 1999-01-07 16:06:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

#include "G4ios.hh"
#include "G4HEPEvtInterface.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4HEPEvtParticle.hh"
#include "G4Event.hh"

G4HEPEvtInterface::G4HEPEvtInterface(char* evfile)
{
  inputFile.open(evfile);
  fileName = evfile;
}

G4HEPEvtInterface::G4HEPEvtInterface(G4String evfile)
{
  const char* fn = evfile.data();
  inputFile.open((char*)fn);
  fileName = evfile;
}

G4HEPEvtInterface::~G4HEPEvtInterface()
{;}

void G4HEPEvtInterface::GeneratePrimaryVertex(G4Event* evt)
{
  int NHEP;  // number of entries
  inputFile >> NHEP;
  if( inputFile.eof() ) 
  {
    G4Exception("End-Of-File : HEPEvt input file");
    return;
  }

  for( int IHEP=0; IHEP<NHEP; IHEP++ )
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

    // create G4PrimaryParticle object
    G4PrimaryParticle* particle 
      = new G4PrimaryParticle( IDHEP, PHEP1*GeV, PHEP2*GeV, PHEP3*GeV );
    particle->SetMass( PHEP5*GeV );

    // create G4HEPEvtParticle object
    G4HEPEvtParticle* hepParticle
      = new G4HEPEvtParticle( particle, ISTHEP, JDAHEP1, JDAHEP2 );

    // Store
    HPlist.insert( hepParticle );
  }

  // check if there is at least one particle
  if( HPlist.entries() == 0 ) return; 

  // make connection between daughter particles decayed from 
  // the same mother
  for( int i=0; i<HPlist.entries(); i++ )
  {
    if( HPlist(i)->GetJDAHEP1() > 0 ) //  it has daughters
    {
      int jda1 = HPlist(i)->GetJDAHEP1()-1; // FORTRAN index starts from 1
      int jda2 = HPlist(i)->GetJDAHEP2()-1; // but C++ starts from 0.
      G4PrimaryParticle* mother = HPlist(i)->GetTheParticle();
      for( int j=jda1; j<=jda2; j++ )
      {
        G4PrimaryParticle* daughter = HPlist(j)->GetTheParticle();
        if(HPlist(j)->GetISTHEP()>0)
        {
          mother->SetDaughter( daughter );
          HPlist(j)->Done();
        }
      }
    }
  }

  // create G4PrimaryVertex object
  G4double x0 = 0.;  // vertex point X
  G4double y0 = 0.;  // vertex point Y
  G4double z0 = 0.;  // vertex point Z
  G4double t0 = 0.;  // vertex time
  G4PrimaryVertex* vertex = new G4PrimaryVertex(x0,y0,z0,t0);

  // put initial particles to the vertex
  for( int ii=0; ii<HPlist.entries(); ii++ )
  {
    if( HPlist(ii)->GetISTHEP() > 0 ) // ISTHEP of daughters had been 
                                       // set to negative
    {
      G4PrimaryParticle* initialParticle = HPlist(ii)->GetTheParticle();
      vertex->SetPrimary( initialParticle );
    }
  }

  // clear G4HEPEvtParticles
  HPlist.clearAndDestroy();

  // Put the vertex to G4Event object
  evt->AddPrimaryVertex( vertex );
}

