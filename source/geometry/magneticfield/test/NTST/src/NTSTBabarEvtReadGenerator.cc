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
// $Id: NTSTBabarEvtReadGenerator.cc,v 1.2 2003-12-09 15:35:28 gunter Exp $
//
// -- Bogus -- BaBar Object-Oriented Geant-based Unified Simulation
//
// NTSTBabarEvtReadGenerator
//
// Description:
//   cloned from BabarEvtReadGenerator, it interpets the information in the
//   stdhep ascii files produced by GenFwkInt and passes it along to G4
//   (both vertex and 4-vector information). 
//
//   This class is used in the Bogus standalone application.
//
// Author List:
//   Bill Lockman
//
// Modification History:
//
//-----------------------------------------------------------------------------

// ====== C/C++ headers ======
// #include <sstream>
// #include <iostream>
// #include <iomanip>

// -- Older C++ headers
#include <strstream.h>
#include <iostream.h>
#include <iomanip.h>

// ====== This class header ======
#include "NTSTBabarEvtReadGenerator.hh"

// ====== Collaborating classes ======
#include "G4ios.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"

NTSTBabarEvtReadGenerator::NTSTBabarEvtReadGenerator(const char* evfile)
  : fileName(evfile){
  inputFile.open(evfile);
  if (inputFile == 0){
    G4cerr << "NTSTBabarEvtReadGenerator:: cannot open file " << fileName << G4endl;
    abort();
  }
}

NTSTBabarEvtReadGenerator::~NTSTBabarEvtReadGenerator()
{
  inputFile.close();
}

void NTSTBabarEvtReadGenerator::GeneratePrimaryVertex(G4Event* anEvent){
  G4int nevhep; // event number
  G4int nhep;  // number of entries
  
  inputFile >> nevhep >> nhep;
  if( inputFile.eof() ) {
    G4cerr << "End-Of-File : BabarEvt input file" 
		  << fileName << G4endl;
    abort();
  }
  
  for( G4int ihep=0; ihep<nhep; ihep++ ){
    G4int isthep;   // status code
    G4int idhep;    // PDG code
    G4int moms[2] = {0,0};        // inital and final mother
    G4int daut[2] = {0,0};        // inital and final daughters
    G4double p[5] = {0,0,0,0,0};  // px, py, pz, E, m
    G4double v[4] = {0,0,0,0};    // x,  y,  z,  t
    
    inputFile >> isthep >> idhep >> moms[0] >> daut[0] >> moms[1] >> daut[1];
    inputFile >> p[0] >> p[1] >> p[2] >> p[3] >> p[4];
    inputFile >> v[0] >> v[1] >> v[2] >> v[3];
    
    if (isthep == 1) { // stable particle
      
      // create primary vertex - distance units are mm
      G4PrimaryVertex* primVtx = 
	new G4PrimaryVertex(v[0]*mm, v[1]*mm, v[2]*mm, v[3]*ns);
      
      // create primary particle
      G4PrimaryParticle* primPart = 
	new G4PrimaryParticle(idhep, p[0]*GeV, p[1]*GeV, p[2]*GeV);
      
      primPart->SetMass(p[4]*GeV);
      
      // no need to set daughter since we are keeping only the stable particles

      // add it to the primary vertex (ownership is passing to the vertex)
      primVtx->SetPrimary( primPart );
      
      // add primary vertex to the event
      anEvent->AddPrimaryVertex( primVtx );
    }
  }
  //  for (G4int ivtx=0; ivtx < anEvent->GetNumberOfPrimaryVertex(); ivtx++){
  //    anEvent->GetPrimaryVertex(ivtx)->Print();
  //  }
}





