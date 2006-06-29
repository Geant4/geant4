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
// $Id: OlapGenerator.cc,v 1.4 2006-06-29 17:22:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapGenerator
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapGenerator.hh"
#include "G4Event.hh"
#include "G4VisExtent.hh"
#include "G4PrimaryParticle.hh"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"

//#define OLAPDEBUG1

OlapGrid::OlapGrid()
 : count(3,0), gsize(3,3), axis(0), eventsPerRun(27)
{
  Reset(); 
  G4cout << "OlapGrid(): " << G4endl
       << "   gsize: " << gsize[0] << " " << gsize[1]
       << " " << gsize[2] << G4endl;    
}

void OlapGrid::Next()
{
    G4int out = (axis+1)%3;  // outer loop
    G4int in  = (axis+2)%3;  // inner loop
    if ( !(count[in] = (count[in]+1)%gsize[in]) )
    {
      if( !(count[out] = (count[out]+1)%gsize[out]) )
      {
        axis = (axis+1)%3;
      }
    }   
}

OlapGenerator::OlapGenerator() : autoinc(false)
{
   if (!ext.size())
   {
     ext.push_back(1);
     ext.push_back(2);
     ext.push_back(3);
   } 
   /*
   ext[0]=1.;
   ext[1]=2.;
   ext[2]=3.;
   */
}


OlapGenerator::~OlapGenerator()
{
}


void OlapGenerator::GeneratePrimaries(G4Event * anEvent)
{
  G4int out = (grid.axis+1)%3;  // outer loop
  G4int in  = (grid.axis+2)%3;  // inner loop

  posAB[grid.axis]=-ext[grid.axis]/2.;
  posAB[in]=-ext[in]/2.+ext[in]/G4double(grid.gsize[in]-1)
                               *G4double(grid.count[in]);
  posAB[out]=-ext[out]/2.+ext[out]/G4double(grid.gsize[out]-1)
                                  *G4double(grid.count[out]);
  
  posBA = posAB;
  posBA[grid.axis]=-posBA[grid.axis];
  
  G4ThreeVector dirAB, dirBA;
  dirAB[grid.axis] =  1;
  dirBA[grid.axis] = -1;
  
  #ifdef OLAPDEBUG1
    G4int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    G4int evtID = anEvent->GetEventID();
    G4cout << "generator: " 
           << "run=" << runID << " evt=" << evtID
           << " axis=" << grid.axis << " out=" << grid.count[out]
           << " in=" << grid.count[in] 
           << " posAB=" << posAB 
	   << " posBA=" << posBA << G4endl;
  #endif
  
  // now generator 2 geantinos flying in opposite direction from A->B and B->A
// ==========================
  G4ParticleDefinition* pd = G4Geantino::GeantinoDefinition();
  
  // create 2 opposite positioned vertices
  G4PrimaryVertex* vertexAB = 
    new G4PrimaryVertex(posAB,0);
  G4PrimaryVertex* vertexBA =
    new G4PrimaryVertex(posBA,0);
    
  // create one geantino per vertex as primary
  G4double mass =  pd->GetPDGMass();
  G4double energy = 1. + mass;
  G4double pmom = std::sqrt(energy*energy-mass*mass);
  G4double pxAB = pmom*dirAB.x();
  G4double pyAB = pmom*dirAB.y();
  G4double pzAB = pmom*dirAB.z();
  G4double pxBA = pmom*dirBA.x();
  G4double pyBA = pmom*dirBA.y();
  G4double pzBA = pmom*dirBA.z();
  
  G4PrimaryParticle* particleAB = new G4PrimaryParticle(pd,pxAB,pyAB,pzAB);
  G4PrimaryParticle* particleBA = new G4PrimaryParticle(pd,pxBA,pyBA,pzBA);
  particleAB->SetMass( mass );  
  particleBA->SetMass( mass );  
  particleAB->SetCharge(0.);
  particleBA->SetCharge(0.);
  vertexAB->SetPrimary( particleAB );
  vertexBA->SetPrimary( particleBA );
  anEvent->AddPrimaryVertex( vertexAB );
  anEvent->AddPrimaryVertex( vertexBA );
    
  // ok, set all counters to be ready for the next position  
  if (autoinc)
    grid.Next();     
  /*
  if ( !(count[in] = (count[in]+1)%grid[in]) ) {
    if( !(count[out] = (count[out]+1)%grid[out]) ) {
      axis = (axis+1)%3;
    }
  } 
  */   

}    


void OlapGenerator::SetExtent(const G4VisExtent & anExt)
{
   ext[0]=(anExt.GetXmax()-anExt.GetXmin());
   ext[1]=(anExt.GetYmax()-anExt.GetYmin());
   ext[2]=(anExt.GetZmax()-anExt.GetZmin());
   for (G4int i=0; i<3; i++)
     ext[i] += ext[i]/2000.;   // enlarge by 1%
   
   Reset();
     
}

void OlapGenerator::SetExtent(G4double e)
{
  ext[0]=e;
  ext[1]=e;
  ext[2]=e;
}


void OlapGenerator::Reset()
{
   grid.Reset();
   /*
   for (G4int i=0; i<3; i++) {
      count[i]=0;  
   }
   
   axis=0;
   */
}


void OlapGenerator::SetGrid(G4int x, G4int y, G4int z)
{
   if (x>2)
     grid.gsize[0]=x;
   else
     G4cerr << "Warning in OlapGenerator::SetGrid(): wrong x-grid parameter: " << x << G4endl;
     

   if (y>2)
     grid.gsize[1]=y;
   else
     G4cerr << "Warning in OlapGenerator::SetGrid(): wrong y-grid parameter: " << y << G4endl;

   if (z>2)
     grid.gsize[2]=z;
   else
     G4cerr << "Warning in OlapGenerator::SetGrid(): wrong z-grid parameter: " << z << G4endl;
   
   grid.eventsPerRun =  x*y + y*z + x*z;  
   Reset();
}
