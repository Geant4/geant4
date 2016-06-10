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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPField.hh"
#include "G4HadronicException.hh"
#include "G4ios.hh"


  G4ParticleHPField::G4ParticleHPField()
  {
    theData = new G4ParticleHPFieldPoint[100]; 
    nPoints=100;
    nEntries=0;
    theData->SetData(0,0,0);
  }
  
  G4ParticleHPField::~G4ParticleHPField(){ delete [] theData;}
  
  G4double G4ParticleHPField::GetY(G4double e, G4int j)
  {
    G4int low   = 0;
    G4int high  = 0;
    G4int i;
    for (i=1; i<nEntries/10; i++)
    {
      if(theData[10*i].GetX()>e) break;
    }
    if(i==(nEntries/10))
    {
      i=10*i;
      while (i<nEntries) // Loop checking, 11.05.2015, T. Koi
      {
        if(theData[i++].GetX()>e) break;
      } 
      if (i==nEntries)
      {
        low  = nEntries-1;
        high = nEntries-2;
      }else{
        low = i-1;
        high = i;
      }
    }else{
      for (G4int jj=0; jj<10; jj++)
      {
        if(theData[i].GetX()<e) break;
        i--;
      }
      low = i;
      high = i+1;
    }
    G4double x1, x2, y1, y2, x, y;
    x = e;
    x1 = theData[low] .GetX();
    x2 = theData[high].GetX();
    y1 = theData[low] .GetY(j);
    y2 = theData[high].GetY(j);
    y = x*(y2-y1)/(x2-x1);
    return y += y2-x2*(y2-y1)/(x2-x1);
  }

  void G4ParticleHPField::Dump()
  {
    G4cout << nEntries<<G4endl;
    for(G4int i=0; i<nEntries; i++)
    {
      G4cout << theData[i].GetX()<<" ";
      for(G4int j=0; j<theData[i].GetDepth(); j++)
      {
        G4cout << theData[i].GetY(j)<<" ";
      }
      G4cout << G4endl;
    }
  }
  
  void G4ParticleHPField::Check(G4int i)
  {
    if(i>nEntries) throw G4HadronicException(__FILE__, __LINE__, "Skipped some index numbers in G4ParticleHPField");
    if(i==nPoints)
    {
      nPoints += 50;
      G4ParticleHPFieldPoint * buff = new G4ParticleHPFieldPoint[nPoints];
//      G4cout << "copying 1"<<G4endl;
      for (G4int j=0; j<nEntries; j++) 
      {
        buff[j] = theData[j];
      }
//      G4cout << "copying 2"<<G4endl;
      delete [] theData;
      theData = buff;
    }
    if(i==nEntries) nEntries=i+1;
  }
