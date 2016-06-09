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
#include "G4NeutronHPList.hh"
#include "G4HadronicException.hh"

  void G4NeutronHPList::Check(G4int i)
  {
    if(i<0) 
    {
      throw G4HadronicException(__FILE__, __LINE__, "G4NeutronHPList::Check(G4int) called with negative index");
    }
    if(i>nEntries) throw G4HadronicException(__FILE__, __LINE__, "Skipped some index numbers in G4NeutronHPList");
    if(i==nPoints)
    {
      nPoints = static_cast<G4int>(1.5*nPoints);
      G4double * buff = new G4double[nPoints];
      for (G4int j=0; j<nEntries; j++) buff[j] = theData[j];
      delete [] theData;
      theData = buff;
    }
    if(i==nEntries) nEntries++;
  }
  
  void G4NeutronHPList::Init(std::ifstream & aDataFile, G4int nPar, G4double unit)
  {
    G4int i;
    G4double y;
    for (i=0; i<nPar; i++)
    {
      aDataFile >> y;
      SetValue(i,y*unit);
    }
  }

  void G4NeutronHPList::Init(std::ifstream & aDataFile, G4double unit)
  {
    G4int total, i;
    aDataFile >> total;
    G4double y;
    for (i=0;i<total;i++)
    {
      aDataFile >>y;
      SetValue(i,y*unit);
    }
  }
  
  G4double G4NeutronHPList::GetValue(G4int i) 
  { 
//    G4cout << "TestList "<<i<<" "<<nEntries<<G4endl;
    if(nEntries<0)
    {
//      G4cout <<nPoints<<" "<<nEntries<<" "<<theData<<G4endl;
//      for(G4int ii=0; ii<2; ii++) G4cout << theData[ii]<<" ";
//      G4cout << G4endl;
    }
    if (i<0) i=0;
    if(i>=GetListLength()) i=GetListLength()-1;
    return theData[i];
  }
  
