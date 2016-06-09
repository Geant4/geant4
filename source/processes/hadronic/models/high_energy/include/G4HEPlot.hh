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
// $Id: G4HEPlot.hh,v 1.9 2006/06/29 20:29:37 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//
// G4 Gheisha friend class G4HEPlot  -- header file
// last modified: H. Fesefeldt 18-November-1996

#include <stdio.h>
#include "globals.hh"

class G4HEPlot
 {
  public:
     G4double Xstart;
     G4double Xbin;
     G4double Weight;
     G4int    Nbin;
     G4int Entries;
     G4int EntriesOverflow;
     G4int EntriesUnderflow;
     G4double WeightOverflow;
     G4double WeightUnderflow;
     G4double* Xvalue;
     G4double* Yvalue;

   G4HEPlot(){ };

  ~G4HEPlot(){ };

   G4double  getWeight()
     { return Weight;}; 

   G4double getWeightUnderflow()
     { return WeightUnderflow; };

   G4double getWeightOverflow()
     { return WeightOverflow;};

   G4int getNumberOfEntries()
     { return Entries;};

   G4int getNumberOfEntriesOverflow()
     { return EntriesOverflow;};

   G4int getNumberOfEntriesUnderflow()
     { return EntriesUnderflow;};

   G4int getNumberOfBins()
     { return Nbin;}; 

   G4double getBinsize()
     { return Xbin;};

   G4double getXstart()
     { return Xstart;}; 

   void Init(G4int nbin, G4double xstart, G4double Xbin); 

   void Add( G4double s1, G4double s2, const G4HEPlot & p1, const G4HEPlot & p2 );

   void Multiply( G4double s1, G4double s2, const G4HEPlot & p1, const G4HEPlot & p2 );

   void Divide( G4double s1, G4double s2, const G4HEPlot & p1, const G4HEPlot & p2 );

   void Scale( G4double s, const G4HEPlot & p);

   void XScale( G4double a, G4double b); 

   void Shift( G4int nshift);

   void Log( G4double s, const G4HEPlot & p);

   void Sqrt( G4double s, const G4HEPlot & p);

   void Reset();

   void LinearFit(G4double& a, G4double& b);

   void Fill(G4double x, G4double w);

   void Print(G4String name, G4int i);

   void DumpToFile(G4int aPlot, G4String aName);

   void GetFromFile(G4int aPlot, G4String aName);
 
};

