// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEPlot.hh,v 1.3 1999-12-15 14:52:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Gheisha friend class G4HEPlot  -- header file
// last modified: H. Fesefeldt 18-November-1996

#include <stdio.h>
#include "globals.hh"

class G4HEPlot
 {
  protected:
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

  public:

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

   void XScale( G4double a, G4double b, const G4HEPlot & p); 

   void Log( G4double s, const G4HEPlot & p);

   void Sqrt( G4double s, const G4HEPlot & p);

   void Reset();

   void Fill(G4double x, G4double w);

   void Print(G4String name, G4int i);

   void DumpToFile(G4int aPlot, G4String aName);

   void GetFromFile(G4int aPlot, G4String aName);
 
};

