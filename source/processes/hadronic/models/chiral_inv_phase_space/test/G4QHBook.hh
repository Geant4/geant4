// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QHBook.hh,v 1.2 2000-09-11 09:03:53 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QHBook_h
#define G4QHBook_h 1

// ----------------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QHBook ----------------
//      by Mikhail Kossov and P. Degtiarenko, Nov 1999.
//      class for booking and filling histograms and ntuples
//      for the main CHIPStest routine
// ----------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "G4QHadronVector.hh"

class G4QHBook
{
public:
  G4QHBook();                                     // Default Constructor
  ~G4QHBook();                                    // Destructor

  // Specific Modifiers
  void  FillEvt(const G4QHadronVector* hadrons);  // Fill Histos & ntuples for the event

private: 
  G4int       nEvnt; // Consecutive number of call to fill the histograms 
  ofstream histNevt; // 1D histogram, id=1, to store number of calls to fill
  ofstream tuplEvtA; // Ntuple, id=25, to be filled once per event (All particles)
  ofstream tuplEvtQ; // Ntuple, id=27, to be filled once per event (Quasmon particles)
  ofstream tuplIncl; // Ntuple, id=20, to be filled once per track
  ofstream tuple3pi; // Ntuple, id=22, to be filled once per event with 3 pions
};

#endif
