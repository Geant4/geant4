// $Id: G4VMCTruthIO.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VMCTruthIO.hh
//
// History:
//   '01.11.18  Youhei Morita  Initial creation

#ifndef V_M_C_TRUTH_I_O_HH
#define V_M_C_TRUTH_I_O_HH 1

#include "G4MCTEvent.hh"

// Class Description:
//   Abstract base class for storing and retrieving MCTruth events.

class G4VMCTruthIO
{
    public: // With description
      G4VMCTruthIO();
      // Constructor

      virtual ~G4VMCTruthIO() {};
      // Destructor

    public: // With description
      virtual G4bool Store(G4MCTEvent*) =0;
      // Pure virtual method for storing MCTruth Event.
      // Each persistency package should implement a concrete method
      // of storing the G4MCTEvent with this signature.

      virtual G4bool Retrieve(G4MCTEvent*&) =0;
      // Pure virtual method for retrieving MCTruth Event.
      // Each persistency package should implement a concrete method
      // of storing the G4MCTEvent with this signature.

      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

    protected:
      G4int m_verbose;

}; // End of class G4VMCTruthIO

#endif

