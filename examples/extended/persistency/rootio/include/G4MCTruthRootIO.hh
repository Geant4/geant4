// $Id: G4MCTruthRootIO.hh,v 1.2 2002-12-04 14:12:26 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4MCTruthRootIO.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef MCTRUTH_ROOT_IO_HH
#define MCTRUTH_ROOT_IO_HH 1

// Class inherited:
#include "G4VMCTruthIO.hh"

// Class Description:
//   Concreate class for storing and retrieving MCTruth events.

class G4MCTruthRootIO
 : public G4VMCTruthIO
{
    public: // With description
      G4MCTruthRootIO() {};
      // Constructor

      virtual ~G4MCTruthRootIO() {};
      // Destructor

    public: // With description
      virtual bool Store(G4MCTEvent* mctevent);
      // Method for storing MCTruth Event.

      virtual bool Retrieve(G4MCTEvent* & mctevent);
      // Method for retrieving MCTruth GenEvent.

      static G4MCTruthRootIO* GetMCTruthRootIO();
      // method #3: GetMCTruthRootIO()

    private:
      static G4MCTruthRootIO* thePointer;

}; // End of class G4MCTruthRootIO

#endif

