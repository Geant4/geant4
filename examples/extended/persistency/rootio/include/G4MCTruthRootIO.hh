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

