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

