// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4Data.hh
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 1 Giugno 1999
//
//      Modifications: 
//   
//
// Class description:
// Utility for Low Energy electromagnetic e/photon processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4DATA_HH
#define G4DATA_HH

#include "G4DataVector.hh"

class G4Data : public G4DataVector{


public:

 ~G4Data();

  G4bool operator == (const G4Data& ) const;

  G4bool operator < (const G4Data&) const;

};

#endif






