// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: tst2SimpleReporter.hh,v 1.1 1999-06-17 04:43:05 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
#ifndef tst2SimpleReporter_h
#define tst2SimpleReporter_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "tst2VParticleReporter.hh"

class tst2SimpleReporter: public tst2VParticleReporter
{
 public:
  //constructors
    tst2SimpleReporter();

  //destructor
    virtual ~tst2SimpleReporter();

 public:
	virtual void Print(const tst2ParticleContainer& container, 
                       const G4String& option="");
};


#endif
