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
//
//
// $Id: G4ParallelTransport.hh,v 1.2 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelTransport
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelTransport_hh
#define G4ParallelTransport_hh

#include "g4std/strstream"

#include "G4VProcess.hh"
#include "globals.hh"

class G4VPGeoDriver;
class G4VParallelStepper;

class G4ParallelTransport : public G4VProcess
{

public:  // with description

  G4ParallelTransport(G4VPGeoDriver &pgeodriver, 
		      G4VParallelStepper &aStepper,
		      const G4String &aName = "ParallelTransport");
  virtual ~G4ParallelTransport();

  // the post step pair of functions may be overwritten
  // in certain derived classes like the importance process

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				       G4double   previousStepSize,
				       G4ForceCondition* condition);
  virtual G4VParticleChange *  PostStepDoIt(const G4Track&, const G4Step&);

  G4VParallelStepper &GetPStepper(){return fPStepper;}

public:  // without description
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt

  G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                        G4double  ,
                                        G4double  ,
                                        G4double& ,
                                        G4GPILSelection*) {return -1.0;}
  G4double AtRestGetPhysicalInteractionLength(const G4Track& ,
                                        G4ForceCondition*) {return -1.0;}
  
  G4VParticleChange*  AtRestDoIt(const G4Track&, const G4Step&) {return 0;}
  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) {return 0;}

protected:

  virtual void Error(const G4String &m);
  virtual void Warning(const G4String &m);

  G4ParticleChange *fParticleChange;


private:

  G4ParallelTransport(const G4ParallelTransport &);
  G4ParallelTransport &operator=(const G4ParallelTransport &);

private:

  G4VPGeoDriver &fPgeodriver;
  G4VParallelStepper &fPStepper;
  G4bool fCrossBoundary;

};

#endif
