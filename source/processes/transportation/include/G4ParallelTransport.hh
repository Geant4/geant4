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
// $Id: G4ParallelTransport.hh,v 1.10 2006/06/29 21:10:14 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4ParallelTransport
//
// Class description:
//
// Used internally by importance sampling and scoring in a "parallel"
// geometry.
// This process "moves" a particle in a "parallel" geometry.
// It should be placed in the post step process do-it vector after the 
// G4Transportation process.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelTransport_hh
#define G4ParallelTransport_hh G4ParallelTransport_hh

#include "G4Types.hh"

#include "G4VProcess.hh"

class G4VPGeoDriver;
class G4VParallelStepper;

class G4ParallelTransport : public G4VProcess
{

public:  // with description

  G4ParallelTransport(G4VPGeoDriver &pgeodriver, 
                      G4VParallelStepper &aStepper,
                      const G4String &aName = "ParallelTransport");
    // create G4ParticleChange

  virtual ~G4ParallelTransport();
    // delete G4ParticleChange

  // the post step pair of functions may be overwritten
  // in certain derived classes like the importance process

  virtual G4double 
  PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                       G4double   previousStepSize,
                                       G4ForceCondition* condition);
    // straight distance to the boundary in the "parallel" geometry
    // in the direction of the track   

  virtual G4VParticleChange *PostStepDoIt(const G4Track&, const G4Step&);
    // move a track in a "parallel" geometry
    // and update a G4PArallelStepper

  virtual void StartTracking(G4Track*);
  virtual void EndTracking();


  const G4VParallelStepper &GetPStepper() const;
    // get the fPStepper

public:  // without description
  
  //  no operation in  AtRestDoIt and  AlongStepDoIt

   virtual G4double 
   AlongStepGetPhysicalInteractionLength(const G4Track&,
                                       G4double  ,
                                       G4double  ,
                                       G4double& ,
                                       G4GPILSelection*);
   virtual G4double 
   AtRestGetPhysicalInteractionLength(const G4Track& ,
                                      G4ForceCondition*);

   virtual G4VParticleChange*  AtRestDoIt(const G4Track&, 
                                          const G4Step&);
   virtual G4VParticleChange* AlongStepDoIt(const G4Track&, 
                                            const G4Step&);
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
  G4bool fInitStep;

};

// ------------------------------------------------------------

inline const G4VParallelStepper &
G4ParallelTransport::GetPStepper() const
{
  return fPStepper;
}


#endif
