// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleChangeForTransport.hh,v 1.2 1999-11-07 16:32:01 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// 
// ------------------------------------------------------------
//   Implemented for the new scheme                 10 May. 1998  H.Kurahige
//
// Class Description
//  This class is a concrete class for ParticleChange for transportation
//        
#ifndef G4ParticleChangeForTransport_h
#define G4ParticleChangeForTransport_h 1

#include "globals.hh"
#include "G4ios.hh"
class G4VTouchable;
#include "G4ParticleChange.hh"

class G4ParticleChangeForTransport: public G4ParticleChange
{ 
  public:
    // default constructor
    G4ParticleChangeForTransport();

    // destructor
    virtual ~G4ParticleChangeForTransport();

  protected:
    // hide copy constructor and assignment operaor as protected
    G4ParticleChangeForTransport(const G4ParticleChangeForTransport &right);
    G4ParticleChangeForTransport & operator=(const G4ParticleChangeForTransport &right);

  public: // with description
    // ----------------------------------------------------
    // --- the following methods are for updating G4Step -----   
    // Return the pointer to the G4Step after updating the Step information
    // by using final state information of the track given by a physics
    // process    
    virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
    virtual G4Step* UpdateStepForAtRest(G4Step* Step);
    virtual G4Step* UpdateStepForPostStep(G4Step* Step);
    // A physics process gives the final state of the particle 
    // based on information of G4Track (or equivalently the PreStepPoint)
 
    virtual void Initialize(const G4Track&);
    // Initialize all propoerties by using G4Track information
           
    // ----------------------------------------------------
    //--- methods to keep information of the final state--
    //  IMPORTANT NOTE: Although the name of the class and methods are
    //   "Change", what it stores (and returns in get) are the "FINAL" 
    //   values of the Position, Momentum, etc.

    G4VTouchable* GetTouchableChange() const;
    void  SetTouchableChange(G4VTouchable* fTouchable);
    //  Get/Set the touchable of the current particle.
    //  Note: Touchable in PostStepPoint will be updated only after PostStepDoIt

    G4bool GetMomentumChanged() const;
    void SetMomentumChanged(G4bool b);

  public:
    virtual void DumpInfo() const;

  protected:
    G4VTouchable* theTouchableChange;
    //  The changed touchable of a given particle.

  private:
    G4bool isMomentumChanged;

};

#include "G4ParticleChangeForTransport.icc"

#endif
















