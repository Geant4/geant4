// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VFastSimulationModel.hh,v 1.3 1999-12-15 14:53:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//
//  G4VFastSimulationModel.hh
//
//  Description:
//    Base class for fast simulation models.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------


#ifndef G4VFastSimulationModel_h
#define G4VFastSimulationModel_h

#include "G4FastTrack.hh"
#include "G4FastStep.hh"

//---------------------------
// For possible future needs:
//---------------------------
typedef G4LogicalVolume G4Envelope;

//-------------------------------------------
//
//        G4VFastSimulationModel class
//
//-------------------------------------------

// Class Description:
//   This is the abstract class for the implementation of parameterisations. 
//   You have to inherit from it to implement your concrete parameterisation 
//   model.
//

 class G4VFastSimulationModel 
{
 public: // With description

  G4VFastSimulationModel(const G4String& aName);
  // aName identifies the parameterisation model.

  G4VFastSimulationModel(const G4String& aName, G4Envelope*, 
			 G4bool IsUnique=FALSE);
  // This constructor allows you to get a quick "getting started".
  // In addition to the model name, this constructor accepts a G4LogicalVolume 
  // pointer. This volume will automatically becomes the envelope, and the 
  // needed G4FastSimulationManager object is constructed if necessary giving 
  // it the G4LogicalVolume pointer and the boolean value. If it already 
  // exists, the model is simply added to this manager. However the 
  // G4VFastSimulationModel object will not keep track of the envelope given 
  // in the constructor.
  // The boolean argument is there for optimization purpose: if you know that 
  // the G4LogicalVolume envelope is placed only once you can turn this 
  // boolean value to "true" (an automated mechanism is foreseen here.)

public: // Without description
  virtual ~G4VFastSimulationModel() {};

public: // With description

  virtual G4bool IsApplicable(const G4ParticleDefinition&) = 0;
  // In your implementation, you have to return "true" when your model is 
  // applicable to the G4ParticleDefinition passed to this method. The 
  // G4ParticleDefinition provides all intrisic particle informations (mass, 
  // charge, spin, name ...).

  virtual G4bool ModelTrigger(const G4FastTrack &) = 0;
  // You have to return "true" when the dynamics conditions to trigger your
  // parameterisation are fulfiled. The G4FastTrack provides you access to 
  // the current G4Track, gives simple access to envelope related features 
  // (G4LogicalVolume, G4VSolid, G4AffineTransform references between the 
  // global and the envelope local coordinates systems) and simple access to 
  // the position, momentum expressed in the envelope coordinate system. 
  // Using those quantities and the G4VSolid methods, you can for example 
  // easily check how far you are from the envelope boundary. 

  virtual void DoIt(const G4FastTrack&, G4FastStep&) = 0;
  // Your parameterisation properly said. The G4FastTrack reference provides 
  // input informations. The final state of the particles after parameterisation
  // has to be returned through the G4FastStep reference. This final state is 
  // described has "requests" the tracking will apply after your 
  // parameterisation has been invoked.

  // ---------------------------
  // -- Idem for AtRest methods:
  // ---------------------------
  // -- A default dummy implementation is provided.

  virtual 
  G4bool AtRestModelTrigger(const G4FastTrack& fastTrack) {return false;}
  // You have to return "true" when the dynamics conditions to trigger your
  // parameterisation are fulfiled. The G4FastTrack provides you access to 
  // the current G4Track, gives simple access to envelope related features 
  // (G4LogicalVolume, G4VSolid, G4AffineTransform references between the 
  // global and the envelope local coordinates systems) and simple access to 
  // the position, momentum expressed in the envelope coordinate system. 
  // Using those quantities and the G4VSolid methods, you can for example 
  // easily check how far you are from the envelope boundary. 

  virtual 
  void   AtRestDoIt  (const G4FastTrack& fastTrack, G4FastStep& fastStep) {}
  // Your parameterisation properly said. The G4FastTrack reference provides 
  // input informations. The final state of the particles after parameterisation
  // has to be returned through the G4FastStep reference. This final state is 
  // described has "requests" the tracking will apply after your 
  // parameterisation has been invoked.
  
public: // Without description

  // Useful public methods :
  const G4String GetName() const;
  G4bool operator == ( const G4VFastSimulationModel&) const;

private:
  //-------------
  // Model Name:
  //-------------
  G4String theModelName;
};

inline const G4String G4VFastSimulationModel::GetName() const 
{
  return theModelName;
}

inline G4bool 
G4VFastSimulationModel::operator == (const G4VFastSimulationModel& fsm) const
{
  return (this==&fsm) ? true : false;
}
#endif
