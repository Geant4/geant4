// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VFastSimulationModel.hh,v 1.1 1999-01-07 16:14:04 gunter Exp $
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
class G4VFastSimulationModel
{
public:
  //------------------------
  // Constructor/Destructor
  //------------------------
  //  constructors require at least the model name.
  G4VFastSimulationModel(const G4String& aName);

  // About the IsUnique flag, by default the envelope can
  // be placed n-Times. If the user is sure that it'll be 
  // placed just one time the IsUnique flag should be set 
  // TRUE to avoid the G4AffineTransform re-calculations each 
  // time we reach the envelope.
  G4VFastSimulationModel(const G4String& aName, G4Envelope*, 
			 G4bool IsUnique=FALSE);

  ~G4VFastSimulationModel() {};

  //
  //----------------------------------------------
  // Interface IsApplicable virtual method for the 
  // G4FastSimulationManager
  //----------------------------------------------
  virtual G4bool IsApplicable(const G4ParticleDefinition&) = 0;

  //
  //----------------------------------------------
  // Interface trigger virtual method for the 
  // G4FastSimulationManager
  //----------------------------------------------
  virtual G4bool ModelTrigger(const G4FastTrack &) = 0;

  // -- Pure virtual DoIt() method : the FastSimulationModel 
  //    properly said. The user has to compute the detector 
  //    response and to manage the particle changes state:
  // --    * change position/momentum or kill the primary
  // --    * create secondaries
  virtual void DoIt(const G4FastTrack&, G4FastStep&) = 0;

  // ---------------------------
  // -- Idem for AtRest methods:
  // ---------------------------
  // -- A default dummy implementation is provided.
  virtual 
  G4bool AtRestModelTrigger(const G4FastTrack& fastTrack) {return false;}
  virtual 
  void   AtRestDoIt  (const G4FastTrack& fastTrack, G4FastStep& fastStep) {}
  
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
