// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VProcess.cc,v 1.1 1999-01-07 16:13:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------
//   removed thePhysicsTable           02 Aug. 1998 H.Kurashige
//   Modified DumpInfo                 15 Aug. 1998 H.Kurashige

#include "G4PhysicsTable.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"
#include "G4VProcess.hh"

//G4VProcess::G4VProcess()
//{
//  G4Exception("G4VProcess:: default constructor is called");
//}

G4VProcess::G4VProcess(const G4String& aName, G4ProcessType   aType )
                  : theProcessName(aName),
		    theProcessType(aType),
		    pParticleChange(NULL),
                    theNumberOfInteractionLengthLeft(-1.0),
                    currentInteractionLength(-1.0),
                    verboseLevel(0)
{
  pParticleChange = &aParticleChange;
}

G4VProcess::~G4VProcess()
{
}

G4VProcess::G4VProcess(G4VProcess& right):
	    pParticleChange(NULL),
            theProcessName(right.theProcessName),
            theProcessType(right.theProcessType),
            theNumberOfInteractionLengthLeft(-1.0),
            currentInteractionLength(-1.0)
{
}


void G4VProcess::SubtractNumberOfInteractionLengthLeft(
                                  G4double previousStepSize )
{
  if (currentInteractionLength>0.0) {
    theNumberOfInteractionLengthLeft -= previousStepSize/currentInteractionLength;
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" <<endl;
      G4cerr << " currentInteractionLength = " << currentInteractionLength/cm << " [cm]";
      G4cerr << " previousStepSize = " << previousStepSize/cm << " [cm]";
      G4cerr << endl;
    }
#endif
    G4Exception("G4VProcess::SubtractNumberOfInteractionLengthLeft()  negative currentInteractionLength" );
  }
}

void G4VProcess::StartTracking()
{
  currentInteractionLength = -1.0;
  theNumberOfInteractionLengthLeft = -1.0;
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cerr << "G4VProcess::StartTracking() [" << theProcessName << "]" <<endl;
  }
#endif
}

void G4VProcess::EndTracking()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cerr << "G4VProcess::EndTracking() [" << theProcessName << "]" <<endl;
  }
#endif
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
}


G4String G4VProcess::GetProcessTypeName(G4ProcessType aType ) 
{
  if (aType ==   fNotDefined) {
    return "NotDefined";
  } else if  (aType ==   fTransportation ) {
    return "Transportation";
  } else if  (aType ==   fElectromagnetic ) {
    return "Electromagnetic";
  } else if  (aType ==   fOptical ) {
    return "Optical";
  } else if  (aType ==   fHadronic ) {
    return "Hadronic";
  } else if  (aType ==   fPhotolepton_hadron ) {
    return "Photolepton_hadron";
  } else if  (aType ==   fDecay ) {
    return "Decay";
  } else if  (aType ==   fGeneral ) {
    return "General";
  } else if  (aType ==   fParameterisation ) {
    return "Parameterisation";
  } else if  (aType ==   fUserDefined ) {
    return "UserDefined";
  } else {
    return "------";   // Do not modify this end mark !!!!       
  }
}

G4VProcess & G4VProcess::operator=(const G4VProcess &)
{
  G4Exception("G4VProcess::assignment operator is called");
  return *this;
}

G4int G4VProcess::operator==(const G4VProcess &right) const
{
  return (this == (G4VProcess *) &right);
}

G4int G4VProcess::operator!=(const G4VProcess &right) const
{
  return (this != (G4VProcess *) &right);
}

void G4VProcess::DumpInfo() const
{
  G4cout << "Process Name " << theProcessName ;   
  G4cout << " : Type[" << GetProcessTypeName(theProcessType) << "]"<< endl;
}












