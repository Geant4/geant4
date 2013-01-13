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
// $Id$
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------
//   removed thePhysicsTable           02 Aug. 1998 H.Kurashige
//   Modified DumpInfo                 15 Aug. 1998 H.Kurashige

#include "G4VProcess.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhysicsTable.hh"
#include "G4MaterialTable.hh"
#include "G4ElementTable.hh"
#include "G4ElementVector.hh"

G4VProcess::G4VProcess(const G4String& aName, G4ProcessType   aType )
                  : aProcessManager(0),
	            pParticleChange(0),
                    theNumberOfInteractionLengthLeft(-1.0),
                    currentInteractionLength(-1.0),
		    theInitialNumberOfInteractionLength(-1.0),
                    theProcessName(aName),
		    theProcessType(aType),
		    theProcessSubType(-1),
                    thePILfactor(1.0),
                    enableAtRestDoIt(true),
                    enableAlongStepDoIt(true),
                    enablePostStepDoIt(true),
                    verboseLevel(0)

{
  pParticleChange = &aParticleChange;
}

G4VProcess::~G4VProcess()
{
}

G4VProcess::G4VProcess(const G4VProcess& right)
          : aProcessManager(0),
	    pParticleChange(0),
            theNumberOfInteractionLengthLeft(-1.0),
            currentInteractionLength(-1.0),
	    theInitialNumberOfInteractionLength(-1.0),
            theProcessName(right.theProcessName),
            theProcessType(right.theProcessType),
	    theProcessSubType(right.theProcessSubType),
            thePILfactor(1.0),
            enableAtRestDoIt(right.enableAtRestDoIt),
            enableAlongStepDoIt(right.enableAlongStepDoIt),
            enablePostStepDoIt(right.enablePostStepDoIt),
            verboseLevel(right.verboseLevel)
{
}


void G4VProcess::ResetNumberOfInteractionLengthLeft()
{
  theNumberOfInteractionLengthLeft =  -std::log( G4UniformRand() );
  theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 
}

void G4VProcess::SubtractNumberOfInteractionLengthLeft(
                                  G4double previousStepSize )
{
  if (currentInteractionLength>0.0) {
    theNumberOfInteractionLengthLeft -= previousStepSize/currentInteractionLength;
    if(theNumberOfInteractionLengthLeft<0.) {
       theNumberOfInteractionLengthLeft=perMillion;
    }          

  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VProcess::SubtractNumberOfInteractionLengthLeft()";
      G4cerr << " [" << theProcessName << "]" <<G4endl;
      G4cerr << " currentInteractionLength = " << currentInteractionLength/cm << " [cm]";
      G4cerr << " previousStepSize = " << previousStepSize/cm << " [cm]";
      G4cerr << G4endl;
    }
#endif
    G4String msg = "Negative currentInteractionLength for ";
    msg += 	theProcessName;
    G4Exception("G4VProcess::SubtractNumberOfInteractionLengthLeft()",
		"ProcMan201",EventMustBeAborted,
		msg);
  }
}

void G4VProcess::StartTracking(G4Track*)
{
  currentInteractionLength = -1.0;
  theNumberOfInteractionLengthLeft = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4VProcess::StartTracking() [" << theProcessName << "]" <<G4endl;
  }
#endif
}

void G4VProcess::EndTracking()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4VProcess::EndTracking() [" << theProcessName << "]" <<G4endl;
  }
#endif
  theNumberOfInteractionLengthLeft = -1.0;
  currentInteractionLength = -1.0;
  theInitialNumberOfInteractionLength=-1.0;
}


const G4String& G4VProcess::GetProcessTypeName(G4ProcessType aType ) 
{
  static __thread G4String *typeNotDefined_G4MT_TLS_ = 0 ; if (!typeNotDefined_G4MT_TLS_) {typeNotDefined_G4MT_TLS_ = new  G4String  ; *typeNotDefined_G4MT_TLS_= "NotDefined" ; }  G4String &typeNotDefined = *typeNotDefined_G4MT_TLS_;
  static __thread G4String *typeTransportation_G4MT_TLS_ = 0 ; if (!typeTransportation_G4MT_TLS_) {typeTransportation_G4MT_TLS_ = new  G4String  ; *typeTransportation_G4MT_TLS_= "Transportation" ; }  G4String &typeTransportation = *typeTransportation_G4MT_TLS_;
  static __thread G4String *typeElectromagnetic_G4MT_TLS_ = 0 ; if (!typeElectromagnetic_G4MT_TLS_) {typeElectromagnetic_G4MT_TLS_ = new  G4String  ; *typeElectromagnetic_G4MT_TLS_= "Electromagnetic" ; }  G4String &typeElectromagnetic = *typeElectromagnetic_G4MT_TLS_;
  static __thread G4String *typeOptical_G4MT_TLS_ = 0 ; if (!typeOptical_G4MT_TLS_) {typeOptical_G4MT_TLS_ = new  G4String  ; *typeOptical_G4MT_TLS_= "Optical" ; }  G4String &typeOptical = *typeOptical_G4MT_TLS_;
  static __thread G4String *typeHadronic_G4MT_TLS_ = 0 ; if (!typeHadronic_G4MT_TLS_) {typeHadronic_G4MT_TLS_ = new  G4String  ; *typeHadronic_G4MT_TLS_= "Hadronic" ; }  G4String &typeHadronic = *typeHadronic_G4MT_TLS_;
  static __thread G4String *typePhotolepton_hadron_G4MT_TLS_ = 0 ; if (!typePhotolepton_hadron_G4MT_TLS_) {typePhotolepton_hadron_G4MT_TLS_ = new  G4String  ; *typePhotolepton_hadron_G4MT_TLS_= "Photolepton_hadron" ; }  G4String &typePhotolepton_hadron = *typePhotolepton_hadron_G4MT_TLS_;
  static __thread G4String *typeDecay_G4MT_TLS_ = 0 ; if (!typeDecay_G4MT_TLS_) {typeDecay_G4MT_TLS_ = new  G4String  ; *typeDecay_G4MT_TLS_= "Decay" ; }  G4String &typeDecay = *typeDecay_G4MT_TLS_;
  static __thread G4String *typeGeneral_G4MT_TLS_ = 0 ; if (!typeGeneral_G4MT_TLS_) {typeGeneral_G4MT_TLS_ = new  G4String  ; *typeGeneral_G4MT_TLS_= "General" ; }  G4String &typeGeneral = *typeGeneral_G4MT_TLS_;
  static __thread G4String *typeParameterisation_G4MT_TLS_ = 0 ; if (!typeParameterisation_G4MT_TLS_) {typeParameterisation_G4MT_TLS_ = new  G4String  ; *typeParameterisation_G4MT_TLS_= "Parameterisation" ; }  G4String &typeParameterisation = *typeParameterisation_G4MT_TLS_;
  static __thread G4String *typeUserDefined_G4MT_TLS_ = 0 ; if (!typeUserDefined_G4MT_TLS_) {typeUserDefined_G4MT_TLS_ = new  G4String  ; *typeUserDefined_G4MT_TLS_= "UserDefined" ; }  G4String &typeUserDefined = *typeUserDefined_G4MT_TLS_;
  static __thread G4String *noType_G4MT_TLS_ = 0 ; if (!noType_G4MT_TLS_) {noType_G4MT_TLS_ = new  G4String  ; *noType_G4MT_TLS_= "------" ; }  G4String &noType = *noType_G4MT_TLS_;   // Do not modify this !!!!

  if (aType ==   fNotDefined) {
    return  typeNotDefined;
  } else if  (aType ==   fTransportation ) {
    return typeTransportation;
  } else if  (aType ==   fElectromagnetic ) {
    return typeElectromagnetic;
  } else if  (aType ==   fOptical ) {
    return typeOptical;
  } else if  (aType ==   fHadronic ) {
    return typeHadronic;
  } else if  (aType ==   fPhotolepton_hadron ) {
    return typePhotolepton_hadron;
  } else if  (aType ==   fDecay ) {
    return typeDecay;
  } else if  (aType ==   fGeneral ) {
    return typeGeneral;
  } else if  (aType ==   fParameterisation ) {
    return typeParameterisation;
  } else if  (aType ==   fUserDefined ) {
    return typeUserDefined;
  } else {
    return noType;  
  }
}

G4VProcess & G4VProcess::operator=(const G4VProcess &)
{
  G4Exception("G4VProcess::operator=","ProcMan101",
	      JustWarning,"Assignment operator is called but NO effect");
  return *this;
}

G4int G4VProcess::operator==(const G4VProcess &right) const
{
  return (this == &right);
}

G4int G4VProcess::operator!=(const G4VProcess &right) const
{
  return (this !=  &right);
}

void G4VProcess::DumpInfo() const
{
  G4cout << "Process Name " << theProcessName ;
  G4cout << " : Type[" << GetProcessTypeName(theProcessType) << "]";
  G4cout << " : SubType[" << theProcessSubType << "]"<< G4endl;
}


const G4String&  G4VProcess::GetPhysicsTableFileName(const G4ParticleDefinition* particle,
						     const G4String& directory,
						     const G4String& tableName,
						     G4bool ascii)
{
  G4String thePhysicsTableFileExt;
  if (ascii) thePhysicsTableFileExt = ".asc";
  else       thePhysicsTableFileExt = ".dat";

  thePhysicsTableFileName = directory + "/";
  thePhysicsTableFileName += tableName + "." +  theProcessName + ".";
  thePhysicsTableFileName += particle->GetParticleName() + thePhysicsTableFileExt;
  
  return thePhysicsTableFileName;
}
