// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhotoClusterModel.hh,v 1.1 2000-11-14 16:06:17 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Class for 'fast' parametrisation model describing X-ray ionisation clusters
// created in some G4Envelope. 
// 
// History:
// 14.07.00 V. Grichine first version 
//


#ifndef G4PhotoClusterModel_h
#define G4PhotoClusterModel_h 1


#include "globals.hh"
#include "templates.hh"
#include "G4PAIonisation.hh"
#include "G4VClusterModel.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include <g4rw/tvordvec.h>


class G4PhotoClusterModel : public G4VClusterModel
{
public:

   G4PhotoClusterModel (G4LogicalVolume* anEnvelope);

  ~G4PhotoClusterModel ();

  // Pure virtual functions from base class

  G4bool IsApplicable(const G4ParticleDefinition&);
 
  G4bool ModelTrigger(const G4FastTrack &);
 
  void DoIt(const G4FastTrack&, G4FastStep&)  ;

protected:

  void     ComputePhotoAbsCof() ;
  G4double GetLinearPhotoAbs(G4double) ;



private:

  G4int      fMatIndex         ;
  G4double** fPhotoAbsCof    ;
  G4int      fIntervalNumber ;
  
};

#endif





