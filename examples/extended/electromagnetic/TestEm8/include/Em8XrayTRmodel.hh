// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8XrayTRmodel.hh,v 1.1 2000-02-09 10:46:42 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// base class for 'fast' parametrisation model describing X-ray transition
// created in some G4Envelope. Anglur distribuiton is very rough !!! (see DoIt
// method
// 
// History:
// 23.01.00 V. Grichine first version 
// 09.02.00 V. Grichine, DoIt was transformed from virtual
//


#ifndef Em8XrayTRmodel_h
#define Em8XrayTRmodel_h 1

#include "G4VFastSimulationModel.hh"
#include "G4ForwardXrayTR.hh"

class Em8XrayTRmodel : public G4VFastSimulationModel, public G4ForwardXrayTR
{
public:

   Em8XrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
   virtual  ~Em8XrayTRmodel ();

  // Pure virtual functions from base class

  G4bool IsApplicable(const G4ParticleDefinition&);
 
  G4bool ModelTrigger(const G4FastTrack &);

  // Pure virtuals must be implemented in inherited particular TR radiators
 
  void DoIt(const G4FastTrack&, G4FastStep&)  ;

  virtual  G4double GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle ) = 0  ;

  G4double OneBoundaryXTRNdensity( G4double energy,
                                   G4double gamma,
				   G4double varAngle ) const ;

  void BuildTable() ;
  void BuildEnergyTable() ;
  void BuildAngleTable() ;

  // for photon energy distribution tables

  G4double XTRNSpectralAngleDensity(G4double varAngle) ;
  G4double XTRNSpectralDensity(G4double energy) ;
  
  // for photon angle distribution tables

  G4double XTRNAngleSpectralDensity(G4double energy) ;
  G4double XTRNAngleDensity(G4double varAngle) ;

  void GetNumberOfPhotons() ;  

  void ExampleDoIt(const G4FastTrack&, G4FastStep&) ;

  // Auxiliary functions for plate/gas material parameters

  G4double GetPlateFormationZone(G4double,G4double,G4double) ;
  void     ComputePlatePhotoAbsCof() ;
  G4double GetPlateLinearPhotoAbs(G4double) ;
  void     GetPlateZmuProduct() ;
  G4double GetPlateZmuProduct(G4double,G4double,G4double) ;

  G4double GetGasFormationZone(G4double,G4double,G4double) ;
  void     ComputeGasPhotoAbsCof() ;
  G4double GetGasLinearPhotoAbs(G4double) ;
  void     GetGasZmuProduct() ;
  G4double GetGasZmuProduct(G4double,G4double,G4double) ;


protected:

  G4int fPlateNumber ;

  G4double** fPlatePhotoAbsCof ;
  G4int      fPlateIntervalNumber ;
  G4double   fPlateThick ;
 
  G4double** fGasPhotoAbsCof ;
  G4int      fGasIntervalNumber ;
  G4double   fGasThick ;     


};

#endif
