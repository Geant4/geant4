// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VXrayTRmodel.hh,v 1.1 2000-05-16 13:46:41 grichine Exp $
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


#ifndef G4VXrayTRmodel_h
#define G4VXrayTRmodel_h 1


#include "globals.hh"
#include "templates.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Gamma.hh"

#include "G4VFastSimulationModel.hh"
// #include "G4ForwardXrayTR.hh"

class G4VXrayTRmodel : public G4VFastSimulationModel
// , public G4ForwardXrayTR
{
public:

   G4VXrayTRmodel (G4LogicalVolume *anEnvelope,G4double,G4double);
   virtual  ~G4VXrayTRmodel ();

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

  G4Gamma* fPtrGamma ;  // pointer to TR photon

  G4double* fGammaCutInKineticEnergy ; // TR photon cut in energy array
  G4double  fGammaTkinCut ;            // Tkin cut of TR photon in current mat.

  G4PhysicsTable* fAngleDistrTable ;
  G4PhysicsTable* fEnergyDistrTable ;

  static G4PhysicsLogVector* fProtonEnergyVector ;


  static G4double fTheMinEnergyTR ;            //  static min TR energy
  static G4double fTheMaxEnergyTR ;            //  static max TR energy
         G4double fMinEnergyTR ;               //  min TR energy in material
         G4double fMaxEnergyTR ;               //  max TR energy in material
  static G4double  fTheMaxAngle ;               //  max theta of TR quanta
  static G4double  fTheMinAngle ;               //  max theta of TR quanta
         G4double  fMaxThetaTR ;               //  max theta of TR quanta
  static G4int          fBinTR ;               //  number of bins in TR vectors

  static G4double fMinProtonTkin ;             // min Tkin of proton in tables
  static G4double fMaxProtonTkin ;             // max Tkin of proton in tables
  static G4int    fTotBin        ;             // number of bins in log scale
         G4double fGamma         ;             // current Lorentz factor
         G4double fEnergy ;                    // energy and
         G4double fVarAngle ;                  // angle squared

  static G4double fPlasmaCof ;               // physical consts for plasma energy
  static G4double fCofTR ;

  G4double fSigma1, fSigma2 ;               // plasma energy Sq of matter1/2

  G4int fMatIndex1, fMatIndex2 ;

  G4int fPlateNumber ;

  G4double** fPlatePhotoAbsCof ;
  G4int      fPlateIntervalNumber ;
  G4double   fPlateThick ;
 
  G4double** fGasPhotoAbsCof ;
  G4int      fGasIntervalNumber ;
  G4double   fGasThick ;     


};

#endif













