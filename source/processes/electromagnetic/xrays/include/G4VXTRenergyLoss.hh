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
// $Id: G4VXTRenergyLoss.hh,v 1.1 2002-01-15 16:48:52 grichine Exp $
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
// 15.01.02 V. Grichine first version 
//


#ifndef G4VXTRenergyLoss_h
#define G4VXTRenergyLoss_h 1


#include "globals.hh"
#include "templates.hh"
#include "g4std/complex"
#include "Randomize.hh"

#include "G4LogicalVolume.hh"

#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Gamma.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VContinuousProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Integrator.hh"


class G4VXTRenergyLoss : public G4VContinuousProcess
{
public:

   G4VXTRenergyLoss (G4LogicalVolume *anEnvelope,G4double,G4double,
                     const G4String & processName = "XTRenergyLoss");
   virtual  ~G4VXTRenergyLoss ();

  // Pure virtuals must be implemented in inherited particular TR radiators
 
  virtual  G4double GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle ) = 0  ;


  G4bool IsApplicable(const G4ParticleDefinition&);

  G4double GetContinuousStepLimit(const G4Track& aTrack,
					G4double  ,
					G4double  ,
                                        G4double& );
        // Returns the continuous step limit defined by the XTR process.

  G4VParticleChange* AlongStepDoIt(const G4Track& aTrack, 
				   const G4Step&  aStep);
 
  void BuildTable() ;
  void BuildEnergyTable() ;
  void BuildAngleTable() ;

  G4complex OneInterfaceXTRdEdx( G4double energy, 
                                G4double gamma,
                                G4double varAngle ) ;

  G4double SpectralAngleXTRdEdx(G4double varAngle) ;

  G4double SpectralXTRdEdx(G4double energy) ;

  G4double AngleSpectralXTRdEdx(G4double energy) ;

  G4double AngleXTRdEdx(G4double varAngle) ;


  /////////////////////////////////////////////////////////////

  G4double OneBoundaryXTRNdensity( G4double energy,
                                   G4double gamma,
				   G4double varAngle ) const ;


  // for photon energy distribution tables

  G4double XTRNSpectralAngleDensity(G4double varAngle) ;
  G4double XTRNSpectralDensity(G4double energy) ;
  
  // for photon angle distribution tables

  G4double XTRNAngleSpectralDensity(G4double energy) ;
  G4double XTRNAngleDensity(G4double varAngle) ;

  void GetNumberOfPhotons() ;  

  // Auxiliary functions for plate/gas material parameters

  G4double  GetPlateFormationZone(G4double,G4double,G4double) ;
  G4complex GetPlateComplexFZ(G4double,G4double,G4double) ;
  void      ComputePlatePhotoAbsCof() ;
  G4double  GetPlateLinearPhotoAbs(G4double) ;
  void      GetPlateZmuProduct() ;
  G4double  GetPlateZmuProduct(G4double,G4double,G4double) ;

  G4double  GetGasFormationZone(G4double,G4double,G4double) ;
  G4complex GetGasComplexFZ(G4double,G4double,G4double) ;
  void      ComputeGasPhotoAbsCof() ;
  G4double  GetGasLinearPhotoAbs(G4double) ;
  void      GetGasZmuProduct() ;
  G4double  GetGasZmuProduct(G4double,G4double,G4double) ;

  G4double GetXTRrandomEnergy( G4double scaledTkin, G4int iTkin ) ;
  G4double GetXTRenergy( G4int iPlace, G4double position, G4int iTransfer  ) ;

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
  static G4double fTheMaxAngle ;               //  max theta of TR quanta
  static G4double fTheMinAngle ;               //  max theta of TR quanta
         G4double fMaxThetaTR ;                //  max theta of TR quanta
  static G4int    fBinTR ;                     //  number of bins in TR vectors

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
