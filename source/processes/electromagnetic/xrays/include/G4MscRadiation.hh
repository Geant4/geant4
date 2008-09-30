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
// $Id: G4MscRadiation.hh,v 1.5 2008-09-30 14:39:07 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Class for description of X-ray radiation produced by multiple scattered
// charged particle in absorbing medium (see DoIt
// method
//
// author: V. Grichine (Vladimir.Grichine@cern.ch)
//
// History:
//
// 06.08.08 V. Grichine first implementation using methods of XTR and diffuse elastic
//

#ifndef G4MscRadiation_h
#define G4MscRadiation_h 1

#include <complex>
#include "globals.hh"
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
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Integrator.hh"
#include "G4ParticleChange.hh"

class G4SandiaTable;
class G4VParticleChange;
class G4PhysicsFreeVector;
class G4ParticleDefinition;
class G4DynamicParticle;

class G4MscRadiation : public G4VDiscreteProcess  // G4VContinuousProcess
{
public:

  G4MscRadiation (G4LogicalVolume *anEnvelope,G4Material*,G4Material*,
                    G4double,G4double,G4int,
                    const G4String & processName = "MscRadiation",
                    G4ProcessType type = fElectromagnetic);

  // test constructor

  G4MscRadiation (  G4Material*,G4double,G4DynamicParticle*,
                    const G4String & processName = "MscRadiation",
                    G4ProcessType type = fElectromagnetic);


  G4MscRadiation (  G4Material*,G4double,
                    const G4String & processName = "MscRadiation",
                    G4ProcessType type = fElectromagnetic);

   ~G4MscRadiation ();

  // These add XTR for thin radiators like films
 
  G4double GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle );

  G4bool IsApplicable(const G4ParticleDefinition&);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
				   const G4Step&  aStep);

  G4double GetMeanFreePath(const G4Track& aTrack,
                           G4double previousStepSize,
                           G4ForceCondition* condition);

  void BuildPhysicsTable(const G4ParticleDefinition&);
  void BuildTable();
  void BuildEnergyTable();
  void BuildAngleTable();
  void BuildGlobalAngleTable();

  G4complex OneInterfaceXTRdEdx( G4double energy, 
                                G4double gamma,
                                G4double varAngle );

  G4double SpectralAngleXTRdEdx(G4double varAngle);

  virtual  G4double SpectralXTRdEdx(G4double energy);

  G4double AngleSpectralXTRdEdx(G4double energy);

  G4double AngleXTRdEdx(G4double varAngle);


  /////////////////////////////////////////////////////////////

  G4double OneBoundaryXTRNdensity( G4double energy,
                                   G4double gamma,
				   G4double varAngle ) const;


  // for photon energy distribution tables

  G4double XTRNSpectralAngleDensity(G4double varAngle);
  G4double XTRNSpectralDensity(G4double energy);
  
  // for photon angle distribution tables

  G4double XTRNAngleSpectralDensity(G4double energy);
  G4double XTRNAngleDensity(G4double varAngle);

  void GetNumberOfPhotons();  

  // Auxiliary functions for plate/gas material parameters

  G4double  GetPlateFormationZone(G4double,G4double,G4double);
  G4complex GetPlateComplexFZ(G4double,G4double,G4double);
  void      ComputePlatePhotoAbsCof();
  G4double  GetPlateLinearPhotoAbs(G4double);
  void      GetPlateZmuProduct();
  G4double  GetPlateZmuProduct(G4double,G4double,G4double);

  G4double  GetGasFormationZone(G4double,G4double,G4double);
  G4complex GetGasComplexFZ(G4double,G4double,G4double);
  void      ComputeGasPhotoAbsCof();
  G4double  GetGasLinearPhotoAbs(G4double);
  void      GetGasZmuProduct();
  G4double  GetGasZmuProduct(G4double,G4double,G4double);

  G4double  GetPlateCompton(G4double);
  G4double  GetGasCompton(G4double);
  G4double  GetComptonPerAtom(G4double,G4double);

  G4double  GetXTRrandomEnergy( G4double scaledTkin, G4int iTkin );
  G4double  GetXTRenergy( G4int iPlace, G4double position, G4int iTransfer  );

  G4double  GetRandomAngle( G4double energyXTR, G4int iTkin );
  G4double  GetAngleXTR(G4int iTR,G4double position,G4int iAngle);

  G4double  GetGamma()   {return fGamma;}; 
  G4double  GetEnergy()  {return fEnergy;};                
  G4double  GetVarAngle(){return fVarAngle;};


  // Test methods for msc radiation models 

  G4double CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    );
  G4double CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 );
  G4double CalculateAm( G4double momentum, G4double n, G4double Z);
  G4double CalculateTransportXsc(  const G4ParticleDefinition* aParticle, 
                                          G4double momentum, G4double Z);
  G4double CalculateMscDiffdNdx(  G4DynamicParticle* dParticle, G4double energy);
               
  void // G4double 

  CalculateReciprocalRadLength(); 

  G4double GetRadLength(){return fRadLength;};

  void // G4complex 

  CalculateCorrectionTMGY(G4double energy);

  G4complex GetCorrectionTMGY(){return fCorTMGY;};

  G4complex CalculateMigdalS(G4double energy);

  G4complex CalculateCorrectionMsc(G4double energy);

  G4double CalculateMscMigdalDiffdNdx(  G4DynamicParticle* dParticle, G4double energy);

  G4double CalculateMscE146DiffdNdx(  G4DynamicParticle* dParticle, G4double energy);


  void SetGamma(G4double gamma)      {fGamma    = gamma;}; 
  void SetEnergy(G4double energy)    {fEnergy   = energy;};                
  void SetVarAngle(G4double varAngle){fVarAngle = varAngle;};               
  void SetAngleRadDistr(G4bool pAngleRadDistr){fAngleRadDistr=pAngleRadDistr;};               
  void SetCompton(G4bool pC){fCompton=pC;};               

  G4PhysicsLogVector* GetProtonVector(){ return fProtonEnergyVector;};
  G4int GetTotBin(){return fTotBin;};           
  G4PhysicsFreeVector* GetAngleVector(G4double energy, G4int n);

protected:

  G4ParticleDefinition* fPtrGamma;    // pointer to TR photon
  G4Material* fPlateMaterial; 
  G4double* fGammaCutInKineticEnergy; // TR photon cut in energy array

  G4double         fGammaTkinCut;     // Tkin cut of TR photon in current mat.
  G4LogicalVolume* fEnvelope;
  G4PhysicsTable*  fAngleDistrTable;
  G4PhysicsTable*  fEnergyDistrTable;

  G4PhysicsLogVector* fProtonEnergyVector;
  G4PhysicsLogVector* fXTREnergyVector;

  G4double fTheMinEnergyTR;            //   min TR energy
  G4double fTheMaxEnergyTR;            //   max TR energy
  G4double fMinEnergyTR;               //  min TR energy in material
  G4double fMaxEnergyTR;               //  max TR energy in material
  G4double fTheMaxAngle;               //  max theta of TR quanta
  G4double fTheMinAngle;               //  max theta of TR quanta
  G4double fMaxThetaTR;                //  max theta of TR quanta
  G4int    fBinTR;                     //  number of bins in TR vectors

  G4double fMinProtonTkin;             // min Tkin of proton in tables
  G4double fMaxProtonTkin;             // max Tkin of proton in tables
  G4int    fTotBin;                    // number of bins in log scale
  G4double fGamma;                     // current Lorentz factor
  G4double fEnergy;                    // energy and
  G4double fVarAngle;                  // angle squared
  G4double fLambda;

  G4double fPlasmaCof;                // physical consts for plasma energy
  G4double fCofTR;  

  G4bool   fExitFlux;
  G4bool   fAngleRadDistr;
  G4bool   fCompton;
  G4double fSigma1; 
  G4double fSigma2;                    // plasma energy Sq of matter1/2

  G4int    fMatIndex1;
  G4int    fMatIndex2;
  G4int    fPlateNumber;

  G4double fTotalDist;
  G4double fPlateThick;
  G4double fGasThick;     
  G4double fAlphaPlate;
  G4double fAlphaGas;

  // test fields for msc radiation

  G4double fBeta;
  G4double fZommerfeld;
  G4double fAm;
  G4double fRadLength;
  G4complex fCorTMGY;



  G4SandiaTable* fPlatePhotoAbsCof;
 
  G4SandiaTable* fGasPhotoAbsCof;

  G4ParticleChange fParticleChange;

  G4PhysicsTable*                    fAngleForEnergyTable;
  std::vector<G4PhysicsTable*>       fAngleBank;

};

// inline test methods

////////////////////////////////////////////////////////////////////
//
// return particle beta

inline  G4double G4MscRadiation::CalculateParticleBeta( const G4ParticleDefinition* particle, 
                                 	G4double momentum    )
{
  G4double mass = particle->GetPDGMass();
  G4double a    = momentum/mass;
  fBeta         = a/std::sqrt(1+a*a);

  return fBeta; 
}

////////////////////////////////////////////////////////////////////
//
// return Zommerfeld parameter for Coulomb scattering

inline  G4double G4MscRadiation::CalculateZommerfeld( G4double beta, G4double Z1, G4double Z2 )
{
  fZommerfeld = fine_structure_const*Z1*Z2/beta;

  return fZommerfeld; 
}

////////////////////////////////////////////////////////////////////
//
// return Wentzel correction for Coulomb scattering

inline  G4double G4MscRadiation::CalculateAm( G4double momentum, G4double n, G4double Z)
{
  G4double k   = momentum/hbarc;
  G4double ch  = 1.13 + 3.76*n*n;
  G4double zn  = 1.77*k*std::pow(Z,-1./3.)*Bohr_radius;
  G4double zn2 = zn*zn;
  fAm          = ch/zn2;

  return fAm;
}
////////////////////////////////////////////////////////////////////
//
// return Wentzel correction for Coulomb scattering

inline  G4double G4MscRadiation::CalculateTransportXsc( const G4ParticleDefinition* aParticle, 
                                                        G4double momentum, G4double Z)
{

  G4double z = aParticle->GetPDGCharge();
  if(verboseLevel) 
  {
    G4cout<<"z = "<<z<<G4endl;
  }

  fBeta       = CalculateParticleBeta( aParticle, momentum);
  if(verboseLevel) 
  {
    G4cout<<"fBeta = "<<fBeta<<G4endl;
  }
  fZommerfeld = CalculateZommerfeld( fBeta, z, Z );
  if(verboseLevel) 
  {
    G4cout<<"fZommerfeld = "<<fZommerfeld<<G4endl;
  }
  fAm         = CalculateAm( momentum, fZommerfeld, Z );
  if(verboseLevel) 
  {
    G4cout<<"fAm = "<<fAm<<G4endl;
  }

  G4double xsc = twopi*z*z*Z*Z;

  xsc *= fine_structure_const*fine_structure_const*hbarc*hbarc;

  xsc /=  momentum*momentum*fBeta*fBeta;

  // G4double lnfAm = std::log(1.+1./fAm) - 1./(1.+fAm);
  G4double lnfAm = 1.;


  if(verboseLevel) 
  {
    // G4cout<<"lnfAm = "<<lnfAm<<G4endl;
  }
  xsc *= lnfAm;

  return xsc;
}

///////////////////////////////////////////////////////////////////
//
// return Msc differential xsc


inline  G4double G4MscRadiation::CalculateMscDiffdNdx(  G4DynamicParticle* dParticle, G4double energy)
{

  const G4ParticleDefinition* aParticle = dParticle->GetDefinition(); 
  G4double momentum = dParticle->GetTotalMomentum();
  G4double Tkin = dParticle->GetKineticEnergy();
  if(Tkin <= energy || Tkin <= 0. || energy <= 0.) return 0.;

  G4double y = energy/Tkin;
  if(verboseLevel) 
  {
    G4cout<<"y = "<<y<<G4endl;
  }

  G4double pMass = dParticle->GetMass();
  fGamma = 1. + Tkin/pMass;
  G4double Z = fPlateMaterial->GetZ();
  if(verboseLevel) 
  {
    G4cout<<"Z = "<<Z<<G4endl;
  }
  G4double lambda = hbarc/energy;
  if(verboseLevel) 
  {
    G4cout<<"lambda = "<<lambda<<" mm"<<G4endl;
  }

  G4double transportXsc  = CalculateTransportXsc( aParticle, momentum,  Z);

  transportXsc  *=   fPlateMaterial->GetTotNbOfAtomsPerVolume()*lambda;
  if(verboseLevel) 
  {
    G4cout<<"Nat = "<<fPlateMaterial->GetTotNbOfAtomsPerVolume()<<" 1/mm3"<<G4endl;
  }

  if(verboseLevel) 
  {
    G4cout<<"transportXsc = "<<transportXsc<<" "<<G4endl;
  }

  G4double xscCompton = GetPlateCompton(energy)*lambda;
  if(verboseLevel) 
  {
    G4cout<<"xscCompton = "<<xscCompton<<" "<<G4endl;
  }
  
  G4double xscAbs = GetPlateLinearPhotoAbs(energy)*lambda;
  if(verboseLevel) 
  {
    G4cout<<"xscAbs = "<<xscAbs<<" "<<G4endl;
  }

  G4double mscDiffXsc = fine_structure_const/(pi*hbarc);

  G4double polarisation = 1. +(fGamma*fGamma)*fSigma1/(energy*energy);
  if(verboseLevel) 
  {
    G4cout<<"polarisation = "<<polarisation<<G4endl;
  }
  mscDiffXsc /=  polarisation;

  G4double lnfAm = std::log(183.*std::pow(Z,-1./3.)*std::sqrt(polarisation));
  if(verboseLevel) 
  {
    G4cout<<"lnfAm = "<<lnfAm<<G4endl;
  }

  G4double image = 2*fBeta*transportXsc*lnfAm;

  // image += xscAbs;
  // image += xscCompton;

  mscDiffXsc *= image;
  mscDiffXsc *= fGamma*fGamma;
  mscDiffXsc *= y*y + 4.*(1.- y)/3.;
  

  if( mscDiffXsc < 0.  || y >= 1. )  mscDiffXsc = 0.;

  return mscDiffXsc;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// methods for Migdal averaging in absorbing medium

inline void  // G4double 
G4MscRadiation::CalculateReciprocalRadLength()
{
  G4double onePerL = 4.*fine_structure_const*classic_electr_radius*classic_electr_radius;
  G4double Z       = fPlateMaterial->GetZ();
  onePerL         *= Z*Z*fPlateMaterial->GetTotNbOfAtomsPerVolume();
  onePerL         *= std::log(183.*std::pow(Z,-1./3.));

  fRadLength = 1/onePerL;
  // return onePerL;
}


inline  void // G4complex 
G4MscRadiation::CalculateCorrectionTMGY(G4double energy)
{
  G4double lambda;

  G4double re = 1. +(fGamma*fGamma)*fSigma1/(energy*energy);

  if(energy > 1.*MeV) lambda = fRadLength;
  else lambda = 1./GetPlateLinearPhotoAbs(energy);

  G4double im = -hbarc*fGamma*fGamma/(energy*lambda);
  fCorTMGY = G4complex(re,im);
}

inline  G4complex G4MscRadiation::CalculateMigdalS(G4double energy)
{
  // G4complex tmgy = CalculateCorrectionTMGY(energy);
  G4complex tmgy = fCorTMGY/(8.*fGamma*fGamma); 
  G4double hq = pi*hbarc/(2.*fine_structure_const*fGamma*fGamma*fRadLength);
  G4double ratio = energy/hq;
  return tmgy*std::sqrt(ratio);
}

inline  G4complex G4MscRadiation::CalculateCorrectionMsc(G4double energy)
{
  G4complex order = 6.*CalculateMigdalS(energy);
  return 1. - std::exp(-order);
}


///////////////////////////////////////////////////////////////////
//
// return Msc Migdal averaged in absorbing medium differential xsc


inline  G4double G4MscRadiation::CalculateMscMigdalDiffdNdx(  G4DynamicParticle* dParticle, 
                                                              G4double energy)
{

  // const G4ParticleDefinition* aParticle = dParticle->GetDefinition(); 
  // G4double momentum = dParticle->GetTotalMomentum();

  G4double Tkin = dParticle->GetKineticEnergy();
  if(Tkin <= energy || Tkin <= 0. || energy <= 0.) return 0.;

  G4double y = energy/Tkin;
  if(verboseLevel) 
  {
    G4cout<<"y = "<<y<<G4endl;
  }

  G4double pMass = dParticle->GetMass();
  fGamma = 1. + Tkin/pMass;

  G4double mscDiffXsc = 1/(energy*fRadLength);

  mscDiffXsc *= y*y + 4.*(1.- y)/3.;

  CalculateCorrectionTMGY(energy);

  G4complex cMsc = CalculateCorrectionMsc(energy); 

  G4double corMedium = real(cMsc/fCorTMGY);

  mscDiffXsc *= corMedium;

  if( mscDiffXsc < 0. || y >= 1.)  mscDiffXsc = 0.;

  return mscDiffXsc;
}


///////////////////////////////////////////////////////////////////
//
// return Msc E146 Anthony(9) averaged in absorbing medium differential xsc


inline  G4double G4MscRadiation::CalculateMscE146DiffdNdx(  G4DynamicParticle* dParticle, 
                                                              G4double energy)
{

  // const G4ParticleDefinition* aParticle = dParticle->GetDefinition(); 
  // G4double momentum = dParticle->GetTotalMomentum();

  G4double Tkin = dParticle->GetKineticEnergy();
  if(Tkin <= energy || Tkin <= 0. || energy <= 0.) return 0.;

  G4double y = energy/Tkin;
  if(verboseLevel) 
  {
    G4cout<<"y = "<<y<<G4endl;
  }
  G4double Z = fPlateMaterial->GetZ();
  G4double pMass = dParticle->GetMass();
  fGamma = 1. + Tkin/pMass;

  G4double mscDiffXsc = Z*Z*std::log(184.*std::pow(Z,-1./3.));

  mscDiffXsc += Z*std::log(1194.*std::pow(Z,-2./3.));

  mscDiffXsc *= y*y + 2.*( 1. + (1.-y)*(1.-y) );

  mscDiffXsc += (1.-y)*(Z*Z+Z)/3.;

  mscDiffXsc *= 4.*fine_structure_const*classic_electr_radius*classic_electr_radius;

  mscDiffXsc *= fPlateMaterial->GetTotNbOfAtomsPerVolume()/(energy*3.); 



  CalculateCorrectionTMGY(energy);

  G4complex cMsc = CalculateCorrectionMsc(energy); 

  G4double corMedium = real(cMsc/fCorTMGY);

  mscDiffXsc *= corMedium;

  if( mscDiffXsc < 0. || y >= 1.)  mscDiffXsc = 0.;

  return mscDiffXsc;
}










#endif
