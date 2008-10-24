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
// $Id: G4MscRadiation.hh,v 1.7 2008-10-24 16:12:21 grichine Exp $
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

  G4complex CalculateMigdalXi(G4complex s);

  G4complex CalculateMigdalPhi(G4complex s);

  G4complex CalculateMigdalPsi(G4complex s);



  G4complex CalculateCorrectionMsc(G4double energy);

  G4double CalculateMscMigdalDiffdNdx(  G4DynamicParticle* dParticle, G4double energy);

  G4double CalculateMscE146DiffdNdx(  G4DynamicParticle* dParticle, G4double energy);

  G4double SupressionFunction(G4double kineticEnergy, G4double gammaEnergy);
  void CalcLPMFunctions(G4double kineticEnergy, G4double gammaEnergy);


  G4double GetMigdalS()  { return fMigdalS; };
  G4double GetMigdalXi() { return fMigdalXi; };
  G4double GetMigdalPsi(){ return fMigdalPsi; };
  G4double GetMigdalPhi(){ return fMigdalPhi; };
  G4double GetMigdalG()  { return fMigdalG; };


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
  G4double fZ;
  G4double fY;
  G4double fZommerfeld;
  G4double fAm;
  G4double fRadLength;
  G4complex fCorTMGY;

  G4double fMigdalS;
  G4double fMigdalXi;
  G4double fMigdalPsi;
  G4double fMigdalPhi;
  G4double fMigdalG;



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
  fZ = Z;
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
    // MGY medium corrections

  G4complex s = CalculateMigdalS(energy);

  G4complex xi  = CalculateMigdalXi(s);  
  if(verboseLevel) 
  {
    G4cout<<"xi = "<< xi  <<G4endl;
  }
  G4complex phi = CalculateMigdalPhi(s);  
  if(verboseLevel) 
  {
    G4cout<<"phi = "<< phi  <<G4endl;
  }
  G4complex psi = CalculateMigdalPsi(s);
  if(verboseLevel) 
  {
    G4cout<<"psi = "<< psi  <<G4endl;
  }

  G4complex G = 3.*psi - 2.*phi;

  // G4double sf = SupressionFunction(energy, Tkin);

  mscDiffXsc *=  real(G)/polarisation;
  // mscDiffXsc *=  fMigdalG/polarisation;

  G4double lnfAm = std::log(183.*std::pow(Z,-1./3.)*std::sqrt(polarisation));
	   lnfAm += std::log(1194.*std::pow(Z,-2./3.)*(polarisation))/Z;
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
  G4double onePerL;
  G4double Z       = fPlateMaterial->GetZ();

  onePerL         = Z*Z*std::log(183.*std::pow(Z,-1./3.));
  onePerL         += Z*std::log(1194.*std::pow(Z,-2./3.));
  onePerL         *= fPlateMaterial->GetTotNbOfAtomsPerVolume();
  onePerL         *= 4.*fine_structure_const*classic_electr_radius*classic_electr_radius;

  fRadLength = 1/onePerL;
  // return onePerL;
}


///////////////////////////////////////////////////////////////////////////
//
// Polarisation correction in absorbing medium

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

/////////////////////////////////////////////////////////////////////////
//
// return s-variable of Migdal. Modify!

inline  G4complex G4MscRadiation::CalculateMigdalS(G4double energy)
{
  // G4complex tmgy = CalculateCorrectionTMGY(energy);
  G4complex tmgy = fCorTMGY/(8.*fGamma*fGamma); 
  G4double hq = pi*hbarc/(2.*fine_structure_const*fGamma*fGamma*fRadLength);

  hq *= 1. - fY;  // high energy correction

  G4double ratio = energy/hq;  

  return tmgy*std::sqrt(ratio);
}

/////////////////////////////////////////////////////////////////////////
//
// return Migdal function xi(s)

inline  G4complex G4MscRadiation::CalculateMigdalXi(G4complex s)
{
  G4double s1 = std::pow(fZ,2./3.)/(184.*184.);

  if ( real(s) < s1)       return G4complex(2.,0.);  
  else if ( real(s) >= 1.) return G4complex(1.,0.); 
  else                     return 1. + std::log(s)/std::log(s1); 
}


/////////////////////////////////////////////////////////////////////////
//
// return Migdal function phi(s)

inline  G4complex G4MscRadiation::CalculateMigdalPhi(G4complex x)
{
  G4complex order = 6.*x;
  // order *= 1. + (3. - pi)*x;
  // order -= x*x*x/(0.623+0.796*x+0.658*x*x); 

  return 1. - std::exp(-order);
}


/////////////////////////////////////////////////////////////////////////
//
// return Migdal function psi(s)

inline  G4complex G4MscRadiation::CalculateMigdalPsi(G4complex x)
{
  G4complex order = 4.*x;
  // order +=  8.*x*x/(1.+ 3.96*x + 4.97*x*x - 0.05*x*x*x + 7.5*x*x*x*x); 
  return 1. - std::exp(-order);
}


////////////////////////////////////////////////////////////////////////
//
// Return G-correction of Migdal msc suppression

inline  G4complex G4MscRadiation::CalculateCorrectionMsc(G4double energy)
{
  G4complex order = CalculateMigdalS(energy);

  return 3.*( 1. - std::exp(-4.*order) ) - 2.*( 1. - std::exp(-6.*order) );
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

  // G4double sf = SupressionFunction(Tkin, energy);

  G4complex cMsc = CalculateCorrectionMsc(energy); 

  // G4double corMedium = fMigdalG/real(fCorTMGY); // RelBRmodel

  G4double corMedium = real(cMsc/fCorTMGY); // MGY model

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

  G4double y = energy/(Tkin+electron_mass_c2);

  fY = y; 

  if(verboseLevel) 
  {
    G4cout<<"y = "<<y<<G4endl;
  }
  fZ = fPlateMaterial->GetZ();
  G4double pMass = dParticle->GetMass();
  fGamma = 1. + Tkin/pMass;

  // TM medium correction

  CalculateCorrectionTMGY(energy);

  // G4double mscDiffXsc = fZ*fZ*std::log(184.*std::pow(fZ,-1./3.)*std::sqrt(real(fCorTMGY)));
  G4double mscDiffXsc = fZ*fZ*std::log(184.*std::pow(fZ,-1./3.));

  mscDiffXsc += fZ*std::log(1194.*std::pow(fZ,-2./3.));

  // MGY medium corrections

  G4complex s = CalculateMigdalS(energy);

  G4complex xi  = CalculateMigdalXi(s);  
  if(verboseLevel) 
  {
    G4cout<<"xi = "<< xi  <<G4endl;
  }
  G4complex phi = CalculateMigdalPhi(s);  
  if(verboseLevel) 
  {
    G4cout<<"phi = "<< phi  <<G4endl;
  }
  G4complex psi = CalculateMigdalPsi(s);
  if(verboseLevel) 
  {
    G4cout<<"psi = "<< psi  <<G4endl;
  }

  G4complex G = 3.*psi - 2.*phi;

  // phi /= fCorTMGY;

  // xi /= fCorTMGY;

  G /= fCorTMGY;

    // G4complex cMsc = CalculateCorrectionMsc(energy); 
  
  // G4double sf = SupressionFunction(Tkin, energy);
  // sf /= real(fCorTMGY);

  CalcLPMFunctions(Tkin, energy);
  mscDiffXsc *= ( y*y*fMigdalG + 2.*( 1. + (1.-y)*(1.-y) )*fMigdalPhi )*fMigdalXi;
  mscDiffXsc /= real(fCorTMGY);

  // mscDiffXsc *= real( xi*( y*y*G + 2.*( 1. + (1.-y)*(1.-y) )*phi ) );
  // mscDiffXsc *= real( G*( y*y + 2.*( 1. + (1.-y)*(1.-y) ) ) );
  // mscDiffXsc *= ( y*y + 2.*( 1. + (1.-y)*(1.-y) ) )*sf ;


  if(verboseLevel) 
  {
    G4cout<<"mscDiffXsc = "<< mscDiffXsc  <<G4endl;
  }

  // mscDiffXsc += (1.-y)*(fZ*fZ + fZ)/3.;

  mscDiffXsc *= 4.*fine_structure_const*classic_electr_radius*classic_electr_radius;

  mscDiffXsc *= fPlateMaterial->GetTotNbOfAtomsPerVolume()/(energy*3.); 

  if(verboseLevel) 
  {
    G4cout<<"mscDiffXsc = "<< mscDiffXsc  <<G4endl;
  }


  // G4double corMedium = real(cMsc/fCorTMGY);

  // mscDiffXsc *= corMedium;

  if( mscDiffXsc < 0. || y >= 1.)  mscDiffXsc = 0.;

  return mscDiffXsc;
}

////////////////////////////////////////////////////////////////////////
//
// From G4eBrHEModel. For testing and tuning

inline G4double G4MscRadiation::SupressionFunction(G4double kineticEnergy, G4double gammaEnergy)
{
  G4bool mscLPMflag = true;
  G4double mscMigdalConstant = classic_electr_radius*electron_Compton_length*
                               electron_Compton_length*4.*pi;
  G4double mscLPMconstant = fine_structure_const*electron_mass_c2*electron_mass_c2/
                            (4.*pi*hbarc);

  G4double Eele = kineticEnergy + electron_mass_c2;

  G4double k = gammaEnergy;

  G4double Elpm = mscLPMconstant*(fPlateMaterial->GetRadlen());

  // WARNING tested for single element materials only!

  G4double Z = fPlateMaterial->GetElectronDensity()/
               fPlateMaterial->GetTotNbOfAtomsPerVolume();

  G4double supr = 1.0;

  if (mscLPMflag) 
  {
    G4double y = k/Eele;
    G4double y2 = y*y;
    G4double yone2 = 2.*(1+sqr(1.-y));

    // numerical safety factor
    // WARNING may need modification for very high energy electrons

    if ( k > 0.9999*Eele )  k = 0.9999*Eele;
    

    // WARNING needs Eele and Z dependence!

    if ( Elpm <= 200.*TeV ) 
    {
      if(verboseLevel) 
      {
        G4cout<<"Elpm = "<<Elpm/TeV<<" TeV"<<G4endl;
      }
      // *** calculate lpm variable s & sprime ***
      // Klein eqs. (78) & (79)

      G4double sprime = sqrt(0.125*k*Elpm/(Eele*(Eele-k)));

      if(verboseLevel) 
      {
        G4cout<<"sprime = "<<sprime<<" "<<G4endl;
      }
      G4double s1 = pow(Z,2./3.)/(184.*184.);

      G4double h  = log(sprime)/log(sqrt(2.)*s1);

      G4double xiprime = 2;

      if (sprime > 1.) xiprime = 1.;
      else if (sprime > sqrt(2.)*s1 )
      {
        xiprime = 1+h-0.08*(1-h)*(1-sqr(1-h))/log(sqrt(2.)*s1);
      }
      G4double s = sprime/sqrt(xiprime);

      if(verboseLevel) 
      {
        G4cout<<"s = "<<s<<" "<<G4endl;
      }

      // merging with density effect 
      // using Ter-Mikaelian eq. (20.9)

      G4double k2 = gammaEnergy*gammaEnergy;
      G4double Eele2 = Eele*Eele;
      G4double eDensity = fPlateMaterial->GetElectronDensity();


      G4double kp2 = mscMigdalConstant*Eele2*eDensity;

      s = s * ( 1 + (kp2/k2)*(1.-y) );

      fMigdalS = s;

      if(verboseLevel) 
      {
        G4cout<<"s = "<<s<<" with TM"<<G4endl;
      }


      // recalculate Xi 
      // Klein eq. (75)

      G4double xi = 1;
      if ( s < s1 )                 xi = 2;
      if ( ( s1 < s) && (s <= 1) ) xi = 1 + log(s)/log(s1);

      // Klein eqs. (77)

      G4double s2=s*s;
      G4double s3=s*s2;
      G4double s4=s2*s2;

      G4double psi = 0., phi = 0., G = 0.;

      if ( s < 0.1 ) 
      {
        // high suppression limit

        phi = 6.*s - 18.84955592153876*s2 + 39.47841760435743*s3
              - 57.69873135166053*s4;
        G = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
      }
      else if (s < 1.9516) 
      {
        // intermediate suppression
        // using eq.77 approxim. valid s<2.
        phi = 1-exp(-6*s*(1+(3-pi)*s)
                    +s3/(0.623+0.795*s+0.658*s2));

        if (s<0.415827397755) 
        {
          // using eq.77 approxim. valid 0.07<s<2
          psi = 1-exp(-4*s-8*s2/(1+3.936*s+4.97*s2-0.05*s3+7.50*s4));
          G = 3*psi-2*phi;
        }
        else 
        {
          // using alternative parametrisiation
          G4double par[] = {-0.16072300849123999, 3.7550300067531581,
                            -1.7981383069010097, 0.67282686077812381, 
                            -0.1207722909879257};
          G4double pre =  par[0] + s*(par[1] + s*(par[2] + s*(par[3] + s*par[4])));
          G = tanh(pre);
        }
      }
      else 
      {
        // low suppression limit valid s>2.
        phi = 1. - 0.0119048/s4;
        G = 1. - 0.0230655/s4;
     
      }
      if(verboseLevel) 
      {
        G4cout<<"xi = "<<xi<<" "<<G4endl;
      }
      if(verboseLevel) 
      {
        G4cout<<"phi = "<<phi<<" "<<G4endl;
      }
      if(verboseLevel) 
      {
        G4cout<<"psi = "<<psi<<" "<<G4endl;
      }
      if(verboseLevel) 
      {
        G4cout<<"G = "<<G<<" "<<G4endl;
      }
      if(verboseLevel) 
      {
        G4cout<<"xiprime = "<<xiprime<<" "<<G4endl;
      }
      fMigdalXi  = xi;
      fMigdalPhi = phi;
      fMigdalPsi = psi;
      fMigdalG   = G;


      // make sure suppression is smaller than 1 
      // caused by Migdal approximation in xi    

      // if ( xiprime*phi >1. || s > 0.57 )  xiprime=1./phi;
      if ( xi*phi >1. || s > 0.57 )  xi = 1./phi;


      // calculate suppression factor
 
      // G4double mainLPM = xiprime*(y2 * G + yone2*phi);
      G4double mainLPM = xi*(y2 * G + yone2*phi);

      G4double main = (y2 + yone2); // needed part of classical BH-xsc

      supr = mainLPM/main;
      // if (supr > 1.) supr = 1.;
      // if (supr < 0.) supr = 0.;
    }
  }
  return supr;
}

//////////////////////////////////////////////////////////////////////////////
//
// From G4eBrRelModel for testing and tuning

inline void  G4MscRadiation::CalcLPMFunctions(G4double kineticEnergy, G4double k)
{
  // PROFILE_HERE;
  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)

  G4double mscMigdalConstant = classic_electr_radius*electron_Compton_length*
                               electron_Compton_length*4.*pi;
  G4double mscLPMconstant = fine_structure_const*electron_mass_c2*electron_mass_c2/
                            (4.*pi*hbarc);

  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double y = k/totalEnergy;
  G4double lpmEnergy = mscLPMconstant*(fPlateMaterial->GetRadlen());

  G4double Eele2       = totalEnergy*totalEnergy;
  G4double eDensity    = fPlateMaterial->GetElectronDensity();
  G4double densityCorr = mscMigdalConstant*Eele2*eDensity;

  G4double facFel = log(184.15);
  G4double facFinel = log(1194.);

  G4double preS1 = 1./(184.15*184.15);
  G4double logTwo = log(2.);

  G4double Z = fPlateMaterial->GetElectronDensity()/
               fPlateMaterial->GetTotNbOfAtomsPerVolume();
  G4double z23 = pow(Z,2./3.);
  G4double lnZ = log(Z);


  G4double sprime = sqrt(0.125*k*lpmEnergy/(totalEnergy*(totalEnergy-k)));

  G4double s1 = preS1*z23;
  G4double logS1 = 2./3.*lnZ-2.*facFel;
  G4double logTS1 = logTwo+logS1;

  fMigdalXi = 2.;

  if (sprime>1)    fMigdalXi = 1.;

  else if (sprime>sqrt(2.)*s1) 
  {
    G4double h  = log(sprime)/logTS1;
    fMigdalXi = 1+h-0.08*(1-h)*(1-sqr(1-h))/logTS1;
  }

  G4double s = sprime/sqrt(fMigdalXi);

  // merging with density effect***  should be only necessary in 
  // region "close to" kp, e.g. k<100*kp
  // using Ter-Mikaelian eq. (20.9)

  G4double k2 = k*k;
  s = s * (1 + (densityCorr/k2) );
  fMigdalS = s;
  // recalculate Xi using modified s above
  // Klein eq. (75)
  fMigdalXi = 1.;
  if (s<=s1) fMigdalXi = 2.;
  else if ( (s1<s) && (s<=1) ) fMigdalXi = 1. + log(s)/logS1;


  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  G4double s2=s*s;
  G4double s3=s*s2;
  G4double s4=s2*s2;

  if ( s < 0.1 ) 
  {
    // high suppression limit

    fMigdalPhi = 6.*s - 18.84955592153876*s2 + 39.47841760435743*s3
      - 57.69873135166053*s4;
    fMigdalG = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if ( s < 1.9516 ) 
  {
    // intermediate suppression
    // using eq.77 approxim. valid s<2.

    fMigdalPhi = 1.-exp(-6.*s*(1.+(3.-pi)*s)
                +s3/(0.623+0.795*s+0.658*s2));

    if ( s < 0.415827397755 ) 
    {
      // using eq.77 approxim. valid 0.07<s<2

      G4double fMigdalPsi = 1-exp(-4*s-8*s2/(1+3.936*s+4.97*s2-0.05*s3+7.50*s4));
      fMigdalG = 3*fMigdalPsi-2*fMigdalPhi;
    }
    else 
    {
      // using alternative parametrisiation

      G4double pre = -0.16072300849123999 + s*3.7550300067531581 + s2*-1.7981383069010097
        + s3*0.67282686077812381 + s4*-0.1207722909879257;
      fMigdalG = tanh(pre);
    }
  }
  else 
  {
    // low suppression limit valid s>2.

    fMigdalPhi = 1. - 0.0119048/s4;
    fMigdalG = 1. - 0.0230655/s4;
  }


  // make sure suppression is smaller than 1 
  // caused by Migdal approximation in xi    

  if ( fMigdalXi*fMigdalPhi > 1. || s > 0.57)  fMigdalXi=1./fMigdalPhi;
}





#endif
