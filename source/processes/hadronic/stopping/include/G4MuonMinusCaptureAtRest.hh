// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MuonMinusCaptureAtRest.hh,v 1.3 2000-04-07 16:06:47 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuonMinusCaptureAtRest physics process ------
//                   by Larry Felawka (TRIUMF)
//                     E-mail: felawka@alph04.triumf.ca
//                   and Art Olin (TRIUMF)
//                     E-mail: olin@triumf.ca
//                            April 1998
// ************************************************************
//-----------------------------------------------------------------------------

#ifndef G4MuonMinusCaptureAtRest_h
#define G4MuonMinusCaptureAtRest_h 1
 
#include "globals.hh"
#include "Randomize.hh" 
#include "G4VRestProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4GHEKinematicsVector.hh"
#include "G4StopElementSelector.hh"
#include "G4MuMinusCaptureCascade.hh"

class G4MuonMinusCaptureAtRest : public G4VRestProcess
 
{ 
  private:
  // hide assignment operator as private 
      G4MuonMinusCaptureAtRest& operator=(const G4MuonMinusCaptureAtRest &right);
      G4MuonMinusCaptureAtRest(const G4MuonMinusCaptureAtRest& );
   
  public:
 
     G4MuonMinusCaptureAtRest(const G4String& processName ="MuonMinusCaptureAtRest");
 
    ~G4MuonMinusCaptureAtRest();

     G4bool IsApplicable(const G4ParticleDefinition&);

  // null physics table
     void BuildPhysicsTable(const G4ParticleDefinition&){}

     G4double AtRestGetPhysicalInteractionLength(const G4Track&,
						 G4ForceCondition*);

     G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*);

     G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&); 

  // return number of secondaries produced
     G4int GetNumberOfSecondaries();

  // pointer to array containg kinematics of secondaries
     G4GHEKinematicsVector* GetSecondaryKinematics();

  private:

     void CascadeCorrection(G4double, G4double);
     G4double CoulombBarrier(G4int, G4double);
     void ExcitationEnergyLevel(G4int, G4int, G4double*, G4double*, G4double*);
     G4double CalculateIsotopicMass(G4double, G4double);
     void EvaporationDeexcitation();
     void FermiMotion(G4int);
     void ResidualNucleusCascade(G4int, G4int, G4double, G4double*,
				 G4double*, G4bool*);
     G4double GetIsotopicMass(G4double, G4double);
     void Erup();
     void InitializeMuCapture();
     G4double LevelDensity(G4int, G4int, G4double*, G4double*);
     void GetCaptureIsotope(const G4Track&);
     void NuclearExcitation(G4double, G4double, G4double, G4double*,
			    G4double*, G4double, G4double);
     void DoMuCapture();
     void MuEvaporation();
     void RanPolarAng(G4double*, G4double*);
     G4double NuclearBindingEnergy(G4double, G4double, G4double, G4double);
     void RanDirCos(G4double*, G4double*, G4double*);
     void RanAzimuthalAng(G4double*, G4double*);

  private:

// relative time-of-flight of stopped hadron
     G4double  tDelay;

// atomic mass of target nucleus
     G4double  targetAtomicMass;

// charge of target nucleus
     G4double  targetCharge;

// effective charge squared of target nucleus
     G4double  zeff2;

// nuclear isotopica data
     G4double  isotopicData[11][250];

// residual nuclear masses
     G4double  residualNuclearMasses[297];

     G4GHEKinematicsVector* Fragments;
     G4GHEKinematicsVector* Secondaries;
     G4GHEKinematicsVector* Evaporates;
     G4GHEKinematicsVector* Gkin;
     G4GHEKinematicsVector* Cascade;

     G4int    nGkine;
     G4int    nCascade;

     G4int    nSecPart;
     G4int    nallFragm;
     G4int    nVariousFragm[6];

     G4double nucleonMass[2];
     G4double nucleonMassSquared[2];
     G4double nucleonBindingEn[2];
     G4double nucleonPotWell[2];
     G4double nucleonMaxFermiEn[2];

     G4double atMassCurrent;
     G4double chargeCurrent;

     G4double massTarget;
     G4int    atMassTarget;
     G4int    chargeTarget;

     G4double pairCorrFlag;

     G4double nucleonFermiEn;
     G4double nucleonFermiMomX;
     G4double nucleonFermiMomY;
     G4double nucleonFermiMomZ;

     G4double atMassEvapNucl;
     G4double chargeEvapNucl;
     G4double iniExcitEnEvapNucl;
     G4double finExcitEnEvapNucl;
     G4double recoilEnEvapNucl;

// residual nucleus kinematical quantities
     G4double totMomResidNucl;
     G4double momXResidNucl;
     G4double momYResidNucl;
     G4double momZResidNucl;
     G4double dirCosXResidNucl;
     G4double dirCosYResidNucl;
     G4double dirCosZResidNucl;
     G4double massResidNucl;
     G4double totEnResidNucl;
     G4double kinEnResidNucl;
     G4double excitEnResidNucl;
     G4double recoilEnResidNucl;
     G4int    atMassResidNucl;
     G4int    chargeResidNucl;

     G4double  massGamma;
     G4double  massNeutrinoE;
     G4double  massProton;
     G4double  massNeutron;
     G4double  massFragm[6];

     G4ParticleDefinition* pdefGamma;
     G4ParticleDefinition* pdefNeutrinoE;
     G4ParticleDefinition* pdefMuonMinus;
     G4ParticleDefinition* pdefNeutron;
     G4ParticleDefinition* pdefFragm[6];

     G4StopElementSelector*   pSelector;
     G4MuMinusCaptureCascade* pEMCascade;

};

#endif
 





