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
// $Id$
//
//
// G4 Gheisha High Energy (GHE) model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.
// It includes all Physics Routines from Geant3,
// relevant for simulation of hadronic processes
// above approx. 20 GeV incident momentum.
// Not included are the Physics Routines for
// stopping particles, the low energy neutron
// slowing down description and the low energy
// nuclear reactions a(A,A')b.

#ifndef G4HEInelastic_h
#define G4HEInelastic_h 1

// Class description:
// Each of the high energy parameterized models (e.g. G4HEProtonInelastic)
// derives from the G4HEInelastic class. This class contains the various
// algorithms needed to implement the interaction. These include 
// lambda-fragmentation, meson and nucleon cluster formation and decay, 
// nuclear cascade, and nuclear de-excitation.  
//
// This class is derived from G4HadronicInteraction.

// Class Description - End

#include "G4HEVector.hh"
#include "G4HadronicInteraction.hh"

class G4HEInelastic : public G4HadronicInteraction
{
 public:  // with description 
   G4HEInelastic(const G4String& modelName = "HEInelastic") 
    : G4HadronicInteraction(modelName)
   { 
     SetParticles();
     verboseLevel = 0;
     MAXPART = 0;
     conserveEnergy = true;
   };

   ~G4HEInelastic(){ };
         
   void SetMaxNumberOfSecondaries(const G4int maxnumber)
            {MAXPART = maxnumber;}      
 
   void SetVerboseLevel(const G4int level) {verboseLevel = level;}

   G4int verboseLevel; 
   G4int MAXPART;
   G4bool conserveEnergy;

   void ForceEnergyConservation(G4bool energyConservation)
              {conserveEnergy = energyConservation;}

   G4bool EnergyConservation(void) {return conserveEnergy;} 

   virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

   G4double Amin(G4double a, G4double b);
   G4double Amax(G4double a, G4double b);
   G4int Imin(G4int a, G4int b);
   G4int Imax(G4int a, G4int b);
 
   void FillParticleChange(G4HEVector pv[], G4int aVecLength);

   G4double pmltpc(G4int np, G4int nm, G4int nz,
                   G4int n, G4double b, G4double c);

   G4int Factorial(G4int n); 

   G4double NuclearInelasticity(G4double incidentKineticEnergy,
                                G4double atomicWeight,
                                G4double atomicNumber);

   G4double NuclearExcitation(G4double  incidentKineticEnergy,
                              G4double  atomicWeight,
                              G4double  atomicNumber,
                              G4double& excitationEnergyCascade,
                              G4double& excitationEnergyEvaporation); 
   
   void HighEnergyCascading(G4bool& successful,
                            G4HEVector pv[],
                            G4int& vecLen,
                            G4double& excitationEnergyGNP,
                            G4double& excitationEnergyDTA, 
                            const G4HEVector& incidentParticle,
                            const G4HEVector& targetParticle,
                            G4double atomicWeight,
                            G4double atomicNumber);

   void HighEnergyClusterProduction(G4bool& successful,
                                    G4HEVector pv[],
                                    G4int& vecLen,
                                    G4double& excitationEnergyGNP,
                                    G4double& excitationEnergyDTA, 
                                    const G4HEVector& incidentParticle,
                                    const G4HEVector& targetParticle,
                                    G4double atomicWeight,
                                    G4double atomicNumber);             

   void TuningOfHighEnergyCascading(G4HEVector pv[],
                                    G4int& vecLen,
                                    const G4HEVector& incidentParticle,
                                    const G4HEVector& targetParticle,
                                    G4double atomicWeight,
                                    G4double atomicNumber);   

   void MediumEnergyCascading(G4bool& successful,
                              G4HEVector pv[],
                              G4int& vecLen,
                              G4double& excitationEnergyGNP,
                              G4double& excitationEnergyDTA,
                              const G4HEVector& incidentParticle,
                              const G4HEVector& targetParticle,
                              G4double atomicWeight,
                              G4double atomicNumber);            

   void MediumEnergyClusterProduction(G4bool& successful,
                                      G4HEVector pv[],
                                      G4int& vecLen,
                                      G4double& excitationEnergyGNP,
                                      G4double& excitationEnergyDTA,
                                      const G4HEVector& incidentParticle,
                                      const G4HEVector& targetParticle,
                                      G4double atomicWeight,
                                      G4double atomicNumber);            

   void QuasiElasticScattering(G4bool& successful,
                               G4HEVector pv[],
                               G4int& vecLen,
                               G4double& excitationEnergyGNP,
                               G4double& excitationEnergyDTA, 
                               const G4HEVector& incidentParticle,
                               const G4HEVector& targetParticle,
                               G4double atomicWeight,
                               G4double atomicNumber);

   void ElasticScattering(G4bool& successful,
                          G4HEVector pv[],
                          G4int& vecLen,                      
                          const G4HEVector& incidentParticle,
                          G4double atomicWeight,
                          G4double atomicNumber); 

   G4int rtmi(G4double *x, G4double xli, G4double xri, G4double eps,
              G4int iend, G4double aa, G4double bb, G4double cc, 
              G4double dd, G4double rr);

   G4double fctcos(G4double t, G4double aa, G4double bb,G4double cc, 
                   G4double dd, G4double rr);     
        
   void StrangeParticlePairProduction(const G4double availableEnergy,
                                      const G4double centerOfMassEnergy,
                                      G4HEVector pv[],
                                      G4int& vecLen,
                                      const G4HEVector& incidentParticle,
                                      const G4HEVector& targetParticle); 

   G4double NBodyPhaseSpace(const G4double totalEnergy,
                            const G4bool   constantCrossSection,
                            G4HEVector pv[],
                            G4int &vecLen);
     
   G4double NBodyPhaseSpace(G4int npart, 
                            G4HEVector pv[],
                            G4double wmax,
                            G4double wfcn,
                            G4int maxtrial,
                            G4int ntrial);

   G4double gpdk(G4double a, G4double b, G4double c);

   void QuickSort(G4double arr[], const G4int lidx, const G4int ridx);

   G4double Alam(G4double a, G4double b, G4double c);
 
   G4double CalculatePhaseSpaceWeight( G4int npart);
       
   G4double normal(void);
   G4double GammaRand(G4double avalue);
   G4double Erlang(G4int mvalue);
   G4int Poisson(G4double x);
   void SetParticles(void);

   G4HEVector PionPlus;
   G4HEVector PionZero;
   G4HEVector PionMinus;             
   G4HEVector KaonPlus;
   G4HEVector KaonZero;
   G4HEVector AntiKaonZero;             
   G4HEVector KaonMinus;
   G4HEVector KaonZeroShort; 
   G4HEVector KaonZeroLong;
   G4HEVector Proton;
   G4HEVector AntiProton;
   G4HEVector Neutron;             
   G4HEVector AntiNeutron;
   G4HEVector Lambda;
   G4HEVector AntiLambda;             
   G4HEVector SigmaPlus;
   G4HEVector SigmaZero; 
   G4HEVector SigmaMinus;
   G4HEVector AntiSigmaPlus;
   G4HEVector AntiSigmaZero; 
   G4HEVector AntiSigmaMinus;
   G4HEVector XiZero;             
   G4HEVector XiMinus;
   G4HEVector AntiXiZero; 
   G4HEVector AntiXiMinus;
   G4HEVector OmegaMinus;
   G4HEVector AntiOmegaMinus; 
   G4HEVector Deuteron;
   G4HEVector Triton;
   G4HEVector Alpha;                         
   G4HEVector Gamma;
};

#endif
