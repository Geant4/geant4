// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEInelastic.hh,v 1.4 2000-07-09 09:48:10 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// All routines pass the compiler and give
// some reasonable numbers as output.
// A statistically significant comparison
// with GEANT3 and experimental data has
// still to be done. 

#ifndef G4HEInelastic_h
#define G4HEInelastic_h 1

#include "G4HEVector.hh"
#include "G4HadronicInteraction.hh"

class G4HEInelastic : public G4HadronicInteraction
{
 public: 
         G4HEInelastic()
            { 
              SetMinEnergy(20*GeV);
              SetMaxEnergy(10*TeV);
              MAXPART = 512;
              verboseLevel = 0;
              SetParticles();
              conserveEnergy = false;
            };
        ~G4HEInelastic(){ };
         
         void       SetMaxNumberOfSecondaries( const G4int maxnumber ) 
                         { MAXPART = maxnumber;}      
 
         void       SetVerboseLevel( const G4int level)
                         { verboseLevel = level;}

         G4int      verboseLevel; 
         G4int      MAXPART;
         G4bool     conserveEnergy;
         void       ForceEnergyConservation(G4bool energyConservation)
                         { conserveEnergy = energyConservation;}
         G4bool     EnergyConservation(void)
                         { return conserveEnergy;} 

         void       FillParticleChange(G4HEVector pv[], G4int aVecLength);

         G4double   pmltpc(G4int np, G4int nm, G4int nz, G4int n, G4double b, G4double c);

         G4int      Factorial(G4int n); 

         G4double   NuclearInelasticity(G4double incidentKineticEnergy,
                                        G4double atomicWeight,
                                        G4double atomicNumber);
         G4double   NuclearExcitation(G4double  incidentKineticEnergy,
                                      G4double  atomicWeight,
                                      G4double  atomicNumber,
                                      G4double& excitationEnergyCascade,
                                      G4double& excitationEnergyEvaporation); 
   
         void       HighEnergyCascading(G4bool &successful,
                                        G4HEVector pv[],
                                        G4int &vecLen,
                                        G4double &excitationEnergyGNP,
                                        G4double &excitationEnergyDTA, 
                                        G4HEVector incidentParticle,
                                        G4HEVector targetParticle,
                                        G4double atomicWeight,
                                        G4double atomicNumber);                        

         void       HighEnergyClusterProduction(G4bool &successful,
                                                G4HEVector pv[],
                                                G4int &vecLen,
                                                G4double &excitationEnergyGNP,
                                                G4double &excitationEnergyDTA, 
                                                G4HEVector incidentParticle,
                                                G4HEVector targetParticle,
                                                G4double atomicWeight,
                                                G4double atomicNumber);             

         void       TuningOfHighEnergyCascading( G4HEVector pv[],
                                                 G4int &vecLen,
                                                 G4HEVector incidentParticle,
                                                 G4HEVector targetParticle,
                                                 G4double atomicWeight,
                                                 G4double atomicNumber);   

         void       MediumEnergyCascading(G4bool &successful,
                                          G4HEVector pv[],
                                          G4int &vecLen,
                                          G4double &excitationEnergyGNP,
                                          G4double &excitationEnergyDTA,
                                          G4HEVector incidentParticle,
                                          G4HEVector targetParticle,
                                          G4double atomicWeight,
                                          G4double atomicNumber);            

         void       MediumEnergyClusterProduction(G4bool &successful,
                                                  G4HEVector pv[],
                                                  G4int &vecLen,
                                                  G4double &excitationEnergyGNP,
                                                  G4double &excitationEnergyDTA,
                                                  G4HEVector incidentParticle,
                                                  G4HEVector targetParticle,
                                                  G4double atomicWeight,
                                                  G4double atomicNumber);            

         void       QuasiElasticScattering(G4bool &successful,
                                           G4HEVector pv[],
                                           G4int &vecLen,
                                           G4double &excitationEnergyGNP,
                                           G4double &excitationEnergyDTA, 
                                           G4HEVector incidentParticle,
                                           G4HEVector targetParticle,
                                           G4double atomicWeight,
                                           G4double atomicNumber);

         void       ElasticScattering(G4bool &successful,
                                      G4HEVector pv[],
                                      G4int &vecLen,                      
                                      G4HEVector incidentParticle,
                                      G4double atomicWeight,
                                      G4double atomicNumber); 

         G4int      rtmi(G4double *x, G4double xli, G4double xri, G4double eps,
                         G4int iend,
                         G4double aa, G4double bb, G4double cc, G4double dd, G4double rr);
         G4double   fctcos(G4double t, G4double aa, G4double bb,G4double cc, G4double dd, 
                           G4double rr);     
        
         void       StrangeParticlePairProduction(const G4double availableEnergy,
                                                  const G4double centerOfMassEnergy,
                                                  G4HEVector pv[],
                                                  G4int &vecLen,
                                                  G4HEVector incidentParticle,
                                                  G4HEVector targetParticle); 

         G4double   NBodyPhaseSpace(const G4double totalEnergy,
                                    const G4bool   constantCrossSection,
                                    G4HEVector pv[],
                                    G4int &vecLen);     
         G4double   NBodyPhaseSpace(G4int npart, 
                                    G4HEVector pv[],
                                    G4double wmax,
                                    G4double wfcn,
                                    G4int maxtrial,
                                    G4int ntrial);

         G4double   gpdk(G4double a, G4double b, G4double c);
         void       QuickSort(G4double arr[], const G4int lidx, const G4int ridx);
         G4double   Alam(G4double a, G4double b, G4double c); 
         G4double   CalculatePhaseSpaceWeight( G4int npart);
       
         G4double   normal(void);
         G4double   GammaRand(G4double avalue);
         G4double   Erlang(G4int mvalue);
         G4int      Poisson(G4double x);
         void       SetParticles(void);

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
                                         

