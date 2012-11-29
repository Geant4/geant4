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
//      ---------------- G4QNucleus ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for the nuclei and nuclear environment of the CHIPS Model
// -----------------------------------------------------------------------
//  Short description: a class describing properties of nuclei, which
//  are necessary for the CHIPS Model.
// -----------------------------------------------------------------------

#ifndef G4QNucleus_h
#define G4QNucleus_h 1

#include <utility>
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4RandomDirection.hh"
#include "G4QCandidateVector.hh"
#include "G4QHadronVector.hh"
#include "G4LorentzRotation.hh"
#include "G4QChipolino.hh"

class G4QNucleus : public G4QHadron
{
public:
  G4QNucleus();                                            // Default Constructor
  G4QNucleus(G4int nucPDG);                                // At Rest PDG-Constructor
  G4QNucleus(G4LorentzVector p, G4int nucPDG);             // Full PDG-Constructor
  G4QNucleus(G4QContent nucQC);                            // At Rest QuarkCont-Constructor
  G4QNucleus(G4QContent nucQC, G4LorentzVector p);         // Full QuarkCont-Constructor
  G4QNucleus(G4int z, G4int n, G4int s=0);                 // At Rest ZNS-Constructor
  G4QNucleus(G4int z, G4int n, G4int s, G4LorentzVector p);// Full ZNS-Constructor
  G4QNucleus(G4QNucleus* right, G4bool cop3D = false);     // Copy Constructor by pointer
  G4QNucleus(const G4QNucleus &right, G4bool cop3D=false); // Copy Constructor by value
  ~G4QNucleus();                                           // Public Destructor
  // Overloaded Operators
  const G4QNucleus& operator=(const G4QNucleus& right);
  G4bool operator==(const G4QNucleus &right) const {return this==&right;}
  G4bool operator!=(const G4QNucleus &right) const {return this!=&right;}
  // Specific Selectors
  G4int      GetPDG()      const {return 90000000+1000*(1000*S+Z)+N;}// PDG Code of Nucleus
  G4int      GetZ()        const {return Z;}				    // Get a#of protons
  G4int      GetN()        const {return N;}				    // Get a#of neutrons
  G4int      GetS()        const {return S;}				    // Get a#of lambdas
  G4int      GetA()        const {return Z+N+S;}    // Get A of the nucleus
  G4int      GetDZ()       const {return dZ;}						 // Get a#of protons in dense region
  G4int      GetDN()       const {return dN;}						 // Get a#of neutrons in dense region
  G4int      GetDS()       const {return dS;}						 // Get a#of lambdas in dense region
  G4int      GetDA()       const {return dZ+dN+dS;} // Get A of the dense part of nucleus
  G4int      GetMaxClust() const {return maxClust;} // Get Max BarNum of Clusters
  G4double   GetProbability(G4int bn=0) const {return probVect[bn];} // clust(BarN)probabil
  G4double   GetMZNS()     const {return GetQPDG().GetNuclMass(Z,N,S);} // not H or Q
  G4double   GetTbIntegral(); // Calculate the integral of T(b)
  G4double   GetGSMass()   const {return GetQPDG().GetMass();}//Nucleus GSMass (not Hadron)
  G4QContent GetQCZNS()    const                    // Get ZNS quark content of Nucleus
  {
    if(S>=0) return G4QContent(Z+N+N+S,Z+Z+N+S,S,0,0,0);
    else     return G4QContent(Z+N+N+S,Z+Z+N+S,0,0,0,-S);
  }
  G4int      GetNDefMesonC() const{return nDefMesonC;}; // max#of predefed mesonCandidates
  G4int      GetNDefBaryonC()const{return nDefBaryonC;};// max#of predefed baryonCandidates
  G4double GetDensity(const G4ThreeVector&aPos) {return rho0*GetRelativeDensity(aPos);}
  G4double GetRho0()                 {return rho0;} // One nucleon prob-density 
  G4double GetRelativeDensity(const G4ThreeVector& aPosition); // Densyty/rho0
  G4double GetRelWSDensity(const G4double& r)       // Wood-Saxon rho/rho0(r)
                                        {return 1./(1.+std::exp((r-radius)/WoodSaxonSurf));}    
  G4double GetRelOMDensity(const G4double& r2){return std::exp(-r2/radius);} // OscModelRelDens
  G4double GetRadius(const G4double maxRelativeDenisty=0.5); // Radius of %ofDensity
  G4double GetOuterRadius();                        // Get radius of the most far nucleon
  G4double GetDeriv(const G4ThreeVector& point);    // Derivitive of density
  G4double GetFermiMomentum(G4double density);      // Returns modul of FermyMomentum(dens)
  G4QHadron* GetNextNucleon()
  {
    //G4cout<<"G4QNucleus::GetNextNucleon: cN="<<currentNucleon<<", A="<<GetA()<<G4endl;
    return (currentNucleon>=0&&currentNucleon<GetA()) ? theNucleons[currentNucleon++] : 0;
  }
  void SubtractNucleon(G4QHadron* pNucleon); // Subtract the nucleon from the 3D Nucleus
  void DeleteNucleons();                     // Deletes all residual nucleons
  G4LorentzVector GetNucleons4Momentum()
  {
    G4LorentzVector sum(0.,0.,0.,0.);
    for(unsigned i=0; i<theNucleons.size(); i++) sum += theNucleons[i]->Get4Momentum();
    sum.setE(std::sqrt(sqr(GetGSMass())+sum.v().mag2())); // Energy is corrected !
    return sum;
  }
  std::vector<G4double> const* GetBThickness() {return &Tb;} // T(b) function, step .1 fm

  // Specific Modifiers
  G4bool     EvaporateBaryon(G4QHadron* h1,G4QHadron* h2); // Evaporate Baryon from Nucleus
  void       EvaporateNucleus(G4QHadron* hA, G4QHadronVector* oHV);// Evaporate Nucleus
  //void DecayBaryon(G4QHadron* dB, G4QHadronVector* oHV); // gamma+N or Delt->N+Pi @@later
  void       DecayDibaryon(G4QHadron* dB, G4QHadronVector* oHV);   // deuteron is kept
  void       DecayAntiDibaryon(G4QHadron* dB, G4QHadronVector* oHV);// antiDeuteron is kept
  void       DecayIsonucleus(G4QHadron* dB, G4QHadronVector* oHV); // nP+(Pi+) or nN+(Pi-)
  void       DecayMultyBaryon(G4QHadron* dB, G4QHadronVector* oHV);// A*p, A*n or A*L
  void       DecayAntiStrange(G4QHadron* dB, G4QHadronVector* oHV);// nuclei with K+/K0
  void       DecayAlphaBar(G4QHadron* dB, G4QHadronVector* oHV);   // alpha+p or alpha+n
  void       DecayAlphaDiN(G4QHadron* dB, G4QHadronVector* oHV);   // alpha+p+p
  void       DecayAlphaAlpha(G4QHadron* dB, G4QHadronVector* oHV); // alpha+alpha
  G4int      SplitBaryon();                         // Is it possible to split baryon/alpha
  G4int      HadrToNucPDG(G4int hPDG);              // Converts hadronic PDGCode to nuclear
  G4int      NucToHadrPDG(G4int nPDG);              // Converts nuclear PDGCode to hadronic
  G4bool     Split2Baryons();                       // Is it possible to split two baryons?
  void       ActivateBThickness();                  // Calculate T(b) for nucleus (db=.1fm)
  G4double   GetBThickness(G4double b);             // Calculates T(b)
  G4double   GetThickness(G4double b);              // Calculates T(b)/rho(0)
  void       InitByPDG(G4int newPDG);               // Init existing nucleus by new PDG
  void       InitByQC(G4QContent newQC)             // Init existing nucleus by new QCont
                                {G4int PDG=G4QPDGCode(newQC).GetPDGCode(); InitByPDG(PDG);}
  void       IncProbability(G4int bn);              // Add one cluster to probability
  void       Increase(G4int PDG, G4LorentzVector LV = G4LorentzVector(0.,0.,0.,0.));
  void       Increase(G4QContent QC, G4LorentzVector LV = G4LorentzVector(0.,0.,0.,0.));
  void       Reduce(G4int PDG);                     // Reduce Nucleus by PDG fragment
  void       CalculateMass() {Set4Momentum(G4LorentzVector(0.,0.,0.,GetGSMass()));}
  void       SetMaxClust(G4int maxC){maxClust=maxC;}// Set Max BarNum of Clusters
  void       InitCandidateVector(G4QCandidateVector& theQCandidates,
                                 G4int nM=45, G4int nB=72, G4int nC=117);
  void       PrepareCandidates(G4QCandidateVector& theQCandidates, G4bool piF=false, G4bool
                               gaF=false, G4LorentzVector LV=G4LorentzVector(0.,0.,0.,0.));
  G4int      UpdateClusters(G4bool din);            // Return a#of clusters & calc.probab's
  G4QNucleus operator+=(const G4QNucleus& rhs);     // Add a cluster to the  nucleus
  G4QNucleus operator-=(const G4QNucleus& rhs);     // Subtract a cluster from a nucleus
  G4QNucleus operator*=(const G4int& rhs);          // Multiplication of the Nucleus
  G4bool StartLoop();                               // returns size of theNucleons (cN=0)
  G4bool ReduceSum(G4ThreeVector* vectors, G4ThreeVector sum);// Reduce zero-sum of vectors
  void SimpleSumReduction(G4ThreeVector* vectors, G4ThreeVector sum); // Reduce zero-V-sum
  void DoLorentzBoost(const G4LorentzVector& theBoost) // Boost nucleons by 4-vector
  {
    theMomentum.boost(theBoost);
    for(unsigned i=0; i<theNucleons.size(); i++) theNucleons[i]->Boost(theBoost);
  }
  void DoLorentzRotation(const G4LorentzRotation& theLoRot) // Lorentz Rotate nucleons
  {
    theMomentum=theLoRot*theMomentum;
    for(unsigned i=0; i<theNucleons.size(); i++) theNucleons[i]->LorentzRotate(theLoRot);
  }
  void DoLorentzBoost(const G4ThreeVector& theBeta)// Boost nucleons by v/c
  {
    theMomentum.boost(theBeta);
    for(unsigned i=0; i<theNucleons.size(); i++) theNucleons[i]->Boost(theBeta);
  }
  void DoLorentzContraction(const G4LorentzVector&B){DoLorentzContraction(B.vect()/B.e());}
  void DoLorentzContraction(const G4ThreeVector& theBeta); // Lorentz Contraction by v/c
  void DoTranslation(const G4ThreeVector& theShift); // Used only in G4QFragmentation

  // Static functions
  static void SetParameters(G4double fN=.1,G4double fD=.05, G4double cP=4., G4double mR=1.,
                            G4double nD=.8*CLHEP::fermi);

  // Specific General Functions
  G4int RandomizeBinom(G4double p,G4int N);         // Randomize according to Binomial Law
  G4double CoulombBarGen(const G4double& rZ, const G4double& rA, const G4double& cZ,
                         const G4double& cA);        // CoulombBarrier in MeV (General)
  G4double CoulombBarrier(const G4double& cZ=1, const G4double& cA=1, G4double dZ=0.,
                          G4double dA=0.);          // CoulombBarrier in MeV
  G4double FissionCoulombBarrier(const G4double& cZ, const G4double& cA, G4double dZ=0.,
                          G4double dA=0.);          // Fission CoulombBarrier in MeV
  G4double BindingEnergy(const G4double& cZ=0, const G4double& cA=0, G4double dZ=0.,
                         G4double dA=0.);
  G4double CoulBarPenProb(const G4double& CB, const G4double& E, const G4int& C,
                          const G4int& B);
  std::pair<G4double, G4double> ChooseImpactXandY(G4double maxImpact); // Randomize bbar
  void ChooseNucleons();                            // Initializes 3D Nucleons
  void ChoosePositions();                           // Initializes positions of 3D nucleons
  void ChooseFermiMomenta();                        // Initializes FermyMoms of 3D nucleons
  void InitDensity();                               // Initializes density distribution
  void Init3D();                                    // automatically starts the LOOP
private:  
  // Specific Encapsulated Functions
  void       SetZNSQC(G4int z, G4int n, G4int s);   // Set QC, using Z,N,S
  G4QNucleus GetThis() const {return G4QNucleus(Z,N,S);} // @@ Check for memory leak

// Body
private:
  // Static Parameters
  static const G4int nDefMesonC =45;
  static const G4int nDefBaryonC=72;
  //
  static G4double freeNuc;        // probability of the quasi-free baryon on surface
  static G4double freeDib;        // probability of the quasi-free dibaryon on surface
  static G4double clustProb;      // clusterization probability in dense region
  static G4double mediRatio;      // relative vacuum hadronization probability
  static G4double nucleonDistance;// Distance between nucleons (0.8 fm)
  static G4double WoodSaxonSurf;  // Surface parameter of Wood-Saxon density (0.545 fm)
  // The basic  
  G4int Z;                        // Z of the Nucleus
  G4int N;                        // N of the Nucleus
  G4int S;                        // S of the Nucleus
  // The secondaries
  G4int dZ;                       // Z of the dense region of the nucleus
  G4int dN;                       // N of the dense region of the nucleus
  G4int dS;                       // S of the dense region of the nucleus
  G4int maxClust;                 // Baryon Number of the last calculated cluster
  G4double probVect[256];         // Cluster probability ("a#of issues" can be real) Vector
  // 3D
  G4QHadronVector theNucleons;    // Vector of nucleons of which the Nucleus consists of
  G4int currentNucleon;           // Current nucleon for the NextNucleon (? M.K.)
  G4double rho0;                  // Normalazation density
  G4double radius;                // Nuclear radius
  std::vector<G4double> Tb;       // T(b) function with step .1 fm (@@ make .1 a parameter)
  G4bool TbActive;                // Flag that the T(b) is activated
  G4bool RhoActive;               // Flag that the Density is activated
};

std::ostream& operator<<(std::ostream& lhs, G4QNucleus& rhs);
std::ostream& operator<<(std::ostream& lhs, const G4QNucleus& rhs);

#endif
