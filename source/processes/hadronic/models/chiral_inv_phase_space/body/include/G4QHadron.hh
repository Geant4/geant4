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
// $Id: G4QHadron.hh,v 1.26 2003-12-09 15:38:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QHadron ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Hadrons generated by the CHIPS Model
// ------------------------------------------------------

#ifndef G4QHadron_h
#define G4QHadron_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "G4QParticle.hh"

class G4QHadron
{
public:
  // Constructors
  G4QHadron();                                      // Default Constructor
  G4QHadron(G4LorentzVector p);                     // Kinematical Constructor
  G4QHadron(G4int PDGcode, G4LorentzVector p=G4LorentzVector(0.,0.,0.,0.));// CHIPS World Hadron
  G4QHadron(G4QPDGCode QPDG, G4LorentzVector p=G4LorentzVector(0.,0.,0.,0.));// CHIPS World Had.
  G4QHadron(G4QContent QC, G4LorentzVector p=G4LorentzVector(0.,0.,0.,0.));// QC defined Hadron
  G4QHadron(G4int PDG, G4double m, G4QContent QC);  // Constructor for a Chipolino or a Quasmon
  G4QHadron(G4QPDGCode QPDG, G4double m, G4QContent QC);  // Con. for a Chipolino or a Quasmon
  G4QHadron(G4int PDG, G4LorentzVector p, G4QContent QC);// Constructor for Chipolino or Quasmon
  G4QHadron(G4QPDGCode QPDG, G4LorentzVector p, G4QContent QC);// Con. for Chipolino or Quasmon
  G4QHadron(G4QParticle* pPart, G4double maxM);     // Constructor for a resonance with RANDOM m
  G4QHadron(const G4QHadron& right);                // Copy constructor by object
  G4QHadron(const G4QHadron* right);                // Copy constructor by pointer
  virtual ~G4QHadron();                             // Destructor
  // Operators
  const G4QHadron& operator=(const G4QHadron& right);
  G4bool           operator==(const G4QHadron& right) const;
  G4bool           operator!=(const G4QHadron& right) const;
  // Selectors
  G4int                 GetPDGCode()      const;    // Get PDG code of the Hadron
  G4int                 GetQCode()        const;    // Get Q code of the Hadron
  G4QPDGCode            GetQPDG()         const;    // Get QPDG of the Hadron
  G4LorentzVector       Get4Momentum()    const;    // Get 4-mom of the Hadron
  G4QContent            GetQC()           const;    // Get private quark content
  G4double              GetMass()         const;    // Get a mass of the Hadron
  G4double              GetMass2()        const;    // Get an m^2 value for the Hadron
  G4double              GetWidth()        const;    // Get Width of Hadron 
  G4int                 GetNFragments()   const;    // Get a#of Fragments of this Hadron
  G4int                 GetCharge()       const;    // Get Charge of the Hadron
  G4int                 GetStrangeness()  const;    // Get Strangeness of the Hadron
  G4int                 GetBaryonNumber() const;    // Get Baryon Number of the Hadron
  // Modifiers
  void SetQPDG(const G4QPDGCode& QPDG);             // Set QPDG of the Hadron
  void Set4Momentum(const G4LorentzVector& aMom);   // Set 4-mom of the Hadron
  void SetQC(const G4QContent& newQC);              // Set new private quark content
  void SetNFragments(const G4int& nf);              // Set a#of Fragments of this Hadron
  void   NegPDGCode();
  void MakeAntiHadron();                            // Make AntiHadron of this Hadron
  // General
  G4double RandomizeMass(G4QParticle* pPart, G4double maxM); // Randomize a mass value
  G4bool TestRealNeutral();
  G4bool DecayIn2(G4LorentzVector& f4Mom, G4LorentzVector& s4Mom);
  G4bool CorMDecayIn2(G4double corM, G4LorentzVector& fr4Mom);//->This(new mass corM)+fr4Mom(correct)
  G4bool CorEDecayIn2(G4double corE, G4LorentzVector& fr4Mom);//->This(E+=corE,P)+fr4Mom(fE-=corE,fP)
  G4bool RelDecayIn2(G4LorentzVector& f4Mom, G4LorentzVector& s4Mom, G4LorentzVector& dir,
                  G4double maxCost = 1., G4double minCost = -1.);
  G4bool DecayIn3(G4LorentzVector& f4Mom, G4LorentzVector& s4Mom, G4LorentzVector& t4Mom);

private:
  // Private methods
  void DefineQC(G4int PDGCode);

private:  
  G4QPDGCode            theQPDG;                    // Instance of QPDG for the Hadron
  G4LorentzVector       theMomentum;                // The 4-mom of Hadron
  G4QContent            valQ;                       // QC ( @@ ?? for Quasmon and Chipolino)
  G4int                 nFragm;                     // =0 - stable, =N - decayed in N particles
};

inline G4bool G4QHadron::operator==(const G4QHadron &rhs) const {return this==&rhs;}
inline G4bool G4QHadron::operator!=(const G4QHadron &rhs) const {return this!=&rhs;}
 
inline G4int           G4QHadron::GetPDGCode()      const   {return theQPDG.GetPDGCode();}
inline G4int           G4QHadron::GetQCode()        const   {return theQPDG.GetQCode();}
inline G4QPDGCode      G4QHadron::GetQPDG()         const   {return theQPDG;}
inline G4QContent      G4QHadron::GetQC()           const   {return valQ;}
inline G4LorentzVector G4QHadron::Get4Momentum()    const   {return theMomentum;}
inline G4int           G4QHadron::GetNFragments()   const   {return nFragm;}
//@@ This is an example how to make other inline selectors for the 4-Momentum of the Hadron
inline G4double        G4QHadron::GetMass()         const   {return theMomentum.m();}
inline G4double        G4QHadron::GetMass2()        const   {return theMomentum.m2();}
//@@ This is an example how to make other inline selectors for the Hadron
inline G4int           G4QHadron::GetCharge()       const   {return valQ.GetCharge();}
inline G4int           G4QHadron::GetStrangeness()  const   {return valQ.GetStrangeness();}
inline G4int           G4QHadron::GetBaryonNumber() const   {return valQ.GetBaryonNumber();}

inline void            G4QHadron::MakeAntiHadron()          {if(TestRealNeutral()) NegPDGCode();}
inline void            G4QHadron::SetQPDG(const G4QPDGCode& newQPDG)
{
  theQPDG  = newQPDG;
  G4int PDG= newQPDG.GetPDGCode();
  G4int Q  = newQPDG.GetQCode();
  //G4cout<<"G4QHadron::SetQPDG is called with PDGCode="<<PDG<<", QCode="<<Q<<G4endl;
  if     (Q>-1) valQ=theQPDG.GetQuarkContent();
  else if(PDG>80000000)DefineQC(PDG);
  else
  {
    G4cerr<<"***G4QHadron::SetQPDG: QPDG="<<newQPDG<<G4endl;
    throw G4QException("***G4QHadron::SetQPDG: Impossible QPDG");
  }
}
inline void   G4QHadron::SetQC(const G4QContent& newQC)              {valQ=newQC;}
inline void   G4QHadron::Set4Momentum(const G4LorentzVector& aMom)   {theMomentum=aMom;}
inline void   G4QHadron::SetNFragments(const G4int& nf)              {nFragm=nf;}
inline void   G4QHadron::NegPDGCode()                   {theQPDG.NegPDGCode(); valQ.Anti();}
inline G4bool G4QHadron::TestRealNeutral()              { return theQPDG.TestRealNeutral();}
inline void   G4QHadron::DefineQC(G4int PDGCode)
{
  //G4cout<<"G4QHadron::DefineQC is called with PDGCode="<<PDGCode<<G4endl;
  G4int szn=PDGCode-90000000;
  G4int ds=0;
  G4int dz=0;
  G4int dn=0;
  if(szn<-100000)
  {
    G4int ns=(-szn)/1000000+1;
    szn+=ns*1000000;
    ds+=ns;
  }
  else if(szn<-100)
  {
    G4int nz=(-szn)/1000+1;
    szn+=nz*1000;
    dz+=nz;
  }
  else if(szn<0)
  {
    G4int nn=-szn;
    szn=0;
    dn+=nn;
  }
  G4int sz =szn/1000;
  G4int n  =szn%1000;
  if(n>700)
  {
    n-=1000;
    dz--;
  }
  G4int z  =sz%1000-dz;
  if(z>700)
  {
    z-=1000;
    ds--;
  }
  G4int Sq =sz/1000-ds;
  G4int zns=z+n+Sq;
  G4int Dq=n+zns;
  G4int Uq=z+zns;
  if      (Dq<0&&Uq<0&&Sq<0)valQ=G4QContent(0 ,0 ,0 ,-Dq,-Uq,-Sq);
  else if (Uq<0&&Sq<0)      valQ=G4QContent(Dq,0 ,0 ,0  ,-Uq,-Sq);
  else if (Dq<0&&Sq<0)      valQ=G4QContent(0 ,Uq,0 ,-Dq,0  ,-Sq);
  else if (Dq<0&&Uq<0)      valQ=G4QContent(0 ,0 ,Sq,-Dq,-Uq,0  );
  else if (Uq<0)            valQ=G4QContent(Dq,0 ,Sq,0  ,-Uq,0  );
  else if (Sq<0)            valQ=G4QContent(Dq,Uq,0 ,0  ,0  ,-Sq);
  else if (Dq<0)            valQ=G4QContent(0 ,Uq,Sq,-Dq,0  ,0  );
  else                      valQ=G4QContent(Dq,Uq,Sq,0  ,0  ,0  );
}
#endif

