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
 // $Id: G4PreEquilibrium.hh,v 1.9 2001/07/11 10:03:54 gunter Exp $
 // GEANT4 tag $Name: geant4-05-00 $
 //
 // Hadronic Process: Pre-equilibrium HETC 
 // Joseph L. Chuma, TRIUMF, 24-Mar-2000

#ifndef G4PREEQUILIBRIUM
#define G4PREEQUILIBRIUM
 
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
#include "g4std/vector"
 
 class G4PreEquilibrium
 {
 private:
   struct EXITON
   {
     G4int protons;
     G4int neutrons;
     G4int hols;     // not sure what these are
   };
   
   struct MASTRUCT
   {
     G4ThreeVector momentum;
     G4double array[6];
   };
   
   struct PARZS
   {
     G4int particleType;
     G4double energy;
     G4double angle1;
     G4double angle2;
     G4double what;
     G4double charge;
   };
   
   struct IPSTRUCT
   {
     G4int proton;
     G4int ip1;
     G4int ip2;
     G4int nuclearZone;
     G4int ip4;          // values from 0 to N
   };
   
#define G4Vector G4std::vector
   typedef G4Vector< G4DynamicParticle* > DPvector;
   typedef G4Vector< G4double > Dvector;
   typedef G4Vector< G4int > Ivector;
   typedef G4Vector< MASTRUCT > MAvector;
   typedef G4Vector< PARZS > PARZvector;
   typedef G4Vector< IPSTRUCT > IPvector;
   
 public:
   G4PreEquilibrium( G4Nucleus *aNucleus ) 
       : theNucleus(aNucleus)
    {}
   
   ~G4PreEquilibrium()
    {}
   
   G4PreEquilibrium operator=( const G4PreEquilibrium & );
   
   G4PreEquilibrium( const G4PreEquilibrium & );
   
 private:
   void mashnk( G4int &, const G4double, const G4double, const G4double,
                const G4double, const G4double, const G4double, const G4double );
   
   G4double cemgeo( const G4int, const G4int );
   
   G4double bf( const G4double, const G4double, const G4double );
   
   G4double quadraticInterpolation( const G4double, const G4double *,
                                    const G4double *, const G4int );
   
   G4double massDefect( const G4double, const G4double ) const;
   
   G4double gamagu( const G4int );
   
   G4double tkin( const G4int, const G4double );
   
   void transitionRates( const G4double, const G4double,
                         G4double &, G4double &, G4double & );
   
   G4double tkinm1( const G4int, const G4double, const G4double );
   
   G4double barfit( const G4double, const G4double, const G4int );
   
   G4bool cascem( G4ThreeVector );
   
   void cascem1( G4ThreeVector, MASTRUCT &, IPSTRUCT &, G4int & );
   
   G4bool pointe( MASTRUCT &, IPSTRUCT &, G4double *, G4int *,
                  G4ThreeVector, G4double, G4double, G4double,
                  G4double, G4double, const G4double );
   
   G4double geometricalParticlePath( const MASTRUCT, const IPSTRUCT );
   
   G4bool refrac( const G4int, MASTRUCT &, IPSTRUCT & ) const;
   
   G4double wopt( const G4double, const G4double, const G4int );
   
   void typint( MASTRUCT &, IPSTRUCT &, G4double *, G4int *, G4ThreeVector,
                G4double, G4double, G4double, G4double, G4double, G4int, G4int & );
   
   G4int elasticScattering( const G4ThreeVector, const G4double, const G4double, const MASTRUCT,
                            const IPSTRUCT, const G4int *, const G4int, const G4int,
                            const G4int, const G4int, const G4int, const G4int );
   
   G4bool inelasticScattering( MASTRUCT &, IPSTRUCT &, G4int *, G4int, G4int, G4int,
                               G4int, G4int, G4ThreeVector, G4double, G4double, G4int, G4int );
   
   G4double wim( const G4int, const MASTRUCT, const IPSTRUCT );
   
   void pauliPrinciple( MASTRUCT &, IPSTRUCT &, G4ThreeVector, G4int, G4int, const G4int );
   
   void absorption( const MASTRUCT, const IPSTRUCT, G4double *, G4int, G4int &, G4ThreeVector );
   
   void chargeInAbsorption( const G4int, const G4int, G4int &, G4int & );
   
   void momentaCalc( const MASTRUCT, const G4ThreeVector, const G4double, G4ThreeVector,
                     G4ThreeVector, const G4double, const G4double, const G4double );
   
   void partnerSelection( const MASTRUCT, const IPSTRUCT, G4double *, G4int * );
   
   void chinel( IPSTRUCT &, const G4int, const G4int, const G4int, const G4int,
                const G4int, const G4int, const G4double, const G4int, G4int * );
   
   void direction( G4ThreeVector, const G4double, const G4double, const G4int,
                   const G4int, const G4int, const MASTRUCT, G4int &, const G4int );
   
   G4int coefficientTypeA( const G4int, const G4int, const G4int );
   
   G4bool vmnsp( const MASTRUCT, const IPSTRUCT, const G4double, const G4int,
                 G4int &, const G4int, const G4int, const G4double );
   
   G4double secondaryParticleMomentum( const G4int, const G4double );
   
   G4int coefficientTypeB( const G4int, const G4int, const G4int );
   
   void statisticalModel( const G4double, const G4ThreeVector, const MASTRUCT,
                          const G4int *, const G4int );
   
   void isocem( const G4double, const G4ThreeVector, const G4double, const MASTRUCT,
                const G4int *, const G4int );
   
   G4ThreeVector cms( const G4ThreeVector, const G4ThreeVector, const G4double ) const;
   
   void bertcem( const G4int, const G4int );
   
   G4double fintfis( const G4double, const G4double, const G4double, const G4double );
   
   G4double fints2( G4double, G4double, G4double, G4double );
   
   G4double erupcem( DPvector &, DPvector &, DPvector &, DPvector &, DPvector &, DPvector & );
   
   void precof();
   
   G4double arfaf( G4double * );
   
   G4double fam( const G4double, const G4double, const G4double );
   
   G4double gameqf( const G4int, const G4double, const G4double, const G4double, const G4double );
   
   G4double poten( const G4int, const IPSTRUCT ) const;
   
   G4double cinema( const G4ThreeVector, const G4ThreeVector, G4ThreeVector &,
                    G4double &, G4double &, G4double &, G4double &, const G4double ) const;
   
   G4ThreeVector rotation( const G4ThreeVector, const G4ThreeVector, const G4ThreeVector ) const;
   
   G4double crossSectionInterp( const G4double, const G4int );
   
   G4double sigmat( const G4int, const G4int, const G4int,
                   const G4int, const G4int, const G4double );
   
   void slqek( G4int &, G4int &, G4int &, G4int &, G4int &,
               const G4int, const G4int, const G4int,
               const G4int, const G4int, const G4int,
               const G4int, const G4int );
   
   G4double costa( const G4int, const G4double );
   
   G4double cosex( const G4int, const G4double, const G4double );
   
   G4double cosel( const G4int, const G4int, const G4int,
                  const G4double, const G4double );
   
   void lpoly( const G4double, const G4int, G4double * );
   
   inline void SetINDI( const G4bool v )
    { indi = v; }
   
   inline G4bool GetINDI() const
    { return indi; }
   
   inline void SetING( const G4int i )
    { ing = i; }
   
   inline G4int GetING() const
    { return ing; }
   
   inline void SetCM0( const G4double v )
    { cm0 = v; }
   
   inline G4double GetCM0() const
    { return cm0; }
   
   inline void SetT0( const G4double v )
    { t0 = v; }
   
   inline G4double GetT0() const
    { return t0; }
   
   inline void SetT1( const G4double v )
    { t1 = v; }
   
   inline G4double GetT1() const
    { return t1; }
   
   inline void SetT2( const G4double v )
    { t2 = v; }
   
   inline G4double GetT2() const
    { return t2; }
   
   inline void SetME0( const G4double v )
    { me0 = v; }
   
   inline G4double GetME0() const
    { return me0; }
   
   inline void SetMQ0( const G4double v )
    { mq0 = v; }
   
   inline G4double GetMQ0() const
    { return mq0; }
   
   inline void SetA( const G4double v )
    { a = v; }
   
   inline G4double GetA() const
    { return a; }
   
   inline void SetZ( const G4double v )
    { z = v; }
   
   inline G4double GetZ() const
    { return z; }
   
   inline G4double GetUP() const
    { return up; }
   
   inline void SetUP( const G4double v )
    { up = v; }

   inline void SetU( const G4double v )
    { u = v; }
   
   inline G4double GetU() const
    { return u; }
   
   inline G4double GetT1Y( const G4int i ) const
    { return t1y[i]; }
   
   inline G4double GetT2XY( const G4int i ) const
    { return t2xy[i]; }
   
   inline G4double shell( const G4double a, const G4double z )
    {
      // Cameron (Can.J.Phys.35(1957)1021) shell and pairing corrections
      //
      return t1y[int(z)-1] + t2xy[int(a-z)-1];
    }
   
   inline void SetRBIG( const G4int i, const G4double v )
    { rbig[i] = v; }
   
   inline G4double GetRBIG( const G4int i ) const
    { return rbig[i]; }
   
   inline void SetR0( const G4double v )
    { r0 = v; }
   
   inline G4double GetR0() const
    { return r0; }
   
   inline G4double SetBN( const G4double v )
    { bn = v; }
   
   inline G4double GetBN() const
    { return bn; }
   
   inline void SetAC( const G4double v )
    { ac = v; }
   
   inline G4double GetAC() const
    { return ac; }
   
   inline G4int GetN0() const
    { return n0; }
   
   inline void SetN0( const G4int v )
    { n0 = v; }
   
   inline G4int GetH0() const
    { return h0; }
   
   inline void SetH0( const G4int v )
    { h0 = v; }
   
   inline G4int GetP0() const
    { return p0; }
   
   inline void SetP0( const G4int v )
    { p0 = v; }
   
   inline G4int GetPZ0() const
    { return pz0; }
   
   inline void SetPZ0( const G4int v )
    { pz0 = v; }
   
   inline G4ThreeVector GetLXYZ() const
    { return lxyz; }
   
   inline void SetLXYZ( const G4ThreeVector v )
    { lxyz = v; }
   
   inline void SetWF( const G4double v )
    { wf = v; }
   
   inline G4double GetWF() const
    { return wf; }
   
   inline void SetANG1( const G4int i, const G4int j, const G4int k, const G4double v )
    { ang1[i][j][k] = v; }
   
   inline G4double GetANG1( const G4int i, const G4int j, const G4int k ) const
    { return ang1[i][j][k]; }
   
   inline void SetANG2( const G4int i, const G4int j, const G4int k, const G4double v )
    { ang2[i][j][k] = v; }
   
   inline G4double GetANG2( const G4int i, const G4int j, const G4int k ) const
    { return ang2[i][j][k]; }
   
   inline void SetEPART( const G4int i, const G4int j, const G4double v )
    { epart[i][j] = v; }
   
   inline G4double GetEPART( const G4int i, const G4int j ) const
    { return epart[i][j]; }
   
   inline void SetHEPART( const G4int i, const G4int j, const G4double v )
    { hepart[i][j] = v; }
   
   inline G4double GetHEPART( const G4int i, const G4int j ) const
    { return hepart[i][j]; }
   
   inline void SetNHIST( const G4int v )
    { nhist = v; }
   
   inline G4double GetNHIST() const
    { return nhist; } 
   
   inline G4double GetZ_FCOMON( const G4int i, const G4int j ) const
    { return z_fcomon[i][j]; }
   
   inline G4double GetA_FCOMON( const G4int i, const G4int j ) const
    { return a_fcomon[i][j]; }
   
 private:
   static G4double t1y[130];
   static G4double t2xy[200];
   
   static G4double oneThird;
   static G4double twoThirds;
   static G4double fourThirds;
   
   EXITON exitons;
   
   G4ThreeVector lxyz;
   
   Dvector rbig;
   Dvector af;
   Dvector ngen;
   Dvector gb;
   Dvector alj;
   Dvector afj;
   Dvector zfj;
   Dvector rj;
   Dvector vj;
   Dvector bj;
   Dvector r0j;
   Dvector rsm;
   Dvector rhon;
   Dvector rhop;
   Dvector tfp;
   Dvector tfn;
   
   G4Vector< Dvector* > spt;
   
   MAvector pmemo;
   PARZvector parz;
   IPvector imemo;
   
   G4Nucleus *theNucleus;
   
   G4bool indi;
   
   G4int n0;
   G4int nhist;
   G4int h0;
   G4int p0;
   G4int pz0;
   G4int KTOT;
   G4int ing;
   
   G4double up; // excitation energy before evaporation
   G4double ap; // atomic weight before evaporation
   G4double zp; // atomic number before evaporation
   G4double u;  // excitation energy after evaporation
   G4double a;  // atomic weight after evaporation
   G4double z;  // atomic number after evaporation
   
   G4double wf;
   G4double ac;
   G4double exn;
   G4double r0;
   G4double cm0;
   G4double t0;
   G4double t1;
   G4double t2;
   G4double mq0;
   G4double me0;
   G4double bn;
   G4double wam;
   G4double aNucl;
   G4double zNucl;
   
   G4double a_fcomon[10][17];
   G4double z_fcomon[10][17];
   
   G4double epart[100][2];
   G4double hepart[100][4];
   
   G4double ang1[100][3][2];
   G4double ang2[100][3][4];
 };
 
#endif
