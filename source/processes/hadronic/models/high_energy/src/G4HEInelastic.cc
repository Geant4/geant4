// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEInelastic.cc,v 1.2 1999-06-16 04:35:57 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#include "globals.hh"
#include "G4ios.hh"

//
// G4 Process: Gheisha High Energy Collision model.
// This includes the high energy cascading model, the two-body-resonance model
// and the low energy two-body model. Not included are the low energy stuff like
// nuclear reactions, nuclear fission without any cascading and all processes for
// particles at rest.  
// H. Fesefeldt, RWTH-Aachen, 23-October-1996
// Last modified: 29-July-1998 
// HPW, fixed bug in getting pdgencoding for nuclei
// Hisaya, fixed HighEnergyCascading
 
#include "G4HEInelastic.hh"
#include "G4HEVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"

void  G4HEInelastic::FillParticleChange(G4HEVector pv[], G4int aVecLength)
{
  theParticleChange.SetNumberOfSecondaries( aVecLength );
  for (G4int i=0; i<aVecLength; i++)
  {
    G4int pdgCode = pv[i].getCode();
    G4ParticleDefinition * aDefinition;
    if(pdgCode == 0)
    {
      G4int bNumber = pv[i].getBaryonNumber();
      if(bNumber==2) aDefinition = G4Deuteron::Deuteron();
      if(bNumber==3) aDefinition = G4Triton::Triton();
      if(bNumber==4) aDefinition = G4Alpha::Alpha();
    }
    else
    {
     aDefinition =  G4ParticleTable::GetParticleTable()->FindParticle(pdgCode);
    }
    G4DynamicParticle * aParticle = new G4DynamicParticle();
    aParticle->SetDefinition(aDefinition);
    aParticle->SetMomentum(pv[i].getMomentum()*GeV);
    theParticleChange.AddSecondary(aParticle);
    G4ParticleDefinition * dummy = G4KaonZero::KaonZero();
    dummy = G4AntiKaonZero::AntiKaonZero();
  }
}

void
G4HEInelastic::SetParticles()
  {
    PionPlus.setDefinition("PionPlus");
    PionZero.setDefinition("PionZero");
    PionMinus.setDefinition("PionMinus");
    KaonPlus.setDefinition("KaonPlus");
    KaonZero.setDefinition("KaonZero");
    AntiKaonZero.setDefinition("AntiKaonZero");
    KaonMinus.setDefinition("KaonMinus"); 
    KaonZeroShort.setDefinition("KaonZeroShort");
    KaonZeroLong.setDefinition("KaonZeroLong"); 
    Proton.setDefinition("Proton");
    AntiProton.setDefinition("AntiProton");
    Neutron.setDefinition("Neutron");
    AntiNeutron.setDefinition("AntiNeutron");
    Lambda.setDefinition("Lambda");
    AntiLambda.setDefinition("AntiLambda");
    SigmaPlus.setDefinition("SigmaPlus"); 
    SigmaZero.setDefinition("SigmaZero");
    SigmaMinus.setDefinition("SigmaMinus");
    AntiSigmaPlus.setDefinition("AntiSigmaPlus");
    AntiSigmaZero.setDefinition("AntiSigmaZero");
    AntiSigmaMinus.setDefinition("AntiSigmaMinus");
    XiZero.setDefinition("XiZero");
    XiMinus.setDefinition("XiMinus"); 
    AntiXiZero.setDefinition("AntiXiZero");
    AntiXiMinus.setDefinition("AntiXiMinus");
    OmegaMinus.setDefinition("OmegaMinus");
    AntiOmegaMinus.setDefinition("AntiOmegaMinus");
    Deuteron.setDefinition("Deuteron"); 
    Triton.setDefinition("Triton"); 
    Alpha.setDefinition("Alpha");
    Gamma.setDefinition("Gamma");
    return;
  }
 
G4double
G4HEInelastic::NuclearInelasticity(G4double incidentKineticEnergy,
                                   G4double atomicWeight,
                                   G4double atomicNumber)
  {
    G4double expu   = log(MAXFLOAT);
    G4double expl   = -expu;
    G4double ala    = log(atomicWeight);
    G4double ale    = log(incidentKineticEnergy);
    G4double sig1   = 0.5;
    G4double sig2   = 0.5;
    G4double em     = min(0.239 + 0.0408*ala*ala, 1.);
    G4double cinem  = min(0.0019*pow(ala,3.), 0.15);
    G4double sig    = (ale > em) ? sig2 : sig1; 
    G4double corr   = min(max(-pow(ale-em,2.)/(2.*sig*sig),expl), expu);
    G4double dum1   = -(incidentKineticEnergy)*cinem;
    G4double dum2   = abs(dum1);
    G4double dum3   = exp(corr);
    G4double cinema = 0.;
    if (dum2 >= 1.)           cinema = dum1*dum3;
    else if (dum3 > 1.e-10) cinema = dum1*dum3;    
    cinema = - max(-incidentKineticEnergy, cinema);
    if(verboseLevel > 1)
      { G4cout << " NuclearInelasticity: " << ala << " " <<  ale << " "  << em << " " << corr
             << " " <<  dum1 << " "  << dum2 << " " <<  dum3 << " " <<  cinema << endl;
      }                 
    return cinema;
  }

G4double
G4HEInelastic::NuclearExcitation(G4double  incidentKineticEnergy,
                                 G4double  atomicWeight,
                                 G4double  atomicNumber,
                                 G4double& excitationEnergyGPN,
                                 G4double& excitationEnergyDTA)
  {
    G4double neutronMass  = Neutron.getMass();
    G4double electronMass = 0.000511;
    G4double exnu         = 0.; 
    excitationEnergyGPN   = 0.;
    excitationEnergyDTA   = 0.;

    if (atomicWeight > (neutronMass + 2.*electronMass))
       {
         G4int    magic = ((G4int)(atomicNumber+0.1) == 82) ? 1 : 0;
         G4double ekin  = min(max(incidentKineticEnergy, 0.1), 4.);
         G4double cfa   = max(0.35 +((0.35 - 0.05)/2.3)*log(ekin), 0.15);
                  exnu  = 7.716*cfa*exp(-cfa);
         G4double atno  = min(atomicWeight, 120.);
                  cfa   = ((atno - 1.)/120.) * exp(-(atno-1.)/120.);
                  exnu  = exnu * cfa;
         G4double fpdiv = max(1.-0.25*ekin*ekin, 0.5);
         G4double gfa   = 2.*((atomicWeight-1.)/70.) 
                            * exp(-(atomicWeight-1.)/70.);

         excitationEnergyGPN = exnu * fpdiv;
         excitationEnergyDTA = exnu - excitationEnergyGPN;  
        
         G4double ran1 = 0., ran2 = 0.;
	 if (!magic)
	    { ran1 = normal();
              ran2 = normal();
            }
         excitationEnergyGPN = max(excitationEnergyGPN*(1.+ran1*gfa),0.);
         excitationEnergyDTA = max(excitationEnergyDTA*(1.+ran2*gfa),0.);
         exnu = excitationEnergyGPN + excitationEnergyDTA;
         if(verboseLevel > 1)
           {
             G4cout << " NuclearExcitation: " << magic << " " <<  ekin << " " << cfa
                  << " " <<  atno << " " << fpdiv << " " <<  gfa << " " << excitationEnergyGPN
                  << " " <<  excitationEnergyDTA << endl;
           } 

         while (exnu >= incidentKineticEnergy)
	    {
              excitationEnergyGPN *= (1. - 0.5*normal());
              excitationEnergyDTA *= (1. - 0.5*normal());
              exnu = excitationEnergyGPN + excitationEnergyDTA;
            }
       }             
    return exnu;
  }     
     
G4double 
G4HEInelastic::pmltpc(G4int np, G4int nm, G4int nz, G4int n, G4double b, G4double c)
 { 
   G4double expxu = log(MAXFLOAT);
   G4double expxl = -expxu;
   G4int i;
   G4double npf = 0.0, nmf = 0.0, nzf = 0.0;
   for(i=2;i<=np;i++) npf += log((G4double)i);
   for(i=2;i<=nm;i++) nmf += log((G4double)i);
   for(i=2;i<=nz;i++) nzf += log((G4double)i);
   G4double r = min(expxu,max(expxl,-(np-nm+nz+b)*(np-nm+nz+b)/(2*c*c*n*n)-npf-nmf-nzf));
   return exp(r);
 }

G4int G4HEInelastic::Factorial(G4int n)
 { 
   G4int result = 1;
   if(n<0) exit(EXIT_FAILURE);
   while(n>1) result *= n--;
   return result;
 } 

G4double G4HEInelastic::normal()
 {
   G4double ran = -6.0;
   for(G4int i=0; i<12; i++) ran += G4UniformRand();
   return ran;
 }

G4int G4HEInelastic::Poisson( G4double x )
 {
   G4int i, iran = 0;
   G4double ran;
   if ( x > 9.9 )
      {
        iran = G4int( max( 0.0, x + normal() * sqrt( x ) ) );
      }
   else
      {
        G4int mm = G4int( 5.0 * x );
        if ( mm <= 0 )
	   {
             G4double p1 = x * exp( -x );
             G4double p2 = x * p1/2.;
             G4double p3 = x * p2/3.;
             ran = G4UniformRand();
             if      ( ran < p3 ) iran = 3;
             else if ( ran < p2 ) iran = 2;
             else if ( ran < p1 ) iran = 1;
           }
        else
	   { G4double r = exp( -x );
             ran = G4UniformRand();
             if (ran > r)
	        {
                  G4double rrr;
                  G4double rr = r;
                  for (i=1; i <= mm; i++)
		      {
                        iran++;
                        if ( i > 5 ) rrr = exp(i*log(x)-(i+0.5)*log((G4double)i)+i-0.9189385);
                        else     rrr = pow(x,i)*Factorial(i);
                        rr += r * rrr;
                        if (ran <= rr) break;
		      }
		}         
	   }
      }
   return iran;
 }
G4double 
G4HEInelastic::GammaRand( G4double avalue )
 {
   G4double ga = avalue -1.;
   G4double la = sqrt(2.*avalue - 1.);
   G4double ep = 1.570796327 + atan(ga/la);
   G4double ro = 1.570796327 - ep;
   G4double y  = 1.;
   G4double xtrial;
   repeat:
   xtrial = ga + la * tan(ep*G4UniformRand() + ro);
   if(xtrial == 0.) goto repeat;
   y = log(1.+sqr((xtrial-ga)/la))+ga*log(xtrial/ga)-xtrial+ga;
   if(log(G4UniformRand()) > y) goto repeat;  
   return xtrial;
 }
G4double 
G4HEInelastic::Erlang( G4int mvalue )
 {
   G4double ran = G4UniformRand();
   G4double xtrial = 0.62666*log((1.+ran)/(1.-ran));
   if(G4UniformRand()<0.5) xtrial = -xtrial;
   return mvalue+xtrial*sqrt(G4double(mvalue));
 }  

void
G4HEInelastic::StrangeParticlePairProduction(
                       const G4double availableEnergy,
                       const G4double centerOfMassEnergy,
                       G4HEVector pv[],
                       G4int &vecLen,
                       G4HEVector incidentParticle,
                       G4HEVector targetParticle )

   // Choose charge combinations K+ K-, K+ K0, K0 K0, K0 K-,
   //                            K+ Y0, K0 Y+,  K0 Y-
   // For antibaryon induced reactions half of the cross sections KB YB
   // pairs are produced.  Charge is not conserved, no experimental data available
   // for exclusive reactions, therefore some average behaviour assumed.
   // The ratio L/SIGMA is taken as 3:1 (from experimental low energy data)

 {
   static G4double avrs[]  = {3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.};
   static G4double avkkb[] = {0.0015,0.0050,0.0120,0.0285,0.0525,0.0750,0.0975,
                              0.1230,0.2800,0.3980,0.4950,0.5730};
   static G4double kkb[]   = {0.2500,0.3750,0.5000,0.5625,0.6250,0.6875,0.7500,
                              0.8750,1.0000};
   static G4double ky[]    = {0.2000,0.3000,0.4000,0.5500,0.6250,0.7000,0.8000,
                              0.8500,0.9000,0.9500,0.9750,1.0000};
   static G4int ipakkb[]   = {10,13,10,11,10,12,11,11,11,12,12,11,12,12,
                              11,13,12,13};
   static G4int ipaky[]    = {18,10,18,11,18,12,20,10,20,11,20,12,21,10,
                              21,11,21,12,22,10,22,11,22,12};
   static G4int ipakyb[]   = {19,13,19,12,19,11,23,13,23,12,23,11,24,13,
                              24,12,24,11,25,13,25,12,25,11};
   static G4double avky[]  = {0.0050,0.0300,0.0640,0.0950,0.1150,0.1300,0.1450,
                              0.1550,0.2000,0.2050,0.2100,0.2120};
   static G4double avnnb[] ={0.00001,0.0001,0.0006,0.0025,0.0100,0.0200,0.0400,
                              0.0500,0.1200,0.1500,0.1800,0.2000};

   G4int i, ibin, i3, i4;       // misc. local variables
   G4double avk, avy, avn, ran;

   G4double protonMass = Proton.getMass();
   G4double sigmaMinusMass = SigmaMinus.getMass();
   G4int antiprotonCode = AntiProton.getCode();
   G4int antineutronCode = AntiNeutron.getCode();
   G4int antilambdaCode = AntiLambda.getCode();    

   G4double incidentMass = incidentParticle.getMass();
   G4double incidentEnergy = incidentParticle.getEnergy();
   G4int incidentCode = incidentParticle.getCode();

   G4double targetMass = targetParticle.getMass();

                     // protection against annihilation processes like pbar p -> pi pi.

   if (vecLen <= 2) return;   

                     // determine the center of mass energy bin

   i = 1;
   while ( (i<12) && (centerOfMassEnergy > avrs[i]) )i++;
   if    ( i == 12 ) ibin = 11;
   else              ibin = i;

                     // the fortran code chooses a random replacement of produced kaons
                     //  but does not take into account charge conservation 

   if( vecLen == 3 ) 
     {                                        // we know that vecLen > 2
        i3 = 2;
        i4 = 3;                               // note that we will be adding a new 
     }                                        // secondary particle in this case only
   else 
     {                                        // otherwise  2 <= i3,i4 <= vecLen
        i4 = i3 = 2 + G4int( (vecLen-2)*G4UniformRand() );
        while ( i3 == i4 ) i4 = 2 + G4int( (vecLen-2)*G4UniformRand() );
     }

                // use linear interpolation or extrapolation by y=centerofmassEnergy*x+b

   avk = (log(avkkb[ibin])-log(avkkb[ibin-1]))*(centerOfMassEnergy-avrs[ibin-1])
          /(avrs[ibin]-avrs[ibin-1]) + log(avkkb[ibin-1]);
   avk = exp(avk);

   avy = (log(avky[ibin])-log(avky[ibin-1]))*(centerOfMassEnergy-avrs[ibin-1])
          /(avrs[ibin]-avrs[ibin-1]) + log(avky[ibin-1]);
   avy = exp(avy);

   avn = (log(avnnb[ibin])-log(avnnb[ibin-1]))*(centerOfMassEnergy-avrs[ibin-1])
          /(avrs[ibin]-avrs[ibin-1]) + log(avnnb[ibin-1]);
   avn = exp(avn);

   if ( avk+avy+avn <= 0.0 ) return;

   if ( incidentMass < protonMass ) avy /= 2.0;
   avy += avk+avn;
   avk += avn;

   ran = G4UniformRand();
   if (  ran < avn ) 
      {                                                // p pbar && n nbar production
         if ( availableEnergy < 2.0) return;
         if ( vecLen == 3 )                 
            {                                          // add a new secondary
              if ( G4UniformRand() < 0.5 ) 
                {
                  pv[i3] = Neutron;;
                  pv[vecLen++] = AntiNeutron;
                } 
              else 
                {
                  pv[i3] = Proton;
                  pv[vecLen++] = AntiProton;
                }
            } 
         else 
            {                                          // replace two secondaries
              if ( G4UniformRand() < 0.5 ) 
                {
                  pv[i3] = Neutron;
                  pv[i4] = AntiNeutron;
                } 
              else 
                {
                  pv[i3] = Proton;
                  pv[i4] = AntiProton;
                }
            }
      } 
   else if ( ran < avk ) 
      {                                                // K Kbar production
        if ( availableEnergy < 1.0) return;
        G4double ran1 = G4UniformRand();
        i = 0;
        while( (i<9) && (ran1>kkb[i]) )i++;
        if ( i == 9 ) return;

               // ipakkb[] = { 10,13, 10,11, 10,12, 11, 11,  11,12, 12,11, 12,12, 11,13, 12,13 };
               // charge       K+ K-  K+ K0S K+ K0L K0S K0S K0S K0L K0LK0S K0LK0L K0S K- K0LK-

        switch( ipakkb[i*2] ) 
          {
            case 10: pv[i3] = KaonPlus;     break;
            case 11: pv[i3] = KaonZeroShort;break;
            case 12: pv[i3] = KaonZeroLong; break;
            case 13: pv[i3] = KaonMinus;    break;
          }

        if( vecLen == 2 ) 
          {                                                // add a secondary
            switch( ipakkb[i*2+1] ) 
              {
                case 10: pv[vecLen++] = KaonPlus;     break;
                case 11: pv[vecLen++] = KaonZeroShort;break;
                case 12: pv[vecLen++] = KaonZeroLong; break;
                case 13: pv[vecLen++] = KaonMinus;    break;
              }
          } 
        else 
          {                                        // replace
            switch( ipakkb[i*2+1] ) 
              {
                case 10: pv[i4] = KaonPlus;     break;
                case 11: pv[i4] = KaonZeroShort;break;
                case 12: pv[i4] = KaonZeroLong; break;
                case 13: pv[i4] = KaonMinus;    break;
              }
          }
      } 
   else if ( ran < avy ) 
      {                                            // Lambda K && Sigma K
        if( availableEnergy < 1.6) return;
        G4double ran1 = G4UniformRand();
        i = 0;
        while( (i<12) && (ran1>ky[i]) )i++;
        if ( i == 12 ) return;
        if ( (incidentMass<protonMass) || (G4UniformRand()<0.5) ) 
           {

                          // ipaky[] = { 18,10, 18,11, 18,12, 20,10, 20,11, 20,12,
                          //             L0 K+  L0 K0S L0 K0L S+ K+  S+ K0S S+ K0L 
                          //
                          //             21,10, 21,11, 21,12, 22,10, 22,11, 22,12 }
                          //             S0 K+  S0 K0S S0 K0L S- K+  S- K0S S- K0L

              switch( ipaky[i*2] ) 
                 {
                   case 18: pv[1] = Lambda;    break;
                   case 20: pv[1] = SigmaPlus; break;
                   case 21: pv[1] = SigmaZero; break;
                   case 22: pv[1] = SigmaMinus;break;
                 }
              switch( ipaky[i*2+1] ) 
                 {
                   case 10: pv[i3] = KaonPlus;     break;
                   case 11: pv[i3] = KaonZeroShort;break;
                   case 12: pv[i3] = KaonZeroLong; break;
                 }
           } 
        else 
           {                                          // Lbar K && Sigmabar K production 

                          // ipakyb[] = { 19,13, 19,12, 19,11,  23,13,  23,12,  23,11,
                          //              Lb K-  Lb K0L Lb K0S  S+b K- S+b K0L S+b K0S 
                          //              24,13, 24,12, 24,11, 25,13, 25,12, 25,11 };
                          //             S0b K-  S0BK0L S0BK0S S-BK- S-B K0L S-BK0S    
      
              if( (incidentCode==antiprotonCode) || (incidentCode==antineutronCode) ||
                  (incidentCode==antilambdaCode) || (incidentMass>sigmaMinusMass) ) 
                {
                  switch( ipakyb[i*2] ) 
                   {
                    case 19:pv[0] = AntiLambda;    break;
                    case 23:pv[0] = AntiSigmaPlus; break;
                    case 24:pv[0] = AntiSigmaZero; break;
                    case 25:pv[0] = AntiSigmaMinus;break;
                   }
                  switch( ipakyb[i*2+1] ) 
                   {
                    case 11:pv[i3] = KaonZeroShort;break;
                    case 12:pv[i3] = KaonZeroLong; break;
                    case 13:pv[i3] = KaonMinus;    break;
                   }
                } 
              else 
                {
                  switch( ipaky[i*2] ) 
                   {
                    case 18:pv[0] = Lambda;    break;
                    case 20:pv[0] = SigmaPlus; break;
                    case 21:pv[0] = SigmaZero; break;
                    case 22:pv[0] = SigmaMinus;break;
                   }
                  switch( ipaky[i*2+1] ) 
                   {
                    case 10: pv[i3] = KaonPlus;     break;
                    case 11: pv[i3] = KaonZeroShort;break;
                    case 12: pv[i3] = KaonZeroLong; break;
                   }
                }
           }
      } 
   else 
      return;

         //  check the available energy
         //   if there is not enough energy for kkb/ky pair production
         //   then reduce the number of secondary particles 
         //  NOTE:
         //        the number of secondaries may have been changed
         //        the incident and/or target particles may have changed
         //        charge conservation is ignored (as well as strangness conservation)

   incidentMass = incidentParticle.getMass();
   targetMass   = targetParticle.getMass();

   G4double energyCheck = centerOfMassEnergy-(incidentMass+targetMass);
   if (verboseLevel > 1) G4cout << "Particles produced: " ;
  
   for ( i=0; i < vecLen; i++ ) 
       {
         energyCheck -= pv[i].getMass(); 
         if (verboseLevel > 1) G4cout << pv[i].getCode() << " " ;
         if( energyCheck < 0.0 ) 
           {
             if( i > 0 ) vecLen = --i;      // chop off the secondary list
             return;
           }
       }
   if (verboseLevel > 1) G4cout << endl;
   return;
 }

void
G4HEInelastic::HighEnergyCascading(G4bool &successful,
                                   G4HEVector pv[],
                                   G4int &vecLen,	
                                   G4double &excitationEnergyGNP,
                                   G4double &excitationEnergyDTA,
                                   G4HEVector incidentParticle,
                                   G4HEVector targetParticle,
                                   G4double atomicWeight,
                                   G4double atomicNumber)
 {   
//
//  The multiplicity of particles produced in the first interaction has been
//  calculated in one of the FirstIntInNuc.... routines. The nuclear
//  cascading particles are parametrized from experimental data.
//  A simple single variable description E D3S/DP3= F(Q) with
//  Q^2 = (M*X)^2 + PT^2 is used. Final state kinematic is produced
//  by an FF-type iterative cascade method.
//  Nuclear evaporation particles are added at the end of the routine.

//  All quantities on the G4HEVector Array pv are in GeV- units.
//  The method is a copy of MediumEnergyCascading with some special tuning
//  for high energy interactions.


   G4int protonCode       = Proton.getCode();
   G4double protonMass    = Proton.getMass();
   G4int neutronCode      = Neutron.getCode();
   G4double neutronMass   = Neutron.getMass();
   G4double kaonPlusMass  = KaonPlus.getMass();
   G4int kaonPlusCode     = KaonPlus.getCode();   
   G4int kaonMinusCode    = KaonMinus.getCode();
   G4int kaonZeroSCode    = KaonZeroShort.getCode(); 
   G4int kaonZeroLCode    = KaonZeroLong.getCode();
   G4int kaonZeroCode     = KaonZero.getCode();
   G4int antiKaonZeroCode = AntiKaonZero.getCode(); 
   G4int pionPlusCode     = PionPlus.getCode();    
   G4int pionZeroCode     = PionZero.getCode();    
   G4int pionMinusCode    = PionMinus.getCode(); 
   G4String mesonType     = "meson";
   G4String baryonType    = "baryon"; 
   
   G4int    targetCode   = targetParticle.getCode();
   G4double targetMass   = targetParticle.getMass();
   G4double targetCharge = targetParticle.getCharge();

   G4int    incidentCode          = incidentParticle.getCode();
   G4double incidentMass          = incidentParticle.getMass();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   G4double incidentCharge        = incidentParticle.getCharge();
   G4double incidentEnergy        = incidentParticle.getEnergy();
   G4double incidentKineticEnergy = incidentParticle.getKineticEnergy();
   G4String incidentType          = incidentParticle.getType();
//   G4double incidentTOF           = incidentParticle.getTOF();   
   G4double incidentTOF           = 0.;
   
                                                              // some local variables

   G4int i, j, k, l;

   if (verboseLevel > 1) G4cout << " G4HEInelastic::HighEnergyCascading " << endl;
   successful = false; 
   if(incidentTotalMomentum < 25. + G4UniformRand()*25.) return;
 
                              // define annihilation channels.
                                 
   G4bool annihilation = false;
   if (incidentCode < 0 && incidentType == baryonType && 
       pv[0].getType() != baryonType &&
       pv[1].getType() != baryonType   )
         { 
           annihilation = true;
         }   
     


   G4double twsup[] = { 1., 1., 0.7, 0.5, 0.3, 0.2, 0.1, 0.0 };

   if( annihilation ) goto start;
   if( vecLen >= 8)   goto start;
   if( incidentKineticEnergy < 1.) return; 
   if(   (   incidentCode == kaonPlusCode  || incidentCode == kaonMinusCode
          || incidentCode == kaonZeroCode  || incidentCode == antiKaonZeroCode
          || incidentCode == kaonZeroSCode || incidentCode == kaonZeroLCode )
      && (   G4UniformRand() < 0.5) ) goto start;
   if( G4UniformRand() > twsup[vecLen-1]) goto start;
   if( incidentKineticEnergy > (G4UniformRand()*200 + 50.) ) goto start;
   return;
   
   start:

   if (annihilation)
      {                        // do some corrections of incident particle kinematic
        G4double ekcor = max( 1., 1./incidentKineticEnergy);
        incidentKineticEnergy = 2*targetMass + incidentKineticEnergy*(1.+ekcor/atomicWeight);
        G4double excitation = NuclearExcitation(incidentKineticEnergy,
                                                atomicWeight,
                                                atomicNumber,
                                                excitationEnergyGNP,
                                                excitationEnergyDTA);
        incidentKineticEnergy -= excitation;
        if (incidentKineticEnergy < excitationEnergyDTA) incidentKineticEnergy = 0.;
        incidentEnergy = incidentKineticEnergy + incidentMass;
        incidentTotalMomentum =
                 sqrt( max(0., incidentEnergy*incidentEnergy - incidentMass*incidentMass));
      }  
    
   G4HEVector pTemp;
   for (i = 2; i<vecLen; i++)
     {
       j = min( vecLen-1, (G4int)(2. + G4UniformRand()*(vecLen - 2)));
       pTemp = pv[j];
       pv[j] = pv[i];
       pv[i] = pTemp;
     }                         
                              // randomize the first two leading particles
                              // for kaon induced reactions only 
                              // (need from experimental data)

   if(   (incidentCode==kaonPlusCode  || incidentCode==kaonMinusCode    ||
          incidentCode==kaonZeroCode  || incidentCode==antiKaonZeroCode ||
          incidentCode==kaonZeroSCode || incidentCode==kaonZeroLCode)   
      && (G4UniformRand() > 0.9) ) 
     {
         pTemp = pv[1];
         pv[1] = pv[0];
         pv[0] = pTemp;
     }
                              // mark leading particles for incident strange particles 
                              // and antibaryons, for all other we assume that the first 
                              // and second particle are the leading particles. 
                              // We need this later for kinematic aspects of strangeness
                              // conservation.
                          
   G4int lead = 0;                   
   G4HEVector leadParticle;
   if( (incidentMass >= kaonPlusMass-0.05) && (incidentCode != protonCode)          
                                           && (incidentCode != neutronCode) ) 
         {       
           G4double pMass = pv[0].getMass();
           G4int    pCode = pv[0].getCode();
           if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                            && (pCode != neutronCode) ) 
                  {       
                    lead = pCode; 
                    leadParticle = pv[0];                           
                  } 
           else   
                  {
                    pMass = pv[1].getMass();
                    pCode = pv[1].getCode();
                    if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                                     && (pCode != neutronCode) ) 
                        {       
                          lead = pCode;
                          leadParticle = pv[1];
                        }
                  }
         }

                              //  distribute particles in forward and backward 
                              //  hemispheres in center of mass system
                              //  incident goes in forward hemisphere
   
   G4HEVector pvI = incidentParticle;  // for the incident particle
   pvI.setSide( 1 );

   G4HEVector pvT = targetParticle;   // for the target particle
   pvT.setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvT.setSide( -1 );
   pvT.setTOF( -1.); 


   G4double centerOfMassEnergy = sqrt( sqr(incidentMass)+sqr(targetMass)
                                      +2.0*targetMass*incidentEnergy );
   G4double availableEnergy    = centerOfMassEnergy - ( targetMass + incidentMass );

   G4double tavai1      = centerOfMassEnergy/2.0 - incidentMass;
   G4double tavai2      = centerOfMassEnergy/2.0 - targetMass;           
   
                              // define G4HEVector- array for kinematic manipulations,
                              // with a one by one correspondence to the pv-Array. 
   
   G4int ntb = 1;
   for( i=0; i < vecLen; i++ ) 
      {
        if      (i == 0) pv[i].setSide(  1 );
        else if (i == 1) pv[i].setSide( -1 );
        else  
           { if( G4UniformRand() < 0.5 ) 
               {
                 pv[i].setSide( -1 );
                 ntb++;
                } 
             else 
                 pv[i].setSide( 1 );
           }
        pv[i].setTOF(    incidentTOF);
      }
   G4double tb = 2. * ntb;
   if (centerOfMassEnergy < (2. + G4UniformRand())) 
       tb = (2. * ntb + vecLen)/2.;     

   if (verboseLevel > 1)
      { G4cout << " pv Vector after Randomization " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      } 

                              //  add particles from intranuclear cascade
                              //  nuclearCascadeCount = number of new secondaries 
                              //  produced by nuclear cascading.
                              //  extraCount = number of nucleons within these new secondaries
   
   G4double s, xtarg, ran;
   s = centerOfMassEnergy*centerOfMassEnergy;
   G4double afc;
   afc = min( 0.5, 0.312 + 0.200 * log(log(s))+ pow(s,1.5)/6000.0 ); 
   xtarg = max(0.01, afc * (pow(atomicWeight, 0.33) - 1.0) * tb);

   G4int nstran = Poisson( 0.03*xtarg);
     
   G4int momentumBin = 0;
   G4double nucsup[] = { 1.00, 0.7, 0.5, 0.4, 0.5, 0.5 };
   G4double psup[]   = {   3.,  6., 20., 50., 100., 1000. };
   while( (momentumBin < 6) && (incidentTotalMomentum > psup[momentumBin])) momentumBin++;
   momentumBin = min(5, momentumBin);
   G4double xpnhmf = max(0.01,xtarg*nucsup[momentumBin]);
   G4double xshhmf = max(0.01,xtarg - xpnhmf);
   G4double rshhmf = 0.25*xshhmf;
   G4double rpnhmf = 0.81*xpnhmf;
   G4double xhmf;
   if(verboseLevel > 1)
     { G4cout << "xtarg= " << xtarg << " xpnhmf = " << xpnhmf << endl;
     }
   G4int nshhmf, npnhmf;
   if (rshhmf > 1.1)
     {
        rshhmf = xshhmf/(rshhmf-1.);
        if (rshhmf <= 20.) 
            xhmf = GammaRand( rshhmf );
        else
            xhmf = Erlang( G4int(rshhmf+0.5) );
        xshhmf *= xhmf/rshhmf;
     }
   nshhmf = Poisson( xshhmf );   
   if(verboseLevel > 1)
     { G4cout << "xshhmf = " << xshhmf << " xhmf = " << xhmf << " rshhmf = " << rshhmf << endl;
     }

   if (rpnhmf > 1.1)
     {
        rpnhmf = xpnhmf/(rpnhmf -1.);
        if (rpnhmf <= 20.)
            xhmf = GammaRand( rpnhmf );
        else
            xhmf = Erlang( G4int(rpnhmf+0.5) );
        xpnhmf *= xhmf/rpnhmf;
     }
   npnhmf = Poisson( xpnhmf );
   if(verboseLevel > 1)
     { G4cout << "nshhmf = " << nshhmf << " npnhmf = " <<  npnhmf << " nstran = " << nstran << endl;
     } 
   G4int ntarg = nshhmf + npnhmf + nstran;          

   G4int targ = 0;
   
   while (npnhmf > 0)  
     {
       if ( G4UniformRand() > (1. - atomicNumber/atomicWeight))
          pv[vecLen] = Proton;
       else
          pv[vecLen] = Neutron;
       targ++;
       pv[vecLen].setSide( -2 );
       pv[vecLen].setFlag( true );
       pv[vecLen].setTOF( incidentTOF );
       vecLen++;
       npnhmf--;
     }
   while (nstran > 0)
     {
       ran = G4UniformRand();
       if (ran < 0.14)      pv[vecLen] = Lambda;
       else if (ran < 0.20) pv[vecLen] = SigmaZero;
       else if (ran < 0.43) pv[vecLen] = KaonPlus;
       else if (ran < 0.66) pv[vecLen] = KaonZero;
       else if (ran < 0.89) pv[vecLen] = AntiKaonZero;
       else                 pv[vecLen] = KaonMinus;
       if (G4UniformRand() > 0.2)
         { 
           pv[vecLen].setSide( -2 );
           pv[vecLen].setFlag( true );
         } 
       else
         {
           pv[vecLen].setSide( 1 );
           pv[vecLen].setFlag( false );
           ntarg--;
         } 
       pv[vecLen].setTOF( incidentTOF );
       vecLen++;         
       nstran--;   
     } 
   while (nshhmf > 0)
     {
       ran = G4UniformRand();
       if( ran < 0.33333 ) 
           pv[vecLen] = PionPlus;
       else if( ran < 0.66667 ) 
           pv[vecLen] = PionZero;
       else 
           pv[vecLen] = PionMinus;
       if (G4UniformRand() > 0.2)
          {
            pv[vecLen].setSide( -2 );                             // backward cascade particles
            pv[vecLen].setFlag( true );                           // true is the same as IPA(i)<0
          }
       else
          {
            pv[vecLen].setSide( 1 );
            pv[vecLen].setFlag( false );
            ntarg--;
          }  
       pv[vecLen].setTOF( incidentTOF );
       vecLen++; 
       nshhmf--;
     }
                                         
                                         //  assume conservation of kinetic energy 
                                         //  in forward & backward hemispheres
   G4int is, iskip, iavai1;
   if(vecLen <= 1) return;

   tavai1 = centerOfMassEnergy/2.;
   iavai1 = 0;
 
   for (i = 0; i < vecLen; i++) 
       { 
         if (pv[i].getSide() > 0)
            { 
               tavai1 -= pv[i].getMass();
               iavai1++;
            }    
       } 
   if ( iavai1 == 0) return;

   while( tavai1 <= 0.0 ) 
        {                                // must eliminate a particle from the forward side
           iskip = G4int(G4UniformRand()*iavai1) + 1; 
           is = 0;  
           for( i=vecLen-1; i>=0; i-- ) 
              {
                if( pv[i].getSide() > 0 ) 
                  {
                    if (++is == iskip) 
                        {
                           tavai1 += pv[i].getMass();
                           iavai1--;            
                           if ( i != vecLen-1)
                              { 
                                 for( j=i; j<vecLen; j++ ) 
                                    {         
                                       pv[j] = pv[j+1];
                                    }
                              }
                           if( --vecLen == 0 ) return;     // all the secondaries except of the 
                           break;            // --+
                        }                    //   |
                  }                          //   v
              }                              // break goes down to here
        }                                    // to the end of the for- loop.
                                       

     tavai2 = (targ+1)*centerOfMassEnergy/2.;
     G4int iavai2 = 0;

     for (i = 0; i < vecLen; i++)
         {
           if (pv[i].getSide() < 0)
              {
                 tavai2 -= pv[i].getMass();
                 iavai2++;
              }
         }
     if (iavai2 == 0) return;

     while( tavai2 <= 0.0 ) 
        {                               // must eliminate a particle from the backward side
           iskip = G4int(G4UniformRand()*iavai2) + 1; 
           is = 0;
           for( i = vecLen-1; i >= 0; i-- ) 
              {
                if( pv[i].getSide() < 0 ) 
                  {
                    if( ++is == iskip ) 
                       {
                         tavai2 += pv[i].getMass();
                         iavai2--;
                         if (pv[i].getSide() == -2) ntarg--;
                         if (i != vecLen-1)
                            {
                              for( j=i; j<vecLen; j++)
                                 { 
                                   pv[j] = pv[j+1];
                                 } 
                            }
                         if (--vecLen == 0) return;
                         break;     
                       }
                  }   
              }
        }

   if (verboseLevel > 1)
      { G4cout << " pv Vector after Energy checks "
             << vecLen << " " << tavai1 << " " << iavai1 << " " << tavai2
             << " " <<  iavai2 << " " << ntarg << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      } 
   
                                 //  define some vectors for Lorentz transformations
   
   G4HEVector* pvmx = new G4HEVector [10];
   
   pvmx[0].setMass( incidentMass );
   pvmx[0].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   pvmx[1].setMass( protonMass);
   pvmx[1].setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvmx[3].setMass( protonMass*(1+targ));
   pvmx[3].setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvmx[4].setZero();
   pvmx[5].setZero();
   pvmx[7].setZero();
   pvmx[8].setZero();
   pvmx[8].setMomentum( 1.0, 0.0 );
   pvmx[2].Add( pvmx[0], pvmx[1] );
   pvmx[3].Add( pvmx[3], pvmx[0] );
   pvmx[0].Lor( pvmx[0], pvmx[2] );
   pvmx[1].Lor( pvmx[1], pvmx[2] );

   if (verboseLevel > 1)
     { G4cout << " General Vectors after Definition " << endl;
       for (i=0; i<10; i++) pvmx[i].Print(i);
     }
                                                 //  main loop for 4-momentum generation
                                                 //  see Pitha-report (Aachen) 
                                                 //  for a detailed description of the method
                                                 // process the secondary particles in reverse order
   G4double dndl[20];
   G4double binl[20];
   G4double pvMass, pvEnergy;
   G4int    pvCode; 
   G4double aspar, pt, phi, et, xval;
   G4double ekin  = 0.;
   G4double ekin1 = 0.;
   G4double ekin2 = 0.;
   G4int npg   = 0;
   G4double rmg0 = 0.;
   G4int targ1 = 0;                           // no fragmentation model for nucleons from the
   for( i=vecLen-1; i>=0; i-- )               // intranuclear cascade. Mark them with -3 and 
      {                                       // leave the loop
        if( i == 1)
          {
            if ( (pv[i].getMass() > neutronMass + 0.05) && (G4UniformRand() < 0.2))
               { 
                 if(++npg < 19)
                   {
                     pv[i].setSide(-3);
                     rmg0 += pv[i].getMass();
                     targ++;
                     continue;
                   }
               }
            else if ( pv[i].getMass() > protonMass - 0.05)
               {
                 if(++npg < 19)
                   {
                     pv[i].setSide(-3);
                     rmg0 += pv[i].getMass();
                     targ++;
                     continue;
                   }
               }  
          }  
        if( pv[i].getSide() == -2) 
          { 
            if ( pv[i].getName() == "Proton" || pv[i].getName() == "Neutron")
               {                                 
                 if( ++npg < 19 ) 
                   {
                     pv[i].setSide( -3 );
                     rmg0 += pv[i].getMass(); 
                     targ1++;
                     continue;                                    // leave the for loop !!
                   }     
               }
          }
                                              //  set pt and phi values, 
                                              //  they are changed somewhat in the iteration loop
                                              //  set mass parameter for lambda fragmentation model
        G4double maspar[] = { 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.20, 0.10};
        G4double     bp[] = { 4.00, 2.50, 2.20, 3.00, 3.00, 1.70, 3.50, 3.50};
        G4double   ptex[] = { 1.70, 1.70, 1.50, 1.70, 1.40, 1.20, 1.70, 1.20};    

                                              // set parameters for lambda simulation 
                                              // pt is the average transverse momentum
                                              // aspar is average transverse mass  
        pvMass = pv[i].getMass();       
        j = 2;                                              
        if (pv[i].getType() == mesonType ) j = 1;
        if ( pv[i].getMass() < 0.4 ) j = 0;
        if ( i <= 1 ) j += 3;
        if (pv[i].getSide() <= -2) j = 6;
        if (j == 6 && pv[i].getType() == baryonType) j = 7;
        pt    = sqrt(pow(-log(1.-G4UniformRand())/bp[j],ptex[j]));
        if(pt<0.05) pt = max(0.001, 0.3*G4UniformRand());
        aspar = maspar[j]; 
        phi   = G4UniformRand()*M_2PI; 
        pv[i].setMomentum( pt*cos(phi), pt*sin(phi) );        // set x- and y-momentum

        for( j=0; j<20; j++ ) binl[j] = j/(19.*pt);           // set the lambda - bins.
     
        if( pv[i].getSide() > 0 )
           et = pvmx[0].getEnergy();
        else
           et = pvmx[1].getEnergy();
     
        dndl[0] = 0.0;     
                                                              // start of outer iteration loop
        G4int outerCounter = 0, innerCounter = 0;             // three times.
        G4bool eliminateThisParticle = true;
        G4bool resetEnergies = true;
        while( ++outerCounter < 3 ) 
             {
               for( l=1; l<20; l++ ) 
                  {
                    xval  = (binl[l]+binl[l-1])/2.;              // x = lambda /GeV 
                    if( xval > 1./pt )
                       dndl[l] = dndl[l-1];
                    else
                       dndl[l] = dndl[l-1] + 
                                 aspar/sqrt( pow((1.+aspar*xval*aspar*xval),3) ) *
                                 (binl[l]-binl[l-1]) * et / 
                                 sqrt( pt*xval*et*pt*xval*et + pt*pt + pvMass*pvMass );
                  }  
       
                                                               // start of inner iteration loop
               innerCounter = 0;                               // try this not more than 7 times. 
               while( ++innerCounter < 7 ) 
                    {
                      l = 1;
                      ran = G4UniformRand()*dndl[19];
                      while( ( ran >= dndl[l] ) && ( l < 20 ) )l++;
                      l = min( 19, l );
                      xval = min( 1.0, pt*(binl[l-1] + G4UniformRand()*(binl[l]-binl[l-1])/2.) );
                      if( pv[i].getSide() < 0 ) xval *= -1.;
                      pv[i].setMomentumAndUpdate( xval*et );                  // set the z-momentum
                      pvEnergy = pv[i].getEnergy();
                      if( pv[i].getSide() > 0 )                               // forward side 
                        {
                          if ( i < 2 )
                             { 
                               ekin = tavai1 - ekin1;
                               if (ekin < 0.) ekin = 0.04*fabs(normal());
                               G4double pp1 = pv[i].Length();
                               if (pp1 >= 1.e-6)
                                  {
                                    G4double pp = sqrt(ekin*(ekin+2*pvMass));
                                    pp = max(0., pp*pp - pt*pt);
                                    pp = sqrt(pp)*pv[i].getSide()/fabs(G4double(pv[i].getSide())); // cast for aCC 
                                    pv[i].setMomentum( pp );
                                  }
                               else
                                  {
                                    pv[i].setMomentum(0.,0.,0.); 
                                    pv[i].setKineticEnergyAndUpdate( ekin);
                                  }
                               pvmx[4].Add( pvmx[4], pv[i]);
                               resetEnergies = false;
                               eliminateThisParticle = false; 
                               continue;
                             }      
                          else if( (ekin1+pvEnergy-pvMass) < 0.95*tavai1 ) 
                            {
                              pvmx[4].Add( pvmx[4], pv[i] );
                              ekin1 += pvEnergy - pvMass;
                              pvmx[6].Add( pvmx[4], pvmx[5] );
                              pvmx[6].setMomentum( 0.0 );
                              phi = pvmx[6].Ang( pvmx[8] );
                              if( pvmx[6].getMomentum().y() < 0.0 ) phi = M_2PI - phi;
                              phi += M_PI + normal()*M_PI/12.0;
                              if( phi > M_2PI )phi -= M_2PI;
                              if( phi < 0.0 )phi = M_2PI - phi;
                              outerCounter = 2;                     // leave outer loop
                              eliminateThisParticle = false;        // don't eliminate this particle
                              resetEnergies = false;
                              break;                                // next particle
                            }
                          if( innerCounter > 5 ) break;             // leave inner loop
                          
                          if( tavai2 >= pvMass ) 
                            {                                       // switch sides
                              pv[i].setSide( -1 );
                              tavai1 += pvMass;
                              tavai2 -= pvMass;
                              iavai2++;
                            }
                        } 
                      else 
                        {                                            // backward side
                          xval = min(0.999,0.95+0.05*targ/20.0);
                          if( (ekin2+pvEnergy-pvMass) < xval*tavai2 ) 
                            {
                              pvmx[5].Add( pvmx[5], pv[i] );
                              ekin2 += pvEnergy - pvMass;
                              pvmx[6].Add( pvmx[4], pvmx[5] );
                              pvmx[6].setMomentum( 0.0 );             // set z-momentum
                              phi = pvmx[6].Ang( pvmx[8] );
                              if( pvmx[6].getMomentum().y() < 0.0 )phi = M_2PI - phi;
                              phi += M_PI + normal() * M_PI / 12.0;
                              if( phi > M_2PI )phi -= M_2PI;
                              if( phi < 0.0 )phi = M_2PI - phi;
                              outerCounter = 2;                    // leave outer iteration
                              eliminateThisParticle = false;       // don't eliminate this particle
                              resetEnergies = false;
                              break;                               // leave inner iteration
                            }
                          if( innerCounter > 5 )break;             // leave inner iteration
                          
                          if( tavai1 >= pvMass ) 
                            {                                      // switch sides
                              pv[i].setSide( 1 );
                              tavai1  -= pvMass;
                              tavai2  += pvMass;
                              iavai2--;
                            }
                        }
                      pv[i].setMomentum( pv[i].getMomentum().x() * 0.9,
                                         pv[i].getMomentum().y() * 0.9);
                      pt *= 0.9;
                      dndl[19] *= 0.9;
                    }                                           // closes inner loop

               if (resetEnergies)
                    {
                      if (verboseLevel > 1) 
                         {
                           G4cout << " Reset energies for index " << i << " " << ekin1 
                                << " " << tavai1 << endl;
                           pv[i].Print(i);
                         }
                      ekin1 = 0.0;
                      ekin2 = 0.0;
                      pvmx[4].setZero();
                      pvmx[5].setZero();

                      for( l=i+1; l < vecLen; l++ ) 
                         {
                           if( (pv[l].getMass() < protonMass) || (pv[l].getSide() > 0) ) 
                             {
                                pvEnergy = pv[l].getMass() + 0.95*pv[l].getKineticEnergy(); 
                                pv[l].setEnergyAndUpdate( pvEnergy );
                                if( pv[l].getSide() > 0) 
                                  {
                                    ekin1 += pv[l].getKineticEnergy();
                                    pvmx[4].Add( pvmx[4], pv[l] );
                                  } 
                                else 
                                  {
                                    ekin2 += pv[l].getKineticEnergy();
                                    pvmx[5].Add( pvmx[5], pv[l] );
                                  }
                             }
                         }
                    }
             }                                                           // closes outer iteration

        if( eliminateThisParticle )                                      // not enough energy, 
          {                                                              // eliminate this particle
            if (verboseLevel > 1)
               {
                  G4cout << " Eliminate particle index " << i << endl;
                  pv[i].Print(i);
               }
            if(i != vecLen-1)
               {
                 for( j=i; j < vecLen-1; j++ ) 
                   {                                                         // shift down
                     pv[j] = pv[j+1];
                   }
               }
            vecLen--;
            if(vecLen < 2) return;
            i++;
            pvmx[6].Add( pvmx[4], pvmx[5] );
            pvmx[6].setMomentum( 0.0 );                                  // set z-momentum
            phi = pvmx[6].Ang( pvmx[8] );
            if( pvmx[6].getMomentum().y() < 0.0 )phi = M_2PI - phi;
            phi += M_PI + normal() * M_PI / 12.0;
            if( phi > M_2PI )phi -= M_2PI;
            if( phi < 0.0 )phi = M_2PI - phi;
          }
      }                                                                  // closes main for loop
   if (verboseLevel > 1)
      { G4cout << " pv Vector after lambda fragmentation " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
        for (i=0; i < 10; i++) pvmx[i].Print(i);
      } 

   
   
                                                                   // backward nucleons produced 
                                                                   // with a cluster model
   G4double gpar[] = {2.6, 2.6, 1.80, 1.30, 1.20};
   G4double cpar[] = {0.6, 0.6, 0.35, 0.15, 0.10};
 
   if (npg > 0)
     {
       G4double rmg = rmg0;
       if (npg > 1)
          {
            G4int npg1 = npg-1;
            if (npg1 >4) npg1 = 4;
            rmg += pow( -log(1.-G4UniformRand()), cpar[npg1])/gpar[npg1];
          }
       G4double ga = 1.2;
       G4double ekit1 = 0.04, ekit2 = 0.6;
       if(incidentKineticEnergy < 5.)
         {
           ekit1 *= sqr(incidentKineticEnergy)/25.;
           ekit2 *= sqr(incidentKineticEnergy)/25.;
         }
       G4double avalue = (1.-ga)/(pow(ekit2,1.-ga)-pow(ekit1,1.-ga));
       for (i = 0; i < vecLen; i++)
          {
            if (pv[i].getSide() == -3)
              {
                G4double ekit = pow(G4UniformRand()*(1.-ga)/avalue + pow(ekit1,1.-ga), 1./(1.-ga) );
                G4double cost = max(-1., min(1., log(2.23*G4UniformRand()+0.383)/0.96));
                G4double sint = sqrt(1. - cost*cost);
                G4double phi  = M_2PI*G4UniformRand();
                G4double pp   = sqrt(ekit*(ekit+2*pvMass));
                pv[i].setMomentum( pp*sint*sin(phi),
                                   pp*sint*cos(phi),
                                   pp*cost          );
                pv[i].Lor( pv[i], pvmx[2] );
                pvmx[5].Add( pvmx[5], pv[i] );
              }
          } 
     }        
   if  (vecLen <= 2)
     {
       successful = false;
       return;
     }  


                                                   //  Lorentz transformation in lab system
   targ = 0;
   for( i=0; i < vecLen; i++ ) 
      {
        if( pv[i].getType() == baryonType )targ++;
        if(verboseLevel > 1) pv[i].Print(i); 
        pv[i].Lor( pv[i], pvmx[1] );
        if(verboseLevel > 1) pv[i].Print(i);
      }
   targ = max( 1, targ );

   G4bool dum;
   if( lead ) 
     {
       for( i=0; i<vecLen; i++ ) 
          {
            if( pv[i].getCode() == lead ) 
              {
                dum = false;
                break;
              }
          }
       if( dum ) 
         {
           i = 0;          
 
           if(   (    (leadParticle.getType() == baryonType)
                   && (pv[1].getType() == baryonType))
              || (    (leadParticle.getType() == mesonType)
                   && (pv[1].getType() == mesonType)))
             {
               i = 1;
             } 
            ekin = pv[i].getKineticEnergy();
            pv[i] = leadParticle;  
            if( pv[i].getFlag() )
                pv[i].setTOF( -1.0 );
            else
                pv[i].setTOF( 1.0 );
            pv[i].setKineticEnergyAndUpdate( ekin );
         }
     }

   pvmx[3].setMass( incidentMass);
   pvmx[3].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   
   G4double ekin0 = pvmx[3].getKineticEnergy();
   
   pvmx[4].setMass( protonMass * targ);
   pvmx[4].setEnergy( protonMass * targ);
   pvmx[4].setKineticEnergy(0.);
   pvmx[4].setMomentum(0., 0., 0.);
   ekin = pvmx[3].getEnergy() + pvmx[4].getEnergy();

   pvmx[5].Add( pvmx[3], pvmx[4] );
   pvmx[3].Lor( pvmx[3], pvmx[5] );
   pvmx[4].Lor( pvmx[4], pvmx[5] );
   
   G4double tecm = pvmx[3].getEnergy() + pvmx[4].getEnergy();

   pvmx[7].setZero();
   
   ekin1 = 0.0;   
   G4double teta, wgt; 
   
   for( i=0; i < vecLen; i++ ) 
      {
        pvmx[7].Add( pvmx[7], pv[i] );
        ekin1 += pv[i].getKineticEnergy();
        ekin  -= pv[i].getMass();
      }
   
   if( vecLen > 1 && vecLen < 19 ) 
     {
       G4bool constantCrossSection = true;
       G4HEVector pw[19];
       for(i=0; i<vecLen; i++) pw[i] = pv[i]; 
       wgt = NBodyPhaseSpace( tecm, constantCrossSection, pw, vecLen );
       ekin = 0.0;
       for( i=0; i < vecLen; i++ ) 
          {
            pvmx[6].setMass( pw[i].getMass());
            pvmx[6].setMomentum( pw[i].getMomentum() );
            pvmx[6].SmulAndUpdate( pvmx[6], 1. );
            pvmx[6].Lor( pvmx[6], pvmx[4] );
            ekin += pvmx[6].getKineticEnergy();
          }
       teta = pvmx[7].Ang( pvmx[3] );
       if (verboseLevel > 1) 
         { 
           G4cout << " vecLen > 1 && vecLen < 19 " << teta << " " << ekin0 << " " << ekin1 << " " << ekin << endl;
         }
     }

   if( ekin1 != 0.0 ) 
     {
       pvmx[6].setZero();
       wgt = ekin/ekin1;
       ekin1 = 0.;
       for( i=0; i < vecLen; i++ ) 
          {
            pvMass = pv[i].getMass();
            ekin   = pv[i].getKineticEnergy() * wgt;
            pv[i].setKineticEnergyAndUpdate( ekin );
            ekin1 += ekin;
            pvmx[6].Add( pvmx[6], pv[i] );
          }
       teta = pvmx[6].Ang( pvmx[3] );
       if (verboseLevel > 1) 
	 {
            G4cout << " ekin1 != 0 " << teta << " " <<  ekin0 << " " <<  ekin1 << endl;
            incidentParticle.Print(0);
            targetParticle.Print(1);
            for(i=0;i<vecLen;i++) pv[i].Print(i);
         }
     }

                          //    Rotate in direction of z-axis, this does disturb in some way our
                          //    inclusive distributions, but it is necessary for momentum conservation

   pvmx[6].setZero();
   for( i=0; i < vecLen; i++ ) pvmx[6].Add( pvmx[6], pv[i] );
   
                          //  Some smearing in transverse direction from Fermi motion
   


                                                  // RANNOR
   G4double ry   = G4UniformRand();
   G4double rz   = G4UniformRand();
   G4double rx   = M_2PI*rz;
   G4double a1   = sqrt(-2.0*log(ry));
   G4double ran1 = a1*sin(rx);
   G4double ran2 = a1*cos(rx);
                                                  
   pvmx[6].setMomentum( pvmx[6].getMomentum().x()+ran1*0.020*targ,
                        pvmx[6].getMomentum().y()+ran2*0.020*targ );
   pvmx[7].Defs( pvmx[3], pvmx[6], pvmx[8], pvmx[9] );
   pvmx[6].setZero();
   for (i = 0; i < vecLen; i++) 
      { 
        pv[i].Trac( pv[i], pvmx[7] , pvmx[8], pvmx[9] );
        pvmx[6].Add( pvmx[6], pv[i] );  
      }
   teta = pvmx[6].Ang( pvmx[3] );   
   if (verboseLevel > 1) G4cout << " After smearing " << teta << endl;

                               //  Rotate in direction of primary particle, subtract binding energies
                               //   and make some further corrections if required

   G4double dekin = 0.0;
   G4int npions = 0;    
   G4double ek1 = 0.0;
   G4double alekw, xxh;
   G4double cfa = 0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.);
   G4double alem[] = {1.40, 2.30, 2.70, 3.00, 3.40, 4.60, 7.00, 10.00};
   G4double val0[] = {0.00, 0.40, 0.48, 0.51, 0.54, 0.60, 0.65,  0.70};
   

   for (i = 0; i < vecLen; i++)
      { 
        pv[i].Defs1( pv[i], pvI );
        if (atomicWeight > 1.5)
           {
             ekin  = max( 1.e-6,pv[i].getKineticEnergy() - cfa*( 1. + 0.5*normal()));
             alekw = log( incidentKineticEnergy );
             xxh   = 1.;
             if(incidentCode == pionPlusCode || incidentCode == pionMinusCode)
	       {
                 if(pv[i].getCode() == pionZeroCode)
                   {
                     if(G4UniformRand() < log(atomicWeight))
                       { 
                         if (alekw > alem[0])
                           {
                              for (j = 1; j < 8; j++)
                                 {
                                    if(alekw < alem[j]) break;
                                 } 
                              xxh =   (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alekw
                                     + val0[j-1] - (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alem[j-1];
                              xxh = 1. - xxh;
                           }
                       }      
                    }
               }  
             dekin += ekin*(1.-xxh);
             ekin *= xxh;
             pv[i].setKineticEnergyAndUpdate( ekin );
             pvCode = pv[i].getCode();
             if ((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode))
                {
                  npions += 1;
                  ek1 += ekin; 
                }
           }
      }
   if( (ek1 > 0.0) && (npions > 0) ) 
      {
        dekin = 1.+dekin/ek1;
        for (i = 0; i < vecLen; i++)
          {
            pvCode = pv[i].getCode();
            if((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode)) 
              {
                ekin = max( 1.0e-6, pv[i].getKineticEnergy() * dekin );
                pv[i].setKineticEnergyAndUpdate( ekin );
              }
          }
      }
   if (verboseLevel > 1)
      { G4cout << " Lab-System " <<  ek1 << " " << npions << endl;
        incidentParticle.Print(0);
        targetParticle.Print(1);
        for (i=0; i<vecLen; i++) pv[i].Print(i);
      }

                                       //  add black track particles
                                       //  the total number of particles produced is restricted to 198
                                       //  this may have influence on very high energies

   if (verboseLevel > 1) 
      G4cout << " Evaporation : " <<  atomicWeight << " " << excitationEnergyGNP << " " 
           <<  excitationEnergyDTA << endl;

   G4double sprob = 0.;
   if (incidentKineticEnergy > 5.)
       sprob = min(1., (0.25-0.02*log(atomicWeight))*log(incidentKineticEnergy-4.) );

   if( atomicWeight > 1.5 && G4UniformRand() > sprob ) 
     {

       G4double ekw, tex, cost, sint, pp, eka, etb;
       G4int spall, nbl;
                                     //  first add protons and neutrons

       if( excitationEnergyGNP >= 0.001 ) 
         {
                                     //  nbl = number of proton/neutron black track particles
                                     //  tex is their total kinetic energy (GeV)
       
           nbl = Poisson( (1.5+1.25*targ)*excitationEnergyGNP/
                                         (excitationEnergyGNP+excitationEnergyDTA));
           if( targ+nbl > atomicWeight ) nbl = atomicWeight - targ;
           if (verboseLevel > 1) 
              G4cout << " evaporation " << targ << " " << nbl << " " << sprob << endl; 
           spall = targ;
           if( nbl > 0) 
             {
               ekin = (excitationEnergyGNP)/nbl;
               ekin2 = 0.0;
               for( i=0; i<nbl; i++ ) 
                  {
                    if( G4UniformRand() < sprob ) 
                      {
                        if(verboseLevel > 1) G4cout << " Particle skipped " << endl;
                        continue;
                      }
                    if( ekin2 > excitationEnergyGNP) break;
                    ran = G4UniformRand();
                    ekin1 = -ekin*log(ran) - cfa*(1.0+0.5*normal());
                    if (ekin1 < 0) ekin1 = -0.010*log(ran);
                    ekin2 += ekin1;
                    if( ekin2 > excitationEnergyGNP)
                    ekin1 = max( 1.0e-6, excitationEnergyGNP-(ekin2-ekin1) );
                    if( G4UniformRand() > (1.0-atomicNumber/(atomicWeight)))
                       pv[vecLen] = Proton;
                    else
                       pv[vecLen] = Neutron;
                    spall++;
                    cost = G4UniformRand() * 2.0 - 1.0;
                    sint = sqrt(fabs(1.0-cost*cost));
                    phi = M_2PI * G4UniformRand();
                    pv[vecLen].setFlag( true );                 // true is the same as IPA(i)<0
                    pv[vecLen].setSide( -4 );
                    pv[vecLen].setTOF( 1.0 );
                    pvMass = pv[vecLen].getMass();
                    pvEnergy = ekin1 + pvMass;
                    pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                    pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                     pp*sint*cos(phi),
                                                     pp*cost );
                    if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                    vecLen++;
                  }
               if( (atomicWeight >= 10.0 ) && (incidentKineticEnergy <= 2.0) ) 
                  {
                    G4int ika, kk = 0;
                    eka = incidentKineticEnergy;
                    if( eka > 1.0 )eka *= eka;
                    eka = max( 0.1, eka );
                    ika = G4int(3.6*exp((atomicNumber*atomicNumber
                                /atomicWeight-35.56)/6.45)/eka);
                    if( ika > 0 ) 
                      {
                        for( i=(vecLen-1); i>=0; i-- ) 
                           {
                             if( (pv[i].getCode() == protonCode) && pv[i].getFlag() ) 
                               {
                                 pTemp = pv[i];
                                 pv[i].setDefinition("Neutron");
                                 pv[i].setMomentumAndUpdate(pTemp.getMomentum());
                                 if (verboseLevel > 1) pv[i].Print(i);
                                 if( ++kk > ika ) break;
                               }
                           }
                      }
                  }
             }
         }

                                              // finished adding proton/neutron black track particles
                                              //  now, try to add deuterons, tritons and alphas
     
     if( excitationEnergyDTA >= 0.001 ) 
       {
         nbl = Poisson( (1.5+1.25*targ)*excitationEnergyDTA
                                      /(excitationEnergyGNP+excitationEnergyDTA));
       
                                      //    nbl is the number of deutrons, tritons, and alphas produced

         if (verboseLevel > 1) 
            G4cout << " evaporation " << targ << " " << nbl << " " << sprob << endl;        
         if( nbl > 0 ) 
           {
             ekin = excitationEnergyDTA/nbl;
             ekin2 = 0.0;
             for( i=0; i<nbl; i++ ) 
                {
                  if( G4UniformRand() < sprob ) 
                    {
                      if(verboseLevel > 1) G4cout << " Particle skipped " << endl;
                      continue;
                    } 
                  if( ekin2 > excitationEnergyDTA) break;
                  ran = G4UniformRand();
                  ekin1 = -ekin*log(ran)-cfa*(1.+0.5*normal());
                  if( ekin1 < 0.0 ) ekin1 = -0.010*log(ran);
                  ekin2 += ekin1;
                  if( ekin2 > excitationEnergyDTA)
                     ekin1 = max( 1.0e-6, excitationEnergyDTA-(ekin2-ekin1));
                  cost = G4UniformRand()*2.0 - 1.0;
                  sint = sqrt(fabs(1.0-cost*cost));
                  phi = M_2PI*G4UniformRand();
                  ran = G4UniformRand();
                  if( ran <= 0.60 ) 
                      pv[vecLen] = Deuteron;
                  else if (ran <= 0.90)
                      pv[vecLen] = Triton;
                  else
                      pv[vecLen] = Alpha;
                  spall += pv[vecLen].getMass() * 1.066;
                  if( spall > atomicWeight ) break;
                  pv[vecLen].setFlag( true );                    // true is the same as IPA(i)<0
                  pv[vecLen].setSide( -4 );
                  pvMass = pv[vecLen].getMass();
                  pv[vecLen].setTOF( 1.0 );
                  pvEnergy = pvMass + ekin1;
                  pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                  pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                   pp*sint*cos(phi),
                                                   pp*cost );
                  if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                  vecLen++;
                }
            }
        }
    }
   if( centerOfMassEnergy <= (4.0+G4UniformRand()) ) 
     {
       for( i=0; i<vecLen; i++ ) 
         {
           G4double etb = pv[i].getKineticEnergy();
           if( etb >= incidentKineticEnergy ) 
              pv[i].setKineticEnergyAndUpdate( incidentKineticEnergy );
         }
     }

   if(verboseLevel > 1)
     { G4cout << "Call TuningOfHighEnergyCacsading vecLen = " << vecLen << endl;
       incidentParticle.Print(0);
       targetParticle.Print(1);
       for (i=0; i<vecLen; i++) pv[i].Print(i);
     } 

   TuningOfHighEnergyCascading( pv, vecLen, 
                                incidentParticle, targetParticle, 
                                atomicWeight, atomicNumber);

   if(verboseLevel > 1)
     { G4cout << " After Tuning: " << endl;
       for(i=0; i<vecLen; i++) pv[i].Print(i);
     }
                                                     //  calculate time delay for nuclear reactions

   G4double tof = incidentTOF;
   if(     (atomicWeight >= 1.5) && (atomicWeight <= 230.0) 
        && (incidentKineticEnergy <= 0.2) )
           tof -= 500.0 * exp(-incidentKineticEnergy /0.04) * log( G4UniformRand() );

   successful = true;
   delete [] pvmx;
   return;
 }

void
G4HEInelastic::TuningOfHighEnergyCascading( G4HEVector pv[],
                                            G4int &vecLen,
                                            G4HEVector incidentParticle,
                                            G4HEVector targetParticle,
                                            G4double atomicWeight,
                                            G4double atomicNumber)
 {
   G4int i,j;
   G4double incidentKineticEnergy   = incidentParticle.getKineticEnergy();
   G4double incidentTotalMomentum   = incidentParticle.getTotalMomentum();
   G4double incidentCharge          = incidentParticle.getCharge(); 
   G4double incidentMass            = incidentParticle.getMass();
   G4int    incidentCode            = incidentParticle.getCode(); 
   G4double targetMass              = targetParticle.getMass();
   G4int    pionPlusCode            = PionPlus.getCode();
   G4int    pionMinusCode           = PionMinus.getCode();
   G4int    pionZeroCode            = PionZero.getCode();  
   G4int    protonCode              = Proton.getCode();
   G4int    neutronCode             = Neutron.getCode();
   G4HEVector *pvmx   = new G4HEVector [10];
   G4double   *reddec = new G4double [7];

   if ( incidentKineticEnergy > (25.+G4UniformRand()*75.) )
     {
       G4double reden = -0.7 + 0.29*log10(incidentKineticEnergy);
       G4double redat =  1.0 - 0.40*log10(atomicWeight);
       G4double pmax   = -200.;
       G4double pmapim = -200.;
       G4double pmapi0 = -200.;
       G4double pmapip = -200.;
       G4int ipmax  = 0;
       G4int maxpim = 0;
       G4int maxpi0 = 0;
       G4int maxpip = 0;
       G4int iphmf; 
       if (   (G4UniformRand() > (atomicWeight/100.-0.28)) 
           && (incidentCharge != 0.) )
         { 
           for (i=0; i < vecLen; i++)
             { 
               iphmf = pv[i].getCode();
               G4double ppp = pv[i].Length();
               if ( ppp > pmax) 
                 { 
                   pmax  = ppp; ipmax = i;
                 }
               if (iphmf == pionPlusCode && ppp > pmapip ) 
                 { 
                   pmapip = ppp; maxpip = i;       
                 }   
               else if (iphmf == pionZeroCode && ppp > pmapi0)
                 { 
                   pmapi0 = ppp; maxpi0 = i;
                 }
               else if (iphmf == pionMinusCode && ppp > pmapim)
                 { 
                   pmapim = ppp; maxpim = i;  
                 }                       
             }
         }
       if(verboseLevel > 1)
         {
           G4cout << "ipmax, pmax   " << ipmax  << " " << pmax   << endl;
           G4cout << "maxpip,pmapip " << maxpip << " " << pmapip << endl;
           G4cout << "maxpi0,pmapi0 " << maxpi0 << " " << pmapi0 << endl;
           G4cout << "maxpim,pmapim " << maxpim << " " << pmapim << endl; 
         }

       if ( vecLen > 2)
         { 
           for (i=2; i<vecLen; i++)
             { 
               iphmf = pv[i].getCode();
               if (    ((iphmf==protonCode)||(iphmf==neutronCode)||(pv[i].getMass()>1.7))   
                    && (pv[i].Length()<1.5)
                    && ((G4UniformRand()<reden)||(G4UniformRand()<redat)))
		  {
                    pv[i].setMomentumAndUpdate( 0., 0., 0.);
                    if(verboseLevel > 1)
                      {
                        G4cout << "zero Momentum for particle " << endl;
                        pv[i].Print(i);
                      }                                  
                  } 
             }  
         }
       if (maxpi0 == ipmax)
         { 
           if (G4UniformRand() < pmapi0/incidentTotalMomentum) 
            { 
              if ((incidentCharge > 0.) && (maxpip != 0))
                { 
                  G4ParticleMomentum mompi0 = pv[maxpi0].getMomentum();
                  pv[maxpi0].setMomentumAndUpdate( pv[maxpip].getMomentum() );
                  pv[maxpip].setMomentumAndUpdate( mompi0);
                  if(verboseLevel > 1)
                    {
                      G4cout << " exchange Momentum for " << maxpi0 << " and " << maxpip << endl;
                    } 
                }
              else if ((incidentCharge < 0.) && (maxpim != 0))
                { 
                  G4ParticleMomentum mompi0 = pv[maxpi0].getMomentum();
                  pv[maxpi0].setMomentumAndUpdate( pv[maxpim].getMomentum() );
                  pv[maxpim].setMomentumAndUpdate( mompi0);
                  if(verboseLevel > 1)
                    {
                      G4cout << " exchange Momentum for " << maxpi0 << " and " << maxpip << endl;
                    }   
                }
            }
         }
       G4double bntot = - incidentParticle.getBaryonNumber() - atomicWeight;
       for (i=0; i<vecLen; i++) bntot += pv[i].getBaryonNumber();
       bntot = ( (atomicWeight < 1.5) || (atomicWeight > (75.+G4UniformRand()*25.))) ? 
               0. : 1.+bntot/atomicWeight;
       if(verboseLevel > 1) 
         {
           G4cout << " Calculated Baryon- Number " << bntot << endl;
         }

       j = 0;
       for (i=0; i<vecLen; i++)
         { 
           G4double ppp = pv[i].Length();
           if ( ppp > 1.e-6)
             { 
               iphmf = pv[i].getCode();
               if(    (bntot > 0.3) 
                   && ((iphmf == protonCode) || (iphmf == neutronCode) || (pv[i].getMass() > 1.7) )
                   && (G4UniformRand() < 0.5) 
                   && (ppp < 1.2)  ) 
                 { 
                   if(verboseLevel > 1)
                     {
                       G4cout << " skip Baryon: " << endl;
                       pv[i].Print(i);
                     }
                   continue; 

                 }   
               if (j != i)
                 { 
                   pv[j] = pv[i];
                 }
               j++;
             }
         }
       vecLen = j;     
     }

   pvmx[0] = incidentParticle;
   pvmx[1] = targetParticle;
   pvmx[8].setZero();
   pvmx[2].Add(pvmx[0], pvmx[1]);
   pvmx[3].Lor(pvmx[0], pvmx[2]);
   pvmx[4].Lor(pvmx[1], pvmx[2]);
   
   if(verboseLevel > 1)
     {
       pvmx[0].Print(0);
       incidentParticle.Print(0);
       pvmx[1].Print(1);
       targetParticle.Print(1);
       pvmx[2].Print(2);
       pvmx[3].Print(3);
       pvmx[4].Print(4);
     }
 
   G4int ledpar = 0;
   G4double redpar = 0.;
   G4double gespar = 0.;
   G4int incidentS =   incidentParticle.getQuarkContent( 3 )
                     - incidentParticle.getAntiQuarkContent( 3 );
   G4int incidentB =   incidentParticle.getBaryonNumber();   

   if(verboseLevel > 1)
     {
       G4cout << " incidentS, incidentB " << incidentS << " " << incidentB << endl;
     }

   for (i=0; i<vecLen; i++)
     { G4int iphmf = pv[i].getCode();
       G4double ppp   = pv[i].Length();
       if (ppp > 1.e-3)
         { 
           pvmx[5].Lor( pv[i], pvmx[2] );
           G4double cost = pvmx[3].CosAng( pvmx[5] );
           G4int particleS =   pv[i].getQuarkContent( 3 )
                             - pv[i].getAntiQuarkContent( 3 );
           G4int particleB =   pv[i].getBaryonNumber();
           G4double hfmass;
           if (cost > 0.)
             { 
               reddec[0] = fabs( incidentMass - pv[i].getMass() );
               reddec[1] = fabs( incidentCharge - pv[i].getCharge());
               reddec[2] = fabs( G4double(incidentS - particleS) ); // cast for aCC
               reddec[3] = fabs( G4double(incidentB - particleB) ); // cast for aCC
               hfmass = incidentMass;  
             }
           else
             { 
               reddec[0] = fabs( targetMass - pv[i].getMass() );
               reddec[1] = fabs( atomicNumber/atomicWeight - pv[i].getCharge());
               reddec[2] = fabs( G4double(particleS) ); // cast for aCC
               reddec[3] = fabs( 1. - particleB );
               hfmass = targetMass;  
             } 
           reddec[5] = reddec[0]+reddec[1]+reddec[2]+reddec[3];
           G4double sbqwgt = reddec[5];
           if ( hfmass < 0.2)
             { 
               sbqwgt = (sbqwgt-2.5)*0.10;
               if(pv[i].getCode() == pionZeroCode) sbqwgt = 0.15; 
             } 
           else if (hfmass < 0.6)
               sbqwgt = (sbqwgt-3.0)*0.10;
           else 
             { 
               sbqwgt = (sbqwgt-2.0)*0.10;
               if(pv[i].getCode() == pionZeroCode) sbqwgt = 0.15;  
             }
           ppp = pvmx[5].Length();
           if ( (sbqwgt>0.) && (ppp>1.e-6))
             { 
               G4double pthmf = ppp*sqrt(1.-cost*cost);
               G4double plhmf = ppp*cost*(1.-sbqwgt);
               pvmx[7].Cross( pvmx[3], pvmx[5] );
               pvmx[7].Cross( pvmx[7], pvmx[3] );
               pvmx[6].SmulAndUpdate( pvmx[3], plhmf/pvmx[3].Length() );
               pvmx[7].SmulAndUpdate( pvmx[7], pthmf/pvmx[7].Length() );
               pvmx[5].Add3(pvmx[6], pvmx[7] );
               pvmx[5].setEnergy( sqrt(sqr(pvmx[5].Length()) + sqr(pv[i].getMass())));
               pv[i].Lor( pvmx[5], pvmx[4] );
               if(verboseLevel > 1) 
                 {
                   G4cout << " Particle Momentum changed to: " << endl;
                   pv[i].Print(i); 
		 }
             }    
           if (iphmf != pionZeroCode)
             { 
               pvmx[9].setMomentum( incidentParticle.getMomentum());
               pvmx[9].Smul( pvmx[9], 1.);
               pvmx[7].Sub3( pvmx[9], pv[i] );
               reddec[4] = pvmx[7].Length()/incidentTotalMomentum;
               reddec[6] = reddec[4]*2./3. + reddec[5]/12.;
               reddec[6] = max(0., 1. - reddec[6]);
               gespar   += reddec[6];
               if ( (reddec[5] <= 3.75) && (reddec[6] > redpar) )
                 { ledpar = i;
                   redpar = reddec[6]; 
                 }   
             } 
         }
       pvmx[8].Add3(pvmx[8], pv[i] );
     } 

  if (ledpar > 0)
     { 
      if(verboseLevel > 1)
        {
          G4cout << " Leading Particle found : " << ledpar << endl;
          pv[ledpar].Print(ledpar);
        }
       pvmx[4].Sub3(pvmx[8], pvmx[9] );
       G4double ptled = pv[ledpar].Length() * sqrt(1. - sqr(pvmx[9].CosAng( pv[ledpar] )));
       pv[ledpar].Sub3( pv[ledpar], pvmx[4] );
       pvmx[5].Smul( pvmx[9], pv[ledpar].Length()*pvmx[9].CosAng(pv[ledpar])/pvmx[9].Length());
       pvmx[6].Cross( pvmx[9], pv[ledpar]);
       pvmx[6].Cross( pvmx[6], pvmx[9]);
       pvmx[6].Smul( pvmx[6], ptled/pvmx[6].Length());
       pv[ledpar].Add3( pvmx[5], pvmx[6] );
       pv[ledpar].setMassAndUpdate( pv[ledpar].getMass());
       if(verboseLevel > 1)
	 {
           pv[ledpar].Print(ledpar);
         }  
     }
  delete [] reddec;
  delete [] pvmx;
  return;
 }     

void
G4HEInelastic::HighEnergyClusterProduction(G4bool &successful,
                                           G4HEVector pv[],
                                           G4int &vecLen,	
                                           G4double &excitationEnergyGNP,
                                           G4double &excitationEnergyDTA,
                                           G4HEVector incidentParticle,
                                           G4HEVector targetParticle,
                                           G4double atomicWeight,
                                           G4double atomicNumber)
 {   
// For low multiplicity in the first intranuclear interaction the cascading process
// as described in G4HEInelastic::MediumEnergyCascading does not work 
// satisfactorily. From experimental data it is strongly suggested to use 
// a two- body resonance model.   
//  
//  All quantities on the G4HEVector Array pv are in GeV- units.

   G4int protonCode       = Proton.getCode();
   G4double protonMass    = Proton.getMass();
   G4int neutronCode      = Neutron.getCode();
   G4double neutronMass   = Neutron.getMass();
   G4double kaonPlusMass  = KaonPlus.getMass();
   G4int kaonPlusCode     = KaonPlus.getCode();   
   G4int kaonMinusCode    = KaonMinus.getCode();
   G4int kaonZeroSCode    = KaonZeroShort.getCode(); 
   G4int kaonZeroLCode    = KaonZeroLong.getCode();
   G4int kaonZeroCode     = KaonZero.getCode();
   G4int antiKaonZeroCode = AntiKaonZero.getCode(); 
   G4int pionPlusCode     = PionPlus.getCode();    
   G4int pionZeroCode     = PionZero.getCode();    
   G4int pionMinusCode    = PionMinus.getCode(); 
   G4String mesonType     = PionPlus.getType();
   G4String baryonType    = Proton.getType(); 
   
   G4int    targetCode   = targetParticle.getCode();
   G4double targetMass   = targetParticle.getMass();
   G4double targetCharge = targetParticle.getCharge();

   G4int    incidentCode          = incidentParticle.getCode();
   G4double incidentMass          = incidentParticle.getMass();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   G4double incidentCharge        = incidentParticle.getCharge();
   G4double incidentEnergy        = incidentParticle.getEnergy();
   G4double incidentKineticEnergy = incidentParticle.getKineticEnergy();
   G4String incidentType          = incidentParticle.getType();
//   G4double incidentTOF           = incidentParticle.getTOF();   
   G4double incidentTOF           = 0.;
   
                                                              // some local variables

   G4int i, j, k, l;
   
   if(verboseLevel > 1) G4cout << " G4HEInelastic::HighEnergyClusterProduction " << endl;

   successful = false; 
   if(incidentTotalMomentum < 25. + G4UniformRand()*25.) return;

   G4double centerOfMassEnergy = sqrt( sqr(incidentMass) + sqr(targetMass)
                                       +2.*targetMass*incidentEnergy);

   G4HEVector pvI = incidentParticle;  // for the incident particle
   pvI.setSide( 1 );

   G4HEVector pvT = targetParticle;   // for the target particle
   pvT.setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvT.setSide( -1 );
   pvT.setTOF( -1.); 
                              // distribute particles in forward and backward
                              // hemispheres. Note that only low multiplicity
                              // events from FirstIntInNuc.... should go into 
                              // this routine. 
   G4int targ  = 0;  
   G4int ifor  = 0; 
   G4int iback = 0;
   G4int pvCode;
   G4double pvMass, pvEnergy; 

   pv[0].setSide(  1 );
   pv[1].setSide( -1 );
   for(i = 0; i < vecLen; i++)
      {
        if (i > 1)
           {
              if( G4UniformRand() < 0.5) 
                {
                  pv[i].setSide( 1 );
                  if (++ifor > 18) 
                     { 
                       pv[i].setSide( -1 );
                       ifor--;
                       iback++;
                     }
                }
              else
                {
                  pv[i].setSide( -1 );
                  if (++iback > 18)
                     { 
                       pv[i].setSide( 1 );
                       ifor++;
                       iback--;
                     }
                }
           }

        pvCode = pv[i].getCode();

        if (   (    (incidentCode == protonCode) || (incidentCode == neutronCode)
                 || (incidentType == mesonType) )
            && (    (pvCode == pionPlusCode) || (pvCode == pionMinusCode) )
            && (    (G4UniformRand() < (10.-incidentTotalMomentum)/6.) )
            && (    (G4UniformRand() < atomicWeight/300.) ) )
           { 
               if (G4UniformRand() > atomicNumber/atomicWeight) 
                  pv[i].setDefinition( "Neutron" );
               else
                  pv[i].setDefinition( "Proton" );
               targ++;
           }    
        pv[i].setTOF( incidentTOF );                 
      }    
   G4double tb = 2. * iback;
   if (centerOfMassEnergy < (2+G4UniformRand())) tb = (2.*iback + vecLen)/2.;
   
   G4double nucsup[] = { 1.0, 0.7, 0.5, 0.4, 0.35, 0.3};
   G4double  psup[]   = { 3. , 6. , 20., 50., 100.,1000.};   
   G4double s = centerOfMassEnergy*centerOfMassEnergy;

   G4double xtarg = max(0.01, min(0.50, 0.312+0.2*log(log(s))+pow(s,1.5)/6000.) 
                             * (pow(atomicWeight,0.33)-1.) * tb);
   G4int momentumBin = 0;
   while( (momentumBin < 6) && (incidentTotalMomentum > psup[momentumBin])) momentumBin++;
   momentumBin = min(5, momentumBin);
   G4double xpnhmf = max(0.01,xtarg*nucsup[momentumBin]);
   G4double xshhmf = max(0.01,xtarg-xpnhmf);
   G4double rshhmf = 0.25*xshhmf;
   G4double rpnhmf = 0.81*xpnhmf;
   G4double xhmf;
   G4int nshhmf, npnhmf;
   if (rshhmf > 1.1)
     {
       rshhmf = xshhmf/(rshhmf-1.);
       if (rshhmf <= 20.)
         xhmf = GammaRand( rshhmf );
       else
         xhmf = Erlang( G4int(rshhmf+0.5) );
       xshhmf *= xhmf/rshhmf;
     }
   nshhmf = Poisson( xshhmf );
   if (rpnhmf > 1.1)
     {
       rpnhmf = xpnhmf/(rpnhmf -1.);
       if (rpnhmf <= 20.)
         xhmf = GammaRand( rpnhmf );
       else
         xhmf = Erlang( G4int(rpnhmf+0.5) );
       xpnhmf *= xhmf/rpnhmf;
     }
   npnhmf = Poisson( xpnhmf );
   G4int ntarg = nshhmf + npnhmf;

   while (npnhmf > 0)
     {
       if ( G4UniformRand() > (1. - atomicNumber/atomicWeight))
          pv[vecLen].setDefinition( "Proton" );
       else
          pv[vecLen].setDefinition( "Neutron" );
       targ++;
       pv[vecLen].setSide( -2 );
       pv[vecLen].setFlag( true );
       pv[vecLen].setTOF( incidentTOF );
       vecLen++;
       npnhmf--;
     }
   while (nshhmf > 0)
     {
       G4double ran = G4UniformRand();
       if (ran < 0.333333 )
          pv[vecLen].setDefinition( "PionPlus" );
       else if (ran < 0.6667)
          pv[vecLen].setDefinition( "PionZero" );
       else
          pv[vecLen].setDefinition( "PionMinus" );
       pv[vecLen].setSide( -2 );
       pv[vecLen].setFlag( true );
       pv[vecLen].setTOF( incidentTOF );
       vecLen++;
       nshhmf--;
     }    
                              // mark leading particles for incident strange particles 
                              // and antibaryons, for all other we assume that the first 
                              // and second particle are the leading particles. 
                              // We need this later for kinematic aspects of strangeness
                              // conservation.
                          
   G4int lead = 0;                   
   G4HEVector leadParticle;
   if( (incidentMass >= kaonPlusMass-0.05) && (incidentCode != protonCode)          
                                           && (incidentCode != neutronCode) ) 
         {       
           G4double pMass = pv[0].getMass();
           G4int    pCode = pv[0].getCode();
           if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                            && (pCode != neutronCode) ) 
                  {       
                    lead = pCode; 
                    leadParticle = pv[0];                           
                  } 
           else   
                  {
                    pMass = pv[1].getMass();
                    pCode = pv[1].getCode();
                    if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                                     && (pCode != neutronCode) ) 
                        {       
                          lead = pCode;
                          leadParticle = pv[1];
                        }
                  }
         }

   if (verboseLevel > 1)
      { G4cout << " pv Vector after initialization " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      }     

   G4double tavai = 0.;
   for(i=0;i<vecLen;i++) if(pv[i].getSide() != -2) tavai  += pv[i].getMass();

   while (tavai > centerOfMassEnergy)
      {
         for (i=vecLen-1; i >= 0; i--)
            {
              if (pv[i].getSide() != -2)
                 {
                    tavai -= pv[i].getMass();
                    if( i != vecLen-1) 
                      {
                        for (j=i; j < vecLen; j++)
                           {
                              pv[j]  = pv[j+1];
                           }
                      }
                    if ( --vecLen < 2)
                      {
                        successful = false;
                        return;
                      }
                    break;
                 }
            }    
      }
                                          // now we produce 3 Clusters:
                                          // 1. forward cluster
                                          // 2. backward meson cluster
                                          // 3. backward nucleon cluster
   G4double rmc0 = 0., rmd0 = 0., rme0 = 0.;
   G4int    ntc  = 0,  ntd  = 0,  nte  = 0;
   
   for (i=0; i < vecLen; i++)
      {
        if(pv[i].getSide() > 0)
           {
             if(ntc < 17)
               { 
                 rmc0 += pv[i].getMass();
                 ntc++;
               }
             else
               {
                 if(ntd < 17)
                   {
                     pv[i].setSide(-1);
                     rmd0 += pv[i].getMass();
                     ntd++;
                   }
                 else
                   {
                     pv[i].setSide(-2);
                     rme0 += pv[i].getMass();
                     nte++;
                   }
               }
           }  
        else if (pv[i].getSide() == -1)
           {
             if(ntd < 17)
                {
                   rmd0 += pv[i].getMass();
                   ntd++;
                }
             else
                {
                   pv[i].setSide(-2); 
                   rme0 += pv[i].getMass();
                   nte++;
                }           
           }
        else
	   {
             rme0 += pv[i].getMass();
             nte++;
           } 
      }         

   G4double cpar[] = {0.6, 0.6, 0.35, 0.15, 0.10};
   G4double gpar[] = {2.6, 2.6, 1.80, 1.30, 1.20};
    
   G4double rmc = rmc0, rmd = rmd0, rme = rme0; 
   G4int ntc1 = min(4,ntc-1);
   G4int ntd1 = min(4,ntd-1);
   G4int nte1 = min(4,nte-1);
   if (ntc > 1) rmc = rmc0 + pow(-log(1.-G4UniformRand()),cpar[ntc1])/gpar[ntc1];
   if (ntd > 1) rmd = rmd0 + pow(-log(1.-G4UniformRand()),cpar[ntd1])/gpar[ntd1];
   if (nte > 1) rme = rme0 + pow(-log(1.-G4UniformRand()),cpar[nte1])/gpar[nte1];
   while( (rmc+rmd) > centerOfMassEnergy)
      {
        if ((rmc == rmc0) && (rmd == rmd0))
          { 
            rmd *= 0.999*centerOfMassEnergy/(rmc+rmd);
            rmc *= 0.999*centerOfMassEnergy/(rmc+rmd);
          }
        else
          {
            rmc = 0.1*rmc0 + 0.9*rmc;
            rmd = 0.1*rmd0 + 0.9*rmd;
          }   
      }             
   if(verboseLevel > 1) G4cout << " Cluster Masses: " << ntc << " " << rmc << " " << ntd
                             << " " << rmd << " " << nte << " " << rme << endl;
 
   
   G4HEVector* pvmx = new G4HEVector[11];

   pvmx[1].setMass( incidentMass);
   pvmx[1].setMomentumAndUpdate( 0., 0., incidentTotalMomentum);
   pvmx[2].setMass( targetMass);
   pvmx[2].setMomentumAndUpdate( 0., 0., 0.);
   pvmx[0].Add( pvmx[1], pvmx[2] );
   pvmx[1].Lor( pvmx[1], pvmx[0] );
   pvmx[2].Lor( pvmx[2], pvmx[0] );

   G4double pf = sqrt(max(0.0001,  sqr(sqr(centerOfMassEnergy) + rmd*rmd -rmc*rmc)
                                 - 4*sqr(centerOfMassEnergy)*rmd*rmd))/(2.*centerOfMassEnergy);
   pvmx[3].setMass( rmc );
   pvmx[4].setMass( rmd );
   pvmx[3].setEnergy( sqrt(pf*pf + rmc*rmc) );
   pvmx[4].setEnergy( sqrt(pf*pf + rmd*rmd) );
   
   G4double tvalue = -MAXFLOAT;
   G4double bvalue = max(0.01, 4.0 + 1.6*log(incidentTotalMomentum));
   if (bvalue != 0.0) tvalue = log(G4UniformRand())/bvalue;
   G4double pin = pvmx[1].Length();
   G4double tacmin = sqr( pvmx[1].getEnergy() - pvmx[3].getEnergy()) - sqr( pin - pf);
   G4double ctet   = max(-1., min(1., 1.+2.*(tvalue-tacmin)/max(1.e-10, 4.*pin*pf)));
   G4double stet   = sqrt(max(0., 1.0 - ctet*ctet));
   G4double phi    = M_2PI * G4UniformRand();
   pvmx[3].setMomentum( pf * stet * sin(phi), 
                        pf * stet * cos(phi),
                        pf * ctet           );
   pvmx[4].Smul( pvmx[3], -1.);
   
   if (nte > 0)
      {
        G4double ekit1 = 0.04;
        G4double ekit2 = 0.6;
        G4double gaval = 1.2;
        if (incidentKineticEnergy <= 5.)
           {
             ekit1 *= sqr(incidentKineticEnergy)/25.;
             ekit2 *= sqr(incidentKineticEnergy)/25.;
           }
        G4double avalue = (1.-gaval)/(pow(ekit2, 1.-gaval)-pow(ekit1, 1.-gaval));
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -2)
                 {
                   G4double ekit = pow(G4UniformRand()*(1.-gaval)/avalue +pow(ekit1, 1.-gaval),
                                       1./(1.-gaval));
                   pv[i].setKineticEnergyAndUpdate( ekit );
                   ctet = max(-1., min(1., log(2.23*G4UniformRand()+0.383)/0.96));
                   stet = sqrt( max( 0.0, 1. - ctet*ctet ));
                   phi  = G4UniformRand()*M_2PI;
                   G4double pp = pv[i].Length();
                   pv[i].setMomentum( pp * stet * sin(phi),
                                      pp * stet * cos(phi),
                                      pp * ctet           );
                   pv[i].Lor( pv[i], pvmx[0] );
                 }              
            }             
      }
   pvmx[1] = pvmx[3];
   pvmx[2] = pvmx[4];
   pvmx[5].SmulAndUpdate( pvmx[3], -1.);
   pvmx[6].SmulAndUpdate( pvmx[4], -1.);

   if (verboseLevel > 1)
      { G4cout << " General vectors before Phase space Generation " << endl;
        for (i=0; i<7; i++) pvmx[i].Print(i);
      }  


   G4HEVector* tempV = new G4HEVector[18];
   G4bool constantCrossSection = true;
   G4double wgt;
   G4int npg;

   if (ntc > 1)
      {
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() > 0)
                 {
                   tempV[npg++] = pv[i];
                 }
            }
        wgt = NBodyPhaseSpace( pvmx[3].getMass(), constantCrossSection, tempV, npg);
                     
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() > 0)
                 {
                   pv[i].setMomentum( tempV[npg++].getMomentum());
                   pv[i].SmulAndUpdate( pv[i], 1. );
                   pv[i].Lor( pv[i], pvmx[5] );
                 }
            }
      }
   else if(ntc == 1)
      {
        for(i=0; i<vecLen; i++)
          {
            if(pv[i].getSide() > 0) pv[i].setMomentumAndUpdate(pvmx[3].getMomentum());
          }
      }
   else
      {
      }
     
   if (ntd > 1)
      {
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -1)
                 {
                   tempV[npg++] = pv[i];
                 }
            }
        wgt = NBodyPhaseSpace( pvmx[4].getMass(), constantCrossSection, tempV, npg);
                     
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -1)
                 {
                   pv[i].setMomentum( tempV[npg++].getMomentum());
                   pv[i].SmulAndUpdate( pv[i], 1.);
                   pv[i].Lor( pv[i], pvmx[6] );
                 }
            }
      }
   else if(ntd == 1)
      {
        for(i=0; i<vecLen; i++)
          {
            if(pv[i].getSide() == -1) pv[i].setMomentumAndUpdate(pvmx[4].getMomentum());
          }
      }
   else
      {
      }

   if(verboseLevel > 1)
     {
       G4cout << " Vectors after PhaseSpace generation " << endl;
       for(i=0; i<vecLen; i++) pv[i].Print(i);
     }
                                                   //  Lorentz transformation in lab system
   targ = 0;
   for( i=0; i < vecLen; i++ ) 
      {
        if( pv[i].getType() == baryonType )targ++;
        pv[i].Lor( pv[i], pvmx[2] );
      }
   targ = max( 1, targ );

   if(verboseLevel > 1)
     {
       G4cout << " Transformation in Lab- System " << endl;
       for(i=0; i<vecLen; i++) pv[i].Print(i);
     }

   G4bool dum;
   G4double ekin, teta;

   if( lead ) 
     {
       for( i=0; i<vecLen; i++ ) 
          {
            if( pv[i].getCode() == lead ) 
              {
                dum = false;
                break;
              }
          }
       if( dum ) 
         {
           i = 0;          
 
           if(   (    (leadParticle.getType() == baryonType)
                   && (pv[1].getType() == baryonType))
              || (    (leadParticle.getType() == mesonType)
                   && (pv[1].getType() == mesonType)))
             {
               i = 1;
             } 

            ekin = pv[i].getKineticEnergy();
            pv[i] = leadParticle;
            if( pv[i].getFlag() )
                pv[i].setTOF( -1.0 );
            else
                pv[i].setTOF( 1.0 );
            pv[i].setKineticEnergyAndUpdate( ekin );
         }
     }

   pvmx[4].setMass( incidentMass);
   pvmx[4].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   
   G4double ekin0 = pvmx[4].getKineticEnergy();
   
   pvmx[5].setMass ( protonMass * targ);
   pvmx[5].setEnergy( protonMass * targ);
   pvmx[5].setKineticEnergy(0.);
   pvmx[5].setMomentum( 0.0, 0.0, 0.0 );

   ekin = pvmx[4].getEnergy() + pvmx[5].getEnergy();

   pvmx[6].Add( pvmx[4], pvmx[5] );
   pvmx[4].Lor( pvmx[4], pvmx[6] );
   pvmx[5].Lor( pvmx[5], pvmx[6] );
   
   G4double tecm = pvmx[4].getEnergy() + pvmx[5].getEnergy();

   pvmx[8].setZero();
   
   G4double ekin1 = 0.0;   
   
   for( i=0; i < vecLen; i++ ) 
      {
        pvmx[8].Add( pvmx[8], pv[i] );
        ekin1 += pv[i].getKineticEnergy();
        ekin  -= pv[i].getMass();
      }
   
   if( vecLen > 1 && vecLen < 19 ) 
     {
       constantCrossSection = true;
       G4HEVector pw[18];
       for(i=0;i<vecLen;i++) pw[i] = pv[i];
       wgt = NBodyPhaseSpace( tecm, constantCrossSection, pw, vecLen );
       ekin = 0.0;
       for( i=0; i < vecLen; i++ ) 
          {
            pvmx[7].setMass( pw[i].getMass());
            pvmx[7].setMomentum( pw[i].getMomentum() );
            pvmx[7].SmulAndUpdate( pvmx[7], 1.);
            pvmx[7].Lor( pvmx[7], pvmx[5] );
            ekin += pvmx[7].getKineticEnergy();
          }
       teta = pvmx[8].Ang( pvmx[4] );
       if (verboseLevel > 1) 
           G4cout << " vecLen > 1 && vecLen < 19 " << teta << " " << ekin0 << " " << ekin1 << " " << ekin << endl;
     }

   if( ekin1 != 0.0 ) 
     {
       pvmx[7].setZero();
       wgt = ekin/ekin1;
       ekin1 = 0.;
       for( i=0; i < vecLen; i++ ) 
          {
            pvMass = pv[i].getMass();
            ekin   = pv[i].getKineticEnergy() * wgt;
            pv[i].setKineticEnergyAndUpdate( ekin );
            ekin1 += ekin;
            pvmx[7].Add( pvmx[7], pv[i] );
          }
       teta = pvmx[7].Ang( pvmx[4] );
       if (verboseLevel > 1) G4cout << " ekin1 != 0 " << teta << " " << ekin0 << " " << ekin1 << endl;
     }

   if(verboseLevel > 1)
     {
       G4cout << " After energy- correction " << endl;
       for(i=0; i<vecLen; i++) pv[i].Print(i);
     }      

                          //    Rotate in direction of z-axis, this does disturb in some way our
                          //    inclusive distributions, but it is necessary for momentum conservation

   pvmx[7].setZero();
   for( i=0; i < vecLen; i++ ) pvmx[7].Add( pvmx[7], pv[i] );
   
                          //  Some smearing in transverse direction from Fermi motion
   


                                                  // RANNOR
   G4double ry   = G4UniformRand();
   G4double rz   = G4UniformRand();
   G4double rx   = M_2PI*rz;
   G4double a1   = sqrt(-2.0*log(ry));
   G4double ran1 = a1*sin(rx);
   G4double ran2 = a1*cos(rx);
                                                  
   pvmx[7].setMomentum( pvmx[7].getMomentum().x()+ran1*0.020*targ,
                        pvmx[7].getMomentum().y()+ran2*0.020*targ );
   pvmx[8].Defs( pvmx[4], pvmx[7], pvmx[9], pvmx[10] );
   pvmx[7].setZero();
   for (i = 0; i < vecLen; i++) 
      { 
        pv[i].Trac( pv[i], pvmx[8] , pvmx[9], pvmx[10] );
        pvmx[7].Add( pvmx[7], pv[i] );  
      }
   teta = pvmx[7].Ang( pvmx[4] );   
   if (verboseLevel > 1) G4cout << " After smearing " << teta << endl;

                               //  Rotate in direction of primary particle, subtract binding energies
                               //   and make some further corrections if required

   G4double dekin = 0.0;
   G4int npions = 0;    
   G4double ek1 = 0.0;
   G4double alekw, xxh;
   G4double cfa = 0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.);
   G4double alem[] = {1.40, 2.30, 2.70, 3.00, 3.40, 4.60, 7.00, 10.0};
   G4double val0[] = {0.00, 0.40, 0.48, 0.51, 0.54, 0.60, 0.65, 0.70};
   

   for (i = 0; i < vecLen; i++)
      { 
        pv[i].Defs1( pv[i], pvI );
        if (atomicWeight > 1.5)
           {
             ekin  = max( 1.e-6,pv[i].getKineticEnergy() - cfa*( 1. + 0.5*normal()));
             alekw = log( incidentKineticEnergy );
             xxh   = 1.;
             xxh   = 1.;
             if(incidentCode == pionPlusCode || incidentCode == pionMinusCode)
	       {
                 if(pv[i].getCode() == pionZeroCode)
                   {
                     if(G4UniformRand() < log(atomicWeight))
                       { 
                         if (alekw > alem[0])
                           {
                              for (j = 1; j < 8; j++)
                                 {
                                    if(alekw < alem[j]) break;
                                 } 
                              xxh =   (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alekw
                                     + val0[j-1] - (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alem[j-1];
                              xxh = 1. - xxh;
                           }
                       }      
                    }
               }  
             dekin += ekin*(1.-xxh);
             ekin *= xxh;
             pv[i].setKineticEnergyAndUpdate( ekin );
             pvCode = pv[i].getCode();
             if ((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode))
                {
                  npions += 1;
                  ek1 += ekin; 
                }
           }
      }
   if( (ek1 > 0.0) && (npions > 0) ) 
      {
        dekin = 1.+dekin/ek1;
        for (i = 0; i < vecLen; i++)
          {
            pvCode = pv[i].getCode();
            if((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode)) 
              {
                ekin = max( 1.0e-6, pv[i].getKineticEnergy() * dekin );
                pv[i].setKineticEnergyAndUpdate( ekin );
              }
          }
      }
   if (verboseLevel > 1)
      { G4cout << " Lab-System " << ek1 << " " << npions << endl;
        for (i=0; i<vecLen; i++) pv[i].Print(i);
      }

                                       //  add black track particles
                                       //  the total number of particles produced is restricted to 198
                                       //  this may have influence on very high energies

   if (verboseLevel > 1) 
       G4cout << " Evaporation " << atomicWeight << " " << excitationEnergyGNP
            << " " << excitationEnergyDTA << endl;

   G4double sprob = 0.;
   if (incidentKineticEnergy > 5. )
       sprob = min( 1., (0.25-0.02*log(atomicWeight))*log(incidentKineticEnergy-4.) );

   if( atomicWeight > 1.5 && G4UniformRand() > sprob) 
     {

       G4double ekw, tex, cost, sint, ekin2, ran, pp, eka;
       G4int spall, nbl;
                                     //  first add protons and neutrons

       if( excitationEnergyGNP >= 0.001 ) 
         {
                                     //  nbl = number of proton/neutron black track particles
                                     //  tex is their total kinetic energy (GeV)
       
           nbl = Poisson( (1.5+1.25*targ)*excitationEnergyGNP/
                                         (excitationEnergyGNP+excitationEnergyDTA));
           if( targ+nbl > atomicWeight ) nbl = atomicWeight - targ;
           if (verboseLevel > 1) G4cout << " evaporation " << targ << " " << nbl << " " << sprob << endl; 
           spall = targ;
           if( nbl > 0) 
             {
               ekin = excitationEnergyGNP/nbl;
               ekin2 = 0.0;
               for( i=0; i<nbl; i++ ) 
                  {
                    if( G4UniformRand() < sprob ) continue;
                    if( ekin2 > excitationEnergyGNP) break;
                    ran = G4UniformRand();
                    ekin1 = -ekin*log(ran) - cfa*(1.0+0.5*normal());
                    if (ekin1 < 0) ekin1 = -0.010*log(ran);
                    ekin2 += ekin1;
                    if( ekin2 > excitationEnergyGNP)
                    ekin1 = max( 1.0e-6, excitationEnergyGNP-(ekin2-ekin1) );
                    if( G4UniformRand() > (1.0-atomicNumber/(atomicWeight)))
                       pv[vecLen].setDefinition( "Proton");
                    else
                       pv[vecLen].setDefinition( "Neutron" );
                    spall++;
                    cost = G4UniformRand() * 2.0 - 1.0;
                    sint = sqrt(fabs(1.0-cost*cost));
                    phi = M_2PI * G4UniformRand();
                    pv[vecLen].setFlag( true );                 // true is the same as IPA(i)<0
                    pv[vecLen].setSide( -4 );
                    pvMass = pv[vecLen].getMass();
                    pv[vecLen].setTOF( 1.0 );
                    pvEnergy = ekin1 + pvMass;
                    pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                    pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                     pp*sint*cos(phi),
                                                     pp*cost );
                    if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                    vecLen++;
                  }
               if( (atomicWeight >= 10.0 ) && (incidentKineticEnergy <= 2.0) ) 
                  {
                    G4int ika, kk = 0;
                    eka = incidentKineticEnergy;
                    if( eka > 1.0 )eka *= eka;
                    eka = max( 0.1, eka );
                    ika = G4int(3.6*exp((atomicNumber*atomicNumber
                                /atomicWeight-35.56)/6.45)/eka);
                    if( ika > 0 ) 
                      {
                        for( i=(vecLen-1); i>=0; i-- ) 
                           {
                             if( (pv[i].getCode() == protonCode) && pv[i].getFlag() ) 
                               {
                                 G4HEVector pTemp = pv[i];
                                 pv[i].setDefinition( "Neutron" );
                                 pv[i].setMomentumAndUpdate(pTemp.getMomentum());
                                 if (verboseLevel > 1) pv[i].Print(i);
                                 if( ++kk > ika ) break;
                               }
                           }
                      }
                  }
             }
         }

                                              // finished adding proton/neutron black track particles
                                              //  now, try to add deuterons, tritons and alphas
     
     if( excitationEnergyDTA >= 0.001 ) 
       {
         nbl = Poisson( (1.5+1.25*targ)*excitationEnergyDTA
                                      /(excitationEnergyGNP+excitationEnergyDTA));
       
                                      //    nbl is the number of deutrons, tritons, and alphas produced
       
         if( nbl > 0 ) 
           {
             ekin = excitationEnergyDTA/nbl;
             ekin2 = 0.0;
             for( i=0; i<nbl; i++ ) 
                {
                  if( G4UniformRand() < sprob ) continue;
                  if( ekin2 > excitationEnergyDTA) break;
                  ran = G4UniformRand();
                  ekin1 = -ekin*log(ran)-cfa*(1.+0.5*normal());
                  if( ekin1 < 0.0 ) ekin1 = -0.010*log(ran);
                  ekin2 += ekin1;
                  if( ekin2 > excitationEnergyDTA )
                     ekin1 = max( 1.0e-6, excitationEnergyDTA-(ekin2-ekin1));
                  cost = G4UniformRand()*2.0 - 1.0;
                  sint = sqrt(fabs(1.0-cost*cost));
                  phi = M_2PI*G4UniformRand();
                  ran = G4UniformRand();
                  if( ran <= 0.60 ) 
                      pv[vecLen].setDefinition( "Deuteron");
                  else if (ran <= 0.90)
                      pv[vecLen].setDefinition( "Triton" );
                  else
                      pv[vecLen].setDefinition( "Alpha" );
                  spall += pv[vecLen].getMass() * 1.066;
                  if( spall > atomicWeight ) break;
                  pv[vecLen].setFlag( true );                    // true is the same as IPA(i)<0
                  pv[vecLen].setSide( -4 );
                  pvMass = pv[vecLen].getMass();
                  pv[vecLen].setTOF( 1.0 );
                  pvEnergy = pvMass + ekin1;
                  pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                  pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                   pp*sint*cos(phi),
                                                   pp*cost );
                  if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                  vecLen++;
                }
            }
        }
    }
   if( centerOfMassEnergy <= (4.0+G4UniformRand()) ) 
     {
       for( i=0; i<vecLen; i++ ) 
         {
           G4double etb = pv[i].getKineticEnergy();
           if( etb >= incidentKineticEnergy ) 
              pv[i].setKineticEnergyAndUpdate( incidentKineticEnergy );
         }
     }

   TuningOfHighEnergyCascading( pv, vecLen,
                                incidentParticle, targetParticle,
                                atomicWeight, atomicNumber);    

                                                     //  calculate time delay for nuclear reactions
   G4double tof = incidentTOF;
   if(     (atomicWeight >= 1.5) && (atomicWeight <= 230.0) 
        && (incidentKineticEnergy <= 0.2) )
           tof -= 500.0 * exp(-incidentKineticEnergy /0.04) * log( G4UniformRand() );
   for ( i=0; i < vecLen; i++)     
     { 
       
       pv[i].setTOF ( tof );
//       vec[i].SetTOF ( tof );
     }

   successful = true;
   delete [] pvmx;
   delete [] tempV;
   return;
 }

void
G4HEInelastic::MediumEnergyCascading(G4bool &successful,
                                     G4HEVector pv[],
                                     G4int &vecLen,	
                                     G4double &excitationEnergyGNP,
                                     G4double &excitationEnergyDTA,
                                     G4HEVector incidentParticle,
                                     G4HEVector targetParticle,
                                     G4double atomicWeight,
                                     G4double atomicNumber)
 {   
//
//  The multiplicity of particles produced in the first interaction has been
//  calculated in one of the FirstIntInNuc.... routines. The nuclear
//  cascading particles are parametrized from experimental data.
//  A simple single variable description E D3S/DP3= F(Q) with
//  Q^2 = (M*X)^2 + PT^2 is used. Final state kinematic is produced
//  by an FF-type iterative cascade method.
//  Nuclear evaporation particles are added at the end of the routine.

//  All quantities on the G4HEVector Array pv are in GeV- units.

   G4int protonCode       = Proton.getCode();
   G4double protonMass    = Proton.getMass();
   G4int neutronCode      = Neutron.getCode();
   G4double neutronMass   = Neutron.getMass();
   G4double kaonPlusMass  = KaonPlus.getMass();
   G4int kaonPlusCode     = KaonPlus.getCode();   
   G4int kaonMinusCode    = KaonMinus.getCode();
   G4int kaonZeroSCode    = KaonZeroShort.getCode(); 
   G4int kaonZeroLCode    = KaonZeroLong.getCode();
   G4int kaonZeroCode     = KaonZero.getCode();
   G4int antiKaonZeroCode = AntiKaonZero.getCode(); 
   G4int pionPlusCode     = PionPlus.getCode();    
   G4int pionZeroCode     = PionZero.getCode();    
   G4int pionMinusCode    = PionMinus.getCode(); 
   G4String mesonType     = PionPlus.getType();
   G4String baryonType    = Proton.getType(); 
   
   G4int    targetCode   = targetParticle.getCode();
   G4double targetMass   = targetParticle.getMass();
   G4double targetCharge = targetParticle.getCharge();

   G4int    incidentCode          = incidentParticle.getCode();
   G4double incidentMass          = incidentParticle.getMass();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   G4double incidentCharge        = incidentParticle.getCharge();
   G4double incidentEnergy        = incidentParticle.getEnergy();
   G4double incidentKineticEnergy = incidentParticle.getKineticEnergy();
   G4String incidentType          = incidentParticle.getType();
//   G4double incidentTOF           = incidentParticle.getTOF();   
   G4double incidentTOF           = 0.;
   
                                                              // some local variables

   G4int i, j, k, l;

   if(verboseLevel > 1) G4cout << " G4HEInelastic::MediumEnergyCascading " << endl;

                              // define annihilation channels.
                                 
   G4bool annihilation = false;
   if (incidentCode < 0 && incidentType == baryonType && 
       pv[0].getType() != baryonType &&
       pv[1].getType() != baryonType   )
         { 
           annihilation = true;
         }   
     
   successful = false; 

   G4double twsup[] = { 1., 1., 0.7, 0.5, 0.3, 0.2, 0.1, 0.0 };
   
   if(annihilation) goto start;
   if(vecLen >= 8) goto start;
   if(incidentKineticEnergy < 1.) return;
   if(    (   incidentCode == kaonPlusCode  || incidentCode == kaonMinusCode
           || incidentCode == kaonZeroCode  || incidentCode == antiKaonZeroCode
           || incidentCode == kaonZeroSCode || incidentCode == kaonZeroLCode )
       && ( G4UniformRand() < 0.5)) goto start;
   if(G4UniformRand() > twsup[vecLen-1]) goto start;
   return;

   start:  
   
   if (annihilation)
      {                        // do some corrections of incident particle kinematic
        G4double ekcor = max( 1., 1./incidentKineticEnergy);
        incidentKineticEnergy = 2*targetMass + incidentKineticEnergy*(1.+ekcor/atomicWeight);
        G4double excitation = NuclearExcitation(incidentKineticEnergy,
                                                atomicWeight,
                                                atomicNumber,
                                                excitationEnergyGNP,
                                                excitationEnergyDTA);
        incidentKineticEnergy -= excitation;
        if (incidentKineticEnergy < excitationEnergyDTA) incidentKineticEnergy = 0.;
        incidentEnergy = incidentKineticEnergy + incidentMass;
        incidentTotalMomentum =
                 sqrt( max(0., incidentEnergy*incidentEnergy - incidentMass*incidentMass));
      }  
                  
   G4HEVector pTemp;
   for(i = 2; i<vecLen; i++)
     {
       j = min(vecLen-1, (G4int)(2. + G4UniformRand()*(vecLen-2)));
       pTemp = pv[j];
       pv[j] = pv[i];
       pv[i] = pTemp;
     }

                              // randomize the first two leading particles
                              // for kaon induced reactions only 
                              // (need from experimental data)

   if(   (incidentCode==kaonPlusCode  || incidentCode==kaonMinusCode    ||
          incidentCode==kaonZeroCode  || incidentCode==antiKaonZeroCode ||
          incidentCode==kaonZeroSCode || incidentCode==kaonZeroLCode)   
      && (G4UniformRand()>0.7) ) 
     {
       pTemp = pv[1];
       pv[1] = pv[0];
       pv[0] = pTemp;
     }
                              // mark leading particles for incident strange particles 
                              // and antibaryons, for all other we assume that the first 
                              // and second particle are the leading particles. 
                              // We need this later for kinematic aspects of strangeness
                              // conservation.
                          
   G4int lead = 0;                   
   G4HEVector leadParticle;
   if( (incidentMass >= kaonPlusMass-0.05) && (incidentCode != protonCode)          
                                           && (incidentCode != neutronCode) ) 
         {       
           G4double pMass = pv[0].getMass();
           G4int    pCode = pv[0].getCode();
           if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                            && (pCode != neutronCode) ) 
                  {       
                    lead = pCode; 
                    leadParticle = pv[0];                           
                  } 
           else   
                  {
                    pMass = pv[1].getMass();
                    pCode = pv[1].getCode();
                    if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                                     && (pCode != neutronCode) ) 
                        {       
                          lead = pCode;
                          leadParticle = pv[1];
                        }
                  }
         }

                              //  distribute particles in forward and backward 
                              //  hemispheres in center of mass system
                              //  incident goes in forward hemisphere
   
   G4HEVector pvI = incidentParticle;  // for the incident particle
   pvI.setSide( 1 );

   G4HEVector pvT = targetParticle;   // for the target particle
   pvT.setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvT.setSide( -1 );
   pvT.setTOF( -1.);  


   G4double centerOfMassEnergy = sqrt( sqr(incidentMass)+sqr(targetMass)
                                      +2.0*targetMass*incidentEnergy );
   G4double availableEnergy    = centerOfMassEnergy - ( targetMass + incidentMass );

   G4double tavai1      = centerOfMassEnergy/2.0 - incidentMass;
   G4double tavai2      = centerOfMassEnergy/2.0 - targetMass;           
   
   G4int ntb = 1;
   for( i=0; i < vecLen; i++ ) 
      {
        if      (i == 0) pv[i].setSide(  1 );
        else if (i == 1) pv[i].setSide( -1 );
        else  
           { if( G4UniformRand() < 0.5 ) 
               {
                 pv[i].setSide( -1 );
                 ntb++;
                } 
             else 
                 pv[i].setSide( 1 );
           }
        pv[i].setTOF(    incidentTOF);
      }
   G4double tb = 2. * ntb;
   if (centerOfMassEnergy < (2. + G4UniformRand())) 
       tb = (2. * ntb + vecLen)/2.;     

   if (verboseLevel > 1)
      { G4cout << " pv Vector after Randomization " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      } 

                              //  add particles from intranuclear cascade
                              //  nuclearCascadeCount = number of new secondaries 
                              //  produced by nuclear cascading.
                              //  extraCount = number of nucleons within these new secondaries
   
   G4double s, xtarg, ran;
   s = centerOfMassEnergy*centerOfMassEnergy;
   xtarg = max( 0.01, min( 0.75, 0.312 + 0.200 * log(log(s)) 
                                       + pow(s,1.5)/6000.0 ) 
                     *(pow(atomicWeight, 0.33) - 1.0) * tb);

   G4int ntarg = Poisson( xtarg );
   G4int targ = 0;
   
   if( ntarg > 0 ) 
     {
       G4double nucsup[] = { 1.00, 0.7, 0.5, 0.4, 0.35,   0.3 };
       G4double psup[]   = {   3.,  6., 20., 50., 100., 1000. };
       G4int momentumBin = 0;
       while( (momentumBin < 6) && (incidentTotalMomentum > psup[momentumBin]) )
             momentumBin++;
       momentumBin = min( 5, momentumBin );
     
                              // NOTE: in GENXPT, these new particles were given negative codes
                              //       here I use  flag = true  instead
     
       for( i=0; i<ntarg; i++ ) 
        {
          if( G4UniformRand() < nucsup[momentumBin] ) 
            {
              if( G4UniformRand() > 1.0-atomicNumber/atomicWeight ) 
                  pv[vecLen].setDefinition( "Proton" );
              else 
                  pv[vecLen].setDefinition( "Neutron" );
              targ++;
            } 
          else 
            {
              ran = G4UniformRand();
              if( ran < 0.33333 ) 
                  pv[vecLen].setDefinition( "PionPlus");
              else if( ran < 0.66667 ) 
                  pv[vecLen].setDefinition( "PionZero");
              else 
                  pv[vecLen].setDefinition( "PionMinus" );
            }
          pv[vecLen].setSide( -2 );                             // backward cascade particles
          pv[vecLen].setFlag( true );                           // true is the same as IPA(i)<0
          pv[vecLen].setTOF( incidentTOF );
          vecLen++; 
        }
     }
                                         
                                         //  assume conservation of kinetic energy 
                                         //  in forward & backward hemispheres
   G4int is, iskip;
   tavai1 = centerOfMassEnergy/2.;
   G4int iavai1 = 0;
 
   for (i = 0; i < vecLen; i++) 
       { 
         if (pv[i].getSide() > 0)
            { 
               tavai1 -= pv[i].getMass();
               iavai1++;
            }    
       } 
   if ( iavai1 == 0) return;

   while( tavai1 <= 0.0 ) 
        {                                // must eliminate a particle from the forward side
           iskip = G4int(G4UniformRand()*iavai1) + 1; 
           is = 0;  
           for( i=vecLen-1; i>=0; i-- ) 
              {
                if( pv[i].getSide() > 0 ) 
                  {
                    if (++is == iskip) 
                        {
                           tavai1 += pv[i].getMass();
                           iavai1--;            
                           if ( i != vecLen-1)
                              { 
                                 for( j=i; j<vecLen; j++ ) 
                                    {         
                                       pv[j] = pv[j+1];
                                    }
                              }
                           if( --vecLen == 0 ) return;     // all the secondaries except of the 
                           break;            // --+
                        }                    //   |
                  }                          //   v
              }                              // break goes down to here
        }                                    // to the end of the for- loop.
                                       

     tavai2 = (targ+1)*centerOfMassEnergy/2.;
     G4int iavai2 = 0;

     for (i = 0; i < vecLen; i++)
         {
           if (pv[i].getSide() < 0)
              {
                 tavai2 -= pv[i].getMass();
                 iavai2++;
              }
         }
     if (iavai2 == 0) return;

     while( tavai2 <= 0.0 ) 
        {                               // must eliminate a particle from the backward side
           iskip = G4int(G4UniformRand()*iavai2) + 1; 
           is = 0;
           for( i = vecLen-1; i >= 0; i-- ) 
              {
                if( pv[i].getSide() < 0 ) 
                  {
                    if( ++is == iskip ) 
                       {
                         tavai2 += pv[i].getMass();
                         iavai2--;
                         if (pv[i].getSide() == -2) ntarg--;
                         if (i != vecLen-1)
                            {
                              for( j=i; j<vecLen; j++)
                                 { 
                                   pv[j] = pv[j+1];
                                 } 
                            }
                         if (--vecLen == 0) return;
                         break;     
                       }
                  }   
              }
        }

   if (verboseLevel > 1)
      { G4cout << " pv Vector after Energy checks " << vecLen << " " << tavai1 << " " 
             << iavai1 << " " << tavai2 << " " << iavai2 << " " << ntarg << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      } 
   
                                 //  define some vectors for Lorentz transformations
   
   G4HEVector* pvmx = new G4HEVector [10];
   
   pvmx[0].setMass( incidentMass );
   pvmx[0].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   pvmx[1].setMass( protonMass);
   pvmx[1].setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvmx[3].setMass( protonMass*(1+targ));
   pvmx[3].setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvmx[4].setZero();
   pvmx[5].setZero();
   pvmx[7].setZero();
   pvmx[8].setZero();
   pvmx[8].setMomentum( 1.0, 0.0 );
   pvmx[2].Add( pvmx[0], pvmx[1] );
   pvmx[3].Add( pvmx[3], pvmx[0] );
   pvmx[0].Lor( pvmx[0], pvmx[2] );
   pvmx[1].Lor( pvmx[1], pvmx[2] );

   if (verboseLevel > 1)
     { G4cout << " General Vectors after Definition " << endl;
       for (i=0; i<10; i++) pvmx[i].Print(i);
     }
                                                 //  main loop for 4-momentum generation
                                                 //  see Pitha-report (Aachen) 
                                                 //  for a detailed description of the method
                                                 // process the secondary particles in reverse order
   G4double dndl[20];
   G4double binl[20];
   G4double pvMass, pvEnergy;
   G4int    pvCode; 
   G4double aspar, pt, phi, et, xval;
   G4double ekin  = 0.;
   G4double ekin1 = 0.;
   G4double ekin2 = 0.;
   G4int npg   = 0;
   G4int targ1 = 0;                           // no fragmentation model for nucleons from the
   for( i=vecLen-1; i>=0; i-- )               // intranuclear cascade. Mark them with -3 and 
      {                                       // leave the loop.
        if( (pv[i].getSide() == -2) || (i == 1) ) 
          { 
            if ( pv[i].getType() == baryonType)
               {                                 
                 if( ++npg < 19 ) 
                   {
                     pv[i].setSide( -3 );
                     targ1++;
                     continue;                                    // leave the for loop !!
                   }     
               }
          }
                                              //  set pt and phi values, 
                                              //  they are changed somewhat in the iteration loop
                                              //  set mass parameter for lambda fragmentation model
        G4double maspar[] = { 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.75, 0.20};
        G4double     bp[] = { 3.50, 3.50, 3.50, 6.00, 5.00, 4.00, 3.50, 3.50};
        G4double   ptex[] = { 1.70, 1.70, 1.50, 1.70, 1.40, 1.20, 1.70, 1.20};    

                                              // set parameters for lambda simulation 
                                              // pt is the average transverse momentum
                                              // aspar is average transverse mass  
        pvMass = pv[i].getMass();       
        j = 2;                                              
        if ( pv[i].getType() == mesonType ) j = 1;
        if ( pv[i].getMass() < 0.4 ) j = 0;
        if ( i <= 1 ) j += 3;
        if (pv[i].getSide() <= -2) j = 6;
        if (j == 6 && pv[i].getType() == baryonType) j = 7;
        pt    = max(0.001, sqrt(pow(-log(1.-G4UniformRand())/bp[j],ptex[j])));
        aspar = maspar[j]; 
        phi   = G4UniformRand()*M_2PI; 
        pv[i].setMomentum( pt*cos(phi), pt*sin(phi) );        // set x- and y-momentum

        for( j=0; j<20; j++ ) binl[j] = j/(19.*pt);           // set the lambda - bins.
     
        if( pv[i].getSide() > 0 )
           et = pvmx[0].getEnergy();
        else
           et = pvmx[1].getEnergy();
     
        dndl[0] = 0.0;
     
                                                              // start of outer iteration loop
        G4int outerCounter = 0, innerCounter = 0;             // three times.
        G4bool eliminateThisParticle = true;
        G4bool resetEnergies = true;
        while( ++outerCounter < 3 ) 
             {
               for( l=1; l<20; l++ ) 
                  {
                    xval  = (binl[l]+binl[l-1])/2.;              // x = lambda /GeV 
                    if( xval > 1./pt )
                       dndl[l] = dndl[l-1];
                    else
                       dndl[l] = dndl[l-1] + 
                                 aspar/sqrt( pow((1.+aspar*xval*aspar*xval),3) ) *
                                 (binl[l]-binl[l-1]) * et / 
                                 sqrt( pt*xval*et*pt*xval*et + pt*pt + pvMass*pvMass );
                  }  
       
                                                               // start of inner iteration loop
               innerCounter = 0;                               // try this not more than 7 times. 
               while( ++innerCounter < 7 ) 
                    {
                      l = 1;
                      ran = G4UniformRand()*dndl[19];
                      while( ( ran >= dndl[l] ) && ( l < 20 ) )l++;
                      l = min( 19, l );
                      xval = min( 1.0, pt*(binl[l-1] + G4UniformRand()*(binl[l]-binl[l-1])/2.) );
                      if( pv[i].getSide() < 0 ) xval *= -1.;
                      pv[i].setMomentumAndUpdate( xval*et );                  // set the z-momentum
                      pvEnergy = pv[i].getEnergy();
                      if( pv[i].getSide() > 0 )                               // forward side 
                        {
                          if ( i < 2 )
                             { 
                               ekin = tavai1 - ekin1;
                               if (ekin < 0.) ekin = 0.04*fabs(normal());
                               pv[i].setKineticEnergyAndUpdate( ekin);
                               pvmx[4].Add( pvmx[4], pv[i]);
                               continue;
                             }      
                          else if( (ekin1+pvEnergy-pvMass) < 0.95*tavai1 ) 
                            {
                              pvmx[4].Add( pvmx[4], pv[i] );
                              ekin1 += pvEnergy - pvMass;
                              pvmx[6].Add( pvmx[4], pvmx[5] );
                              pvmx[6].setMomentum( 0.0 );
                              phi = pvmx[6].Ang( pvmx[8] );
                              if( pvmx[6].getMomentum().y() < 0.0 ) phi = M_2PI - phi;
                              phi += M_PI + normal()*M_PI/12.0;
                              if( phi > M_2PI )phi -= M_2PI;
                              if( phi < 0.0 )phi = M_2PI - phi;
                              outerCounter = 2;                     // leave outer loop
                              eliminateThisParticle = false;        // don't eliminate this particle
                              resetEnergies = false;
                              break;                                // next particle
                            }
                          if( innerCounter > 5 ) break;             // leave inner loop
                          
                          if( tavai2 >= pvMass ) 
                            {                                       // switch sides
                              pv[i].setSide( -1 );
                              tavai1 += pvMass;
                              tavai2 -= pvMass;
                              iavai2++;
                            }
                        } 
                      else 
                        {                                            // backward side
                          xval = min(0.999,0.95+0.05*targ/20.0);
                          if( (ekin2+pvEnergy-pvMass) < xval*tavai2 ) 
                            {
                              pvmx[5].Add( pvmx[5], pv[i] );
                              ekin2 += pvEnergy - pvMass;
                              pvmx[6].Add( pvmx[4], pvmx[5] );
                              pvmx[6].setMomentum( 0.0 );             // set z-momentum
                              phi = pvmx[6].Ang( pvmx[8] );
                              if( pvmx[6].getMomentum().y() < 0.0 )phi = M_2PI - phi;
                              phi += M_PI + normal() * M_PI / 12.0;
                              if( phi > M_2PI )phi -= M_2PI;
                              if( phi < 0.0 )phi = M_2PI - phi;
                              outerCounter = 2;                    // leave outer iteration
                              eliminateThisParticle = false;       // don't eliminate this particle
                              resetEnergies = false;
                              break;                               // leave inner iteration
                            }
                          if( innerCounter > 5 )break;             // leave inner iteration
                          
                          if( tavai1 >= pvMass ) 
                            {                                      // switch sides
                              pv[i].setSide( 1 );
                              tavai1  -= pvMass;
                              tavai2  += pvMass;
                              iavai2--;
                            }
                        }
                      pv[i].setMomentum( pv[i].getMomentum().x() * 0.9,
                                         pv[i].getMomentum().y() * 0.9);
                      pt *= 0.9;
                      dndl[19] *= 0.9;
                    }                                           // closes inner loop

               if (resetEnergies)
                    {
                      ekin1 = 0.0;
                      ekin2 = 0.0;
                      pvmx[4].setZero();
                      pvmx[5].setZero();
                      if (verboseLevel > 1) G4cout << " Reset energies for index " << i << endl;
                      for( l=i+1; l < vecLen; l++ ) 
                         {
                           if( (pv[l].getMass() < protonMass) || (pv[l].getSide() > 0) ) 
                             {
                                pvEnergy = max( pv[l].getMass(),   0.95*pv[l].getEnergy() 
                                                                 + 0.05*pv[l].getMass() );
                                pv[l].setEnergyAndUpdate( pvEnergy );
                                if( pv[l].getSide() > 0) 
                                  {
                                    ekin1 += pv[l].getKineticEnergy();
                                    pvmx[4].Add( pvmx[4], pv[l] );
                                  } 
                                else 
                                  {
                                    ekin2 += pv[l].getKineticEnergy();
                                    pvmx[5].Add( pvmx[5], pv[l] );
                                  }
                             }
                         }
                    }
             }                                                           // closes outer iteration

        if( eliminateThisParticle )                                      // not enough energy, 
          {                                                              // eliminate this particle
            if (verboseLevel > 1)
               {
                  G4cout << " Eliminate particle with index " << i << endl;
                  pv[i].Print(i);
               }
            for( j=i; j < vecLen; j++ ) 
               {                                                         // shift down
                  pv[j] = pv[j+1];
               }
            vecLen--;
            if(vecLen < 2) return;
            i++;
            pvmx[6].Add( pvmx[4], pvmx[5] );
            pvmx[6].setMomentum( 0.0 );                                  // set z-momentum
            phi = pvmx[6].Ang( pvmx[8] );
            if( pvmx[6].getMomentum().y() < 0.0 )phi = M_2PI - phi;
            phi += M_PI + normal() * M_PI / 12.0;
            if( phi > M_2PI )phi -= M_2PI;
            if( phi < 0.0 )phi = M_2PI - phi;
          }
      }                                                                  // closes main for loop
   if (verboseLevel > 1)
      { G4cout << " pv Vector after lambda fragmentation " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
        for (i=0; i < 10; i++) pvmx[i].Print(i);
      } 

   
   
                                                                   // backward nucleons produced 
                                                                   // with a cluster model
   
   pvmx[6].Lor( pvmx[3], pvmx[2] );
   pvmx[6].Sub( pvmx[6], pvmx[4] );
   pvmx[6].Sub( pvmx[6], pvmx[5] );
   if (verboseLevel > 1) pvmx[6].Print(6);

   npg = 0;
   G4double rmb0 = 0.;
   G4double rmb;
   G4double wgt;
   G4bool constantCrossSection = true;
   for (i = 0; i < vecLen; i++)
     {
       if(pv[i].getSide() == -3) 
          { 
            npg++;
            rmb0 += pv[i].getMass();
          }
     }  
   if( targ1 == 1 || npg < 2) 
     {                                                  // target particle is the only backward nucleon
       ekin = min( tavai2-ekin2, centerOfMassEnergy/2.0-protonMass );
       if( ekin < 0.04 ) ekin = 0.04 * fabs( normal() );
       pv[1].setKineticEnergyAndUpdate( ekin );
       pvmx[5].Add( pvmx[5], pv[1] );
     } 
   else 
     {
       G4double cpar[] = { 0.6, 0.6, 0.35, 0.15, 0.10 };
       G4double gpar[] = { 2.6, 2.6, 1.80, 1.30, 1.20 };

       G4int tempCount = min( 5, targ1 ) - 1;

       rmb = rmb0 + pow(-log(1.0-G4UniformRand()), cpar[tempCount])/gpar[tempCount];
       pvEnergy = pvmx[6].getEnergy();
       if ( rmb > pvEnergy ) rmb = pvEnergy; 
       pvmx[6].setMass( rmb );
       pvmx[6].setEnergyAndUpdate( pvEnergy );
       pvmx[6].Smul( pvmx[6], -1. );
       if (verboseLevel > 1)
          { G4cout << " General Vectors before input to NBodyPhaseSpace "
                 << targ1 << " " << tempCount << " " << rmb0 << " " << rmb << " " 
                 << pvEnergy << endl;
            for (i=0; i<10; i++) pvmx[i].Print(i);
          }
                                                             //  tempV contains the backward nucleons

     
       G4HEVector* tempV = new G4HEVector[18];
       G4double totalMass = 0.0;
       npg = 0;
       for( i=0; i < vecLen; i++ )  
         {
           if( pv[i].getSide() == -3 ) tempV[npg++] = pv[i]; 
         }

       wgt = NBodyPhaseSpace( pvmx[6].getMass(), constantCrossSection, tempV, npg );

       npg = 0;
       for( i=0; i < vecLen; i++ ) 
         {
           if( pv[i].getSide() == -3 ) 
             {
               pv[i].setMomentum( tempV[npg++].getMomentum());
               pv[i].SmulAndUpdate( pv[i], 1.);
               pv[i].Lor( pv[i], pvmx[6] );
               pvmx[5].Add( pvmx[5], pv[i] );
             }
         }
     }
   if( vecLen <= 2 ) 
     {
       successful = false;
       return;
     }
                                                   //  Lorentz transformation in lab system
   targ = 0;
   for( i=0; i < vecLen; i++ ) 
      {
        if( pv[i].getType() == baryonType )targ++;
        pv[i].Lor( pv[i], pvmx[1] );
      }
   targ = max( 1, targ );

   G4bool dum;
   if( lead ) 
     {
       for( i=0; i<vecLen; i++ ) 
          {
            if( pv[i].getCode() == lead ) 
              {
                dum = false;
                break;
              }
          }
       if( dum ) 
         {
           i = 0;          
 
           if(   (    (leadParticle.getType() == baryonType)
                   && (pv[1].getType() == baryonType))
              || (    (leadParticle.getType() == mesonType)
                   && (pv[1].getType() == mesonType)))
             {
               i = 1;
             } 
            ekin = pv[i].getKineticEnergy();
            pv[i] = leadParticle;
            if( pv[i].getFlag() )
                pv[i].setTOF( -1.0 );
            else
                pv[i].setTOF( 1.0 );
            pv[i].setKineticEnergyAndUpdate( ekin );
         }
     }

   pvmx[3].setMass( incidentMass);
   pvmx[3].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   
   G4double ekin0 = pvmx[3].getKineticEnergy();
   
   pvmx[4].setMass ( protonMass * targ);
   pvmx[4].setEnergy( protonMass * targ);
   pvmx[4].setMomentum(0.,0.,0.);
   pvmx[4].setKineticEnergy(0.);

   ekin = pvmx[3].getEnergy() + pvmx[4].getEnergy();

   pvmx[5].Add( pvmx[3], pvmx[4] );
   pvmx[3].Lor( pvmx[3], pvmx[5] );
   pvmx[4].Lor( pvmx[4], pvmx[5] );
   
   G4double tecm = pvmx[3].getEnergy() + pvmx[4].getEnergy();

   pvmx[7].setZero();
   
   ekin1 = 0.0;   
   G4double teta; 
   
   for( i=0; i < vecLen; i++ ) 
      {
        pvmx[7].Add( pvmx[7], pv[i] );
        ekin1 += pv[i].getKineticEnergy();
        ekin  -= pv[i].getMass();
      }
   
   if( vecLen > 1 && vecLen < 19 ) 
     {
       constantCrossSection = true;
       G4HEVector pw[18];
       for(i=0;i<vecLen;i++) pw[i] = pv[i];
       wgt = NBodyPhaseSpace( tecm, constantCrossSection, pw, vecLen );
       ekin = 0.0;
       for( i=0; i < vecLen; i++ ) 
          {
            pvmx[6].setMass( pw[i].getMass());
            pvmx[6].setMomentum( pw[i].getMomentum() );
            pvmx[6].SmulAndUpdate( pvmx[6], 1.);
            pvmx[6].Lor( pvmx[6], pvmx[4] );
            ekin += pvmx[6].getKineticEnergy();
          }
       teta = pvmx[7].Ang( pvmx[3] );
       if (verboseLevel > 1) G4cout << " vecLen > 1 && vecLen < 19 " << teta << " " << ekin0
                                  << " " << ekin1 << " " << ekin << endl;
     }

   if( ekin1 != 0.0 ) 
     {
       pvmx[6].setZero();
       wgt = ekin/ekin1;
       ekin1 = 0.;
       for( i=0; i < vecLen; i++ ) 
          {
            pvMass = pv[i].getMass();
            ekin   = pv[i].getKineticEnergy() * wgt;
            pv[i].setKineticEnergyAndUpdate( ekin );
            ekin1 += ekin;
            pvmx[6].Add( pvmx[6], pv[i] );
          }
       teta = pvmx[6].Ang( pvmx[3] );
       if (verboseLevel > 1) G4cout << " ekin1 != 0 " << teta << " " << ekin0 << " " 
                                  << ekin1 << endl;
     }

                          //    Rotate in direction of z-axis, this does disturb in some way our
                          //    inclusive distributions, but it is necessary for momentum conservation

   pvmx[6].setZero();
   for( i=0; i < vecLen; i++ ) pvmx[6].Add( pvmx[6], pv[i] );
   
                          //  Some smearing in transverse direction from Fermi motion
   


                                                  // RANNOR
   G4double ry   = G4UniformRand();
   G4double rz   = G4UniformRand();
   G4double rx   = M_2PI*rz;
   G4double a1   = sqrt(-2.0*log(ry));
   G4double ran1 = a1*sin(rx);
   G4double ran2 = a1*cos(rx);
                                                  
   pvmx[6].setMomentum( pvmx[6].getMomentum().x()+ran1*0.020*targ,
                        pvmx[6].getMomentum().y()+ran2*0.020*targ );
   pvmx[7].Defs( pvmx[3], pvmx[6], pvmx[8], pvmx[9] );
   pvmx[6].setZero();
   for (i = 0; i < vecLen; i++) 
      { 
        pv[i].Trac( pv[i], pvmx[7] , pvmx[8], pvmx[9] );
        pvmx[6].Add( pvmx[6], pv[i] );  
      }
   teta = pvmx[6].Ang( pvmx[3] );   
   if (verboseLevel > 1) G4cout << " After smearing " << teta << endl;

                               //  Rotate in direction of primary particle, subtract binding energies
                               //   and make some further corrections if required

   G4double dekin = 0.0;
   G4int npions = 0;    
   G4double ek1 = 0.0;
   G4double alekw, xxh;
   G4double cfa = 0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.);
   G4double alem[] = {1.40, 2.30, 2.70, 3.00, 3.40, 4.60, 7.00, 10.00};
   G4double val0[] = {0.00, 0.40, 0.48, 0.51, 0.54, 0.60, 0.65,  0.70};
   

   for (i = 0; i < vecLen; i++)
      { 
        pv[i].Defs1( pv[i], pvI );
        if (atomicWeight > 1.5)
           {
             ekin  = max( 1.e-6,pv[i].getKineticEnergy() - cfa*( 1. + 0.5*normal()));
             alekw = log( incidentKineticEnergy );
             xxh   = 1.;
             if(incidentCode == pionPlusCode || incidentCode == pionMinusCode)
               {
                 if(pv[i].getCode() == pionZeroCode)
                   {
                     if(G4UniformRand() < log(atomicWeight))
                       {            
                         if (alekw > alem[0])
                           {
                             for (j = 1; j < 8; j++)
                                {
                                  if(alekw < alem[j]) break;
                                }
                             xxh =   (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alekw
                                    + val0[j-1] - (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alem[j-1];
                             xxh = 1. - xxh;
                           }
                       }
                    }
               }  
             dekin += ekin*(1.-xxh);
             ekin *= xxh;
             pv[i].setKineticEnergyAndUpdate( ekin );
             pvCode = pv[i].getCode();
             if ((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode))
                {
                  npions += 1;
                  ek1 += ekin; 
                }
           }
      }
   if( (ek1 > 0.0) && (npions > 0) ) 
      {
        dekin = 1.+dekin/ek1;
        for (i = 0; i < vecLen; i++)
          {
            pvCode = pv[i].getCode();
            if((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode)) 
              {
                ekin = max( 1.0e-6, pv[i].getKineticEnergy() * dekin );
                pv[i].setKineticEnergyAndUpdate( ekin );
              }
          }
      }
   if (verboseLevel > 1)
      { G4cout << " Lab-System " << ek1 << " " << npions << endl;
        for (i=0; i<vecLen; i++) pv[i].Print(i);
      }

                                       //  add black track particles
                                       //  the total number of particles produced is restricted to 198
                                       //  this may have influence on very high energies

   if (verboseLevel > 1) G4cout << " Evaporation " << atomicWeight << " " <<
                     excitationEnergyGNP << " " << excitationEnergyDTA << endl;

   if( atomicWeight > 1.5 ) 
     {

       G4double ekw, tex, sprob, cost, sint, pp, eka, etb;
       G4int spall, nbl;
                                     //  sprob is the probability of self-absorption in heavy molecules

       if( incidentKineticEnergy < 5.0 )
         sprob = 0.0;
       else
         sprob = min( 1.0, 0.6*log(incidentKineticEnergy-4.0) );

                                     //  first add protons and neutrons

       if( excitationEnergyGNP >= 0.001 ) 
         {
                                     //  nbl = number of proton/neutron black track particles
                                     //  tex is their total kinetic energy (GeV)
       
           nbl = Poisson( (1.5+1.25*targ)*excitationEnergyGNP/
                                         (excitationEnergyGNP+excitationEnergyDTA));
           if( targ+nbl > atomicWeight ) nbl = atomicWeight - targ;
           if (verboseLevel > 1) G4cout << " evaporation " << targ << " " << nbl << " " << sprob << endl; 
           spall = targ;
           if( nbl > 0) 
             {
               ekin = excitationEnergyGNP/nbl;
               ekin2 = 0.0;
               for( i=0; i<nbl; i++ ) 
                  {
                    if( G4UniformRand() < sprob ) continue;
                    if( ekin2 > excitationEnergyGNP) break;
                    ran = G4UniformRand();
                    ekin1 = -ekin*log(ran) - cfa*(1.0+0.5*normal());
                    if (ekin1 < 0) ekin1 = -0.010*log(ran);
                    ekin2 += ekin1;
                    if( ekin2 > excitationEnergyGNP)
                    ekin1 = max( 1.0e-6, excitationEnergyGNP-(ekin2-ekin1) );
                    if( G4UniformRand() > (1.0-atomicNumber/(atomicWeight)))
                       pv[vecLen].setDefinition( "Proton");
                    else
                       pv[vecLen].setDefinition( "Neutron");
                    spall++;
                    cost = G4UniformRand() * 2.0 - 1.0;
                    sint = sqrt(fabs(1.0-cost*cost));
                    phi = M_2PI * G4UniformRand();
                    pv[vecLen].setFlag( true );                 // true is the same as IPA(i)<0
                    pv[vecLen].setSide( -4 );
                    pvMass = pv[vecLen].getMass();
                    pv[vecLen].setTOF( 1.0 );
                    pvEnergy = ekin1 + pvMass;
                    pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                    pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                     pp*sint*cos(phi),
                                                     pp*cost );
                    if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                    vecLen++;
                  }
               if( (atomicWeight >= 10.0 ) && (incidentKineticEnergy <= 2.0) ) 
                  {
                    G4int ika, kk = 0;
                    eka = incidentKineticEnergy;
                    if( eka > 1.0 )eka *= eka;
                    eka = max( 0.1, eka );
                    ika = G4int(3.6*exp((atomicNumber*atomicNumber
                                /atomicWeight-35.56)/6.45)/eka);
                    if( ika > 0 ) 
                      {
                        for( i=(vecLen-1); i>=0; i-- ) 
                           {
                             if( (pv[i].getCode() == protonCode) && pv[i].getFlag() ) 
                               {
                                 pTemp = pv[i];
                                 pv[i].setDefinition( "Neutron");
                                 pv[i].setMomentumAndUpdate(pTemp.getMomentum());
                                 if (verboseLevel > 1) pv[i].Print(i);
                                 if( ++kk > ika ) break;
                               }
                           }
                      }
                  }
             }
         }

                                              // finished adding proton/neutron black track particles
                                              //  now, try to add deuterons, tritons and alphas
     
     if( excitationEnergyDTA >= 0.001 ) 
       {
         nbl = Poisson( (1.5+1.25*targ)*excitationEnergyDTA
                                      /(excitationEnergyGNP+excitationEnergyDTA));
       
                                      //    nbl is the number of deutrons, tritons, and alphas produced
       
         if( nbl > 0 ) 
           {
             ekin = excitationEnergyDTA/nbl;
             ekin2 = 0.0;
             for( i=0; i<nbl; i++ ) 
                {
                  if( G4UniformRand() < sprob ) continue;
                  if( ekin2 > excitationEnergyDTA) break;
                  ran = G4UniformRand();
                  ekin1 = -ekin*log(ran)-cfa*(1.+0.5*normal());
                  if( ekin1 < 0.0 ) ekin1 = -0.010*log(ran);
                  ekin2 += ekin1;
                  if( ekin2 > excitationEnergyDTA)
                     ekin1 = max( 1.0e-6, excitationEnergyDTA-(ekin2-ekin1));
                  cost = G4UniformRand()*2.0 - 1.0;
                  sint = sqrt(fabs(1.0-cost*cost));
                  phi = M_2PI*G4UniformRand();
                  ran = G4UniformRand();
                  if( ran <= 0.60 ) 
                      pv[vecLen].setDefinition( "Deuteron");
                  else if (ran <= 0.90)
                      pv[vecLen].setDefinition( "Triton");
                  else
                      pv[vecLen].setDefinition( "Alpha");
                  spall += pv[vecLen].getMass() * 1.066;
                  if( spall > atomicWeight ) break;
                  pv[vecLen].setFlag( true );                    // true is the same as IPA(i)<0
                  pv[vecLen].setSide( -4 );
                  pvMass = pv[vecLen].getMass();
                  pv[vecLen].setSide( pv[vecLen].getCode());
                  pv[vecLen].setTOF( 1.0 );
                  pvEnergy = pvMass + ekin1;
                  pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                  pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                   pp*sint*cos(phi),
                                                   pp*cost );
                  if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                  vecLen++;
                }
            }
        }
    }
   if( centerOfMassEnergy <= (4.0+G4UniformRand()) ) 
     {
       for( i=0; i<vecLen; i++ ) 
         {
           G4double etb = pv[i].getKineticEnergy();
           if( etb >= incidentKineticEnergy ) 
              pv[i].setKineticEnergyAndUpdate( incidentKineticEnergy );
         }
     }
                                                     //  calculate time delay for nuclear reactions

   G4double tof = incidentTOF;
   if(     (atomicWeight >= 1.5) && (atomicWeight <= 230.0) 
        && (incidentKineticEnergy <= 0.2) )
           tof -= 500.0 * exp(-incidentKineticEnergy /0.04) * log( G4UniformRand() );
   for ( i=0; i < vecLen; i++)     
     { 
       
       pv[i].setTOF ( tof );
//       vec[i].SetTOF ( tof );
     }
   successful = true;
   delete [] pvmx;
   return;
 }

void
G4HEInelastic::MediumEnergyClusterProduction(G4bool &successful,
                                             G4HEVector pv[],
                                             G4int &vecLen,	
                                             G4double &excitationEnergyGNP,
                                             G4double &excitationEnergyDTA,
                                             G4HEVector incidentParticle,
                                             G4HEVector targetParticle,
                                             G4double atomicWeight,
                                             G4double atomicNumber)
 {   
// For low multiplicity in the first intranuclear interaction the cascading process
// as described in G4HEInelastic::MediumEnergyCascading does not work 
// satisfactorily. From experimental data it is strongly suggested to use 
// a two- body resonance model.   
//  
//  All quantities on the G4HEVector Array pv are in GeV- units.

   G4int protonCode       = Proton.getCode();
   G4double protonMass    = Proton.getMass();
   G4int neutronCode      = Neutron.getCode();
   G4double neutronMass   = Neutron.getMass();
   G4double kaonPlusMass  = KaonPlus.getMass();
   G4int kaonPlusCode     = KaonPlus.getCode();   
   G4int kaonMinusCode    = KaonMinus.getCode();
   G4int kaonZeroSCode    = KaonZeroShort.getCode(); 
   G4int kaonZeroLCode    = KaonZeroLong.getCode();
   G4int kaonZeroCode     = KaonZero.getCode();
   G4int antiKaonZeroCode = AntiKaonZero.getCode(); 
   G4int pionPlusCode     = PionPlus.getCode();    
   G4int pionZeroCode     = PionZero.getCode();    
   G4int pionMinusCode    = PionMinus.getCode(); 
   G4String mesonType     = PionPlus.getType();
   G4String baryonType    = Proton.getType(); 
   
   G4int    targetCode   = targetParticle.getCode();
   G4double targetMass   = targetParticle.getMass();
   G4double targetCharge = targetParticle.getCharge();

   G4int    incidentCode          = incidentParticle.getCode();
   G4double incidentMass          = incidentParticle.getMass();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   G4double incidentCharge        = incidentParticle.getCharge();
   G4double incidentEnergy        = incidentParticle.getEnergy();
   G4double incidentKineticEnergy = incidentParticle.getKineticEnergy();
   G4String incidentType          = incidentParticle.getType();
//   G4double incidentTOF           = incidentParticle.getTOF();   
   G4double incidentTOF           = 0.;
   
                                                              // some local variables

   G4int i, j, k, l;
   
   if(verboseLevel > 1) G4cout << " G4HEInelastic::MediumEnergyClusterProduction " << endl;

   if (incidentTotalMomentum < 0.01)
      {
        successful = false;
        return;
      }
   G4double centerOfMassEnergy = sqrt( sqr(incidentMass) + sqr(targetMass)
                                       +2.*targetMass*incidentEnergy);

   G4HEVector pvI = incidentParticle;  // for the incident particle
   pvI.setSide( 1 );

   G4HEVector pvT = targetParticle;   // for the target particle
   pvT.setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvT.setSide( -1 );
   pvT.setTOF( -1.); 
                              // distribute particles in forward and backward
                              // hemispheres. Note that only low multiplicity
                              // events from FirstIntInNuc.... should go into 
                              // this routine. 
   G4int targ  = 0;  
   G4int ifor  = 0; 
   G4int iback = 0;
   G4int pvCode;
   G4double pvMass, pvEnergy; 

   pv[0].setSide(  1 );
   pv[1].setSide( -1 );
   for(i = 0; i < vecLen; i++)
      {
        if (i > 1)
           {
              if( G4UniformRand() < 0.5) 
                {
                  pv[i].setSide( 1 );
                  if (++ifor > 18) 
                     { 
                       pv[i].setSide( -1 );
                       ifor--;
                       iback++;
                     }
                }
              else
                {
                  pv[i].setSide( -1 );
                  if (++iback > 18)
                     { 
                       pv[i].setSide( 1 );
                       ifor++;
                       iback--;
                     }
                }
           }

        pvCode = pv[i].getCode();

        if (   (    (incidentCode == protonCode) || (incidentCode == neutronCode)
                 || (incidentType == mesonType) )
            && (    (pvCode == pionPlusCode) || (pvCode == pionMinusCode) )
            && (    (G4UniformRand() < (10.-incidentTotalMomentum)/6.) )
            && (    (G4UniformRand() < atomicWeight/300.) ) )
           { 
               if (G4UniformRand() > atomicNumber/atomicWeight) 
                  pv[i].setDefinition( "Neutron");
               else
                  pv[i].setDefinition( "Proton");
               targ++;
           }    
        pv[i].setTOF( incidentTOF );                 
      }    
   G4double tb = 2. * iback;
   if (centerOfMassEnergy < (2+G4UniformRand())) tb = (2.*iback + vecLen)/2.;
   
   G4double nucsup[] = { 1.0, 0.8, 0.6, 0.5, 0.4}; 

   G4double xtarg = max(0.01,  (0.312+0.2*log(log(centerOfMassEnergy*centerOfMassEnergy)))
                             * (pow(atomicWeight,0.33)-1.) * tb);
   G4int ntarg = Poisson(xtarg);
   if (ntarg > 0)
      {
        G4int ipx = min(4, (G4int)(incidentTotalMomentum/3.));
        for (i=0; i < ntarg; i++)
            { 
               if (G4UniformRand() < nucsup[ipx] )
                  {
                     if (G4UniformRand() < (1.- atomicNumber/atomicWeight))
                        pv[vecLen].setDefinition( "Neutron");
                     else
                        pv[vecLen].setDefinition( "Proton");
                     targ++;
                  }   
               else
                  {
                     G4double ran = G4UniformRand();
                     if (ran < 0.3333 ) 
                        pv[vecLen].setDefinition( "PionPlus");
                     else if (ran < 0.6666)
                        pv[vecLen].setDefinition( "PionZero");
                     else
                        pv[vecLen].setDefinition( "PionMinus");
                  }
               pv[vecLen].setSide( -2 );        
               pv[vecLen].setFlag( true );
               pv[vecLen].setTOF( incidentTOF );
               vecLen++;
            }
      }
                              // mark leading particles for incident strange particles 
                              // and antibaryons, for all other we assume that the first 
                              // and second particle are the leading particles. 
                              // We need this later for kinematic aspects of strangeness
                              // conservation.
                          
   G4int lead = 0;                   
   G4HEVector leadParticle;
   if( (incidentMass >= kaonPlusMass-0.05) && (incidentCode != protonCode)          
                                           && (incidentCode != neutronCode) ) 
         {       
           G4double pMass = pv[0].getMass();
           G4int    pCode = pv[0].getCode();
           if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                            && (pCode != neutronCode) ) 
                  {       
                    lead = pCode; 
                    leadParticle = pv[0];                           
                  } 
           else   
                  {
                    pMass = pv[1].getMass();
                    pCode = pv[1].getCode();
                    if( (pMass >= kaonPlusMass-0.05) && (pCode != protonCode) 
                                                     && (pCode != neutronCode) ) 
                        {       
                          lead = pCode;
                          leadParticle = pv[1];
                        }
                  }
         }

   if (verboseLevel > 1)
      { G4cout << " pv Vector after initialization " << vecLen << endl;
        pvI.Print(-1);
        pvT.Print(-1);
        for (i=0; i < vecLen ; i++) pv[i].Print(i);
      }     

   G4double tavai = 0.;
   for(i=0;i<vecLen;i++) if(pv[i].getSide() != -2) tavai  += pv[i].getMass();

   while (tavai > centerOfMassEnergy)
      {
         for (i=vecLen-1; i >= 0; i--)
            {
              if (pv[i].getSide() != -2)
                 {
                    tavai -= pv[i].getMass();
                    if( i != vecLen-1) 
                      {
                        for (j=i; j < vecLen; j++)
                           {
                              pv[j]  = pv[j+1];
                           }
                      }
                    if ( --vecLen < 2)
                      {
                        successful = false;
                        return;
                      }
                    break;
                 }
            }    
      }
                                          // now we produce 3 Clusters:
                                          // 1. forward cluster
                                          // 2. backward meson cluster
                                          // 3. backward nucleon cluster
   G4double rmc0 = 0., rmd0 = 0., rme0 = 0.;
   G4int    ntc  = 0,  ntd  = 0,  nte  = 0;
   
   for (i=0; i < vecLen; i++)
      {
        if(pv[i].getSide() > 0)
           {
             if(ntc < 17)
               { 
                 rmc0 += pv[i].getMass();
                 ntc++;
               }
             else
               {
                 if(ntd < 17)
                   {
                     pv[i].setSide(-1);
                     rmd0 += pv[i].getMass();
                     ntd++;
                   }
                 else
                   {
                     pv[i].setSide(-2);
                     rme0 += pv[i].getMass();
                     nte++;
                   }
               }
           }  
        else if (pv[i].getSide() == -1)
           {
             if(ntd < 17)
                {
                   rmd0 += pv[i].getMass();
                   ntd++;
                }
             else
                {
                   pv[i].setSide(-2); 
                   rme0 += pv[i].getMass();
                   nte++;
                }           
           }
        else
	   {
             rme0 += pv[i].getMass();
             nte++;
           } 
      }         

   G4double cpar[] = {0.6, 0.6, 0.35, 0.15, 0.10};
   G4double gpar[] = {2.6, 2.6, 1.80, 1.30, 1.20};
    
   G4double rmc = rmc0, rmd = rmd0, rme = rme0; 
   G4int ntc1 = min(4,ntc-1);
   G4int ntd1 = min(4,ntd-1);
   G4int nte1 = min(4,nte-1);
   if (ntc > 1) rmc = rmc0 + pow(-log(1.-G4UniformRand()),cpar[ntc1])/gpar[ntc1];
   if (ntd > 1) rmd = rmd0 + pow(-log(1.-G4UniformRand()),cpar[ntd1])/gpar[ntd1];
   if (nte > 1) rme = rme0 + pow(-log(1.-G4UniformRand()),cpar[nte1])/gpar[nte1];
   while( (rmc+rmd) > centerOfMassEnergy)
      {
        if ((rmc == rmc0) && (rmd == rmd0))
          { 
            rmd *= 0.999*centerOfMassEnergy/(rmc+rmd);
            rmc *= 0.999*centerOfMassEnergy/(rmc+rmd);
          }
        else
          {
            rmc = 0.1*rmc0 + 0.9*rmc;
            rmd = 0.1*rmd0 + 0.9*rmd;
          }   
      }             
   if(verboseLevel > 1) G4cout << " Cluster Masses: " << ntc << " " << rmc << " " << ntd << " "
                             << rmd << " " << nte << " " << rme << endl;
 
   
   G4HEVector* pvmx = new G4HEVector[11];

   pvmx[1].setMass( incidentMass);
   pvmx[1].setMomentumAndUpdate( 0., 0., incidentTotalMomentum);
   pvmx[2].setMass( targetMass);
   pvmx[2].setMomentumAndUpdate( 0., 0., 0.);
   pvmx[0].Add( pvmx[1], pvmx[2] );
   pvmx[1].Lor( pvmx[1], pvmx[0] );
   pvmx[2].Lor( pvmx[2], pvmx[0] );

   G4double pf = sqrt(max(0.0001,  sqr(sqr(centerOfMassEnergy) + rmd*rmd -rmc*rmc)
                                 - 4*sqr(centerOfMassEnergy)*rmd*rmd))/(2.*centerOfMassEnergy);
   pvmx[3].setMass( rmc );
   pvmx[4].setMass( rmd );
   pvmx[3].setEnergy( sqrt(pf*pf + rmc*rmc) );
   pvmx[4].setEnergy( sqrt(pf*pf + rmd*rmd) );
   
   G4double tvalue = -MAXFLOAT;
   G4double bvalue = max(0.01, 4.0 + 1.6*log(incidentTotalMomentum));
   if (bvalue != 0.0) tvalue = log(G4UniformRand())/bvalue;
   G4double pin = pvmx[1].Length();
   G4double tacmin = sqr( pvmx[1].getEnergy() - pvmx[3].getEnergy()) - sqr( pin - pf);
   G4double ctet   = max(-1., min(1., 1.+2.*(tvalue-tacmin)/max(1.e-10, 4.*pin*pf)));
   G4double stet   = sqrt(max(0., 1.0 - ctet*ctet));
   G4double phi    = M_2PI * G4UniformRand();
   pvmx[3].setMomentum( pf * stet * sin(phi), 
                      pf * stet * cos(phi),
                      pf * ctet           );
   pvmx[4].Smul( pvmx[3], -1.);
   
   if (nte > 0)
      {
        G4double ekit1 = 0.04;
        G4double ekit2 = 0.6;
        G4double gaval = 1.2;
        if (incidentKineticEnergy <= 5.)
           {
             ekit1 *= sqr(incidentKineticEnergy)/25.;
             ekit2 *= sqr(incidentKineticEnergy)/25.;
           }
        G4double avalue = (1.-gaval)/(pow(ekit2, 1.-gaval)-pow(ekit1, 1.-gaval));
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -2)
                 {
                   G4double ekit = pow(G4UniformRand()*(1.-gaval)/avalue +pow(ekit1, 1.-gaval),
                                       1./(1.-gaval));
                   pv[i].setKineticEnergyAndUpdate( ekit );
                   ctet = max(-1., min(1., log(2.23*G4UniformRand()+0.383)/0.96));
                   stet = sqrt( max( 0.0, 1. - ctet*ctet ));
                   phi  = G4UniformRand()*M_2PI;
                   G4double pp = pv[i].Length();
                   pv[i].setMomentum( pp * stet * sin(phi),
                                      pp * stet * cos(phi),
                                      pp * ctet           );
                   pv[i].Lor( pv[i], pvmx[0] );
                 }              
            }             
      }
   pvmx[1] = pvmx[3];
   pvmx[2] = pvmx[4];
   pvmx[5].SmulAndUpdate( pvmx[3], -1.);
   pvmx[6].SmulAndUpdate( pvmx[4], -1.);

   if (verboseLevel > 1)
      { G4cout << " General vectors before Phase space Generation " << endl;
        for (i=0; i<7; i++) pvmx[i].Print(i);
      }  


   G4HEVector* tempV = new G4HEVector[18];
   G4bool constantCrossSection = true;
   G4double wgt;
   G4int npg;

   if (ntc > 1)
      {
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() > 0)
                 {
                    tempV[npg++] = pv[i];
                    if(verboseLevel > 1) pv[i].Print(i);
                 }
            }
        wgt = NBodyPhaseSpace( pvmx[3].getMass(), constantCrossSection, tempV, npg);
                     
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() > 0)
                 {
                   pv[i].setMomentum( tempV[npg++].getMomentum());
                   pv[i].SmulAndUpdate( pv[i], 1. );
                   pv[i].Lor( pv[i], pvmx[5] );
                   if(verboseLevel > 1) pv[i].Print(i);
                 }
            }
      }
   else if(ntc == 1)
      {
        for(i=0; i<vecLen; i++)
          {
            if(pv[i].getSide() > 0) pv[i].setMomentumAndUpdate(pvmx[3].getMomentum());
            if(verboseLevel > 1) pv[i].Print(i); 
          }
      }
   else
      {
      }
     
   if (ntd > 1)
      {
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -1)
                 {
                    tempV[npg++] = pv[i];
                    if(verboseLevel > 1) pv[i].Print(i);
                 }
            }
        wgt = NBodyPhaseSpace( pvmx[4].getMass(), constantCrossSection, tempV, npg);
                     
        npg = 0;
        for (i=0; i < vecLen; i++)
            {
              if (pv[i].getSide() == -1)
                 {
                   pv[i].setMomentum( tempV[npg++].getMomentum());
                   pv[i].SmulAndUpdate( pv[i], 1.);
                   pv[i].Lor( pv[i], pvmx[6] );
                   if(verboseLevel > 1) pv[i].Print(i);
                 }
            }
      }
   else if(ntd == 1)
      {
        for(i=0; i<vecLen; i++)
          {
            if(pv[i].getSide() == -1) pv[i].setMomentumAndUpdate(pvmx[4].getMomentum());
            if(verboseLevel > 1) pv[i].Print(i);
          }
      }
   else
      {
      }

   if(verboseLevel > 1)
     {
       G4cout << " Vectors after PhaseSpace generation " << endl;
       for(i=0;i<vecLen; i++) pv[i].Print(i);
     } 

                                                   //  Lorentz transformation in lab system
   targ = 0;
   for( i=0; i < vecLen; i++ ) 
      {
        if( pv[i].getType() == baryonType )targ++;
        pv[i].Lor( pv[i], pvmx[2] );
      }

   targ = max( 1, targ );

   if(verboseLevel > 1)
     {
       G4cout << " Transformation in Lab- System " << endl;
       for(i=0; i<vecLen; i++) pv[i].Print(i);
     }

   G4bool dum;
   G4double ekin, teta;

   if( lead ) 
     {
       for( i=0; i<vecLen; i++ ) 
          {
            if( pv[i].getCode() == lead ) 
              {
                dum = false;
                break;
              }
          }
       if( dum ) 
         {
           i = 0;          
 
           if(   (    (leadParticle.getType() == baryonType)
                   && (pv[1].getType() == baryonType))
              || (    (leadParticle.getType() == mesonType)
                   && (pv[1].getType() == mesonType)))
             {
               i = 1;
             } 

            ekin = pv[i].getKineticEnergy();
            pv[i] = leadParticle;
            if( pv[i].getFlag() )
                pv[i].setTOF( -1.0 );
            else
                pv[i].setTOF( 1.0 );
            pv[i].setKineticEnergyAndUpdate( ekin );
         }
     }

   pvmx[4].setMass( incidentMass);
   pvmx[4].setMomentumAndUpdate( 0.0, 0.0, incidentTotalMomentum );
   
   G4double ekin0 = pvmx[4].getKineticEnergy();
   
   pvmx[5].setMass ( protonMass * targ);
   pvmx[5].setMomentumAndUpdate( 0.0, 0.0, 0.0 );

   ekin = pvmx[4].getEnergy() + pvmx[5].getEnergy();

   pvmx[6].Add( pvmx[4], pvmx[5] );
   pvmx[4].Lor( pvmx[4], pvmx[6] );
   pvmx[5].Lor( pvmx[5], pvmx[6] );
   
   G4double tecm = pvmx[4].getEnergy() + pvmx[5].getEnergy();

   pvmx[8].setZero();
   
   G4double ekin1 = 0.0;   
   
   for( i=0; i < vecLen; i++ ) 
      {
        pvmx[8].Add( pvmx[8], pv[i] );
        ekin1 += pv[i].getKineticEnergy();
        ekin  -= pv[i].getMass();
      }
   
   if( vecLen > 1 && vecLen < 19 ) 
     {
       constantCrossSection = true;
       G4HEVector pw[18];
       for(i=0;i<vecLen;i++) pw[i] = pv[i];
       wgt = NBodyPhaseSpace( tecm, constantCrossSection, pw, vecLen );
       ekin = 0.0;
       for( i=0; i < vecLen; i++ ) 
          {
            pvmx[7].setMass( pw[i].getMass());
            pvmx[7].setMomentum( pw[i].getMomentum() );
            pvmx[7].SmulAndUpdate( pvmx[7], 1.);
            pvmx[7].Lor( pvmx[7], pvmx[5] );
            ekin += pvmx[7].getKineticEnergy();
          }
       teta = pvmx[8].Ang( pvmx[4] );
       if (verboseLevel > 1) G4cout << " vecLen > 1 && vecLen < 19 " << teta << " " << ekin0 
                                  << " " << ekin1 << " " << ekin << endl;
     }

   if( ekin1 != 0.0 ) 
     {
       pvmx[7].setZero();
       wgt = ekin/ekin1;
       ekin1 = 0.;
       for( i=0; i < vecLen; i++ ) 
          {
            pvMass = pv[i].getMass();
            ekin   = pv[i].getKineticEnergy() * wgt;
            pv[i].setKineticEnergyAndUpdate( ekin );
            ekin1 += ekin;
            pvmx[7].Add( pvmx[7], pv[i] );
          }
       teta = pvmx[7].Ang( pvmx[4] );
       if (verboseLevel > 1) G4cout << " ekin1 != 0 " << teta << " " << ekin0 << " " 
                                  << ekin1 << endl;
     }

                          //    Rotate in direction of z-axis, this does disturb in some way our
                          //    inclusive distributions, but it is necessary for momentum conservation

   pvmx[7].setZero();
   for( i=0; i < vecLen; i++ ) pvmx[7].Add( pvmx[7], pv[i] );
   
                          //  Some smearing in transverse direction from Fermi motion
   


                                                  // RANNOR
   G4double ry   = G4UniformRand();
   G4double rz   = G4UniformRand();
   G4double rx   = M_2PI*rz;
   G4double a1   = sqrt(-2.0*log(ry));
   G4double ran1 = a1*sin(rx);
   G4double ran2 = a1*cos(rx);
                                                  
   pvmx[7].setMomentum( pvmx[7].getMomentum().x()+ran1*0.020*targ,
                        pvmx[7].getMomentum().y()+ran2*0.020*targ );
   pvmx[8].Defs( pvmx[4], pvmx[7], pvmx[9], pvmx[10] );
   pvmx[7].setZero();
   for (i = 0; i < vecLen; i++) 
      { 
        pv[i].Trac( pv[i], pvmx[8] , pvmx[9], pvmx[10] );
        pvmx[7].Add( pvmx[7], pv[i] );  
      }
   teta = pvmx[7].Ang( pvmx[4] );   
   if (verboseLevel > 1) G4cout << " After smearing " << teta << endl;

                               //  Rotate in direction of primary particle, subtract binding energies
                               //   and make some further corrections if required

   G4double dekin = 0.0;
   G4int npions = 0;    
   G4double ek1 = 0.0;
   G4double alekw, xxh;
   G4double cfa = 0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.);
   G4double alem[] = {1.40, 2.30, 2.70, 3.00, 3.40, 4.60, 7.00};
   G4double val0[] = {0.00, 0.40, 0.48, 0.51, 0.54, 0.60, 0.65};
   

   for (i = 0; i < vecLen; i++)
      { 
        pv[i].Defs1( pv[i], pvI );
        if (atomicWeight > 1.5)
           {
             ekin  = max( 1.e-6,pv[i].getKineticEnergy() - cfa*( 1. + 0.5*normal()));
             alekw = log( incidentKineticEnergy );
             xxh   = 1.;
             xxh   = 1.;
             if(incidentCode == pionPlusCode || incidentCode == pionMinusCode)
	       {
                 if(pv[i].getCode() == pionZeroCode)
                   {
                     if(G4UniformRand() < log(atomicWeight))
                       { 
                         if (alekw > alem[0])
                           {
                              for (j = 1; j < 8; j++)
                                 {
                                    if(alekw < alem[j]) break;
                                 } 
                              xxh =   (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alekw
                                     + val0[j-1] - (val0[j]-val0[j-1])/(alem[j]-alem[j-1])*alem[j-1];
                              xxh = 1. - xxh;
                           }
                       }      
                    }
               }  
             dekin += ekin*(1.-xxh);
             ekin *= xxh;
             pv[i].setKineticEnergyAndUpdate( ekin );
             pvCode = pv[i].getCode();
             if ((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode))
                {
                  npions += 1;
                  ek1 += ekin; 
                }
           }
      }
   if( (ek1 > 0.0) && (npions > 0) ) 
      {
        dekin = 1.+dekin/ek1;
        for (i = 0; i < vecLen; i++)
          {
            pvCode = pv[i].getCode();
            if((pvCode == pionPlusCode) || (pvCode == pionMinusCode) || (pvCode == pionZeroCode)) 
              {
                ekin = max( 1.0e-6, pv[i].getKineticEnergy() * dekin );
                pv[i].setKineticEnergyAndUpdate( ekin );
              }
          }
      }
   if (verboseLevel > 1)
      { G4cout << " Lab-System " << ek1 << " " << npions << endl;
        for (i=0; i<vecLen; i++) pv[i].Print(i);
      }

                                       //  add black track particles
                                       //  the total number of particles produced is restricted to 198
                                       //  this may have influence on very high energies

   if (verboseLevel > 1) G4cout << " Evaporation " <<  atomicWeight << " " <<
                     excitationEnergyGNP << " " << excitationEnergyDTA << endl;

   if( atomicWeight > 1.5 ) 
     {

       G4double ekw, tex, sprob, cost, sint, ekin2, ran, pp, eka;
       G4int spall, nbl;
                                     //  sprob is the probability of self-absorption in heavy molecules

       if( incidentKineticEnergy < 5.0 )
         sprob = 0.0;
       else
         sprob = min( 1.0, 0.6*log(incidentKineticEnergy-4.0) );

                                     //  first add protons and neutrons

       if( excitationEnergyGNP >= 0.001 ) 
         {
                                     //  nbl = number of proton/neutron black track particles
                                     //  tex is their total kinetic energy (GeV)
       
           nbl = Poisson( (1.5+1.25*targ)*excitationEnergyGNP/
                                         (excitationEnergyGNP+excitationEnergyDTA));
           if( targ+nbl > atomicWeight ) nbl = atomicWeight - targ;
           if (verboseLevel > 1) G4cout << " evaporation " << targ << " " << nbl << " " 
                                      << sprob << endl; 
           spall = targ;
           if( nbl > 0) 
             {
               ekin = excitationEnergyGNP/nbl;
               ekin2 = 0.0;
               for( i=0; i<nbl; i++ ) 
                  {
                    if( G4UniformRand() < sprob ) continue;
                    if( ekin2 > excitationEnergyGNP) break;
                    ran = G4UniformRand();
                    ekin1 = -ekin*log(ran) - cfa*(1.0+0.5*normal());
                    if (ekin1 < 0) ekin1 = -0.010*log(ran);
                    ekin2 += ekin1;
                    if( ekin2 > excitationEnergyGNP )
                    ekin1 = max( 1.0e-6, excitationEnergyGNP-(ekin2-ekin1) );
                    if( G4UniformRand() > (1.0-atomicNumber/(atomicWeight)))
                       pv[vecLen].setDefinition( "Proton");
                    else
                       pv[vecLen].setDefinition( "Neutron");
                    spall++;
                    cost = G4UniformRand() * 2.0 - 1.0;
                    sint = sqrt(fabs(1.0-cost*cost));
                    phi = M_2PI * G4UniformRand();
                    pv[vecLen].setFlag( true );                 // true is the same as IPA(i)<0
                    pv[vecLen].setSide( -4 );
                    pvMass = pv[vecLen].getMass();
                    pv[vecLen].setTOF( 1.0 );
                    pvEnergy = ekin1 + pvMass;
                    pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                    pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                     pp*sint*cos(phi),
                                                     pp*cost );
                    if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                    vecLen++;
                  }
               if( (atomicWeight >= 10.0 ) && (incidentKineticEnergy <= 2.0) ) 
                  {
                    G4int ika, kk = 0;
                    eka = incidentKineticEnergy;
                    if( eka > 1.0 )eka *= eka;
                    eka = max( 0.1, eka );
                    ika = G4int(3.6*exp((atomicNumber*atomicNumber
                                /atomicWeight-35.56)/6.45)/eka);
                    if( ika > 0 ) 
                      {
                        for( i=(vecLen-1); i>=0; i-- ) 
                           {
                             if( (pv[i].getCode() == protonCode) && pv[i].getFlag() ) 
                               {
                                 G4HEVector pTemp = pv[i];
                                 pv[i].setDefinition( "Neutron");
                                 pv[i].setMomentumAndUpdate(pTemp.getMomentum());
                                 if (verboseLevel > 1) pv[i].Print(i);
                                 if( ++kk > ika ) break;
                               }
                           }
                      }
                  }
             }
         }

                                              // finished adding proton/neutron black track particles
                                              //  now, try to add deuterons, tritons and alphas
     
     if( excitationEnergyDTA >= 0.001 ) 
       {
         nbl = Poisson( (1.5+1.25*targ)*excitationEnergyDTA
                                      /(excitationEnergyGNP+excitationEnergyDTA));
       
                                      //    nbl is the number of deutrons, tritons, and alphas produced
       
         if( nbl > 0 ) 
           {
             ekin = excitationEnergyDTA/nbl;
             ekin2 = 0.0;
             for( i=0; i<nbl; i++ ) 
                {
                  if( G4UniformRand() < sprob ) continue;
                  if( ekin2 > excitationEnergyDTA) break;
                  ran = G4UniformRand();
                  ekin1 = -ekin*log(ran)-cfa*(1.+0.5*normal());
                  if( ekin1 < 0.0 ) ekin1 = -0.010*log(ran);
                  ekin2 += ekin1;
                  if( ekin2 > excitationEnergyDTA)
                     ekin1 = max( 1.0e-6, excitationEnergyDTA-(ekin2-ekin1));
                  cost = G4UniformRand()*2.0 - 1.0;
                  sint = sqrt(fabs(1.0-cost*cost));
                  phi = M_2PI*G4UniformRand();
                  ran = G4UniformRand();
                  if( ran <= 0.60 ) 
                      pv[vecLen].setDefinition( "Deuteron");
                  else if (ran <= 0.90)
                      pv[vecLen].setDefinition( "Triton");
                  else
                      pv[vecLen].setDefinition( "Alpha");
                  spall += pv[vecLen].getMass() * 1.066;
                  if( spall > atomicWeight ) break;
                  pv[vecLen].setFlag( true );                    // true is the same as IPA(i)<0
                  pv[vecLen].setSide( -4 );
                  pvMass = pv[vecLen].getMass();
                  pv[vecLen].setTOF( 1.0 );
                  pvEnergy = pvMass + ekin1;
                  pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                  pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                   pp*sint*cos(phi),
                                                   pp*cost );
                  if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                  vecLen++;
                }
            }
        }
    }
   if( centerOfMassEnergy <= (4.0+G4UniformRand()) ) 
     {
       for( i=0; i<vecLen; i++ ) 
         {
           G4double etb = pv[i].getKineticEnergy();
           if( etb >= incidentKineticEnergy ) 
              pv[i].setKineticEnergyAndUpdate( incidentKineticEnergy );
         }
     }
                                                     //  calculate time delay for nuclear reactions

   G4double tof = incidentTOF;
   if(     (atomicWeight >= 1.5) && (atomicWeight <= 230.0) 
        && (incidentKineticEnergy <= 0.2) )
           tof -= 500.0 * exp(-incidentKineticEnergy /0.04) * log( G4UniformRand() );
   for ( i=0; i < vecLen; i++)     
     { 
       
       pv[i].setTOF ( tof );
//       vec[i].SetTOF ( tof );
     }
   successful = true;
   delete [] pvmx;
   delete [] tempV;
   return;
 }

void
G4HEInelastic::QuasiElasticScattering(G4bool &successful,
                                            G4HEVector pv[],
                                            G4int &vecLen,	
                                            G4double &excitationEnergyGNP,
                                            G4double &excitationEnergyDTA,
                                            G4HEVector incidentParticle,
                                            G4HEVector targetParticle,
                                            G4double atomicWeight,
                                            G4double atomicNumber )
 {   
// if the Cascading or Resonance - model fails, we try this,
// QuasiElasticScattering. 
//    
//  All quantities on the G4HEVector Array pv are in GeV- units.

   G4int protonCode       = Proton.getCode();
   G4double protonMass    = Proton.getMass();
   G4int neutronCode      = Neutron.getCode();
   G4double neutronMass   = Neutron.getMass();
   G4double kaonPlusMass  = KaonPlus.getMass();
   G4int kaonPlusCode     = KaonPlus.getCode();   
   G4int kaonMinusCode    = KaonMinus.getCode();
   G4int kaonZeroSCode    = KaonZeroShort.getCode(); 
   G4int kaonZeroLCode    = KaonZeroLong.getCode();
   G4int kaonZeroCode     = KaonZero.getCode();
   G4int antiKaonZeroCode = AntiKaonZero.getCode(); 
   G4int pionPlusCode     = PionPlus.getCode();    
   G4int pionZeroCode     = PionZero.getCode();    
   G4int pionMinusCode    = PionMinus.getCode(); 
   G4String mesonType     = PionPlus.getType();
   G4String baryonType    = Proton.getType(); 
   
   G4int    targetCode   = targetParticle.getCode();
   G4double targetMass   = targetParticle.getMass();
   G4double targetCharge = targetParticle.getCharge();

   G4int    incidentCode          = incidentParticle.getCode();
   G4double incidentMass          = incidentParticle.getMass();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   G4double incidentCharge        = incidentParticle.getCharge();
   G4double incidentEnergy        = incidentParticle.getEnergy();
   G4double incidentKineticEnergy = incidentParticle.getKineticEnergy();
   G4String incidentType          = incidentParticle.getType();
//   G4double incidentTOF           = incidentParticle.getTOF();   
   G4double incidentTOF           = 0.;
   
                                                              // some local variables

   G4int i;
   
   if(verboseLevel > 1) G4cout << " G4HEInelastic::QuasiElasticScattering " << endl;

   if (incidentTotalMomentum < 0.01 || vecLen < 2 )
      {
        successful = false;
        return;
      }
   G4double centerOfMassEnergy = sqrt( sqr(incidentMass) + sqr(targetMass)
                                       +2.*targetMass*incidentEnergy);

   G4HEVector pvI = incidentParticle;  // for the incident particle
   pvI.setSide( 1 );

   G4HEVector pvT = targetParticle;   // for the target particle
   pvT.setMomentumAndUpdate( 0.0, 0.0, 0.0 );
   pvT.setSide( -1 );
   pvT.setTOF( -1.); 

   G4HEVector* pvmx = new G4HEVector[3];

   if (atomicWeight > 1.5)                   // for the following case better use ElasticScattering.
      {
        if (    (pvI.getCode() == pv[0].getCode() )
             && (pvT.getCode() == pv[1].getCode() )
             && (excitationEnergyGNP < 0.001)
             && (excitationEnergyDTA < 0.001)                            )
           {
             successful = false;
             return;
           }
      }

   pv[0].setSide( 1 );
   pv[0].setFlag( false );
   pv[0].setTOF( incidentTOF);
   pv[0].setMomentumAndUpdate( incidentParticle.getMomentum() );
   pv[1].setSide( -1 );
   pv[1].setFlag( false );
   pv[1].setTOF( incidentTOF);
   pv[1].setMomentumAndUpdate(targetParticle.getMomentum() );

   if ( (incidentTotalMomentum > 0.1) && (centerOfMassEnergy > 0.01) )
      {
        if ( pv[1].getType() == mesonType )
           {
             if (G4UniformRand() < 0.5)
                 pv[1].setDefinition( "Proton"); 
             else
                 pv[1].setDefinition( "Neutron");
           }
        pvmx[0].Add( pvI, pvT );
        pvmx[1].Lor( pvI, pvmx[0] );
        pvmx[2].Lor( pvT, pvmx[0] );
        G4double pin = pvmx[1].Length();
        G4double bvalue = max(0.01 , 4.225+1.795*log(incidentTotalMomentum));
        G4double pf =   sqr( sqr(centerOfMassEnergy) + sqr(pv[1].getMass()) - sqr(pv[0].getMass()))
                      - 4 * sqr(centerOfMassEnergy) * sqr(pv[1].getMass());
        if ( pf < 0.001)
           {
             successful = false;
             return;
           }
        pf = sqrt(pf)/(2.*centerOfMassEnergy);
        G4double btrang = 4. * bvalue * pin * pf;
        G4double exindt = -1.;
        if (btrang < 46.) exindt += exp(-btrang);
        G4double tdn = log(1. + G4UniformRand()*exindt)/btrang;
        G4double ctet = max( -1., min(1., 1. + 2.*tdn));
        G4double stet = sqrt((1.-ctet)*(1.+ctet));
        G4double phi  = M_2PI * G4UniformRand();
        pv[0].setMomentumAndUpdate( pf*stet*sin(phi),
                                    pf*stet*cos(phi),
                                    pf*ctet         );
        pv[1].SmulAndUpdate( pv[0], -1.); 
        
        for (i = 0; i < 2; i++)
           {
             pv[i].Lor( pv[i], pvmx[4] );
             pv[i].Defs1( pv[i], pvI );
             if (atomicWeight > 1.5)
               {
                 G4double ekin =    pv[i].getKineticEnergy() 
                                 -  0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.)
                                   *(1. + 0.5*normal()); 
                 ekin = max(0.0001, ekin);
                 pv[i].setKineticEnergyAndUpdate( ekin );
               }
           }
      }
   vecLen = 2;
                                       //  add black track particles
                                       //  the total number of particles produced is restricted to 198
                                       //  this may have influence on very high energies

   if (verboseLevel > 1) G4cout << " Evaporation " << atomicWeight << " " << 
                     excitationEnergyGNP << " " <<  excitationEnergyDTA << endl;

   if( atomicWeight > 1.5 ) 
     {

       G4double ekw, tex, sprob, cost, sint, ekin2, ran, pp, eka;
       G4double ekin, cfa, ekin1, phi, pvMass, pvEnergy;
       G4int spall, nbl;
                                     //  sprob is the probability of self-absorption in heavy molecules

       sprob = 0.;
       cfa   = 0.025*((atomicWeight-1.)/120.)*exp(-(atomicWeight-1.)/120.);
                                     //  first add protons and neutrons

       if( excitationEnergyGNP >= 0.001 ) 
         {
                                     //  nbl = number of proton/neutron black track particles
                                     //  tex is their total kinetic energy (GeV)
       
           nbl = Poisson( excitationEnergyGNP/0.02);
           if( nbl > atomicWeight ) nbl = atomicWeight;
           if (verboseLevel > 1) G4cout << " evaporation " << nbl << " " << sprob << endl; 
           spall = 0;
           if( nbl > 0) 
             {
               ekin = excitationEnergyGNP/nbl;
               ekin2 = 0.0;
               for( i=0; i<nbl; i++ ) 
                  {
                    if( G4UniformRand() < sprob ) continue;
                    if( ekin2 > excitationEnergyGNP) break;
                    ran = G4UniformRand();
                    ekin1 = -ekin*log(ran) - cfa*(1.0+0.5*normal());
                    if (ekin1 < 0) ekin1 = -0.010*log(ran);
                    ekin2 += ekin1;
                    if( ekin2 > excitationEnergyGNP)
                    ekin1 = max( 1.0e-6, excitationEnergyGNP-(ekin2-ekin1) );
                    if( G4UniformRand() > (1.0-atomicNumber/(atomicWeight)))
                       pv[vecLen].setDefinition( "Proton");
                    else
                       pv[vecLen].setDefinition( "Neutron");
                    spall++;
                    cost = G4UniformRand() * 2.0 - 1.0;
                    sint = sqrt(fabs(1.0-cost*cost));
                    phi = M_2PI * G4UniformRand();
                    pv[vecLen].setFlag( true );                 // true is the same as IPA(i)<0
                    pv[vecLen].setSide( -4 );
                    pvMass = pv[vecLen].getMass();
                    pv[vecLen].setTOF( 1.0 );
                    pvEnergy = ekin1 + pvMass;
                    pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                    pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                     pp*sint*cos(phi),
                                                     pp*cost );
                    if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                    vecLen++;
                  }
               if( (atomicWeight >= 10.0 ) && (incidentKineticEnergy <= 2.0) ) 
                  {
                    G4int ika, kk = 0;
                    eka = incidentKineticEnergy;
                    if( eka > 1.0 )eka *= eka;
                    eka = max( 0.1, eka );
                    ika = G4int(3.6*exp((atomicNumber*atomicNumber
                                /atomicWeight-35.56)/6.45)/eka);
                    if( ika > 0 ) 
                      {
                        for( i=(vecLen-1); i>=0; i-- ) 
                           {
                             if( (pv[i].getCode() == protonCode) && pv[i].getFlag() ) 
                               {
                                 pv[i].setDefinition( "Neutron" );
                                 if (verboseLevel > 1) pv[i].Print(i);
                                 if( ++kk > ika ) break;
                               }
                           }
                      }
                  }
             }
         }

                                              // finished adding proton/neutron black track particles
                                              //  now, try to add deuterons, tritons and alphas
     
     if( excitationEnergyDTA >= 0.001 ) 
       {
         nbl = (G4int)(2.*log(atomicWeight));
       
                                      //    nbl is the number of deutrons, tritons, and alphas produced
       
         if( nbl > 0 ) 
           {
             ekin = excitationEnergyDTA/nbl;
             ekin2 = 0.0;
             for( i=0; i<nbl; i++ ) 
                {
                  if( G4UniformRand() < sprob ) continue;
                  if( ekin2 > excitationEnergyDTA) break;
                  ran = G4UniformRand();
                  ekin1 = -ekin*log(ran)-cfa*(1.+0.5*normal());
                  if( ekin1 < 0.0 ) ekin1 = -0.010*log(ran);
                  ekin2 += ekin1;
                  if( ekin2 > excitationEnergyDTA)
                     ekin1 = max( 1.0e-6, excitationEnergyDTA-(ekin2-ekin1));
                  cost = G4UniformRand()*2.0 - 1.0;
                  sint = sqrt(fabs(1.0-cost*cost));
                  phi = M_2PI*G4UniformRand();
                  ran = G4UniformRand();
                  if( ran <= 0.60 ) 
                      pv[vecLen].setDefinition( "Deuteron");
                  else if (ran <= 0.90)
                      pv[vecLen].setDefinition( "Triton");
                  else
                      pv[vecLen].setDefinition( "Alpha");
                  spall += pv[vecLen].getMass() * 1.066;
                  if( spall > atomicWeight ) break;
                  pv[vecLen].setFlag( true );                    // true is the same as IPA(i)<0
                  pv[vecLen].setSide( -4 );
                  pvMass = pv[vecLen].getMass();
                  pv[vecLen].setTOF( 1.0 );
                  pvEnergy = pvMass + ekin1;
                  pp = sqrt( fabs( pvEnergy*pvEnergy - pvMass*pvMass ) );
                  pv[vecLen].setMomentumAndUpdate( pp*sint*sin(phi),
                                                   pp*sint*cos(phi),
                                                   pp*cost );
                  if (verboseLevel > 1) pv[vecLen].Print(vecLen);
                  vecLen++;
                }
            }
        }
    }
                                                     //  calculate time delay for nuclear reactions

   G4double tof = incidentTOF;
   if(     (atomicWeight >= 1.5) && (atomicWeight <= 230.0) 
        && (incidentKineticEnergy <= 0.2) )
           tof -= 500.0 * exp(-incidentKineticEnergy /0.04) * log( G4UniformRand() );
   for ( i=0; i < vecLen; i++)     
     { 
       
       pv[i].setTOF ( tof );
//       vec[i].SetTOF ( tof );
     }
   successful = true;
   delete [] pvmx;
   return;
 }

void
G4HEInelastic::ElasticScattering(G4bool &successful,
                                 G4HEVector pv[],
                                 G4int &vecLen,  
                                 G4HEVector incidentParticle,
                                 G4double atomicWeight,
                                 G4double atomicNumber)
 {
   if(verboseLevel > 1) G4cout << " G4HEInelastic::ElasticScattering " << endl;

   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
   if (verboseLevel > 1)
      G4cout << "DoIt: Incident particle momentum=" << incidentTotalMomentum << " GeV" << endl;
   if (incidentTotalMomentum < 0.01) 
      { 
        successful = false;
        return;
      }
   if (atomicWeight < 0.5) 
      { 
        successful = false;
        return;
      }
   pv[0] = incidentParticle;
   vecLen = 1;

   G4double aa, bb, cc, dd, rr;
   if (atomicWeight <= 62.) 
     {
       aa = pow(atomicWeight, 1.63);
       bb = 14.5*pow(atomicWeight, 0.66);
       cc = 1.4*pow(atomicWeight, 0.33);
       dd = 10.;
     }
   else 
     {
       aa = pow(atomicWeight, 1.33);
       bb = 60.*pow(atomicWeight, 0.33);
       cc = 0.4*pow(atomicWeight, 0.40);
       dd = 10.;
     }
   aa = aa/bb;
   cc = cc/dd;
   G4double ran = G4UniformRand();
   rr = (aa + cc)*ran;
   if (verboseLevel > 1) 
     {
       G4cout << "ElasticScattering: aa,bb,cc,dd,rr" << endl;
       G4cout << aa << " " << bb << " " << cc << " " << dd << " " << rr << endl;
     }
   G4double t1 = -log(ran)/bb;
   G4double t2 = -log(ran)/dd;
   if (verboseLevel > 1) 
     {
       G4cout << "t1,fctcos " << t1 << " " << fctcos(t1, aa, bb, cc, dd, rr) << 
               endl;
       G4cout << "t2,fctcos " << t2 << " " << fctcos(t2, aa, bb, cc, dd, rr) << 
               endl;
     }
   G4double eps = 0.001;
   G4int ind1 = 10;
   G4double t, val;
   G4int ier1;
   ier1 = rtmi(&t, t1, t2, eps, ind1, aa, bb, cc, dd, rr);
   if (verboseLevel > 1) 
     {
       G4cout << "From rtmi, ier1=" << ier1 << endl;
       G4cout << "t, fctcos " << t << " " << fctcos(t, aa, bb, cc, dd, rr) << 
               endl;
     }
   if (ier1 != 0) t = 0.25*(3.*t1 + t2);
   if (verboseLevel > 1) 
     {
       G4cout << "t, fctcos " << t << " " << fctcos(t, aa, bb, cc, dd, rr) << 
               endl;
     }
   G4double phi = G4UniformRand()*2.*M_PI;
   rr = 0.5*t/sqr(incidentTotalMomentum);
   if (rr > 1.) rr = 0.;
   if (verboseLevel > 1)
      G4cout << "rr=" << rr << endl;
   G4double cost = 1. - rr;
   G4double sint = sqrt(max(rr*(2. - rr), 0.));
   if (verboseLevel > 1)
      G4cout << "cos(t)=" << cost << "  sin(t)=" << sint << endl;
                                         // Scattered particle referred to axis of incident particle
   G4HEVector pv0;
   G4HEVector pvI;
   pvI.setMass( incidentParticle.getMass() );
   pvI.setMomentum( incidentParticle.getMomentum() );
   pvI.SmulAndUpdate( pvI, 1. );    
   pv0.setMass( pvI.getMass() );
   
   pv0.setMomentumAndUpdate( incidentTotalMomentum * sint * sin(phi),
                             incidentTotalMomentum * sint * cos(phi),
                             incidentTotalMomentum * cost           );    
   pv0.Defs1( pv0, pvI );
      
   successful = true;
   return;
 }


G4int
G4HEInelastic::rtmi(G4double *x, G4double xli, G4double xri, G4double eps, 
                          G4int iend, 
                          G4double aa, G4double bb, G4double cc, G4double dd, 
                          G4double rr)
 {
   G4int ier = 0;
   G4double xl = xli;
   G4double xr = xri;
   *x = xl;
   G4double tol = *x;
   G4double f = fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;
   G4double fl, fr;
   fl = f;
   *x = xr;
   tol = *x;
   f = fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;
   fr = f;

                                     // Error return in case of wrong input data
   if (fl*fr >= 0.) 
     {
       ier = 2;
       return ier;
     }

                                     // Basic assumption fl*fr less than 0 is satisfied.
                                     // Generate tolerance for function values.
   G4int i = 0;
   G4double tolf = 100.*eps;

                                     // Start iteration loop
   label4:   // <-------------
   i++;

                                     // Start bisection loop
   for (G4int k = 1; k <= iend; k++) 
     {
       *x = 0.5*(xl + xr);
       tol = *x;
       f = fctcos(tol, aa, bb, cc, dd, rr);
       if (f == 0.) return 0;
       if (f*fr < 0.) 
         {                          // Interchange xl and xr in order to get the
           tol = xl;                // same sign in f and fr
           xl = xr;
           xr = tol;
           tol = fl;
           fl = fr;
           fr = tol;
         }
       tol = f - fl;
       G4double a = f*tol;
       a = a + a;
       if (a < fr*(fr - fl) && i <= iend) goto label17;
       xr = *x;
       fr = f;

                                    // Test on satisfactory accuracy in bisection loop
       tol = eps;
       a = fabs(xr);
       if (a > 1.) tol = tol*a;
       if (fabs(xr - xl) <= tol && fabs(fr - fl) <= tolf) goto label14;
     }
                                    // End of bisection loop

                                    // No convergence after iend iteration steps followed by iend
                                    // successive steps of bisection or steadily increasing function
                                    // values at right bounds.  Error return.
   ier = 1;

   label14:  // <---------------
   if (fabs(fr) > fabs(fl)) 
     {
       *x = xl;
       f = fl;
     }
   return ier;

                                     // Computation of iterated x-value by inverse parabolic interp
   label17:  // <---------------
   G4double a = fr - f;
   G4double dx = (*x - xl)*fl*(1. + f*(a - tol)/(a*(fr - fl)))/tol;
   G4double xm = *x;
   G4double fm = f;
   *x = xl - dx;
   tol = *x;
   f = fctcos(tol, aa, bb, cc, dd, rr);
   if (f == 0.) return ier;

                                      // Test on satisfactory accuracy in iteration loop
   tol = eps;
   a = fabs(*x);
   if (a > 1) tol = tol*a;
   if (fabs(dx) <= tol && fabs(f) <= tolf) return ier;

                                       // Preparation of next bisection loop
   if (f*fl < 0.) 
     {
       xr = *x;
       fr = f;
     }
   else 
     {
       xl = *x;
       fl = f;
       xr = xm;
       fr = fm;
     }
   goto label4;
 }


// Test function for root-finder

G4double
G4HEInelastic::fctcos(G4double t, 
                            G4double aa, G4double bb, G4double cc, G4double dd, 
                            G4double rr)
 {
   const G4double expxl = -82.;
   const G4double expxu = 82.;

   G4double test1 = -bb*t;
   if (test1 > expxu) test1 = expxu;
   if (test1 < expxl) test1 = expxl;

   G4double test2 = -dd*t;
   if (test2 > expxu) test2 = expxu;
   if (test2 < expxl) test2 = expxl;

   return aa*exp(test1) + cc*exp(test2) - rr;
 }

 G4double G4HEInelastic::NBodyPhaseSpace
                               ( const G4double totalEnergy,                // MeV
                                 const G4bool constantCrossSection,
                                 G4HEVector  vec[],
                                 G4int& vecLen )
  {
    // derived from original FORTRAN code PHASP by H. Fesefeldt (02-Dec-1986)
    // Returns the weight of the event

    G4int i;
    
    const G4double expxu =  log(FLT_MAX);  // upper bound for arg. of exp
    const G4double expxl = -expxu;         // lower bound for arg. of exp
    
    if( vecLen < 2 ) {
      G4cerr << "*** Error in G4HEInelastic::GenerateNBodyEvent" << endl;
      G4cerr << "    number of particles < 2" << endl;
      G4cerr << "totalEnergy = " << totalEnergy << ", vecLen = " << vecLen << endl;
      return -1.0;
    }
    
    G4double* mass = new G4double [vecLen];    // mass of each particle
    G4double* energy = new G4double [vecLen];  // total energy of each particle
    G4double** pcm;           // pcm is an array with 3 rows and vecLen columns
    pcm = new G4double* [3];
    for( i=0; i<3; ++i )pcm[i] = new G4double [vecLen];
    
    G4double totalMass = 0.0;
    G4double* sm = new G4double [vecLen];
    
    for( i=0; i<vecLen; ++i ) {
      mass[i] = vec[i].getMass();
      vec[i].setMomentum( 0.0, 0.0, 0.0 );
      pcm[0][i] = 0.0;      // x-momentum of i-th particle
      pcm[1][i] = 0.0;      // y-momentum of i-th particle
      pcm[2][i] = 0.0;      // z-momentum of i-th particle
      energy[i] = mass[i];  // total energy of i-th particle
      totalMass += mass[i];
      sm[i] = totalMass;
    }
    if( totalMass >= totalEnergy ) {
      G4cerr << "*** Error in G4HEInelastic::GenerateNBodyEvent" << endl;
      G4cerr << "    total mass (" << totalMass << ") >= total energy ("
           << totalEnergy << ")" << endl;
      delete [] mass;
      delete [] energy;
      for( i=0; i<3; ++i )delete [] pcm[i];
      delete [] pcm;
      delete [] sm;
      return -1.0;
    }
    G4double kineticEnergy = totalEnergy - totalMass;
    G4double* emm = new G4double [vecLen];
    emm[0] = mass[0];
    if( vecLen > 3 ) {          // the random numbers are sorted
      G4double* ran = new G4double [vecLen];
      for( i=0; i<vecLen; ++i )ran[i] = G4UniformRand();
      for( i=0; i<vecLen-1; ++i ) {
        for( G4int j=vecLen-1; j > i; --j ) {
          if( ran[i] > ran[j] ) {
            G4double temp = ran[i];
            ran[i] = ran[j];
            ran[j] = temp;
          }
        }
      }
      for( i=1; i<vecLen; ++i )emm[i] = ran[i-1]*kineticEnergy + sm[i];
      delete [] ran;
    } else {
      emm[1] = G4UniformRand()*kineticEnergy + sm[1];
    }
    emm[vecLen-1] = totalEnergy;
    
    //   Weight is the sum of logarithms of terms instead of the product of terms
    
    G4bool lzero = true;    
    G4double wtmax = 0.0;
    if( constantCrossSection ) {     // this is KGENEV=1 in PHASP
      G4double emmax = kineticEnergy + mass[0];
      G4double emmin = 0.0;
      for( i=1; i<vecLen; ++i ) {
        emmin += mass[i-1];
        emmax += mass[i];
        G4double wtfc = 0.0;
        if( emmax*emmax > 0.0 ) {
          G4double arg = emmax*emmax
            + (emmin*emmin-mass[i]*mass[i])*(emmin*emmin-mass[i]*mass[i])/(emmax*emmax)
            - 2.0*(emmin*emmin+mass[i]*mass[i]);
          if( arg > 0.0 )wtfc = 0.5*sqrt( arg );
        }
        if( wtfc == 0.0 ) {
          lzero = false;
          break;
        }
        wtmax += log( wtfc );
      }
      if( lzero )
        wtmax = -wtmax;
      else
        wtmax = expxu;
    } else {
      wtmax = log( pow( kineticEnergy, vecLen-2 ) *
                   pi * pow( twopi, vecLen-2 ) / Factorial(vecLen-2) );
    }
    lzero = true;
    G4double* pd = new G4double [vecLen-1];
    for( i=0; i<vecLen-1; ++i ) {
      pd[i] = 0.0;
      if( emm[i+1]*emm[i+1] > 0.0 ) {
        G4double arg = emm[i+1]*emm[i+1]
          + (emm[i]*emm[i]-mass[i+1]*mass[i+1])*(emm[i]*emm[i]-mass[i+1]*mass[i+1])
            /(emm[i+1]*emm[i+1])
          - 2.0*(emm[i]*emm[i]+mass[i+1]*mass[i+1]);
        if( arg > 0.0 )pd[i] = 0.5*sqrt( arg );
      }
      if( pd[i] == 0.0 )
        lzero = false;
      else
        wtmax += log( pd[i] );
    }
    G4double weight = 0.0;           // weight is returned by GenerateNBodyEvent
    if( lzero )weight = exp( max(min(wtmax,expxu),expxl) );
    
    G4double bang, cb, sb, s0, s1, s2, c, s, esys, a, b, gama, beta;
    pcm[0][0] = 0.0;
    pcm[1][0] = pd[0];
    pcm[2][0] = 0.0;
    for( i=1; i<vecLen; ++i ) {
      pcm[0][i] = 0.0;
      pcm[1][i] = -pd[i-1];
      pcm[2][i] = 0.0;
      bang = twopi*G4UniformRand();
      cb = cos(bang);
      sb = sin(bang);
      c = 2.0*G4UniformRand() - 1.0;
      s = sqrt( fabs( 1.0-c*c ) );
      if( i < vecLen-1 ) {
        esys = sqrt(pd[i]*pd[i] + emm[i]*emm[i]);
        beta = pd[i]/esys;
        gama = esys/emm[i];
        for( G4int j=0; j<=i; ++j ) {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*s;                           //  rotation
          pcm[1][j] = s0*s + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
          pcm[1][j] = gama*(pcm[1][j] + beta*energy[j]);
        }
      } else {
        for( G4int j=0; j<=i; ++j ) {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*s;                           //  rotation
          pcm[1][j] = s0*s + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
        }
      }
    }
    G4double pModule; 
    for( i=0; i<vecLen; ++i ) {
      kineticEnergy = energy[i] - mass[i];
      pModule = sqrt( sqr(kineticEnergy) + 2*kineticEnergy*mass[i] );    
      vec[i].setMomentum( pcm[0][i]/pModule, 
                          pcm[1][i]/pModule, 
                          pcm[2][i]/pModule );
      vec[i].setKineticEnergyAndUpdate( kineticEnergy );
    }
    delete [] mass;
    delete [] energy;
    for( i=0; i<3; ++i )delete [] pcm[i];
    delete [] pcm;
    delete [] emm;
    delete [] sm;
    delete [] pd;
    return weight;
  }
 
G4double
G4HEInelastic::gpdk( G4double a, G4double b, G4double c )
 {
   if( a == 0.0 ) 
     {
       return 0.0;
      } 
   else 
      {
        G4double arg = a*a+(b*b-c*c)*(b*b-c*c)/(a*a)-2.0*(b*b+c*c);
        if( arg <= 0.0 ) 
          {
            return 0.0;
          } 
        else 
          {
            return 0.5*sqrt(fabs(arg));
          }
      }
 }


G4double
G4HEInelastic::NBodyPhaseSpace(G4int npart, G4HEVector pv[], 
                                     G4double wmax, G4double wfcn, 
                                     G4int maxtrial, G4int ntrial)
 { ntrial = 0;
   G4double wps;
   while ( ntrial < maxtrial)
     { ntrial++;
       G4int i, j;
       G4int nrn = 3*(npart-2)-4;
       G4double *ranarr = new G4double[nrn];
       for (i=0;i<nrn;i++) ranarr[i]=G4UniformRand();
       G4int nrnp = npart-4;
       if(nrnp > 1) QuickSort( ranarr, 0 , nrnp-1 );
       G4HEVector pvcms;
       pvcms.Add(pv[0],pv[1]);
       pvcms.Smul( pvcms, -1.);
       G4double rm = 0.;
       for (i=2;i<npart;i++) rm += pv[i].getMass();
       G4double rm1 = pvcms.getMass() - rm;
       rm -= pv[2].getMass();
       wps = (npart-3)*pow(rm1/sqr(M_2PI), npart-4)/(4*M_PI*pvcms.getMass());
       for (i=3; i=npart-1;i++) wps /= i-2;
       G4double xxx = rm1/sqr(M_2PI);
       for (i=1; i=npart-4; i++) wps /= xxx/i;
       wps /= (4*M_PI*pvcms.getMass());
       G4double p2,cost,sint,phi;
       j = 1;
       while (j)
         { j++;
           rm -= pv[j+1].getMass();
           if(j == npart-2) break;
           G4double rmass = rm + rm1*ranarr[npart-j-1];
           p2 = Alam(sqr(pvcms.getMass()), sqr(pv[j].getMass()),
                     sqr(rmass))/(4.*sqr(pvcms.getMass()));
           cost = 1. - 2.*ranarr[npart+2*j-9];
           sint = sqrt(1.-cost*cost);
           phi  = M_2PI*ranarr[npart+2*j-8];
           p2   = sqrt( max(0., p2));
           wps *= p2;
           pv[j].setMomentumAndUpdate( p2*sint*sin(phi), p2*sint*cos(phi),p2*cost);
           pv[j].Lor(pv[j], pvcms);
           pvcms.Add3( pvcms, pv[j] );
           pvcms.setEnergy(pvcms.getEnergy()-pv[j].getEnergy());
           pvcms.setMass( sqrt(sqr(pvcms.getEnergy()) - sqr(pvcms.Length())));
         }        
       p2 = Alam(sqr(pvcms.getMass()), sqr(pv[j].getMass()),
                 sqr(rm))/(4.*sqr(pvcms.getMass()));
       cost = 1. - 2.*ranarr[npart+2*j-9];
       sint = sqrt(1.-cost*cost);
       phi  = M_2PI*ranarr[npart+2*j-8];
       p2   = sqrt( max(0. , p2));
       wps *= p2;
       pv[j].setMomentumAndUpdate( p2*sint*sin(phi), p2*sint*cos(phi), p2*cost);
       pv[j+1].setMomentumAndUpdate( -p2*sint*sin(phi), -p2*sint*cos(phi), -p2*cost);
       pv[j].Lor( pv[j], pvcms );
       pv[j+1].Lor( pv[j+1], pvcms );
       wfcn = CalculatePhaseSpaceWeight( npart );
       G4double wt = wps * wfcn;
       if (wt > wmax)
         { wmax = wt;
           G4cout << "maximum weight changed to " << wmax << endl;
         }
       wt = wt/wmax;
       if (G4UniformRand() < wt) break; 
     }
   return wps;
 }
      

void 
G4HEInelastic::QuickSort(G4double arr[], const G4int lidx, const G4int ridx)
 {                                   // sorts the Array arr[] in ascending order
   G4double buffer;
   G4int k, e, mid;
   if(lidx>=ridx) return;
   mid = (lidx+ridx)/2.;
   buffer   = arr[lidx];
   arr[lidx]= arr[mid];
   arr[mid] = buffer;
   e = lidx;
   for (k=lidx+1;k<=ridx;k++)
     if (arr[k] < arr[lidx])
       { e++;
         buffer = arr[e];
         arr[e] = arr[k];
         arr[k] = buffer;
       }
   buffer   = arr[lidx];
   arr[lidx]= arr[e];
   arr[e]   = buffer;
   QuickSort(arr, lidx, e-1);
   QuickSort(arr, e+1 , ridx);
   return;
 }

G4double
G4HEInelastic::Alam( G4double a, G4double b, G4double c)
 { return a*a + b*b + c*c - 2.*a*b - 2.*a*c -2.*b*c;
 }    

G4double 
G4HEInelastic::CalculatePhaseSpaceWeight( G4int npart ) 
 { G4double wfcn = 1.;
   return wfcn;
 }      







