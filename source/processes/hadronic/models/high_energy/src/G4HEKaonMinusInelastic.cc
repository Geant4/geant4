// G4 Process: Gheisha High Energy Collision model.
// This includes the high energy cascading model, the two-body-resonance model
// and the low energy two-body model. Not included are the low energy stuff like
// nuclear reactions, nuclear fission without any cascading and all processes for
// particles at rest.  
// First work done by J.L.Chuma and F.W.Jones, TRIUMF, June 96.  
// H. Fesefeldt, RWTH-Aachen, 23-October-1996
// Last modified: 29-July-1998 
 
#include "G4HEKaonMinusInelastic.hh"

G4VParticleChange *  G4HEKaonMinusInelastic::
ApplyYourself( const G4Track &aTrack, G4Nucleus &targetNucleus )
  {
    G4HEVector * pv = new G4HEVector[MAXPART];
    theParticleChange.Initialize( aTrack );
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    const G4double A = targetNucleus.GetN();
    const G4double Z = targetNucleus.GetZ();
    G4HEVector incidentParticle(aParticle);
     
    G4double atomicNumber = Z;
    G4double atomicWeight = A;

    G4int    incidentCode          = incidentParticle.getCode();
    G4double incidentMass          = incidentParticle.getMass();
    G4double incidentTotalEnergy   = incidentParticle.getEnergy();
    G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
    G4double incidentKineticEnergy = incidentTotalEnergy - incidentMass;

    if(incidentKineticEnergy < 1.)
      { 
        cout << "GHEKaonMinusInelastic: incident energy < 1 GeV" << endl;
      }
    if(verboseLevel > 1)
      {
        cout << "G4HEKaonMinusInelastic::ApplyYourself" << endl;
        cout << "incident particle " << incidentParticle.getName()
             << "mass "              << incidentMass
             << "kinetic energy "    << incidentKineticEnergy
             << endl;
        cout << "target material with (A,Z) = (" 
             << atomicWeight << "," << atomicNumber << ")" << endl;
      }

    G4double inelasticity  = NuclearInelasticity(incidentKineticEnergy, 
                                                 atomicWeight, atomicNumber);
    if(verboseLevel > 1)
        cout << "nuclear inelasticity = " << inelasticity << endl;

    incidentKineticEnergy -= inelasticity;
    
    G4double excitationEnergyGNP = 0.;
    G4double excitationEnergyDTA = 0.; 

    G4double excitation    = NuclearExcitation(incidentKineticEnergy,
                                               atomicWeight, atomicNumber,
                                               excitationEnergyGNP,
                                               excitationEnergyDTA);
    if(verboseLevel > 1)
      cout << "nuclear excitation = " << excitation << excitationEnergyGNP 
           << excitationEnergyDTA << endl;             


    incidentKineticEnergy -= excitation;
    incidentTotalEnergy    = incidentKineticEnergy + incidentMass;
    incidentTotalMomentum  = sqrt( (incidentTotalEnergy-incidentMass)                    
                                  *(incidentTotalEnergy+incidentMass));


    G4HEVector targetParticle;
    if(G4UniformRand() < atomicNumber/atomicWeight)
      { 
        targetParticle.setDefinition("Proton");
      }
    else
      { 
        targetParticle.setDefinition("Neutron");
      }

    G4double targetMass         = targetParticle.getMass();
    G4double centerOfMassEnergy = sqrt( incidentMass*incidentMass + targetMass*targetMass
                                       + 2.0*targetMass*incidentTotalEnergy);
    G4double availableEnergy    = centerOfMassEnergy - targetMass - incidentMass;

                                                                // this was the meaning of inElastic in the
                                                                // original Gheisha stand-alone version. 
//    G4bool   inElastic          = InElasticCrossSectionInFirstInt
//                                    (availableEnergy, incidentCode, incidentTotalMomentum);  
                                                                // by unknown reasons, it has been replaced
                                                                // to the following code in Geant???
    G4bool inElastic = true;
//    if (G4UniformRand() < elasticCrossSection/totalCrossSection) inElastic = false;    


    vecLength  = 0;           
        
    if(verboseLevel > 1)
      cout << "ApplyYourself: CallFirstIntInCascade for particle "
           << incidentCode << endl;

    G4bool successful = false; 
    
    if(inElastic || (!inElastic && atomicWeight < 1.5))
      { 
        FirstIntInCasKaonMinus(inElastic, availableEnergy, pv, vecLength,
                               incidentParticle, targetParticle);

        if(verboseLevel > 1)
	   cout << "ApplyYourself::StrangeParticlePairProduction" << endl;  


        if ((vecLength > 0) && (availableEnergy > 1.)) 
                   StrangeParticlePairProduction( availableEnergy, centerOfMassEnergy,
                                                  pv, vecLength,
                                                  incidentParticle, targetParticle);
            HighEnergyCascading( successful, pv, vecLength,
                                 excitationEnergyGNP, excitationEnergyDTA,
                                 incidentParticle, targetParticle,
                                 atomicWeight, atomicNumber);
        if (!successful)
            HighEnergyClusterProduction( successful, pv, vecLength,
                                         excitationEnergyGNP, excitationEnergyDTA,
                                         incidentParticle, targetParticle,
                                         atomicWeight, atomicNumber);
        if (!successful) 
            MediumEnergyCascading( successful, pv, vecLength, 
                                   excitationEnergyGNP, excitationEnergyDTA, 
                                   incidentParticle, targetParticle,
                                   atomicWeight, atomicNumber);

        if (!successful)
            MediumEnergyClusterProduction( successful, pv, vecLength,
                                           excitationEnergyGNP, excitationEnergyDTA,       
                                           incidentParticle, targetParticle,
                                           atomicWeight, atomicNumber);
        if (!successful)
            QuasiElasticScattering( successful, pv, vecLength,
                                    excitationEnergyGNP, excitationEnergyDTA,
                                    incidentParticle, targetParticle, 
                                    atomicWeight, atomicNumber);
      }
    if (!successful)
      { 
            ElasticScattering( successful, pv, vecLength,
                               incidentParticle,    
                               atomicWeight, atomicNumber);
      }

    if (!successful)
      { 
        cout << "GHEInelasticInteraction::ApplyYourself fails to produce final state particles" << endl;
      }
      FillParticleChange(pv,  vecLength);
      delete [] pv;
      theParticleChange.SetStatusChange(fStopAndKill);
      return & theParticleChange;
  }

void
  G4HEKaonMinusInelastic::FirstIntInCasKaonMinus( G4bool &inElastic,
                                                  const G4double availableEnergy,
                                                  G4HEVector pv[],
                                                  G4int &vecLen,
                                                  G4HEVector incidentParticle,
                                                  G4HEVector targetParticle)

// Kaon- undergoes interaction with nucleon within a nucleus.  Check if it is
// energetically possible to produce pions/kaons.  In not, assume nuclear excitation
// occurs and input particle is degraded in energy. No other particles are produced.
// If reaction is possible, find the correct number of pions/protons/neutrons
// produced using an interpolation to multiplicity data.  Replace some pions or
// protons/neutrons by kaons or strange baryons according to the average
// multiplicity per inelastic reaction.

 {
   static const G4double expxu =  log(MAXFLOAT); // upper bound for arg. of exp
   static const G4double expxl = -expxu;         // lower bound for arg. of exp

   static const G4double protb = 0.7;
   static const G4double neutb = 0.7;
   static const G4double     c = 1.25;

   static const G4int   numMul = 1200;
   static const G4int   numSec = 60;

   G4int              neutronCode = Neutron.getCode();
   G4int               protonCode = Proton.getCode();
   G4int            kaonMinusCode = KaonMinus.getCode();
   G4int             kaonZeroCode = KaonZero.getCode();
   G4int         antiKaonZeroCode = AntiKaonZero.getCode();  

   G4int               targetCode = targetParticle.getCode();
   G4double          incidentMass = incidentParticle.getMass();
   G4double        incidentEnergy = incidentParticle.getEnergy();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();

   static G4bool first = true;
   static G4double protmul[numMul], protnorm[numSec];  // proton constants
   static G4double neutmul[numMul], neutnorm[numSec];  // neutron constants

//                                misc. local variables
//                                np = number of pi+,  nm = number of pi-,  nz = number of pi0

   G4int i, counter, nt, np, nm, nz;

   if( first ) 
     {                         // compute normalization constants, this will only be done once
       first = false;
       for( i=0; i<numMul; i++ )protmul[i]  = 0.0;
       for( i=0; i<numSec; i++ )protnorm[i] = 0.0;
       counter = -1;
       for( np=0; np<(numSec/3); np++ ) 
          {
            for( nm=max(0,np-1); nm<=np+1; nm++ ) 
               {
                 for( nz=0; nz<numSec/3; nz++ ) 
                    {
                      if( ++counter < numMul ) 
                        {
                          nt = np+nm+nz;
                          if( (nt>0) && (nt<=numSec) ) 
                            {
                              protmul[counter] =
                                    pmltpc(np,nm,nz,nt,protb,c) ;
                              protnorm[nt-1] += protmul[counter];
                            }
                        }
                    }
               }
          }
       for( i=0; i<numMul; i++ )neutmul[i]  = 0.0;
       for( i=0; i<numSec; i++ )neutnorm[i] = 0.0;
       counter = -1;
       for( np=0; np<numSec/3; np++ ) 
          {
            for( nm=np; nm<=(np+2); nm++ ) 
               {
                 for( nz=0; nz<numSec/3; nz++ ) 
                    {
                      if( ++counter < numMul ) 
                        {
                          nt = np+nm+nz;
                          if( (nt>0) && (nt<=numSec) ) 
                            {
                               neutmul[counter] =
                                      pmltpc(np,nm,nz,nt,neutb,c);
                               neutnorm[nt-1] += neutmul[counter];
                            }
                        }
                    }
               }
          }
       for( i=0; i<numSec; i++ ) 
          {
            if( protnorm[i] > 0.0 )protnorm[i] = 1.0/protnorm[i];
            if( neutnorm[i] > 0.0 )neutnorm[i] = 1.0/neutnorm[i];
          }
     }                                          // end of initialization

         
                                              // initialize the first two places
                                              // the same as beam and target                                    
   pv[0] = incidentParticle;
   pv[1] = targetParticle;
   vecLen = 2;

   if (!inElastic || (availableEnergy <= PionPlus.getMass())) 
      return;

                                        
//                                            inelastic scattering

   np = 0, nm = 0, nz = 0;
   G4double cech[] = { 1., 1., 1., 0.70, 0.60, 0.55, 0.35, 0.25, 0.18, 0.15};
   G4int iplab = G4int( incidentTotalMomentum*5.);
   if( (iplab < 10) && (G4UniformRand() < cech[iplab])) 
     {
       G4int     iplab = min(19, G4int( incidentTotalMomentum*5.));
       G4double cnk0[] = {0.17, 0.18, 0.17, 0.24, 0.26, 0.20, 0.22, 0.21, 0.34, 0.45,
                          0.58, 0.55, 0.36, 0.29, 0.29, 0.32, 0.32, 0.33, 0.33, 0.33};
       if( G4UniformRand() < cnk0[iplab] )
         {
           if( targetCode == protonCode )
             {
               pv[0] = AntiKaonZero;
               pv[1] = Neutron;
               return;
             }
           else
             {
               return;
             }
         }
       G4double ran = G4UniformRand();
       if( targetCode == protonCode )                    // target is a proton 
         {
           if( ran < 0.25 )
             { 
               pv[0] = PionMinus;
               pv[1] = SigmaPlus;
             } 
           else if (ran < 0.50)
             {
               pv[0] = PionZero;
               pv[1] = SigmaZero;
             }
           else if (ran < 0.75)
             { 
               pv[0] = PionPlus;
               pv[1] = SigmaMinus;
             }
           else
             {
               pv[0] = PionZero;
               pv[1] = Lambda;
             }
         } 
       else 
         {                                               // target is a neutron
           if( ran < 0.25 )
             { 
             } 
           else if (ran < 0.50)
             {
               pv[0] = PionMinus;
               pv[1] = SigmaZero;
             }
           else if (ran < 0.75)
             { 
               pv[0] = PionZero;
               pv[1] = SigmaMinus;
             }
           else
             {
               pv[0] = PionMinus;
               pv[1] = Lambda;
             }
         }
       return;
     }
   else
     {
//                    number of total particles vs. centre of mass Energy - 2*proton mass
   
       G4double aleab = log(availableEnergy);
       G4double n     = 3.62567+aleab*(0.665843+aleab*(0.336514
                    + aleab*(0.117712+0.0136912*aleab))) - 2.0;
   
//                    normalization constant for kno-distribution.
//                    calculate first the sum of all constants, check for numerical problems.   
       G4double test, dum, anpn = 0.0;

       for( nt=1; nt<=numSec; nt++ ) 
         {
           test = exp( min( expxu, max( expxl, -(M_PI/4.0)*(nt*nt)/(n*n) ) ) );
           dum = M_PI*nt/(2.0*n*n);
           if( fabs(dum) < 1.0 ) 
             if( test >= 1.0e-10 )anpn += dum*test;
           else 
             anpn += dum*test;
         }
   
       G4double ran = G4UniformRand();
       G4double excs = 0.0;
       if( targetCode == protonCode ) 
         {
           counter = -1;
           for( np=0; np<numSec/3; np++ ) 
              {
                for( nm=max(0,np-1); nm<=np+1; nm++ ) 
                   {
                     for( nz=0; nz<numSec/3; nz++ ) 
                        {
                          if( ++counter < numMul ) 
                            {
                              nt = np+nm+nz;
                              if( (nt>0) && (nt<=numSec) ) 
                                {
                                  test = exp( min( expxu, max( expxl, -(M_PI/4.0)*(nt*nt)/(n*n) ) ) );
                                  dum = (M_PI/anpn)*nt*protmul[counter]*protnorm[nt-1]/(2.0*n*n);
                                  if( fabs(dum) < 1.0 ) 
                                        if( test >= 1.0e-10 )excs += dum*test;
                                   else 
                                        excs += dum*test;
                                   if (ran < excs) goto outOfLoop;      //----------------------->
                                }   
                            }    
                        }     
                   }                                                                                  
              }
       
                                              // 3 previous loops continued to the end
           inElastic = false;                 // quasi-elastic scattering   
           return;
         }
       else   
         {                                         // target must be a neutron
           counter = -1;
           for( np=0; np<numSec/3; np++ ) 
              {
                for( nm=np; nm<=(np+2); nm++ ) 
                   {
                     for( nz=0; nz<numSec/3; nz++ ) 
                        {
                          if( ++counter < numMul ) 
                            {
                              nt = np+nm+nz;
                              if( (nt>=1) && (nt<=numSec) ) 
                                {
                                  test = exp( min( expxu, max( expxl, -(M_PI/4.0)*(nt*nt)/(n*n) ) ) );
                                  dum = (M_PI/anpn)*nt*neutmul[counter]*neutnorm[nt-1]/(2.0*n*n);
                                  if( fabs(dum) < 1.0 ) 
                                      if( test >= 1.0e-10 )excs += dum*test;
                                  else 
                                      excs += dum*test;
                                  if (ran < excs) goto outOfLoop;       // -------------------------->
                                }
                            }
                        }
                   }
              }
                                                  // 3 previous loops continued to the end
           inElastic = false;                     // quasi-elastic scattering.
           return;
         }
     } 
   outOfLoop:           //  <------------------------------------------------------------------------   
    
   if( targetCode == protonCode)
     {
       if( np == (1+nm))
         {
           pv[1] = Neutron;
         }
       else if (np == nm)
         {
           if( G4UniformRand() < 0.75)
             {
             }
           else
             {
               pv[0] = AntiKaonZero;
               pv[1] = Neutron;
             }
         }
       else      
         {
           pv[0] = AntiKaonZero;
         } 
     }  
   else
     {
       if( np == (nm-1))
         {
           if( G4UniformRand() < 0.5)
             {
               pv[1] = Proton;
             }
           else
             {
               pv[0] = AntiKaonZero;
             }
         } 
       else if ( np == nm)
         {
         }
       else
         {
           pv[0] = AntiKaonZero;
           pv[1] = Proton;
         }
     }      

   if( G4UniformRand() < 0.5 )
     {
       if(    (    (pv[0].getCode() == kaonMinusCode)
                && (pv[1].getCode() == neutronCode)  )
           || (    (pv[0].getCode() == kaonZeroCode)
                && (pv[1].getCode() == protonCode)   )
           || (    (pv[0].getCode() == antiKaonZeroCode)
                && (pv[1].getCode() == protonCode)   )   )
         {
           G4double ran = G4UniformRand();
           if( pv[1].getCode() == protonCode)
             { 
               if(ran < 0.68)
                 {
                   pv[0] = PionPlus;
                   pv[1] = Lambda;
                 }
               else if (ran < 0.84)
                 {
                   pv[0] = PionZero;
                   pv[1] = SigmaPlus;
                 }
               else
                 {
                   pv[0] = PionPlus;
                   pv[1] = SigmaZero;
                 }
             }
           else
             {
               if(ran < 0.68)
                 {
                   pv[0] = PionMinus;
                   pv[1] = Lambda;
                 }
               else if (ran < 0.84)
                 {
                   pv[0] = PionMinus;
                   pv[1] = SigmaZero;
                 }
               else
                 {
                   pv[0] = PionZero;
                   pv[1] = SigmaMinus;
                 }
             }
         } 
       else
         {
           G4double ran = G4UniformRand();
           if (ran < 0.67)
              {
                pv[0] = PionZero;
                pv[1] = Lambda;
              }
           else if (ran < 0.78)
              {
                pv[0] = PionMinus;
                pv[1] = SigmaPlus;
              }
           else if (ran < 0.89)
              {
                pv[0] = PionZero;
                pv[1] = SigmaZero;
              }
           else
              {
                pv[0] = PionPlus;
                pv[1] = SigmaMinus;
              }
         }
     }
               

   nt = np + nm + nz;
   while ( nt > 0)
       {
         G4double ran = G4UniformRand();
         if ( ran < (G4double)np/nt)
            { 
              if( np > 0 ) 
                { pv[vecLen++] = PionPlus;
                  np--;
                }
            }
         else if ( ran < (G4double)(np+nm)/nt)
            {   
              if( nm > 0 )
                { 
                  pv[vecLen++] = PionMinus;
                  nm--;
                }
            }
         else
            {
              if( nz > 0 )
                { 
                  pv[vecLen++] = PionZero;
                  nz--;
                }
            }
         nt = np + nm + nz;
       } 
   if (verboseLevel > 1)
      {
        cout << "Particles produced: " ;
        cout << pv[0].getName() << " " ;
        cout << pv[1].getName() << " " ;
        for (i=2; i < vecLen; i++)   
            { 
              cout << pv[i].getName() << " " ;
            }
         cout << endl;
      }
   return;
 }









