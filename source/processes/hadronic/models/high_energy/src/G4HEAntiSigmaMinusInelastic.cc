// G4 Process: Gheisha High Energy Collision model.
// This includes the high energy cascading model, the two-body-resonance model
// and the low energy two-body model. Not included are the low energy stuff like
// nuclear reactions, nuclear fission without any cascading and all processes for
// particles at rest.  
// First work done by J.L.Chuma and F.W.Jones, TRIUMF, June 96.  
// H. Fesefeldt, RWTH-Aachen, 23-October-1996
// Last modified: 29-July-1998 
 
#include "G4HEAntiSigmaMinusInelastic.hh"

G4VParticleChange *  G4HEAntiSigmaMinusInelastic::
ApplyYourself( const G4Track &aTrack, G4Nucleus &targetNucleus )
  {
    G4HEVector * pv = new G4HEVector[MAXPART];
    theParticleChange.Initialize( aTrack );
    const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
    G4DynamicParticle *originalTarget = targetNucleus.ReturnTargetParticle();
    const G4double atomicWeight = targetNucleus.GetN();
    const G4double atomicNumber = targetNucleus.GetZ();
    G4HEVector incidentParticle(aParticle);

    G4int    incidentCode          = incidentParticle.getCode();
    G4double incidentMass          = incidentParticle.getMass();
    G4double incidentTotalEnergy   = incidentParticle.getEnergy();
    G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();
    G4double incidentKineticEnergy = incidentTotalEnergy - incidentMass;

    if(incidentKineticEnergy < 1.)
      { 
        cout << "GHEAntiSigmaMinusInelastic: incident energy < 1 GeV" << endl;
      }
    if(verboseLevel > 1)
      {
        cout << "G4HEAntiSigmaMinusInelastic::ApplyYourself" << endl;
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
        FirstIntInCasAntiSigmaMinus(inElastic, availableEnergy, pv, vecLength,
                                    incidentParticle, targetParticle, atomicWeight);

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
G4HEAntiSigmaMinusInelastic::FirstIntInCasAntiSigmaMinus( G4bool &inElastic,
                                                          const G4double availableEnergy,
                                                          G4HEVector pv[],
                                                          G4int &vecLen,
                                                          G4HEVector incidentParticle,
                                                          G4HEVector targetParticle,
                                                          const G4double atomicWeight)

// AntiSigma- undergoes interaction with nucleon within a nucleus.  Check if it is
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

   static const G4int   numMul   = 1200;
   static const G4int   numMulAn = 400;
   static const G4int   numSec   = 60;

   G4int              neutronCode = Neutron.getCode();
   G4int              protonCode  = Proton.getCode();

   G4int               targetCode = targetParticle.getCode();
   G4double          incidentMass = incidentParticle.getMass();
   G4double        incidentEnergy = incidentParticle.getEnergy();
   G4double incidentTotalMomentum = incidentParticle.getTotalMomentum();

   static G4bool first = true;
   static G4double protmul[numMul],  protnorm[numSec];   // proton constants
   static G4double protmulAn[numMulAn],protnormAn[numSec]; 
   static G4double neutmul[numMul],  neutnorm[numSec];   // neutron constants
   static G4double neutmulAn[numMulAn],neutnormAn[numSec];

                              //  misc. local variables
                              //  np = number of pi+,  nm = number of pi-,  nz = number of pi0

   G4int i, counter, nt, np, nm, nz;

   if( first ) 
     {                         // compute normalization constants, this will only be done once
       first = false;
       for( i=0; i<numMul  ; i++ ) protmul[i]  = 0.0;
       for( i=0; i<numSec  ; i++ ) protnorm[i] = 0.0;
       counter = -1;
       for( np=0; np<(numSec/3); np++ ) 
          {
            for( nm=max(0,np-2); nm<=np; nm++ ) 
               {
                 for( nz=0; nz<numSec/3; nz++ ) 
                    {
                      if( ++counter < numMul ) 
                        {
                          nt = np+nm+nz;
                          if( (nt>0) && (nt<=numSec) ) 
                            {
                              protmul[counter] = pmltpc(np,nm,nz,nt,protb,c);
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
            for( nm=max(0,np-1); nm<=(np+1); nm++ ) 
               {
                 for( nz=0; nz<numSec/3; nz++ ) 
                    {
                      if( ++counter < numMul ) 
                        {
                          nt = np+nm+nz;
                          if( (nt>0) && (nt<=numSec) ) 
                            {
                               neutmul[counter] = pmltpc(np,nm,nz,nt,neutb,c);
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
                                                                   // annihilation
       for( i=0; i<numMulAn  ; i++ ) protmulAn[i]  = 0.0;
       for( i=0; i<numSec    ; i++ ) protnormAn[i] = 0.0;
       counter = -1;
       for( np=1; np<(numSec/3); np++ ) 
          {
            nm = max(0,np-2); 
            for( nz=0; nz<numSec/3; nz++ ) 
               {
                 if( ++counter < numMulAn ) 
                   {
                     nt = np+nm+nz;
                     if( (nt>1) && (nt<=numSec) ) 
                       {
                         protmulAn[counter] = pmltpc(np,nm,nz,nt,protb,c);
                         protnormAn[nt-1] += protmulAn[counter];
                       }
                   }
               }
          }
       for( i=0; i<numMulAn; i++ ) neutmulAn[i]  = 0.0;
       for( i=0; i<numSec;   i++ ) neutnormAn[i] = 0.0;
       counter = -1;
       for( np=0; np<numSec/3; np++ ) 
          {
            nm = np-1; 
            for( nz=0; nz<numSec/3; nz++ ) 
               {
                 if( ++counter < numMulAn ) 
                   {
                     nt = np+nm+nz;
                     if( (nt>1) && (nt<=numSec) ) 
                       {
                          neutmulAn[counter] = pmltpc(np,nm,nz,nt,neutb,c);
                          neutnormAn[nt-1] += neutmulAn[counter];
                       }
                   }
               }
          }
       for( i=0; i<numSec; i++ ) 
          {
            if( protnormAn[i] > 0.0 )protnormAn[i] = 1.0/protnormAn[i];
            if( neutnormAn[i] > 0.0 )neutnormAn[i] = 1.0/neutnormAn[i];
          }
     }                                          // end of initialization

         
                                              // initialize the first two places
                                              // the same as beam and target                                    
   pv[0] = incidentParticle;
   pv[1] = targetParticle;
   vecLen = 2;

   if( !inElastic ) 
     {                                        // some two-body reactions 
       G4double cech[] = {0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.06, 0.04, 0.005, 0.};

       G4int iplab = min(9, G4int( incidentTotalMomentum*2.5 ));
       if( G4UniformRand() < cech[iplab]/pow(atomicWeight,0.42) ) 
         {           
           G4double ran = G4UniformRand();

           if ( targetCode == neutronCode)
             {
               if(ran < 0.2)
                 {
                   pv[0] = AntiSigmaZero;
                   pv[1] = Proton;
                 }
               else if (ran < 0.4)
                 {
                   pv[0] = AntiLambda;
                   pv[1] = Proton;
                 }
               else if (ran < 0.6)
                 {
                   pv[0] = Proton;
                   pv[1] = AntiLambda;
                 }
               else if (ran < 0.8)
                 {
                   pv[0] = Proton;
                   pv[1] = AntiSigmaZero;
                 }
               else
                 {
                   pv[0] = Neutron;
                   pv[1] = AntiSigmaMinus;
                 }     
             }
           else
             {
               pv[0] = Proton;
               pv[1] = AntiSigmaMinus;
             }  
         }   
       return;
     }
   else if (availableEnergy <= PionPlus.getMass())
       return;

                                                  //   inelastic scattering

   np = 0; nm = 0; nz = 0;
   G4double anhl[] = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.97, 0.88, 
                      0.85, 0.81, 0.75, 0.64, 0.64, 0.55, 0.55, 0.45, 0.47, 0.40, 
                      0.39, 0.36, 0.33, 0.10, 0.01};
   G4int            iplab =      G4int( incidentTotalMomentum*10.);
   if ( iplab >  9) iplab = 10 + G4int( (incidentTotalMomentum  -1.)*5. );          
   if ( iplab > 14) iplab = 15 + G4int(  incidentTotalMomentum  -2.     );
   if ( iplab > 22) iplab = 23 + G4int( (incidentTotalMomentum -10.)/10.); 
                    iplab = min(24, iplab);

   if ( G4UniformRand() > anhl[iplab] )
     {                                           // non- annihilation channels

                         //  number of total particles vs. centre of mass Energy - 2*proton mass
   
           G4double aleab = log(availableEnergy);
           G4double n     = 3.62567+aleab*(0.665843+aleab*(0.336514
                            + aleab*(0.117712+0.0136912*aleab))) - 2.0;
   
                         // normalization constant for kno-distribution.
                         // calculate first the sum of all constants, check for numerical problems.   
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
                    for( nm=max(0,np-2); nm<=np; nm++ ) 
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
                    for( nm=max(0,np-1); nm<=(np+1); nm++ ) 
                       {
                         for( nz=0; nz<numSec/3; nz++ ) 
                            {
                              if( ++counter < numMul ) 
                                {
                                  nt = np+nm+nz;
                                  if( (nt>0) && (nt<=numSec) ) 
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
          
       outOfLoop:           //  <------------------------------------------------------------------------   
    
       ran = G4UniformRand();

       if( targetCode == protonCode)
         {
           if( np == nm)
             {
             }
           else if (np == (nm+1))
             {
               if( ran < 0.50)
                 {
                   pv[1] = Neutron;
                 }
               else if (ran < 0.75)
                 {
                   pv[0] = AntiSigmaZero;
                 }
               else
                 {
                   pv[0] = AntiLambda;
                 }
             }
           else      
             {
               if (ran < 0.5)
                 {
                   pv[0] = AntiSigmaZero;
                   pv[1] = Neutron;
                 }
               else
                 {
                   pv[0] = AntiLambda;
                   pv[1] = Neutron;
                 } 
             } 
         }  
       else
         {
           if( np == nm)
             {
               if (ran < 0.5)
                 {
                 }
               else if(ran < 0.75)
                 {
                   pv[0] = AntiSigmaZero;
                   pv[1] = Proton;
                 }
               else
                 {
                   pv[0] = AntiLambda;
                   pv[1] = Proton;
                 }
             } 
           else if ( np == (nm+1))
             {
               if (ran < 0.5)
                 {
                   pv[0] = AntiSigmaZero;
                 }
               else 
                 {
                   pv[0] = AntiLambda;
                 }
             }
           else 
             {
               pv[1] = Proton;
             }
         }      
     
     }
   else                                                               // annihilation
     {  
       if ( availableEnergy > 2. * PionPlus.getMass() )
         {

           G4double aleab = log(availableEnergy);
           G4double n     = 3.62567+aleab*(0.665843+aleab*(0.336514
                            + aleab*(0.117712+0.0136912*aleab))) - 2.0;
   
                      //   normalization constant for kno-distribution.
                      //   calculate first the sum of all constants, check for numerical problems.   
           G4double test, dum, anpn = 0.0;

           for( nt=2; nt<=numSec; nt++ ) 
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
               for( np=2; np<numSec/3; np++ ) 
                  {
                    nm = np-2; 
                    for( nz=0; nz<numSec/3; nz++ ) 
                      {
                        if( ++counter < numMulAn ) 
                          {
                            nt = np+nm+nz;
                            if( (nt>1) && (nt<=numSec) ) 
                              {
                                test = exp( min( expxu, max( expxl, -(M_PI/4.0)*(nt*nt)/(n*n) ) ) );
                                dum = (M_PI/anpn)*nt*protmulAn[counter]*protnormAn[nt-1]/(2.0*n*n);
                                if( fabs(dum) < 1.0 ) 
                                     if( test >= 1.0e-10 )excs += dum*test;
                                else 
                                     excs += dum*test;
                                if (ran < excs) goto outOfLoopAn;      //----------------------->
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
               for( np=1; np<numSec/3; np++ ) 
                 { 
                   nm = np-1; 
                   for( nz=0; nz<numSec/3; nz++ ) 
                      {
                        if( ++counter < numMulAn ) 
                          {
                            nt = np+nm+nz;
                            if( (nt>1) && (nt<=numSec) ) 
                              {
                                test = exp( min( expxu, max( expxl, -(M_PI/4.0)*(nt*nt)/(n*n) ) ) );
                                dum = (M_PI/anpn)*nt*neutmulAn[counter]*neutnormAn[nt-1]/(2.0*n*n);
                                if( fabs(dum) < 1.0 ) 
                                    if( test >= 1.0e-10 )excs += dum*test;
                                else 
                                excs += dum*test;
                                if (ran < excs) goto outOfLoopAn;       // -------------------------->
                              }
                          }
                      }
                 }
               inElastic = false;                     // quasi-elastic scattering.
               return;
             }
       outOfLoopAn:           //  <------------------------------------------------------------------   
       vecLen = 0;
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










