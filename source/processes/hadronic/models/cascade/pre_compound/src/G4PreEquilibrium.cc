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
 // $Id: G4PreEquilibrium.cc,v 1.9 2001/07/11 10:03:55 gunter Exp $
 // GEANT4 tag $Name: geant4-05-02 $
 //
 // Hadronic Process: Pre-equilibrium HETC 
 // Joseph L. Chuma, TRIUMF, 24-Mar-2000

#include "G4PreEquilibrium.hh"
#include "Randomize.hh"

 G4double oneThird = 1./3.;
 G4double twoThirds = 2./3.;
 G4double fourThirds = 4./3.;
 
 // S(Z)+P(Z) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
 // or Delta M(Z) from Tab. 97 of book [1]
 G4double G4PreEquilibrium::t1y[130] =
 { 20.80,  15.80,  21.00,  16.80,  19.80,
   16.50,  18.80,  16.50,  18.50,  17.20,
   18.26,  15.05,  16.01,  12.04,  13.27,
   11.09,  12.17,  10.26,  11.04,   8.41,
    9.79,   7.36,   8.15,   5.63,   5.88,
    3.17,   3.32,   0.82,   1.83,   0.97,
    2.33,   1.27,   2.92,   1.61,   2.91,
    1.35,   2.40,   0.89,   1.74,   0.36,
    0.95,  -0.65,  -0.04,  -1.73,  -0.96,
   -2.87,  -2.05,  -4.05,  -3.40,  -5.72,
   -3.75,  -4.13,  -2.42,  -2.85,  -1.01,
   -1.33,   0.54,  -0.02,   1.74,   0.75,
    2.24,   1.00,   1.98,   0.79,   1.54,
    0.39,   1.08,   0.00,   0.78,  -0.35,
    0.58,  -0.55,   0.59,  -0.61,   0.59,
   -0.35,   0.32,  -0.96,  -0.52,  -2.08,
   -2.46,  -3.64,  -1.55,  -0.96,   0.97,
    0.88,   2.37,   1.75,   2.72,   1.90,
    2.55,   1.46,   1.93,   0.86,   1.17,
    0.08,   0.39,  -0.76,  -0.39,  -1.51,
   -1.17,  -2.36,  -1.95,  -3.06,  -2.62,
   -3.55,  -2.95,  -3.75,  -3.07,  -3.79,
   -3.06,  -3.77,  -3.05,  -3.78,  -3.12,
   -3.90,  -3.35,  -4.24,  -3.86,  -4.92,
   -5.06,  -6.77,  -7.41,  -9.18, -10.16,
  -11.12,  -9.76,  -9.23,  -7.96,  -7.65 };
 
 // S(N)+P(N) from Tab. 1 from A.G.W. Cameron, Canad. J. Phys., 35(1957)1021
 // or Delta M(N) from Tab. 97 of book [1]
 G4double G4PreEquilibrium::t2xy[] = 
 {   -8.40,  -12.90,   -8.00,   11.90,   -9.20,
    -12.50,  -10.80,  -13.60,  -11.20,  -12.20,
    -12.81,  -15.40,  -13.07,  -15.80,  -13.81,
    -14.98,  -12.63,  -13.76,  -11.37,  -12.38,
     -9.23,   -9.65,   -7.64,   -9.17,   -8.05,
     -9.72,   -8.87,  -10.76,   -8.64,   -8.89,
     -6.60,   -7.13,   -4.77,   -5.33,   -3.06,
     -3.79,   -1.72,   -2.79,   -0.93,   -2.19,
     -0.52,   -1.90,   -0.45,   -2.20,   -1.22,
     -3.07,   -2.42,   -4.37,   -3.94,   -6.08,
     -4.49,   -4.50,   -3.14,   -2.93,   -1.04,
     -1.36,    0.69,    0.21,    2.11,    1.33,
      3.29,    2.46,    4.30,    3.32,    4.79,
      3.62,    4.97,    3.64,    4.63,    3.07,
      4.06,    2.49,    3.30,    1.46,    2.06,
      0.51,    0.74,   -1.18,   -1.26,   -3.54,
     -3.97,   -5.26,   -4.18,   -3.71,   -2.10,
     -1.70,   -0.08,   -0.18,    0.94,    0.27,
      1.13,    0.08,    0.91,   -0.31,    0.49,
     -0.78,    0.08,   -1.15,   -0.23,   -1.41,
     -0.42,   -1.55,   -0.55,   -1.66,   -0.66,
     -1.73,   -0.75,   -1.74,   -0.78,   -1.69,
     -0.78,   -1.60,   -0.75,   -1.46,   -0.67,
     -1.26,   -0.51,   -1.04,   -0.53,   -1.84,
     -2.42,   -4.52,   -4.76,   -6.33,   -6.76,
     -7.81,   -5.80,   -5.37,   -3.63,   -3.35,
     -1.75,   -1.88,   -0.61,   -0.90,    0.09,
     -0.32,    0.55,   -0.13,    0.70,   -0.06,
      0.49,   -0.20,    0.40,   -0.22,    0.36,
     -0.09,    0.58,    0.12,    0.75,    0.15,
      0.70,    0.17,    1.11,    0.89,    1.85,
      1.62,    2.54,    2.29,    3.20,    2.91,
      3.84,    3.53,    4.48,    4.15,    5.12,
      4.78,    5.75,    5.39,    6.31,    5.91,
      6.87,    6.33,    7.13,    6.61,    7.30,
      6.31,    6.27,    4.83,    4.49,    2.85,
      2.32,    0.58,   -0.11,   -0.98,    0.81,
      1.77,    3.37,    4.13,    5.60,    6.15,
      7.29,    7.35,    7.95,    7.67,    8.16,
      7.83,    8.31,    8.01,    8.53,    8.27 
 };
 
 G4double G4PreEquilibrium::cemgeo( const G4int nuclideTypeNumber, const G4int mediumNumber )
  {
    // this routine was called by input
    bertcem( nuclideTypeNumber, mediumNumber );
    return GetA_FCOMON(nuclideTypeNumber,mediumNumber);
  }
 
 G4double G4PreEquilibrium::erupcem( DPvector &neutrons, DPvector &protons, DPvector &deuterons,
                                     DPvector &tritons, DPvector &he3s, DPvector &he4s )
  {
    precof();
    
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aDeuteron = G4Deuteron::Deuteron();
    G4ParticleDefinition *aTriton = G4Triton::Triton();
    G4ParticleDefinition *aHe3 = G4He3::He3();
    //G4ParticleDefinition *aHe4 = G4He4::He4();
    
    G4double protonMass = aProton->GetPDGMass()/MeV;
    G4double neutronMass = aNeutron->GetPDGMass()/MeV;
    G4double deuteronMass = aDeuteron->GetPDGMass()/MeV;
    G4double tritonMass = aTriton->GetPDGMass()/MeV;
    G4double he3Mass = aHe3->GetPDGMass()/MeV;
    //G4double he4Mass = aHe4->GetPDGMass()/MeV;
    
    G4int numberOfProtons = G4int(theNucleus->GetZ()+0.5);
    G4int numberOfNucleons = G4int(theNucleus->GetN()+0.5);
    G4int numberOfNeutrons = numberOfNucleons - numberOfProtons;
    G4double nuclearMass = protonMass*numberOfProtons + neutronMass*numberOfNeutrons;
    G4double P = theNucleus->GetMomentum().mag()/MeV;
    
    G4double totalEnergy = sqrt( P*P + nuclearMass*nuclearMass );
    
    G4double kineticEnergyAfterEvaporation = totalEnergy - nuclearMass;
    
    for( G4int k = 0; k < KTOT; ++k )
    {
      Dvector sptK( *(spt[k]) );
      G4double sint = sptK[0];
      G4double cost = sptK[1];
      G4double sinf = sin( parz[k].angle2 );
      G4double cosf = cos( parz[k].angle2 );
      G4double alc = sint*cosf;
      G4double bec = sint*sinf;
      G4double gac = cost;
      G4double energy = parz[k].energy*GeV/MeV;
      G4DynamicParticle *dp = new G4DynamicParticle();
      G4ThreeVector m;
      switch ( parz[k].particleType )
      {
       case 1:                                       // neutrons
         dp->SetDefinition( aNeutron );
         dp->SetKineticEnergy( energy-neutronMass );
         m.setX( alc ); m.setY( bec ); m.setZ( gac );
         dp->SetMomentum( m );
         neutrons.push_back( dp );
         break;
       case 2:                                       // protons
         dp->SetDefinition( aProton );
         dp->SetKineticEnergy( energy-protonMass );
         m.setX( alc ); m.setY( bec ); m.setZ( gac );
         dp->SetMomentum( m );
         protons.push_back( dp );
         break;
       case 3:                                       // deuterons
         dp->SetDefinition( aDeuteron );
         dp->SetKineticEnergy( energy-deuteronMass );
         m.setX( alc ); m.setY( bec ); m.setZ( gac );
         dp->SetMomentum( m );
         deuterons.push_back( dp );
         break;
       case 4:                                       // tritons
         dp->SetDefinition( aTriton );
         dp->SetKineticEnergy( energy-tritonMass );
         m.setX( alc ); m.setY( bec ); m.setZ( gac );
         dp->SetMomentum( m );
         tritons.push_back( dp );
         break;
       case 5:                                       // he-3
         dp->SetDefinition( aHe3 );
         dp->SetKineticEnergy( energy-he3Mass );
         m.setX( alc ); m.setY( bec ); m.setZ( gac );
         dp->SetMomentum( m );
         he3s.push_back( dp );
         break;
     //case 6:                                       // he-4
         //dp->SetDefinition( aHe4 );
         //dp->SetKineticEnergy( energy-he4Mass );
         //m.setX( alc ); m.setY( bec ); m.setZ( gac );
         //dp->SetMomentum( m );
         //dp->SetMomentum( alc, bec, gac );
         //he4s.push_back( dp );
         //break;
      }
    }
    return kineticEnergyAfterEvaporation;
  }
 
 void G4PreEquilibrium::precof()
  {
    const G4double RM = 1.5;
    
    G4double pevapj[6], gj[7], ami[6];
    G4int i, j;
    SetR0( RM );
    G4double wt = 1.0;
    
    G4double protonMass = G4Proton::Proton()->GetPDGMass()/GeV;
    G4double neutronMass = G4Neutron::Neutron()->GetPDGMass()/GeV;
    G4int numberOfProtons = G4int(theNucleus->GetZ()+0.5);
    G4int numberOfNucleons = G4int(theNucleus->GetN()+0.5);
    G4int numberOfNeutrons = numberOfNucleons - numberOfProtons;
    G4double nuclearMass = protonMass*numberOfProtons + neutronMass*numberOfNeutrons;
    G4double P = theNucleus->GetMomentum().mag()/GeV;
    
    G4double totalEnergy = sqrt( P*P + nuclearMass*nuclearMass );
    G4double kineticEnergy = (totalEnergy-nuclearMass)*GeV/MeV;

    G4ThreeVector momentum = theNucleus->GetMomentum() * (1/totalEnergy*GeV);

    G4ThreeVector uglmom( GetLXYZ() );

    KTOT = 0;
    G4double pz = GetPZ0();
    G4int n = GetN0();
    G4double p = GetP0();
    G4double h = GetH0();
    theNucleus->AddExcitationEnergy( -kineticEnergy );
    //
    if( theNucleus->GetEnergyDeposit() < 3 ||
        theNucleus->GetN() <= 4 ||
        theNucleus->GetZ() <= 2 )return;
    //
    G4double dl = massDefect( theNucleus->GetN(), theNucleus->GetZ() );
    //
    // auxiliary code for nuclear data extraction
    //
    const G4double z1[5] = { 10.00, 20.00, 30.00, 50.00, 70.00 }; // Z:Tab.96
    const G4double a1[5] = {  0.42,  0.58,  0.68,  0.77,  0.80 }; // Kp -//-
    const G4double c1[5] = {  0.50,  0.28,  0.20,  0.15,  0.10 }; // Cp -//-
    const G4double a2[5] = {  0.68,  0.82,  0.91,  0.97,  0.98 }; // Ka -//-
    const G4double c2[5] = {  0.10,  0.10,  0.10,  0.08,  0.06 }; // Ca -//-
    //
    //           particle b:   n      p       d       t      He-3   He-4
    const G4double  aj[6] = { 1.0,   1.0,    2.0,    3.0,    3.0,   4.0   }; // A(b)
    const G4double  zj[6] = { 0.0,   1.0,    1.0,    1.0,    2.0,   2.0   }; // Z(b)
    const G4double dlm[6] = { 8.368, 7.569, 13.835, 15.835, 15.817, 3.607 }; // q(b)
    //
    G4double cc[6], vk[6];
    alj.reserve(6);
    gb.reserve(6);
    afj.reserve(7);
    zfj.reserve(6);
    rj.reserve(7);
    vj.reserve(6);
    bj.reserve(6);
    r0j.reserve(6);
    for( G4int k = 0; k < 100; ++k )
    {
      cc[0] = 0;
      cc[1] = quadraticInterpolation( theNucleus->GetZ(), z1, c1, 4 );
      cc[5] = quadraticInterpolation( theNucleus->GetZ(), z1, c2, 4 );
      cc[2] = cc[1]/2;
      cc[3] = cc[1]/3;
      cc[4] = cc[5]*fourThirds;
      //
      vk[0] = 0;
      vk[1] = quadraticInterpolation( theNucleus->GetZ(), z1, a1, 4 );
      vk[5] = quadraticInterpolation( theNucleus->GetZ(), z1, a2, 4 );
      vk[2] = vk[1] + 0.06;
      vk[3] = vk[1] + 0.12;
      vk[4] = vk[5] - 0.06;
      //
      G4double um2 = uglmom.mag2();
      G4double un = theNucleus->GetN() - theNucleus->GetZ();
      G4double ue = theNucleus->GetEnergyDeposit()-
          12*((1-theNucleus->GetZ()+2*(G4int(theNucleus->GetZ()+0.5)/2))+
              (1-un+2*(G4int(un)/2)))/sqrt(theNucleus->GetN());
      if( ue <= 0.1 )return;
      G4double am = fam( theNucleus->GetN(), theNucleus->GetZ(), ue );
      afj.erase( afj.begin(), afj.end() );
      zfj.erase( zfj.begin(), zfj.end() );
      rj.erase( rj.begin(), rj.end() );
      bj.erase( bj.begin(), bj.end() );
      vj.erase( vj.begin(), vj.end() );
      for( i = 0; i < 6; ++i )
      {
        G4double naj = theNucleus->GetN() - aj[i];
        G4double zzj = theNucleus->GetZ() - zj[i];
        afj.push_back( naj );
        zfj.push_back( zzj );
        G4double pairj = (1-zzj+2*(G4int(zzj)/2)) + (1-(naj-zzj)+2*(G4int(naj-zzj)/2));
        pevapj[i] = 12 * pairj / sqrt(naj);
        G4double uej = theNucleus->GetEnergyDeposit() - pevapj[i];
        //
        if( uej <= 1 )rj.push_back( 0.0 );
        else
        {
          ami[i] = fam( naj, zzj, uej );
          //
          // calculation of coulomb energy
          //
          G4double tmp = vk[i]*1.44/RM*zj[i]*zzj/
              (pow(aj[i],oneThird)+pow(naj,oneThird))*
              (1-theNucleus->GetEnergyDeposit()/(81*theNucleus->GetN()*am));
          tmp < 0 ? vj.push_back( 0.0 ) : vj.push_back( tmp );
          //
          bj.push_back( massDefect(naj,zzj) - (dl-dlm[i]) );
          rj.push_back( uej-bj[i]-vj[i] );
        }
      }
      afj.push_back( theNucleus->GetN() );
      //
      r0j.erase( r0j.begin(), r0j.end() );
      r0j.push_back( 0.76+2.2/pow(afj[1],oneThird) );
      SetBN( (2.12/pow(afj[1],twoThirds)-0.05)/r0j[0] );
      for( i = 1; i < 6; ++i )r0j.push_back( 1+cc[i] );
      G4double ep1, ep2, ep3;
      G4bool ind;
      G4ThreeVector p12;
      G4int lm;
    LABEL:
      if( n < sqrt(1.19*am*theNucleus->GetN()*theNucleus->GetEnergyDeposit()+0.5) )
      {
        // pre-equilibrium emission
        //
        while ( p < 1 )
        {
          ++p;
          ++h;
          n += 2;
        }
        exn = p + h;    
        alj.erase( alj.begin(), alj.end() );
        alj.push_back( p*(exn-1) );
        alj.push_back( alj[0] );
        alj.push_back( alj[1]*(p-1)*(exn-2)/2 );
        alj.push_back( alj[2]*(p-2)*(exn-3)/6 );
        alj.push_back( alj[3] );
        alj.push_back( alj[4]*(p-3)*(exn-4)/12 );
        //
        G4double nucleusN = theNucleus->GetN();
        gb.erase( gb.begin(), gb.end() );
        gb.push_back( 1. );
        gb.push_back( 1. );
        gb.push_back( 16/nucleusN );
        gb.push_back( 243/nucleusN/nucleusN );
        gb.push_back( gb[3] );
        gb.push_back( 4096/nucleusN/nucleusN/nucleusN );
        for( i = 0; i < 6; ++i )
        {
          if( p <= aj[i] - 0.01 )gj[i] = 0;
          else
          {
            if( p+h <= aj[i] + 0.01 )gj[i] = 0;
            else
            {
              if( pz <= zj[i] - 0.01 )gj[i] = 0;
              else
              {
                if( rj[i] )gj[i] = 0;
                else
                {
                  SetAC( 0.585*ami[i] );
                  gj[i] = gamagu(i);
                }
              }
            }
          }
        }
        G4double g = 0;
        for( i = 0; i < 6; ++i )g += gj[i];
        if( n <= 0 )
        {
          // write(16,20)
          // 20 format(20x,20hnumber of exitons = 0)
          return;
        }
        SetAC( 0.595*am );
        G4double c1, c2, c3;
        transitionRates( p, h, c1, c2, c3 );
        G4double c = c1+c2+c3;
        G4double b1 = G4UniformRand();
        G4double p2 = 250*pow(G4UniformRand(),oneThird);
        G4double ct2 = 1-2*G4UniformRand();
        G4double st2 = sqrt(1-ct2*ct2);
        G4double fi2 = 2*pi*G4UniformRand();
        p12.setX( p2*st2*cos(fi2)+theNucleus->GetMomentum().x()*GeV/p );
        p12.setY( p2*st2*sin(fi2)+theNucleus->GetMomentum().y()*GeV/p );
        p12.setZ( p2*ct2+theNucleus->GetMomentum().z()*GeV/p );
        if( b1 > g/(c+g) )
        {
          G4double ran = G4UniformRand();
          if( ran <= c1/c )
          {
            ++p; ++h; ++(++n);
            if( G4UniformRand() <= theNucleus->GetZ()/theNucleus->GetN() )++pz;
          }
          else
          {
            if( ran <= (c2+c1)/c )
            {
              --(--n); --p; --h;
              if( pz > 0 && G4UniformRand() <= theNucleus->GetZ()/theNucleus->GetN() )--pz;
            }
          }
          goto LABEL;
        }
        for( j = 1; j < 6; ++j )gj[j] += gj[j-1];
        G4double b = g*G4UniformRand();
        for( j = 0; j < 6; ++j )
        {
          if( b <= gj[j] )
          {
            lm = j;
            break;
          }
        }
        ep1 = tkinm1( lm, p, h );
        ep2 = zj[lm];
        ep3 = 940*aj[lm];
        ind = true;
      }
      else
      {
        // equilibrium emission
        //
        rj.push_back( 0.0 );
        G4double per = arfaf( ami );
        for( i = 0; i < 7; ++i )
        {
          rj[i] <= 0 ? gj[i] = 0 : gj[i] = gameqf( i, cc[i], per, ami[i], RM );
        }
        G4double gt = 0;
        for( i = 0; i < 7; ++i )gt += gj[i];
        if( gt <= 0 )return;
        wt *= 1-gj[7]/gt;
        SetWF( 1-wt );
        for( j = 1; j < 7; ++j )gj[j] += gj[j-1];
        G4double b = G4UniformRand()*gt;
        for( j = 1; j <= 7; ++j )
        {
          if( b < gj[j] )
          {
            lm = j;
            break;
          }
        }
        if( lm == 7 )return;
        ep1 = tkin( lm, am );
        ep2 = zj[lm];
        ep3 = 940*aj[lm];
        ind = false;
      }
      ++KTOT;
      p -= aj[lm];
      n -= G4int(aj[lm]);
      pz -= ep2;
      
      parz[k].charge = ep2;
      parz[k].particleType = lm;
      if( !ind )parz[k].what = 1000.0;
      else      parz[k].what = 100.0;
      
      theNucleus->SetParameters( afj[lm], zfj[lm] );
      theNucleus->AddExcitationEnergy(
       -theNucleus->GetEnergyDeposit()+rj[lm]-ep1+vj[lm]+pevapj[lm] );
      G4double tl, cot, st, cf, sf;
      G4ThreeVector pl;
      Dvector *vtmp = new Dvector( 5 );
      vtmp->insert( vtmp->begin()+3, ep2 );
      vtmp->insert( vtmp->begin()+4, ep3/GeV );
      if( ind && (momentum.mag2()/GeV > 0.0000001) )
      {
        G4double cx = sqrt(G4UniformRand());
        G4double sx = sqrt(1-cx*cx);
        G4double fi = 2*pi*G4UniformRand();
        G4double pm = sqrt(ep1*(ep1+2*ep3));
        pl.setX( pm*sx*cos(fi) );
        pl.setY( pm*sx*sin(fi) );
        pl.setZ( pm*cx );
        G4ThreeVector ps = rotation( p12, momentum, pl );
        tl = cinema( ps, momentum, pl, cot, st, cf, sf, ep3 );
        vtmp->insert( vtmp->begin(), st );
        vtmp->insert( vtmp->begin()+1, cot );
        vtmp->insert( vtmp->begin()+2, tl/GeV );
      }
      else
      {
        G4double angl[4];
        //
        // choose isotropic distributed angle for particle emitted
        //
        angl[0] = 1-2*G4UniformRand();
        angl[3] = sqrt(1-angl[0]*angl[0]);
        G4double f = 2*pi*G4UniformRand();
        angl[1] = sin(f);
        angl[2] = cos(f);
        //
        G4double pm = sqrt(ep1*(ep1+2*ep3));
        G4ThreeVector ps( pm*angl[3]*angl[2], pm*angl[3]*angl[1], pm*angl[0] );
        if( momentum.mag2() <= 0.0001 )
        {
          vtmp->insert( vtmp->begin(), angl[3] );
          vtmp->insert( vtmp->begin()+1, angl[0] );
          vtmp->insert( vtmp->begin()+2, ep1/GeV );
          cot = angl[0];
          st = angl[2];
          cf = angl[1];
          sf = angl[0];
          pl = ps;
        }
        else
        {
          tl = cinema( ps, momentum, pl, cot, st, cf, sf, ep3 );
          vtmp->insert( vtmp->begin(), st );
          vtmp->insert( vtmp->begin()+1, cot );
          vtmp->insert( vtmp->begin()+2, tl/GeV );
        }
      }
      spt.push_back( vtmp );
      
      theNucleus->GetMomentum().setX( theNucleus->GetMomentum().x()-pl.x()*GeV );
      theNucleus->GetMomentum().setY( theNucleus->GetMomentum().y()-pl.y()*GeV );
      theNucleus->GetMomentum().setZ( theNucleus->GetMomentum().z()-pl.z()*GeV );
      
      numberOfProtons = G4int(theNucleus->GetZ()+0.5);
      numberOfNucleons = G4int(theNucleus->GetN()+0.5);
      numberOfNeutrons = numberOfNucleons - numberOfProtons;
      nuclearMass = protonMass*numberOfProtons + neutronMass*numberOfNeutrons;
      P = theNucleus->GetMomentum().mag()/GeV;
      totalEnergy = sqrt( P*P + nuclearMass*nuclearMass );
      momentum = theNucleus->GetMomentum() * (1/totalEnergy/GeV);

      parz[k].angle1 = atan2(st,cot);
      G4double fi1 = atan2(sf,cf);
      if( fi1 < 0 )fi1 = pi-fi1;
      parz[k].angle2 = fi1;
      parz[k].energy = (*vtmp)[2];
      G4double almax = 0.219327*RM*(pow(theNucleus->GetN(),oneThird)+pow(aj[lm],oneThird))*
                      sqrt(theNucleus->GetN()*aj[lm]*(ep1-vj[lm])/(theNucleus->GetN()+aj[lm]));
      G4double alp = almax * sqrt(G4UniformRand());
      uglmom.setX( uglmom.x() - alp*st*cf );
      uglmom.setY( uglmom.y() - alp*st*sf );
      uglmom.setZ( uglmom.z() - alp*cot );
    }
    // write(16,49)u,a,z
    // format('massiv spt exceeded after evaporation',40x,'u=',f10.5,' a=',f5.1,' z=',f4.1)
    return;
  }
 
 void G4PreEquilibrium::mashnk( G4int &icem, const G4double anucc, const G4double znucc,
                                const G4double ecno, const G4double tipno,
                                const G4double ehipi, const G4double ehin,
                                const G4double ehicut )
  {
    // this routine is called by cascade
    //
    // ehicem is the value of max energy for which siginb is calculated for ifircl > 0
    // siginb is inelastic cross-section for either bert or bertcem
    //
    // icem is trigger for use of bertcem of s.mashik cem code.
    // if imash == 0, bertcem is not used and icem = 0
    // if imash == 1, ehicem is max energy for which icem = 1 and bertcem is called
    //
    icem = 0;
    G4double ehicem;     // ehicem is in common/callcm/ in cascad.f
    G4int mat, icm;      // mat is in COMON.F, icm is in common/callcm/ in cascad.f
    G4double ehic[7];    // ehic is in common/callcm/ in cascad.f
    G4int imash = 0, ihecc = 0;
    
    if( imash == 0 )     // imash is in common/callcm/ in cascad.f
    {
      if( ihecc == 0 )   // ihecc in in common/hecc1/ in cascad.f
      {
        icem = 0;
        ehicem = ehicut;
        icm = icem;
        ehic[mat-1] = ehicem;
        return;
      }
      else if( ihecc == 1 )
      {
        icem = 0;
        ehicem = 10000;
        icm = icem; 
        ehic[mat-1] = ehicem;
        return;
      }
    }
    if( anucc < 14 )     ehicem = 0;
    else if( anucc < 17 )ehicem = 1000;
    else
    {
      // icem = 0 for oxygen and nitrogen if ecno > ehicem
      // ehicem can be set by trial and error for this case
      //
      // next statement is tentative
      //
      tipno <= 1 ? ehicem = ehin : ehicem = ehipi;
      //
      // icem = 0 for energies > ehicem, which can be set at will
      // may want to set it lower if ihecc == 1
      //
      ecno <= ehicem ? icem = 1 : icem = 0; // if icem == 1, bertcem is called
    }
    icm = icem; 
    ehic[mat-1] = ehicem;
    return;
  }
 
 G4double G4PreEquilibrium::bf( const G4double a, const G4double z, const G4double e )
  {
    const G4double xx1[51] =
    { 0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,
      0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
      0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
      0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
      0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.0 };
    //
    // b(x) from Krappe, Nix & Sierk 1979 for a/r0 = 0.5508
    //
    const G4double yy2[51] =
    { 0.1495, 0.1475, 0.1450, 0.1425, 0.1400, 0.1375, 0.1350, 0.1328, 0.1298, 0.1270,
      0.1245, 0.1215, 0.1185, 0.1150, 0.1118, 0.1083, 0.1040, 0.0997, 0.0949, 0.0900,
      0.0849, 0.0795, 0.0739, 0.0685, 0.0630, 0.0578, 0.0520, 0.0460, 0.0408, 0.0342,
      0.0289, 0.0236, 0.0188, 0.0146, 0.0111, 0.0082, 0.0059, 0.0043, 0.0029, 0.0020,
      0.0017, 0.0012, 0.0010, 0.0007, 0.0004, 0.0002, 0.0001, 0.0000, 0.0000, 0.0000, 0.0 };
    //
    const G4double r0m = 1.2;
    //
    // Yukawa-plus-exponential macroscopic model of Krappe, Nix & Sierk (1979)
    //
    const G4double gamma = 3;
    const G4double a2 = 21.7;
    const G4double a3 = 0.7322;
    G4double a3rt = pow(a,oneThird);
    G4double sufnuc = a2*(1-gamma*(1-2*z/a)*(1-2*z/a))*a3rt*a3rt;
    G4double x = 0.5*a3*z*z/a3rt/sufnuc;
    if ( x >= 1 )return 0.0;
    //
    // Cameron (Can.J.Phys.35(1957)1021) shell and pairing corr. for g.s. mass
    // calculate saddle-point shell and/or pairing corrections
    //
    return sufnuc*quadraticInterpolation( x, xx1, yy2, 50 ) -
        GetT1Y(G4int(z)-1) - GetT2XY(G4int(a-z)-1);
  }

 G4double G4PreEquilibrium::quadraticInterpolation( const G4double u, const G4double *e,
                                                    const G4double *f, const G4int n )
  {
    G4double x1, x2, x3, y1, y2, y3;
    if( u <= e[0] )
    {
      x1 = e[0];
      x2 = e[1];
      x3 = e[2];
      y1 = f[0];
      y2 = f[1];
      y3 = f[2];
    }
    else
    {
      if( u >= e[n-1] )
      {
        if( u > e[n] )return f[n];
        x1 = e[n-2];
        x2 = e[n-1];
        x3 = e[n];
        y1 = f[n-2];
        y2 = f[n-1];
        y3 = f[n];
      }
      else
      {
        for( G4int j = 0; j < n; ++j )
        {
          if( u < e[j] )
          {
            x1 = e[j-1];
            x2 = e[j];
            x3 = e[j+1];
            y1 = f[j-1];
            y2 = f[j];
            y3 = f[j+1];
            break;
          }
        }
      }
    }
    return y1*(u-x2)*(u-x3)/((x1-x2)*(x1-x3)) + y2*(u-x1)*(u-x3)/((x2-x1)*(x2-x3)) + 
           y3*(u-x1)*(u-x2)/((x3-x1)*(x3-x2));
  }
 
 G4double G4PreEquilibrium::massDefect( const G4double x, const G4double y ) const
  {
    // calculation of mass defect
    //
    G4double c = pow(x,oneThird);
    G4double d = pow(x,twoThirds);
    G4double e = pow(x,fourThirds);
    G4double f = pow(y,fourThirds);
    G4double a = 1 - 0.62025/d;
    G4double b = (x-2*y)/x;
    G4double es = (25.8357-44.2355*b*b)*a*a*d;
    G4double ec = 0.779*y*(y-1)/c*( 1 - 1.5849/d + 1.2273/x + 1.5772/e );
    G4double ealfa = -0.4323*f/c*( 1 + 0.49597/x - 0.57811/c - 0.14518/d );
    G4double edob = 8.367*x + 31.4506*x*b*b - 0.783*y - 17.0354*x;
    G4int i = G4int(x);
    G4int j = G4int(y);
    return es + ec + ealfa + edob + GetT1Y(j-1) + GetT2XY(i-j-1);
  }
 
 G4double G4PreEquilibrium::gamagu( const G4int j )
  {
    //           particle b:   n      p       d       t      He-3   He-4
    const G4double  aj[6] = { 1.0,   1.0,    2.0,    3.0,    3.0,   4.0   }; // A(b)
    const G4double w[8] = { 0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
                            0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363 };
    const G4double fiks[8] = {  0.9602898565,  0.7966664774,  0.5255324099,  0.1834346425,
                               -0.1834346425, -0.5255324099, -0.7966664774, -0.9602898565 };
    //
    G4double a = 0;
    if( j > 0 )a = vj[j];
    G4double b = theNucleus->GetEnergyDeposit() - bj[j];
    G4double y = 0;
    for( G4int k = 0; k < 8; ++k )
    {
      G4double e = 0.5*( (b-a)*fiks[k] + b + a );
      G4double result;
      if( j > 1 )
        result = gb[j]*r0j[j]*0.104/
            (GetR0()*pow(afj[j],oneThird)*sqrt(aj[j]*theNucleus->GetEnergyDeposit()))*
            alj[j]*((e-vj[j])/theNucleus->GetEnergyDeposit())*
            pow(abs((e+bj[j])/theNucleus->GetEnergyDeposit()),aj[j]-1.5)*
            pow(abs(1-(e+bj[j])/theNucleus->GetEnergyDeposit()),exn-aj[j]-1);
      else if( j == 1 )
        result = 0.000234*GetR0()*GetR0()*pow(afj[1],twoThirds)*r0j[1]*
                 alj[1]/(GetAC()*theNucleus->GetEnergyDeposit()*afj[1])*
                 pow(1-(e+bj[1])/theNucleus->GetEnergyDeposit(),exn-2)*(e-vj[1]);
      else // j must be 0
        result = 0.000234*GetR0()*GetR0()*pow(afj[0],twoThirds)*r0j[0]*
                 alj[0]/(GetAC()*theNucleus->GetEnergyDeposit()*afj[0])*
                 pow(1-(e+bj[0])/theNucleus->GetEnergyDeposit(),exn-2)*(e+GetBN());
      y += w[k]*(b-a)*result/2;
    }
    return y;
  }
 
 G4double G4PreEquilibrium::tkin( const G4int i, const G4double am )
  {
    //
    // kinetic energy for particles in equilibrium decay
    //
    G4double rk, frk;
    G4double rb = 4*am*afj[i]*rj[i];
    do
    {
      G4double b1 = G4UniformRand();
      rk = 1 + 1/sqrt(rb)*log(b1+(1-b1)*exp(-sqrt(rb)));
      if( i==0 )
      {
        G4double q1 = 1 + (2.12/pow(afj[0],twoThirds)-0.05)/
            (0.76+2.2/pow(afj[0],oneThird))/rj[0];
        frk = 3*sqrt(3.)/2/(q1*sqrt(q1))*(q1*rk-rk*rk*rk);
      }
      else
      {
        frk = 3*sqrt(3.)/2*(rk-rk*rk*rk);
      }
    } while( G4UniformRand() > frk );
    return rj[i]*(1-rk*rk) + vj[i];
  }
 
 void G4PreEquilibrium::transitionRates( const G4double p, const G4double h,
                                         G4double &c1, G4double &c2, G4double &c3 )
  {
    // calculation of transition rates (was trncem)
    //
    G4double est = 1.6*45+theNucleus->GetEnergyDeposit()/(p+h);
    G4double b = sqrt(2*est/940);
    G4double sf = (10.63/(b*b)-29.93/b+42.9+34.10/(b*b)-82.20/b+82.2)/2;
    G4double t;
    if( 45/est <= 0.5 )
      t = 1-7*45/est/5;
    else
      t = 1-7*45/est/5+0.4*45/est*pow(2-1/45*est,2.5);
    c1 = 0.00332*sf*t*sqrt(est)/pow(1.2+1/(4.7*b),3);
    c2 = (c1*p*h*(p+h+1)*(p+h-2))/(GetAC()*theNucleus->GetN()*theNucleus->GetEnergyDeposit()*
                                   GetAC()*theNucleus->GetN()*theNucleus->GetEnergyDeposit());
    c3 = c1*(p+h+1)*(p*(p+1)+4*p*h+h*(h-1))/((p+h)*GetAC()*theNucleus->GetN()*
                                             theNucleus->GetEnergyDeposit());
    if( c2 < 0 )c2 = 0;
    return;
  }
 
 G4double G4PreEquilibrium::tkinm1( const G4int j, const G4double p, const G4double h )
  {
    // kinetic energy for particles in pre-equilibrium decay
    //
    //           particle b:   n      p       d       t      He-3   He-4
    const G4double  aj[6] = { 1.0,   1.0,    2.0,    3.0,    3.0,   4.0   }; // A(b)
    //
    G4double e, e1, dj, t3;
    j == 0 ?
        dj = (2.12/pow(afj[0],twoThirds)-0.05)/(0.76+2.2/pow(afj[0],oneThird)) :
        dj = -vj[j];
    G4double t = p + h - aj[j] - 1;
    G4double r2 = rj[j];
    G4double r1 = r2 + vj[j];
    if( j < 2 )
    {
      if( t <= -0.01 )return r1;
      if( t <= 0.1 )
      {
        G4double result;
        j <= 1 ?
            result = -dj+sqrt(dj*dj+(G4UniformRand()*(r2*r2+2*dj*r2))) :
            result = sqrt(G4UniformRand())*r2+vj[j];
        return result;
      }
      e1 = (r1-dj*t)/(t+1);
      do
      {
        e = vj[j] + G4UniformRand()*r2;
        t3 = (e+dj)/(e1+dj)*pow(abs((r1-e)/(r1-e1)),t);
      }
      while ( G4UniformRand() > t3 );
      return e;
    }
    if( t <= -0.1 )return r1;
    if( t <= 0.1 )
    {
      do
      {
        e = vj[j] + G4UniformRand()*r2;
        t3 = pow( abs((theNucleus->GetEnergyDeposit()-r1+e)/
                      (theNucleus->GetEnergyDeposit()-r1+r1)), aj[j]-1.5 ) *
            (e+dj)/(r1+dj);
      }
      while ( G4UniformRand() > t3 );
      return e;
    }
    G4double es = theNucleus->GetEnergyDeposit()*(aj[j]-0.5) +
        (theNucleus->GetEnergyDeposit()-r2)*(p+h-2.5);
    e1 = (es+sqrt(es*es-(theNucleus->GetEnergyDeposit()-r2)*(aj[j]-1.5)*
                  (p+h-1.5)*4.*theNucleus->GetEnergyDeposit()))/
        ((p+h-1.5)*2)-theNucleus->GetEnergyDeposit()+r1;
    do
    {
      e = vj[j] + G4UniformRand()*r2;
      if( e1+0.001 > r1 )break;
      t3 = pow( abs((theNucleus->GetEnergyDeposit()-r1+e)/
                    (theNucleus->GetEnergyDeposit()-r1+e1)), aj[j]-1.5 ) *
           (e+dj)/(e1+dj) * pow( abs((r1-e)/(r1-e1)), t );
    }
    while ( G4UniformRand() > t3 );
    return e;
  }
 
 void G4PreEquilibrium::lpoly( const G4double x, const G4int n, G4double *pl )
  {
    // This subroutine calculates the ordinary Legendre polynomials of
    // order 0 to n-1 of argument  x  and stores them in the vector pl.
    // They are calculated by recursion relation from the first two polynomials.
    //
    // written by A. J. Sierk   LANL  T-9  February,1984
    //
    pl[0] = 1;
    pl[1] = x;
    for( G4int i = 2; i < n; ++i )pl[i] = ((2*i-3)*x*pl[i-1]-(i-2)*pl[i-2])/(i-1);
    return;
  }
 
 G4double G4PreEquilibrium::barfit( const G4double a, const G4double z, const G4int il )
  {
    // This subroutine returns the barrier height, bfis.
    // Arguments:  z, the atomic number
    //             a, the atomic mass number
    //            il, the angular momentum in units of h-bar,
    //                where h-bar is Plancks constant divided by 2*pi
    //
    // The fission barrier for il = 0 is calculated from a 7th order
    // fit in two variables to 638 calculated fission barriers for z values
    // from 20 to 110.  These 638 barriers are fit with an rms deviation of
    // 0.10 MeV by this 49-parameter function.
    // If barfit is called with (z,a) values outside the range of the fit,
    // the barrier height is set to 0
    //
    // For il values not equal to zero, the values of il at which the barrier
    // is 80% and 20% of the il = 0 value are respectively fit to 20-parameter
    // functions of z and a, over a more restricted range of a values, than is
    // the case for il = 0.  The value of il where the barrier disappears, lmax,
    // for 61 nuclei, is fit to a 35-parameter function of z and a, with the
    // same range of z and a values as il-80  and  il-20.
    // Once again, if a (z,a) pair is outside of the range of validity of the
    // fit, the barrier value is set to 0.  These three values
    // (bfis(il=0), il-80, and l-20) and the constraints of bfis = 0 and
    // d(bfis)/dl = 0 at il = lmax and il = 0 lead to a fifth-order fit to
    // bfis(il) for il > l-20.  the first three constraints lead to a third-order
    // fit for the region il < il-20.
    //
    // The ground-state energies are calculated from a 175-parameter
    // fit in z, a, and il to 329 ground-state energies for 36 different
    // z and a values.
    // (the range of z and a is the same as for l-80, l-20, and l-max)
    //
    // The calculated barriers from which the fits were made were calculated
    // in 1983-1985 by A. J. Sierk of Los Alamos National Laboratory, group t-9,
    // using Yukawa-plus-exponential double folded nuclear energy, exact Couloub
    // diffuseness corrections, and diffuse-matter moments of inertia.  The
    // parameters of the model are those derived by Moller and Nix in 1979:
    // r-0 = 1.16 fm, as = 21.13 MeV, kappa-s = 2.3  a = 0.68 fm.
    // The diffuseness of the matter and charge distributions used
    // corresponds to a surface diffuseness parameter (defined by Myers)
    // of 0.99 fm.  The calculated barriers for il = 0 are
    // accurate to a little less than 0.1 MeV; the output from this
    // routine is a little less accurate.  Worst errors may be as large
    // as 0.5 MeV; characteristic uncertainty is in the range of 0.1-0.2
    // MeV.  The values of egs are generally approximated to within
    // about 0.1-0.2 MeV;  the largest deviation is about 0.5 MeV,
    // near il-i for light nuclei.
    //
    // The rms deviation of lmax from the 61 input values is 0.31
    // h-bar.  the approximate value is nearly always within
    // 0.5 h-bar of the calculated one.
    //
    // Below is a table of test values to check implementation of the program
    //  z, a,  l    egnd st  fiss bar      moments of inertia     lmax
    //
    //  28, 58, 0    0.00     33.14        0.816 3.603 3.608      46.1
    //     ,25   21.36     19.50        0.778 3.662 3.662      46.1
    //     ,40   49.66      2.97        0.724 3.648 3.650      46.1
    //     ,46.1 59.14      0.00        0.746 3.160 3.160      46.1
    //  65,153, 0    0.00     28.88        0.621 3.698 3.698      82.3
    //     ,50   19.00     16.16        0.615 3.639 3.639      82.3
    //     ,80   45.24      0.26        0.616 2.765 2.788      82.3
    //     ,82.3 47.04      0.00        0.682 2.231 2.276      82.3
    //  93,229, 0    0.00      3.76        0.715 1.747 1.747      68.1
    //     ,45    8.21      1.26        0.765 1.578 1.578      68.1
    //     ,68.1 17.96      0.00        1.053 1.053 1.236      68.1
    //
    // written by A. J. Sierk,  LANL  t-9
    // version 1.0   February, 1984
    // version 1.1   January, 1985  improved coefficients in egs and lmax
    // version 1.2   September, 1985  improved lmax, egs coefficients
    // version 1.21  June, 1986   minor changes made
    //
    //     Copyright, 1986,  the regents of the University of California.
    //     This software was produced under a U. s. government contract
    //     (w-7405-eng-36) by the Los Alamos National Laboratory, which is
    //     operated by the University of California for the U. S. Department
    //     of Energy.  The U. S. government is licensed to use, reproduce,
    //     and distribute this software.  Permission is granted to the public
    //     to copy and use this software without charge, provided that this
    //     notice and any statement of authorship are reproduced on all
    //     copies.  Neither the government nor the University makes any
    //     warranty, expressed or implied, or assumes any liability
    //     or responsibility for the use of this software.
    //
    G4double pa[7], pz[7], pl[10];
    //
    const G4double emncof[4][5] =
    { { -901.100, -1408.18,  2770.00, -706.695,  889.867 },
      {  13535.5, -20384.7,  10938.4, -4862.97, -618.603 },
      { -3263.67,  1624.47,  1368.56,  1317.31,  153.372 },
      {  7488.63, -12158.1,  5502.81, -1336.30,  0.0505367 } };
    const G4double elmcof[4][5] = 
    { {  1845.42, -5640.02,  5667.30, -3151.50,  954.160 },
      { -2245.77,  8561.33, -9673.48,  5817.44, -1869.97 },
      {  2797.72, -8730.73,  9197.06, -4919.00,  1372.83 },
      { -30.1866,  1411.61, -2859.19,  2130.16, -649.072 } };
    const G4double emxcof[5][7] = {
     {-4106527.32, 10006494.7,-10953375.1, 7847972.52,-3785749.26, 1122379.45,-177561.170},
     {10876333.0,-26375824.5, 28547240.0,-20110746.7, 9483736.41,-2734385.28, 413247.256},
     {-8765309.03, 21425051.3,-23579959.5, 17016134.7,-8237381.90, 2424479.57,-365427.239},
     {6302589.54,-15299900.4, 16564020.0,-11669577.6, 5473691.53,-1549863.42, 215409.246},
     {-1455398.91, 3649618.35,-4212674.23, 3243125.55,-1679279.04, 523795.062,-76657.6599}};
    const G4double elzcof[7][7] = {
     {511819.909,-1303031.86, 1901198.70,-1206282.42, 568208.488, 54834.6483,-24588.3052},
     {-1132694.53, 2977645.90,-4543263.26, 3004648.70,-1449892.74,-102026.610, 62795.9815},
     {1375433.04,-3658089.88, 5477989.99,-3781092.83, 1841317.65, 15366.9695,-69681.7834},
     {-856559.835, 2488722.66,-4073491.28, 3128358.99,-1623940.90, 119797.378, 42573.7058},
     {328723.311,-1098921.75, 2039972.69,-1771857.18, 996051.545,-153305.699,-11298.2954},
     {41585.0238, 72965.3408,-493776.346, 601254.680,-401308.292, 96596.8391,-349960.27},
     {-182751.044, 391386.300,-303639.248, 115782.417,-424992.80,-611772.47, 36682.647}};
    G4double bfis = 0;
    if( G4int(z) < 19 || G4int(z) > 111 )
    {
      // write(16,1000)
      // 1000 format(/10x,'*  *  *  *  barfit called with  z  less than 19 or ',
      //  ' greater than 111.  bfis is set to 0.0.  *  *  *  *')
      return bfis;
    }
    if( G4int(z) > 102 && il > 0 )
    {
      // write(16,1010)
      // 1010 format(/10x,'*  *  *  *  barfit called with  z  greater than 102',
      // ' and  l  not equal to zero.  bfis is set to 0.0.  *  *  *  *')
      return bfis;
    }
    if( a < 1.2*z+0.01*z*z || a > 5.8*z-0.024*z*z )
    {
      // write(16,1020)a
      // 1020 format(/10x,'*  *  *  *  barfit called with  a  = ',i3,', outside ',
      // 'the allowed values for z = ',i3,' *  *  *  *')
      return bfis;
    }
    lpoly( 0.01*z, 7, pz );
    lpoly( 0.0025*a, 7, pa );
    G4int i, j;
    for( i = 0; i < 7; ++i )
    {
      for( j = 0; j < 7; ++j )bfis += elzcof[j][i]*pz[j]*pa[i];
    }
    if( (a < 1.4*z+0.009*z*z-5 || a > 20+3*z+10) && il > 0 )
    {
      // write(16,1030)a,il
      // 1030 format(/10x,'*  *  *  *  barfit called with  a   = ',i3,', outside',
      // ' the allowed values for z = ',i3/26x,'for nonzero  l  = ',i3,
      // '  *  *  *  *')
      bfis = 0;
      return bfis;
    }
    G4double el80 = 0;
    G4double el20 = 0;
    G4double elmax = 0;
    for( i = 0; i < 4; ++i )
    {
      for( j = 0; j < 5; ++j )
      {
        el80 += elmcof[j][i]*pz[j]*pa[i];
        el20 += emncof[j][i]*pz[j]*pa[i];
      }
    }
    for( i = 0; i < 5; ++i )
    {
      for( j = 0; j < 7; ++j )elmax += emxcof[j][i]*pz[j]*pa[i];
    }
    if( il < 1 )return bfis;
    G4double x = el20/elmax;
    G4double y = el80/elmax;
    if( il <= el20 )
    {
      G4double q = 0.2/(el20*el20*el80*el80*(el20-el80));
      bfis *=
          1+q*(4*el80*el80*el80-el20*el20*el20)*il*il-q*(4*el80*el80-el20*el20)*il*il*il;
    }
    else
    {
      G4double aj = (-20*pow(x,5)+25*pow(x,4)-4)*(y-1)*(y-1)*y*y;
      G4double ak = (-20*pow(y,5)+25*pow(y,4)-1)*(x-1)*(x-1)*x*x;
      G4double q = 0.2/((y-x)*((1-x)*(1-y)*x*y)*((1-x)*(1-y)*x*y));
      G4double z = il/elmax;
      bfis *= 4*pow(z,5)-5*pow(z,4)+1+
          (z-1)*(q*(aj*y-ak*x)*(2*z+1)-q*(aj*(2*y+1)-ak*(2*x+1))*z)*z*z*(z-1);
    }
    if( bfis <= 0 || il > elmax )bfis = 0;
    return bfis;
  }
 
 G4bool G4PreEquilibrium::cascem( G4ThreeVector amnucl )
  {
    const G4int N = 7;
    const G4double EPS = 0.007;
    const G4double VPI = 0.025;
    const G4double CM0 = 0.94;
    const G4int ME0 = 1;
    const G4int MQ0 = 1;
    
    G4int mv = 0;
    G4ThreeVector v;
    G4double tin1, sigp, sign, sigabs;
    
    SetINDI( false );
    SetING( 0 );
    exitons.protons = 0;
    exitons.neutrons = 0;
    exitons.hols = 0;
    //
    // calculation of entering point of particle in nucleus
    //
    SetT1( 0 );
    SetT2( 0 );
    G4double temp1 = G4UniformRand();
    G4double temp2 = 2*pi*G4UniformRand();
    //
    G4ThreeVector tmp3vec(  rbig[N] * sqrt(temp1) * cos(temp2),
                            rbig[N] * sqrt(temp1) * sin(temp2),
                           -rbig[N] * sqrt(1-temp1) );
    MASTRUCT p = { tmp3vec, { 0.0, 1.0, 0.0, 1.0, GetT0(), CM0 } };
    //
    IPSTRUCT ip = { ME0, 0, 0, MQ0, N };
    //
    // the following appears to be a bug in the original cem.f
    // T1 is set to 0 just above
    //
    if( GetT0() == 0 )p.array[2] = GetT1()*pow(GetT2()/GetT1(),G4UniformRand());
    const G4double p0 = sqrt(p.array[2]*(p.array[2]+2*p.array[5]));  // not the data member p0
    //
    G4ThreeVector am0( p0*(p.momentum.z()*p.array[0]*p.array[2]-p.momentum.y()*p.array[1]),
                       p0*(p.momentum.x()*p.array[1]-p.momentum.z()*p.array[0]*p.array[3]),
                       p0*(p.momentum.y()*p.array[0]*p.array[3]-p.momentum.x()*p.array[0]*p.array[2]) );
    //
    const G4double t3 = p.array[4];
    p.array[4] += poten( N, ip );
    theNucleus->AddExcitationEnergy( -theNucleus->GetEnergyDeposit() + t3 + MQ0*EPS );
    G4double atwght = aNucl+MQ0;
    G4double charge = zNucl+ME0;
    amnucl = am0;
    const G4double obr = af[6]/2;
 LABEL:
    G4double cutof1 = ip.nuclearZone*( poten(ip.ip4,ip) + ip.proton*obr ) + 0.001;
    if( p.array[4] <= cutof1 )
    {
      exitons.protons += ip.proton;
      exitons.neutrons += 1-ip.proton;
      cascem1( amnucl, p, ip, mv );
      goto LABEL;
    }
    G4int ipe[5];
    G4double pe[9];
    if( pointe( p, ip, pe, ipe, v, u, tin1, sigp, sign, sigabs, t3 ) )
    {
      G4int np;
      if( ip.nuclearZone != 0 )
      {
        G4double cutof3 = cutof1-ip.nuclearZone*obr+0.05;
        if( p.array[4] <= cutof3 )
        {
          G4double r = rsm[N]*p.momentum.mag();
          if( rsm[N] >= r+0.001 )
          {
            G4double wm = 0;
            for( G4int i = 0; i < 50; ++i )
              wm += wim( 1, p, ip ) + wim( 0, p, ip );
            G4double w3 = wopt( p.array[4]-cutof3+0.05, r, ip.proton );
            if( abs((wm/100-w3)/w3) >= 0.3 )
            {
              typint( p, ip, pe, ipe, v, u, tin1, sigp, sign, sigabs, mv, np );
              if( np <= 0 )
              {
                exitons.protons += ip.proton;
                exitons.neutrons += 1-ip.proton;
                cascem1( amnucl, p, ip, mv );
              }
              else
              {
                pauliPrinciple( p, ip, v, mv, np, 0 );
                exitons.protons += ip.proton;
                exitons.neutrons += 1-ip.proton;
                cascem1( amnucl, p, ip, mv );
              }
              goto LABEL;
            }
          }
        }
      }
      typint( p, ip, pe, ipe, v, u, tin1, sigp, sign, sigabs, mv, np );
      if( np <= 0 )
      {
        exitons.protons += ip.proton;
        exitons.neutrons += 1-ip.proton;
        cascem1( amnucl, p, ip, mv );
      }
      else
        pauliPrinciple( p, ip, v, mv, np, 1 );
      goto LABEL;
    }
    ip.proton <= 0 ? temp1 = 0 : temp1 = ip.proton;
    G4double cutof2 = poten( ip.ip4, ip ) + temp1*obr + 0.001;
    if( p.array[4] <= cutof2 )
    {
      exitons.protons += ip.proton;
      exitons.neutrons += 1-ip.proton;
      cascem1( amnucl, p, ip, mv );
      goto LABEL;
    }
    theNucleus->AddExcitationEnergy(
     -(p.array[4]+ip.nuclearZone*EPS+(1-ip.nuclearZone)*(0.14+VPI)*(1-ip.ip1)) );
    atwght -= ip.nuclearZone;
    charge -= ip.proton;
    temp1 = sqrt(p.array[4]*(p.array[4]+2*p.array[5]));
    amnucl.setX( amnucl.x()-temp1*(p.momentum.z()*p.array[0]*p.array[2]-p.momentum.y()*p.array[1]) );
    amnucl.setY( amnucl.y()-temp1*(p.momentum.x()*p.array[1]-p.momentum.z()*p.array[0]*p.array[3]) );
    amnucl.setZ(
     amnucl.z()-temp1*(p.momentum.y()*p.array[0]*p.array[3]-p.momentum.x()*p.array[0]*p.array[2]) );
    if( theNucleus->GetEnergyDeposit() <= 0.0001 &&
        atwght == aNucl && charge == zNucl )return true;
    if( theNucleus->GetEnergyDeposit() < 0 )return true;
    if( KTOT > 100 )
    {
      //write(16,27)
      //27 format (40x,32hmassiv spt exceeded after cascad)
      return false;
    }
    Dvector vtmp( 5 );
    vtmp[0] = p.array[0];
    vtmp[1] = p.array[1];
    vtmp[2] = p.array[4];
    vtmp[3] = ip.proton;
    vtmp[4] = p.array[5];
    spt[KTOT]->swap( vtmp );

    PARZS *pztmp = new PARZS;
    pztmp->charge = ip.proton;
    pztmp->what = GetING();
    if( p.array[5] <= 0.5 )
    {
      if( ip.proton < 0 )
        pztmp->particleType = 7;
      else if( ip.proton == 0 )
        pztmp->particleType = 8;
      else
        pztmp->particleType = 9;
    }
    else
    {
      if( ip.proton <= 0.099 )pztmp->particleType = 1;  // neutron
      else                    pztmp->particleType = 2;  // proton
      if( GetINDI() )pztmp->what = -pztmp->what;
    }
    pztmp->angle1 = atan2( p.array[0], p.array[1] );
    temp1 = atan2( p.array[2], p.array[3] );
    if( temp1 < 0 )temp1 = pi - temp1;
    pztmp->angle2 = temp1;
    pztmp->energy = p.array[4];
    parz.push_back( *pztmp );

    ++KTOT;
    cascem1( amnucl, p, ip, mv );
    goto LABEL;
  }
 
 void G4PreEquilibrium::cascem1( G4ThreeVector a, MASTRUCT &p, IPSTRUCT &ip, G4int &m )
  {
    const G4int N = 7;
    if( m <= 0 )
    {
      G4double temp = 5.07*rsm[N];
      a *= temp;
    }
    else
    {
      G4int i;
      p.momentum = pmemo[m].momentum;
      for( i = 0; i < 6; ++i )p.array[i] = pmemo[m].array[i];
      ip.proton = imemo[m].proton;
      ip.ip1 = imemo[m].ip1;
      ip.ip2 = imemo[m].ip2;
      ip.nuclearZone = imemo[m].nuclearZone;
      ip.ip4 = imemo[m].ip4;
      SetING( ngen[m--] );
    }
    return;
  }
 
 G4bool G4PreEquilibrium::pointe( MASTRUCT &p, IPSTRUCT &ip,
                                  G4double *pe, G4int *ipe,
                                  G4ThreeVector v, G4double u,
                                  G4double tin1, G4double sigp,
                                  G4double sign, G4double sigabs, const G4double t3 )
  {
    // determine interaction point inside nucleus
    //
    const G4int N = 7;
    
    G4bool nout;
    const G4double sk = 3;
    G4double temp1 = log(G4UniformRand());
    do
    {
      G4double s = geometricalParticlePath( p, ip );
      do
      {
        partnerSelection( p, ip, pe, ipe );
        //
        // calculation of t, v, u
        //
        G4double pin = sqrt( p.array[4]*(p.array[4]+2*p.array[5]) );
        G4double pn = sqrt( pe[7]*(pe[7]+2*pe[8]) );
        G4double denom = p.array[4]+p.array[5]+pe[7]+pe[8];
        v.setX( (pin*p.array[0]*p.array[3]+pn*pe[3]*pe[6])/denom );
        v.setY( (pin*p.array[0]*p.array[2]+pn*pe[3]*pe[5])/denom );
        v.setZ( (pin*p.array[1]+pn*pe[4])/denom );
        u = sqrt( 1-v.mag2() )*denom;
        tin1 = (u*u-(p.array[5]+pe[8])*(p.array[5]+pe[8]))/(2*pe[8]);
        //
        G4int i, ms, mq, ksin, ksip, me;
        slqek( i, ms, mq, ksin, me, ip.ip1, ip.ip2, ip.nuclearZone,
               ip.proton, ipe[1], ipe[2], ipe[3], 0 );
        sign = sigmat( i, ms, mq, ksin, 0, tin1 );
        slqek( i, ms, mq, ksip, me, ip.ip1, ip.ip2, ip.nuclearZone,
               ip.proton, ipe[1], ipe[2], ipe[3], 1 );
        sigp = sigmat( i, ms, mq, ksip, 0, tin1 );
        sigabs = sigmat( i, ms, mq, ksip, 3, p.array[4] );
        G4double plambi = 10/((rhon[ip.ip4]*sign+rhop[ip.ip4]*sigp+
                               rhop[ip.ip4]*sigabs)*rsm[N-1]);
        if( -plambi*temp1-6*rhon[0]/rhon[N-1] > 0 )
        {
          p.array[4] = t3;
          return false;
        }
        G4double temp2 = sk;
        G4double deltsi = plambi/temp2;
        G4double temp5 = min( s, deltsi );
        G4double piks = temp1+temp5/plambi;
        if( piks > 0 )
        {
          G4ThreeVector tmp3vec( plambi*temp1*p.array[0]*p.array[3],
                                 plambi*temp1*p.array[0]*p.array[2],
                                 plambi*temp1*p.array[1] );
          p.momentum -= tmp3vec;
          pe[0] = p.momentum.x();
          pe[1] = p.momentum.y();
          pe[2] = p.momentum.z();
          return nout;
        }
        temp1 = piks;
        s -= temp5;
        G4ThreeVector tmp3vec( temp5*p.array[0]*p.array[3],
                               temp5*p.array[0]*p.array[2],
                               temp5*p.array[1] );
        p.momentum += tmp3vec;
      } while ( s != 0 );
      nout = refrac( N, p, ip );
    } while ( ip.ip4 < N+1 );
    G4double temp3 = geometricalParticlePath( p, ip );
    G4ThreeVector tmp3vec( temp3*p.array[0]*p.array[3],
                           temp3*p.array[0]*p.array[2],
                           temp3*p.array[1] );
    p.momentum += tmp3vec;
    nout = refrac( N, p, ip );
    return false;
  }
 
 G4double G4PreEquilibrium::geometricalParticlePath( const MASTRUCT p, const IPSTRUCT ip )
  {
    G4double rin = p.momentum.mag();
    G4double c = (p.momentum.x()*p.array[0]*p.array[3] +
                  p.momentum.y()*p.array[0]*p.array[2] +
                  p.momentum.z()*p.array[1])/rin;
    G4double g, r;
    if( c < 0 )
    {
      g = -1;
      ip.ip4 < 2 ? r = 0 : r = rbig[ip.ip4-2];
    }
    else
    {
      g = 1;
      r = rbig[ip.ip4-1];
    }
    G4double d = r*r - rin*rin*(1-c*c);
    if( d < 0 )
    {
      g = 1;
      r = rbig[ip.ip4-1];
      d = r*r - rin*rin*(1-c*c);
    }
    return g*sqrt(d) - rin*c;
  }
 
 G4bool G4PreEquilibrium::refrac( const G4int N, MASTRUCT &p, IPSTRUCT &ip ) const
  {
    // calculation of energy and direction change of
    // particle while entering another zone of the nucleus
    //
    G4int j;
    G4bool result;
    G4double tmp = p.momentum.x()*p.array[0]*p.array[3] +
                   p.momentum.y()*p.array[0]*p.array[2] + 
                   p.momentum.z()*p.array[1];
    tmp < 0 ? j = ip.ip4-1 : j = ip.ip4+1;
    p.array[4] += poten( j, ip ) - poten( ip.ip4, ip );
    ip.ip4 = j;
    j >= N+2 ? result = false : result = true;
    return result;
  }
 
 G4double G4PreEquilibrium::wopt( const G4double e1, const G4double r, const G4int i )
  {
    const G4double a = aNucl;
    const G4double z = zNucl;
    G4double x = (a-2*z)/a;
    G4double e = e1*GeV;
    G4double rm, am, wv, wsf;  // this rm is not the data member rm
    if( i <= 0 )
    {
      if( e <= 25 )
      {
        rm = 1.26*pow(a,oneThird);
        am = 0.58;
        wv = 0.22*e - 1.56;
        wsf = 13 - 0.25*e - 12*x;
      }
      else
      {
        rm = 1.21*pow(a,oneThird);
        am = 0.6448;
        wv = 0.459 + 0.111*e;
        wsf = 4.28 - 0.0414*e;
      }
    }
    else
    {
      if( e <= 25 )
      {
        rm = 1.32*pow(a,oneThird);
        am = 0.51 + 0.7*x;
        wv = 0.22*e - 2.7;
        wsf = 11.8 - 0.25*e + 12*x;
      }
      else
      {
        rm = 1.37*pow(a,oneThird);
        am = 0.74 - 0.008*e + x;
        wv = 1.2 + 0.09*e;
        wsf = 4.2 - 0.05*e + 15.5*x;
      }
    }
    if( wv < 0 )wv = 0;
    if( wsf < 0 )wsf = 0;
    G4double f = 1/(1+exp((r-rm)/am));
    return f*(wv+4*wsf*(1-f));
  }
 
 void G4PreEquilibrium::typint( MASTRUCT &p, IPSTRUCT &ip,
                                G4double *pe, G4int *ipe, G4ThreeVector v,
                                G4double u, G4double tin1, G4double sigp,
                                G4double sign, G4double sigabs, G4int mv, G4int &np )
  {
    // block of determining of interaction type and calculation
    // of secondary particles characteristics
    //
    G4bool nin = true;
    while ( nin )
    {
      G4int mtemp = ip.ip4-1;
      G4double betabs = rhop[mtemp]*sigabs/
          (rhop[mtemp]*sigp+rhon[mtemp]*sign+rhop[mtemp]*sigabs);
      if( G4UniformRand() <= betabs )
      {
        absorption( p, ip, pe, mv, np, v );
        return;
      }
      G4int i, ms, mq, ksi, me;
      slqek( i, ms, mq, ksi, me, ip.ip1, ip.ip2, ip.nuclearZone,
             ip.proton, ipe[1], ipe[2], ipe[3], ipe[0] );
      G4double betael = (sigmat( i, ms, mq, ksi, 1, tin1 )+
                        sigmat( i, ms, mq, ksi, 2, tin1 ))/
                       sigmat( i, ms, mq, ksi, 0, tin1 );
      if (G4UniformRand() <= betael )
      {
        np = elasticScattering( v, u, tin1, p, ip, ipe, mv, i, ms, mq, ksi, me );
        return;
      }
      nin = inelasticScattering( p, ip, ipe, i, ms, mq, ksi, me, v, u, tin1, mv, np );
    }
    return;
  }
 
 G4int G4PreEquilibrium::elasticScattering( const G4ThreeVector v, const G4double u,
                                            const G4double tin, const MASTRUCT p,
                                            const IPSTRUCT ip, const G4int *ipatne,
                                            const G4int mv, const G4int i,
                                            const G4int ms, const G4int mq,
                                            const G4int ksi, const G4int me )
  {
    // calculate particle characteristics in elastic and charge exchange scattering
    //
    G4double cmi, cmn;
    if( ip.ip1 != 0 )
    {
      cmi = 0.14;
      cmn = 0.94;
    }
    else
    {
      cmi = p.array[5];
      cmn = 0.94;
    }
    G4double tmp1 = sigmat( i, ms, mq, ksi, 2, tin );
    G4double tmp2 = sigmat( i, ms, mq, ksi, 1, tin );
    G4int ie, ne;
    G4double ctsti;
    if( G4UniformRand() >= tmp1/(tmp1+tmp2) )
    {
      ie = ip.proton;
      ne = ipatne[0];
      ctsti = cosel( i, mq, ksi, tin, p.array[5] );
    }
    else
    {
      if( ip.proton != 0 )
      {
        ie = 0;
        ne = me;
      }
      else
      {
        ne = 1 - ipatne[0];
        ie = me - ne;
      }
      ctsti = cosex( i, tin, p.array[5] );
    }
    G4ThreeVector pist, pnst;
    momentaCalc( p, v, u, pist, pnst, ctsti, cmi, cmn );
    //
    MASTRUCT p1 = { pist, { 0.0, 0.0, 0.0, 0.0, 0.0, cmi } };
    pmemo[mv+2] = p1;
    imemo[mv+2].proton = ie;
    imemo[mv+2].ip1 = 0;
    imemo[mv+2].ip2 = ip.ip2;
    imemo[mv+2].nuclearZone = ip.nuclearZone;
    MASTRUCT p2 = { pnst, { 0.0, 0.0, 0.0, 0.0, 0.0, cmn } };
    pmemo[mv] = p2;
    imemo[mv].proton = ne;
    imemo[mv].ip1 = 0;
    imemo[mv].ip2 = 0;
    imemo[mv].nuclearZone = 1;
    return 2;
  }
 
 G4bool G4PreEquilibrium::inelasticScattering( MASTRUCT &p, IPSTRUCT &ip, G4int *ipatne,
                                               G4int l, G4int ms, G4int mq, G4int ksi,
                                               G4int me, G4ThreeVector v, G4double u,
                                               G4double tin1, G4int mv, G4int np )
  {
    if( ip.ip1 > 0 )
    {
      G4double betais = crossSectionInterp(tin1,25)/crossSectionInterp(tin1,26);
      //
      // there is a bug in the original cem.f routine 
      // isocem uses TINI as 3rd parameter
      // it should be TIN1
      //
      if( G4UniformRand() <= betais )
        isocem( u, v, tin1, p, ipatne, mv );
      else
        statisticalModel( u, v, p, ipatne, mv );
      return true;
    }
    G4double betath;
    if( tin1 >= 4 )
      betath = 0;
    else
      betath = sigmat(l,ms,mq,ksi,7,tin1)/
          (sigmat(l,ms,mq,ksi,0,tin1)-sigmat(l,ms,mq,ksi,1,tin1)-sigmat(l,ms,mq,ksi,2,tin1));
    G4int ith;
    G4double th;
    if( G4UniformRand() < betath )
    {
      ith = 1;
      th = 1;
    }
    else
    {
      ith = 0;
      th = 2;
    }
    if( u > p.array[5]+0.14*th+0.96 )
    {
      G4int ik = 0;
      G4int kp = 1;
      while ( kp != 0 )
      {
        G4bool lp = vmnsp( p, ip, u, mv, np, ith, mq, tin1 );
        if( np == 0 )return false;
        if( lp )return true;
        direction( v, u, tin1, mq, mv, np, p, kp, ith );
        if( ++ik >= 50 )return true;
      }
      chinel( ip, l, ms, mq, ksi, np, mv, tin1, me, ipatne );
    }
    return true;
  }
 
 G4double G4PreEquilibrium::wim( const G4int i, const MASTRUCT p, const IPSTRUCT ip )
  {
    const G4double EPS = 0.007;
    const G4int N = 7;
    
    G4double r, r1, a1, fac, tf, rho;
    r = rsm[N]*p.momentum.mag();
    if( i <= 0 )
    {
      r1 = 1.26*pow(aNucl,oneThird);
      a1 = 0.58;
      fac = 1/(1+exp((r-r1)/a1));
      tf = tfn[ip.ip4-1];
      rho = rhon[0]*fac;
    }
    else
    {
      r1 = 1.32*pow(aNucl,oneThird);
      a1 = 0.51+0.7*(aNucl-2*zNucl)/aNucl;
      fac = 1/(1+exp((r-r1)/a1));
      tf = tfp[ip.ip4-1];
      rho = rhop[0]*fac;
    }
    G4double tp = tf*pow(G4UniformRand(),twoThirds);
    SetT0( p.array[5] );
    G4double z = abs(GetT0()-tf-EPS);
    G4double y = sqrt(z*(z+2*p.array[5]));
    G4double ct = 1.0 - 2*G4UniformRand();
    G4double fi = 2*pi*G4UniformRand();
    G4double pp = sqrt(tp*(tp+1.88));
    G4double et = GetT0()+p.array[5]+tp+0.94;
    G4double vx = pp*sqrt(1-ct*ct)*cos(fi)/et;
    G4double vy = pp*sqrt(1-ct*ct)*sin(fi)/et;
    G4double vz = (pp*ct+sqrt(GetT0()*(GetT0()+2*p.array[5])))/et;
    G4double t = (et*et*(1-vx*vx-vy*vy-vz*vz)-(p.array[5]+0.94)*(p.array[5]+0.94))/1.88;
    G4double b = sqrt(t*(t+2*p.array[5]))/(t+p.array[5]);
    G4double s;
    ip.proton == i ? s = 10.63/(b*b)-29.92/b+42.9 : s = 34.10/(b*b)-82.20/b+82.2;
    G4double bs = tf/GetT0();
    G4double dz = 1 - 1.4*bs;
    if( bs > 0.5 )dz += 0.4*bs*pow(2-1/bs,2.5);
    return 10*y/(z+p.array[5])*dz*rho*s;
  }
 
 void G4PreEquilibrium::pauliPrinciple( MASTRUCT &p, IPSTRUCT &ip, G4ThreeVector v,
                                        G4int mv, G4int np, const G4int ipa )
  {
    G4int mtemp, i, j = ip.ip4, cntr = 1;
 L10:
    if( np == 2 )
    {
      if( cntr == 2 )
      {
        ++cntr;
        goto L10;
      }
      else mtemp = mv+cntr-1;
    }
    else mtemp = mv+cntr-1;
    //
    G4ThreeVector pstar( pmemo[mtemp].momentum );
    //
    G4double ct, st, cfi, sfi, t;
    G4ThreeVector ptmp( p.momentum );
    t = cinema( pstar, v, ptmp, ct, st, cfi, sfi, pmemo[mtemp].array[5] );
    
    MASTRUCT tmp = { ptmp, { st, ct, sfi, cfi, t, 0.0 } };
    pmemo.insert( pmemo.begin()+mtemp, tmp );
    
    imemo[mtemp].ip4 = ip.ip4;
    ngen.insert( ngen.begin()+mtemp, GetING()+1 );
    if( imemo[mtemp].nuclearZone == 1 &&
        imemo[mtemp].ip2 == 0 &&
        pmemo[mtemp].array[4]-tfp[j]*imemo[mtemp].proton-tfn[j]*(1-imemo[mtemp].proton) <= 0 )
      return;
    if( cntr < np )
    {
      ++cntr;
      goto L10;
    }
    if( ipa == 0 )return;
    if( p.array[5] == 0.140 && np == 2 && 
        pmemo[mv].array[5] == 0.94 && pmemo[mv+2].array[5] == 0.94 )SetINDI( true );
    for( i = 0; i < 6; ++i )p.array[i] = pmemo[mv+2].array[i];
    ip.proton = imemo[mv+2].proton;
    ip.ip1 = imemo[mv+2].ip2;
    ip.ip2 = imemo[mv+2].ip2;
    ip.nuclearZone = imemo[mv+2].nuclearZone;
    SetING( ngen[mv+2] );
    ++(exitons.hols);
    if( np > 2 )
    {
      G4int ntemp = mv+np;
      for( i = 0; i < 6; ++i )pmemo[mv+2].array[i] = pmemo[ntemp].array[i];
      imemo[mv+2].proton = imemo[ntemp].proton;
      imemo[mv+2].ip1 = imemo[ntemp].ip1;
      imemo[mv+2].ip2 = imemo[ntemp].ip2;
      imemo[mv+2].nuclearZone = imemo[ntemp].nuclearZone;
    }
    mv += np-1;
    return;
  }
 
 void G4PreEquilibrium::absorption( const MASTRUCT p, const IPSTRUCT ip,
                                    G4double *pe, G4int mv, G4int &np, G4ThreeVector v )
  {
    // calculation of particle characteristics in absorption
    //
    G4double p1[9];
    G4int ip1[5];
    //
    partnerSelection( p, ip, p1, ip1 );
    G4double pn1 = sqrt( p1[7]*(p1[7]+1.88) );
    G4double pn = sqrt( pe[7]*(pe[7]+1.88) );
    G4ThreeVector paf( pn*pe[3]*pe[6]+pn1*p1[3]*p1[6],
                       pn*pe[3]*pe[5]+pn1*p1[3]*p1[5],
                       pn*pe[4]+pn1*p1[4] );
    G4double pafm = paf.mag();
    G4double ctn1, stn1, sfn1, cfn1;
    if( paf.z()*paf.z() >= pafm*pafm )
    {
      ctn1 = 1;
      stn1 = 0;
      sfn1 = 0;
      cfn1 = 1;
    }
    else
    {
      ctn1 = paf.z()/pafm;
      stn1 = sqrt(1-paf.z()/pafm*paf.z()/pafm);
      sfn1 = paf.y()/(pafm*stn1);
      cfn1 = paf.x()/(pafm*stn1);
    }
    G4double taf = sqrt(pafm*pafm+1.88*1.88)-1.88;
    //
    // calculation of t, v, u
    //
    G4double pin = sqrt( p.array[4]*(p.array[4]+2*p.array[5]) );
    pn = sqrt( taf*(taf+2*1.88) );
    G4double denom = p.array[4]+p.array[5]+taf+1.88;
    v.setX( (pin*p.array[0]*p.array[3]+pn*stn1*cfn1)/denom );
    v.setY( (pin*p.array[0]*p.array[2]+pn*stn1*sfn1)/denom );
    v.setZ( (pin*p.array[1]+pn*ctn1)/denom );
    G4double u = sqrt( 1-v.mag2() )*denom;
    G4double tin1 = (u*u-(p.array[5]+1.88)*(p.array[5]+1.88))/(2*1.88);
    //
    G4int ne, ie;
    chargeInAbsorption( ip.ip1, ip.proton, ne, ie );
    G4double ctst;
    if( ip.ip1 == 0 )
      ctst = 1-2*G4UniformRand();
    else
      tin1 <= 0.455 ? ctst = costa( 18, tin1 ) : ctst = 1-2*G4UniformRand();
    G4ThreeVector pist, pnst;
    momentaCalc( p, v, u, pist, pnst, ctst, 0.94, 0.94 );
    if( mv <= 97 )
    {
      MASTRUCT pm1 = { pist, { 0.0, 0.0, 0.0, 0.0, 0.0, p1[8] } };
      pmemo.insert( pmemo.begin()+mv+2, pm1 );
      imemo[mv+2].proton = ie;
      imemo[mv+2].ip1 = 0;
      imemo[mv+2].ip2 = 0;
      imemo[mv+2].nuclearZone = 1;
      MASTRUCT pm2 = { (-1)*pist, { 0.0, 0.0, 0.0, 0.0, 0.0, pe[8] } };
      pmemo.insert( pmemo.begin()+mv, pm2 );
      imemo[mv+2].proton = ne;
      imemo[mv+2].ip1 = 0;
      imemo[mv+2].ip2 = 0;
      imemo[mv+2].nuclearZone = 1;
      np = 2;
      return;
    }
    np = 0;
    //write(16,19)
    //19 format(45x,29h memory is exceeded in cascad)
    return;
  }
 
 void G4PreEquilibrium::chargeInAbsorption( const G4int m, const G4int ine,
                                            G4int &ne1, G4int &ne2 )
  {
    // determine charge in absorption
    //
    const G4double a = aNucl;
    const G4double z = zNucl;
    if( ine < 0 )
    {
      if( m != 0 )
      {
        ne1 = 0;
        ne2 = 0;
      }
      else
      {
        G4double temp2 = (z*(z-1))/(2*z*(a-z)+z*(z-1));
        if( G4UniformRand() <= temp2 )
        {
          ne1 = 0;
          ne2 = 1;
        }
        else
        {
          ne1 = 0;
          ne2 = 0;
        }
      }
    }
    else if( ine == 0 )
    {
      if( m != 0 )
      {
        ne1 = 0;
        ne2 = 1;
      }
      else
      {
        G4double temp3 = (2*z*(a-z))/(a*(a-1));
        if( G4UniformRand() <= temp3 )
        {
          ne1 = 0;
          ne2 = 1;
        }
        else
        {
          G4double temp4 = (z*(z-1))/(a*(a-1)) + temp3;
          if( G4UniformRand() <= temp4 )
          {
            ne1 = 1;
            ne2 = 1;
          }
          else
          {
            ne1 = 0;
            ne2 = 0;
          }
        }
      }
    }
    else
    {
      if( m != 0 )
      {
        ne1 = 1;
        ne2 = 1;
      }
      else
      {
        G4double temp1 = (z*(a-z))/(z*(a-z)+(a-z)*(a-z-1)/2);
        if( G4UniformRand() <= temp1 )
        {
          ne1 = 1;
          ne2 = 1;
        }
        else
        {
          ne1 = 0;
          ne2 = 1;
        }
      }
    }
    return;
  }
 
 void G4PreEquilibrium::momentaCalc( const MASTRUCT p, const G4ThreeVector v,
                                     const G4double u, G4ThreeVector pist,
                                     G4ThreeVector pnst, const G4double cti,
                                     const G4double cmi, const G4double cmn )
  {
    // calculate momenta of secondary particles in centre of mass system
    // for absorption and elastic scattering
    //
    G4double angle = 2*pi*G4UniformRand();
    G4double pim = sqrt((u*u+cmi*cmi-cmn*cmn)/(2*u)*(u*u+cmi*cmi-cmn*cmn)/(2*u)-cmi*cmi);
    G4double v2 = v.mag2();
    G4double temp3 = sqrt(p.array[4]*(p.array[4]+2*p.array[5]));
    G4double temp1 = (temp3*p.array[0]*p.array[3]*v.x()+temp3*p.array[1]*v.y()+
                      temp3*p.array[0]*p.array[2]*v.z())*(1/sqrt(1-v2)-1)/v2;
    G4double temp2 = (p.array[4]+p.array[5])/sqrt(1-v2);

    G4ThreeVector pins( temp3*p.array[0]*p.array[3]+v.x()*(temp1-temp2),
                        temp3*p.array[0]*p.array[2]+v.y()*(temp1-temp2),
                        temp3*p.array[1]+v.z()*(temp1-temp2) );
    G4ThreeVector pii( pim*sqrt(1-cti*cti)*cos(angle),
                       pim*sqrt(1-cti*cti)*sin(angle),
                       pim*cti );
    pist = rotation( pins, v, pii );
    pnst = pist * (-1);
    return;
  }
 
 void G4PreEquilibrium::partnerSelection( const MASTRUCT p, const IPSTRUCT ip,
                                          G4double *pe, G4int *ipe )
  {
    G4double t = 0;
    if( G4UniformRand() >= (aNucl-zNucl)/aNucl )
    {
      while( t <= 0 )
        t = tfp[ip.ip4-1]*pow(G4UniformRand(),twoThirds);
      ipe[0] = 1;
    }
    else
    {
      while( t <= 0 )
        t = tfn[ip.ip4-1]*pow(G4UniformRand(),twoThirds);
      ipe[0] = 0;
    }
    pe[4] = 1-2*G4UniformRand();
    pe[3] = sqrt(1-pe[4]*pe[4]);
    G4double phin = 2*pi*G4UniformRand();
    pe[6] = cos(phin);
    pe[5] = sin(phin);
    pe[0] = p.momentum.x();
    pe[1] = p.momentum.y();
    pe[2] = p.momentum.z();
    pe[8] = 0.940;
    pe[7] = t;
    ipe[4] = ip.ip4;
    ipe[1] = 0;
    ipe[2] = 0;
    ipe[3] = 1;
    return;
  }
 
 void G4PreEquilibrium::chinel( IPSTRUCT &ipatin, const G4int l, const G4int ms,
                                const G4int mq, const G4int ksi, const G4int np,
                                const G4int mv, const G4double tin1, const G4int me,
                                G4int *ipatne )
  {
    // determine secondary particles charges in inelastic scattering
    //
    if( np == 3 )
    {
      G4double spi0 = sigmat( l, ms, mq, ksi, 4, tin1 );
      G4double sth = sigmat( l, ms, mq, ksi, 7, tin1 );
      G4double ran = G4UniformRand();
      if( ran >= spi0/sth )
      {
        if( ran >= (spi0+sigmat(l,ms,mq,ksi,5,tin1))/sth )
        {
          imemo[mv].proton = ipatne[1] - (ipatin.ip4-1)*ipatin.ip1;
          imemo[mv+2].proton = (ipatin.ip4-1)*ipatin.ip1*ipatin.ip1 - ipatin.ip4*ipatin.ip1 + 1;
          imemo[mv+1].proton = me - imemo[mv].proton - imemo[mv+2].proton;
        }
        else
        {
          imemo[mv].proton = 1 - ipatne[1];
          imemo[mv+2].proton = ipatin.ip1;
          imemo[mv+1].proton = me - imemo[mv].proton - imemo[mv+2].proton;
        }
      }
      else
      {
        imemo[mv].proton = ipatne[1];
        imemo[mv+1].proton = 0;
        imemo[mv+2].proton = ipatin.ip1;
      }
      return;
    }
    G4int migq;
    do
    {
      G4UniformRand() <= 0.5 ? imemo[mv].proton = 1 : imemo[mv].proton = 0;
      if( mq >= 1 )
      {
        G4UniformRand() <= 0.5 ? imemo[mv+2].proton = 1 : imemo[mv+2].proton = 0;
      }
      G4int lambda = 0;
      do
      {
        ++lambda;
        if( mq <= 1 || lambda != 2 )
        {
          G4double ran = G4UniformRand();
          if( ran >= oneThird )
          {
            ran < twoThirds ? imemo[mv+lambda].proton = 0 : imemo[mv+lambda].proton = -1;
          }
          else imemo[mv+lambda].proton = 1;
        }
      }
      while ( lambda < np );
      migq = 0;
      for( G4int i = 0; i < np; ++i )migq += imemo[mv+i].proton;
    } while( me != migq );
    return;
  }
 
 void G4PreEquilibrium::direction( G4ThreeVector v, const G4double u, const G4double tin1,
                                   const G4int mq, const G4int mv, const G4int np,
                                   const MASTRUCT partin, G4int &kp, const G4int ith )
  {
    // determine direction of secondary particle motion
    //
    G4int nd = 0;
    G4int m1, m2;
    kp = 0;
    if( mq > 1 )
    {
      if( G4UniformRand() >= 0.5 )
      {
        if( G4UniformRand() >= 0.5 )
        {
          m1 = 0;
          m2 = 1;
        }
        else
        {
          m1 = 1;
          m2 = 2;
        }
      }
      else
      {
        m1 = 0;
        m2 = 2;
      }
    }
    else
    {
      if( G4UniformRand() >= 0.5 )
      {
        m1 = 0;
        if( G4UniformRand() >= 0.5 )m2 = 1;
        else                        m2 = 2;
      }
      else
      {
        m1 = 1;
        m2 = 2;
      }
    }
    G4ThreeVector pl, pakv;
    G4int lambda, m1temp, m2temp;
  LABEL:
    pakv.setX( 0.0 );
    pakv.setY( 0.0 );
    pakv.setZ( 0.0 );
    m1temp = mv + m1;
    m2temp = mv + m2;
    lambda = 0;
    do
    {
      if( lambda != m1 && lambda != m2 )
      {
        G4int ja = coefficientTypeA( ith, mq, lambda );
        G4double ctl = costa( ja, tin1 );
        G4double fl = 2*pi*G4UniformRand();
        G4double stl = sqrt(1-ctl*ctl);
        G4double temp1 = cos(fl);
        G4double temp2 = sin(fl);
        G4double temp3 = pmemo[mv+lambda].array[0];
        G4ThreeVector tmp3vec( temp3*stl*temp1, temp3*stl*temp2, temp3*ctl );
        MASTRUCT tmp = { tmp3vec, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } };
        pmemo.insert( pmemo.begin()+mv+lambda, tmp );
        pakv += pmemo[mv+lambda].momentum;
      }
    } while( lambda++ < np );
    G4double pakvm = pakv.mag();
    if( np != 3 && pmemo[m1temp].momentum.z() >= pakvm+pmemo[m2temp].array[0] &&
                   pmemo[m1temp].momentum.z() <= abs(pakvm-pmemo[m2temp].array[0]) )
    {
      if( ++nd >= 100 )
      {
        kp = 2;
        return;
      }
      goto LABEL;
    }
    G4double v2 = v.mag2();
    G4double temp4 = sqrt(partin.array[4]*(partin.array[4]+2*partin.array[5]));
    G4double spv = temp4*partin.array[0]*partin.array[3]*v.x() + temp4*partin.array[1]*v.z() +
                  temp4*partin.array[0]*partin.array[2]*v.y();
    G4double temp5 = spv*(1/sqrt(1-v2)-1)/v2;
    G4double temp6 = (partin.array[4]+partin.array[5])/sqrt(1-v2);
    G4ThreeVector pin( temp4*partin.array[0]*partin.array[3] + v.x()*temp5 - v.x()*temp6,
                       temp4*partin.array[0]*partin.array[2] + v.y()*temp5 - v.y()*temp6,
                       temp4*partin.array[1] + v.z()*temp5 - v.z()*temp6 );
    //
    lambda = 1;
    if( lambda == m1 && lambda < np )++lambda;
    G4int ltemp = mv + lambda;
    G4ThreeVector plst;
    if( lambda == m1 )
    {
      G4ThreeVector pakst;
      pakst = rotation( pin, v, pakv );
      G4double ctm1 = (pmemo[m2temp].array[0]*pmemo[m2temp].array[0]-
                      pmemo[m1temp].array[0]*pmemo[m1temp].array[0]-pakvm*pakvm)/
          (2*pakvm*pmemo[m1temp].array[0]);
      G4double ctm2 = (pmemo[m1temp].array[0]*pmemo[m1temp].array[0]-
                      pmemo[m2temp].array[0]*pmemo[m2temp].array[0]-pakvm*pakvm)/
          (2*pakvm*pmemo[m2temp].array[0]);
      G4double fm1 = 2*pi*G4UniformRand();
      G4double fm2 = pi+fm1;
      G4double stm1 = sqrt(1-ctm1*ctm1);
      G4double stm2 = sqrt(1-ctm2*ctm2);
      G4double cfm1 = cos(fm1);
      G4double sfm1 = sin(fm1);
      G4double cfm2 = cos(fm2);
      G4double sfm2 = sin(fm2);
      pl.setX( pmemo[m1temp].array[0]*stm1*cfm1 );
      pl.setY( pmemo[m1temp].array[0]*stm1*sfm1 );
      pl.setZ( pmemo[m1temp].array[1]*ctm1 );
      plst = rotation( pakst, v, pl );
      pmemo[m1temp].momentum = plst;
      pl.setX( pmemo[m2temp].array[0]*stm2*cfm2 );
      pl.setY( pmemo[m2temp].array[0]*stm2*sfm2 );
      pl.setZ( pmemo[m2temp].array[0]*ctm2 );
      plst = rotation( pakst, v, pl );
      pmemo[m2temp].momentum = plst;
    }
    else if( lambda != m2 )
    {
      pl = pmemo[ltemp].momentum;
      plst = rotation( pin, v, pl );
      pmemo[ltemp].momentum = plst;
    }
    return;
  }
 
 G4int G4PreEquilibrium::coefficientTypeA( const G4int ith, const G4int mq,
                                           const G4int lamb )
  {
    // determine type of coefficients a(n,k)
    //
    G4int result;
    if( ith != 0 )
    {
      if( mq <= 1 )
      {
        lamb > 1 ? result = 25 : result = 24;
      }
      else
      {
        if( lamb > 1 )
        {
          lamb == 3 ? result = 16 : result = 21;
        }
        else result = 20;
      }
    }
    else
    {
      if( mq <= 1 )
      {
        lamb <= 1 ? result = 26 : result = 27;
      }
      else
      {
        if( lamb > 1 )
        {
          lamb == 3 ? result = 22 : result = 23;
        }
        else result = 22;
      }
    }
    return result;
  }
 
 G4bool G4PreEquilibrium::vmnsp( const MASTRUCT partin, const IPSTRUCT ipatin, const G4double u,
                                 const G4int mv, G4int &np, const G4int ith,
                                 const G4int mq, const G4double tin1 )
  {
    // calculate secondary particle number and
    // determine absolute values of momenta in inelastic interaction
    //
    G4int counter = 0;
  LABEL1:
    G4double u1 = u;
    G4int lambda = 0;
  LABEL2:
    if( mv+lambda >= 100 )
    {
      np = 0;
      //write(16,39)
      //39 format(45x,29h memory is exceeded in cascad)
      return false;
    }
    if( lambda == 0 )
    {
      G4ThreeVector tmp3vec( 0.0, 0.0, 0.0 );
      MASTRUCT matmp = { tmp3vec, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.94 } };
      pmemo.insert( pmemo.begin()+mv+lambda, matmp );
      IPSTRUCT iptmp = { 0, 0, 0, 1, 0 };
      imemo.insert( imemo.begin()+mv+lambda, iptmp );
    }
    else if( lambda == 2 )
    {
      G4ThreeVector tmp3vec( 0.0, 0.0, 0.0 );
      MASTRUCT matmp = { tmp3vec, { 0.0, 0.0, 0.0, 0.0, 0.0, partin.array[5] } };
      pmemo.insert( pmemo.begin()+mv+lambda, matmp );
      IPSTRUCT iptmp = { 0, 0, ipatin.ip2, ipatin.nuclearZone, 0 };
      imemo.insert( imemo.begin()+mv+lambda, iptmp );
    }
    else
    {
      G4ThreeVector tmp3vec( 0.0, 0.0, 0.0 );
      MASTRUCT matmp = { tmp3vec, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.14 } };
      pmemo.insert( pmemo.begin()+mv+lambda, matmp );
      IPSTRUCT iptmp = { 0, 0, 0, 0, 0 };
      imemo.insert( imemo.begin()+mv+lambda, iptmp );
    }
    G4int jb = coefficientTypeB( ith, mq, lambda );
    pmemo[mv+lambda].array[0] = secondaryParticleMomentum(jb,tin1);
    G4double el = sqrt(pmemo[mv+lambda].array[0]*pmemo[mv+lambda].array[0]+
                      pmemo[mv+lambda].array[5]*pmemo[mv+lambda].array[5]);
    G4double deltu = u1-el;
    if( lambda == 0 )
    {
      if( deltu <= partin.array[5] )
      {
        if( ++counter == 100 )return true;
        goto LABEL1;
      }
      if( ith == 0 )
      {
        u1 = deltu;
        ++lambda;
        goto LABEL2;
      }
      pmemo[mv+2].array[0] = sqrt(deltu*deltu-partin.array[5]*partin.array[5]);
      pmemo[mv+2].array[5] = partin.array[5];
      imemo[mv+2].ip1 = 0;
      imemo[mv+2].ip2 = ipatin.ip2;
      imemo[mv+2].nuclearZone = ipatin.nuclearZone;
      if( pmemo[mv].array[0] > pmemo[mv+1].array[0]+pmemo[mv+2].array[0] ||
          pmemo[mv].array[0] <= abs(pmemo[mv+1].array[0]-pmemo[mv+2].array[0]) )
      {
        if( ++counter == 100 )return true;
        goto LABEL1;
      }
      np = 3;
      return true;
    }
    if( deltu > 0.14 )
    {
      u1 = deltu;
      ++lambda;
      goto LABEL2;
    }
    if( lambda <= 1 || lambda == 3 )
    {
      if( ++counter == 100 )return true;
      goto LABEL1;
    }
    el += deltu;
    pmemo[mv+lambda].array[0] = sqrt(el*el-pmemo[mv+lambda].array[5]*pmemo[mv+lambda].array[5]);
    np = lambda;
    G4double c = pmemo[mv].array[0];
    G4int i;
    for( i = 0; i < np; ++i )
    {
      if( c < pmemo[mv+i].array[0] )c = pmemo[mv+i].array[0];
    }
    G4double pmax = c;
    G4double sigma = 0;
    for( i = 0; i < np; ++i )sigma += pmemo[mv+i].array[0];
    if( 2*pmax >= sigma )
    {
      if( ++counter == 100 )return true;
      goto LABEL1;
    }
    return false;
  }
 
 G4double G4PreEquilibrium::secondaryParticleMomentum( const G4int j, const G4double t )
  {
    // calculation of secondary particle momentum
    //
    const G4double bnkj[8][4][4] = 
    { { { 0.50278,  3.1442, -7.8172,  8.1667 },
        { 0.93482, -10.590,  29.227, -34.550 },
        {-0.096685,  4.7335, -14.298,  17.685 },
        {-0.025041, -0.62478,  2.0282, -2.5895 } },
      { { 1.1965, -0.82889,  1.0426, -1.9090 },
        { 0.28703, -4.9065,  16.264, -19.904 },
        {-0.24492,  2.9191, -9.5776,  11.938 },
        { 0.037297, -0.42200,  1.3883, -1.7476 } },
      { { 1.3508, -4.3139,  12.291, -15.288 },
        {-0.20086,  1.3641, -3.4030,  3.8559 },
        { 0.012583, -0.083492,  0.18600, -0.20043 },
        {-0.00023628,  0.0013514, -0.0024324,  0.0021906 } },
      { { 1.2419, -4.3633,  13.743, -18.592 },
        {-0.24404,  1.3158, -3.5691,  4.3867 },
        { 0.015693, -0.082579,  0.21427, -0.25846 },
        {-0.00029386,  0.0014060, -0.0033835,  0.0038664 } },
      { { 0.63054, -3.7333,  13.464, -18.594 },
        { 2.1801,  1.5163, -16.380,  27.944 },
        {-1.2886, -2.4570,  15.129, -23.295 },
        { 0.20915,  0.52279, -2.8687,  4.2688 } },
      { { 0.93363, -1.8181,  5.5157, -8.5216 },
        { 1.7811, -8.2927,  20.607, -20.827 },
        {-1.5264,  6.8433, -16.067,  16.845 },
        { 0.27128, -1.1944,  2.7495, -2.9045 } },
      { { 1.9439, -4.6268,  9.7879, -9.6074 },
        {-0.34640,  1.1093, -1.9313,  1.7064 },
        { 0.027054, -0.11638,  0.26969, -0.31853 },
        {-0.00066092,  0.0050728, -0.014995,  0.019605 } },
      { { 1.8693, -5.5678,  14.795, -16.903 },
        {-0.49965,  1.7874, -4.1330,  3.8393 },
        { 0.046194, -0.18536,  0.45315, -0.46273 },
        {-0.0013341,  0.0057710, -0.014554,  0.015554 } } };
    const G4double ckj[8][3] = {
        0.14509,   0.46520,  -0.033005,  0.15376,
        0.27436,  -0.014604,  0.62959,   0.17866,
       -0.0026216, 0.83810,   0.0086137, 0.0032946,
        0.092852,  0.53886,  -0.054493,  0.13032,
        0.40709,  -0.028782,  0.14909,   0.38502,
       -0.012775,  0.18024,   0.33022,  -0.0094491 };
    G4double bnk[4][4];
    G4int n, k;
    for( k = 0; k < 4; ++k )
    {
      for( n = 0; n < 4; ++n )bnk[n][k] = bnkj[n][k][j];
    }
    G4double s1 = 0;
    G4double r1 = G4UniformRand();
    G4double s2 = 0;
    G4double s3 = 0;
    for( n = 0; n < 4; ++n )
    {
      for( k = 0; k < 4; ++k )s1 += bnk[n][k]*pow(t,k)*pow(r1,n);
    }
    for( n = 0; n < 4; ++n )
    {
      for( k = 0; k < 4; ++k )s2 += bnk[n][k]*pow(t,k);
    }
    for( k = 0; k < 3; ++k )s3 += ckj[k][j]*pow(t,k);
    return s3*sqrt(r1)*(s1+(1-s2)*r1*r1*r1*r1);
  }
 
 G4int G4PreEquilibrium::coefficientTypeB( const G4int i, const G4int mq,
                                           const G4int lamb )
  {
    // determine type of coefficients b(n,k)
    //
    G4int result;
    if( i != 0 )
    {
      if( mq <= 1 )
      {
        lamb > 1 ? result = 5 : result = 4;
      }
      else
      {
        if( lamb > 1 )
        {
          lamb == 3 ? result = 0 : result = 1;
        }
        else result = 0;
      }
    }
    else
    {
      if( mq <= 1 )
      {
        lamb <= 1 ? result = 6 : result = 7;
      }
      else
      {
        if( lamb > 1 )
        {
          lamb == 3 ? result = 2 : result = 3;
        }
        else result = 2;
      }
    }
    return result;
  }
 
 void G4PreEquilibrium::statisticalModel( const G4double u, const G4ThreeVector v,
                                          const MASTRUCT partin, const G4int *ipatne,
                                          const G4int mv )
  {
    // determine secondary particle characteristics for
    // gamma-n interaction with statistical model
    //
    G4double tpim = (u*u+0.0196-1.08*1.08)/(2*u) - 0.14;
  LABEL:
    G4double t1 = G4UniformRand()*tpim;
    G4double t2 = G4UniformRand()*tpim;
    G4double e1 = t1 + 0.14;
    G4double e2 = t2 + 0.14;
    //
    if( G4UniformRand() >= 27*e1*e2*(u-e1-e2)/(u*u*u) )goto LABEL;
    //
    G4double t3 = u - e1 - e2 - 0.94;
    if( t3 <= 0 )goto LABEL;
    //
    G4double p1 = sqrt(t1*(t1+0.28));
    G4double p2 = sqrt(t2*(t2+0.28));
    G4double p3 = sqrt(t3*(t3+1.88));
    //
    if( (p1+p2-p3)*(p1-p2+p3)*(p2+p3-p1) <= 0 )goto LABEL;
    //
    G4double ct3 = 1-2*G4UniformRand();
    G4double fi3 = 2*pi*G4UniformRand();
    G4double temp2 = sqrt(1-ct3*ct3);
    //
    G4ThreeVector pv3( p3*temp2*cos(fi3),
                       p3*temp2*sin(fi3),
                       p3*ct3 );
    G4double temp3 = sqrt(partin.array[4]*(partin.array[4]+2*partin.array[5]));
    //
    G4ThreeVector pin( temp3*partin.array[0]*partin.array[3],
                       temp3*partin.array[0]*partin.array[2],
                       temp3*partin.array[1] );
    //
    G4ThreeVector pinst = cms( pin, v, partin.array[4]+partin.array[5] );
    //
    G4ThreeVector ps3 = rotation( pinst, v, pv3 );
    G4double ct1 = (p2*p2-p1*p1-p3*p3)/(2*p3*p1);
    G4double ct2 = (p1*p1-p2*p2-p3*p3)/(2*p3*p2);
    G4double fi1 = 2*pi*G4UniformRand();
    G4double fi2 = pi + fi1;
    G4double st1 = sqrt(1-ct1*ct1);
    G4double st2 = sqrt(1-ct2*ct2);
    //
    G4ThreeVector pv1( p1*st1*cos(fi1),
                       p1*st1*sin(fi1),
                       p1*ct1 );
    G4ThreeVector ps1 = rotation( ps3, v, pv1 );
    //
    G4ThreeVector pv2( p2*st2*cos(fi2),
                       p2*st2*sin(fi2),
                       p2*ct2 );
    // 
    G4ThreeVector ps2 = rotation( ps3, v, pv2 );
    pmemo[mv].momentum = ps1;
    pmemo[mv].array[5] = 0.14;
    pmemo[mv+1].momentum = ps2;
    pmemo[mv+1].array[5] = 0.14;
    pmemo[mv+2].momentum = ps3;
    pmemo[mv+2].array[5] = 0.94;
    if( ipatne[1] <= 0 )
    {
      G4double temp4 = G4UniformRand();
      if( temp4 > oneThird)
      {
        if( temp4 <= twoThirds )
        {
          imemo[mv].proton = 0;
          imemo[mv].ip1 = 0;
          imemo[mv].ip2 = 0;
          imemo[mv].nuclearZone = 0;
          imemo[mv+1].proton = -1;
          imemo[mv+1].ip1 = 0;
          imemo[mv+1].ip2 = 0;
          imemo[mv+1].nuclearZone = 0;
          imemo[mv+2].proton = 1;
          imemo[mv+2].ip1 = 0;
          imemo[mv+2].ip2 = 0;
          imemo[mv+2].nuclearZone = 1;
        }
        else
        {
          imemo[mv].proton = -1;
          imemo[mv].ip1 = 0;
          imemo[mv].ip2 = 0;
          imemo[mv].nuclearZone = 0;
          imemo[mv+1].proton = 1;
          imemo[mv+1].ip1 = 0;
          imemo[mv+1].ip2 = 0;
          imemo[mv+1].nuclearZone = 0;
          imemo[mv+2].proton = 0;
          imemo[mv+2].ip1 = 0;
          imemo[mv+2].ip2 = 0;
          imemo[mv+2].nuclearZone = 1;
        }
      }
      else
      {
        imemo[mv].proton = 0;
        imemo[mv].ip1 = 0;
        imemo[mv].ip2 = 0;
        imemo[mv].nuclearZone = 0;
        imemo[mv+1].proton = 0;
        imemo[mv+1].ip1 = 0;
        imemo[mv+1].ip2 = 0;
        imemo[mv+1].nuclearZone = 0;
        imemo[mv+2].proton = 0;
        imemo[mv+2].ip1 = 0;
        imemo[mv+2].ip2 = 0;
        imemo[mv+2].nuclearZone = 1;
      }
    }
    else
    {
      G4double temp4 = G4UniformRand();
      if( temp4 > oneThird )
      {
        if( temp4 <= twoThirds )
        {
          imemo[mv].proton = 0;
          imemo[mv].ip1 = 0;
          imemo[mv].ip2 = 0;
          imemo[mv].nuclearZone = 0;
          imemo[mv+1].proton = 1;
          imemo[mv+1].ip1 = 0;
          imemo[mv+1].ip2 = 0;
          imemo[mv+1].nuclearZone = 0;
          imemo[mv+2].proton = 0;
          imemo[mv+2].ip1 = 0;
          imemo[mv+2].ip2 = 0;
          imemo[mv+2].nuclearZone = 1;
        }
        else
        {
          imemo[mv].proton = 1;
          imemo[mv].ip1 = 0;
          imemo[mv].ip2 = 0;
          imemo[mv].nuclearZone = 0;
          imemo[mv+1].proton = -1;
          imemo[mv+1].ip1 = 0;
          imemo[mv+1].ip2 = 0;
          imemo[mv+1].nuclearZone = 0;
          imemo[mv+2].proton = 1;
          imemo[mv+2].ip1 = 0;
          imemo[mv+2].ip2 = 0;
          imemo[mv+2].nuclearZone = 1;
        }
      }
      else
      {
        imemo[mv].proton = 0;
        imemo[mv].ip1 = 0;
        imemo[mv].ip2 = 0;
        imemo[mv].nuclearZone = 0;
        imemo[mv+1].proton = 0;
        imemo[mv+1].ip1 = 0;
        imemo[mv+1].ip2 = 0;
        imemo[mv+1].nuclearZone = 0;
        imemo[mv+2].proton = 1;
        imemo[mv+2].ip1 = 0;
        imemo[mv+2].ip2 = 0;
        imemo[mv+2].nuclearZone = 1;
      }
    }
    return;
  }
 
 void G4PreEquilibrium::isocem( const G4double u, const G4ThreeVector v,
                                const G4double tin1, const MASTRUCT partin,
                                const G4int *ipatne, const G4int mv )
  {
    // interaction with (3/2,3/2) isobar production
    //
    //
    G4double a1 = (u*u+0.0196-1.08*1.08)/(2*u);
    G4double a2 = sqrt(a1*a1-0.0196);
    G4double a3 = u-a1;
    //
    G4double epim, bms, pim, p;
    do
    {
      bms = G4UniformRand()*(u-1.22)+1.08;
      epim = (u*u+0.0196-bms*bms)/(2*u);
      pim = sqrt(epim*epim-0.0196);
      p = ((pim*epim*(u-epim))/u)*
          crossSectionInterp((bms*bms-1.08*1.08)/1.88,8)/(200*a1*a2*a3/u);
    }
    while ( G4UniformRand() >= p );
    //
    G4double ctpi;
    tin1 < 1 ? ctpi = costa( 27, tin1 ) : ctpi = costa( 28, tin1 );
    G4double fipi = 2*pi*G4UniformRand();
    G4double epit = (bms*bms+0.0196-0.94*0.94)/(2*bms);
    G4double temp1 = sqrt(1-ctpi*ctpi);
    
    G4ThreeVector ppim( pim*temp1*cos(fipi),
                        pim*temp1*sin(fipi),
                        pim*ctpi );
    
    G4ThreeVector vt = ppim * (1/(epim-u));
    
    G4double ctilpi = 1-2*G4UniformRand();
    G4double ftilpi = 2*pi*G4UniformRand();
    temp1 = sqrt(1-ctilpi*ctilpi);
    G4double pmt = sqrt(epit*epit-0.0196);

    G4ThreeVector ppit( pmt*temp1*cos(ftilpi),
                        pmt*temp1*sin(ftilpi),
                        pmt*ctilpi );
    
    vt *= -1.0;
    
    G4ThreeVector ppi = cms( ppit, vt, epit );

    G4ThreeVector ppt = ppit * (-1.0);
    
    G4ThreeVector pp = cms( ppt, vt, bms-epit );

    G4double temp2 = sqrt(partin.array[4]*(partin.array[4]+2*partin.array[5]));
    G4ThreeVector pin( temp2*partin.array[0]*partin.array[3],
                       temp2*partin.array[0]*partin.array[2],
                       temp2*partin.array[1] );
    
    G4ThreeVector pinst = cms( pin, v, partin.array[4]+partin.array[5] );

    G4ThreeVector ppimst;
    ppimst = rotation( pinst, v, ppim );
    pmemo[mv].momentum = ppimst;
    pmemo[mv].array[5] = 0.14;

    G4ThreeVector ppist;
    ppist = rotation( pinst, v, ppi );
    pmemo[mv+1].momentum = ppist;
    pmemo[mv+1].array[5] = 0.14;

    G4ThreeVector ppst;
    ppst = rotation( pinst, v, pp );
    pmemo[mv+2].momentum = ppst;
    pmemo[mv+2].array[5] = 0.94;
    if( ipatne[1] > 0 )
    {
      imemo[mv].proton = -1;
      imemo[mv].ip1 = 0;
      imemo[mv].ip2 = 0;
      imemo[mv].nuclearZone = 0;
      imemo[mv+1].proton = 1;
      imemo[mv+1].ip1 = 0;
      imemo[mv+1].ip2 = 0;
      imemo[mv+1].nuclearZone = 0;
      imemo[mv+2].proton = 1;
      imemo[mv+2].ip1 = 0;
      imemo[mv+2].ip2 = 0;
      imemo[mv+2].nuclearZone = 1;
    }
    else
    {
      imemo[mv].proton = 1;
      imemo[mv].ip1 = 0;
      imemo[mv].ip2 = 0;
      imemo[mv].nuclearZone = 0;
      imemo[mv+1].proton = -1;
      imemo[mv+1].ip1 = 0;
      imemo[mv+1].ip2 = 0;
      imemo[mv+1].nuclearZone = 0;
      imemo[mv+2].proton = 0;
      imemo[mv+2].ip1 = 0;
      imemo[mv+2].ip2 = 0;
      imemo[mv+2].nuclearZone = 1;
    }
    return;
  }
 
 G4ThreeVector G4PreEquilibrium::cms( const G4ThreeVector p, const G4ThreeVector v,
                                      const G4double tcm ) const
  {
    // momentum calculation in system which has velocity v relative to given one
    //
    G4double temp1 = sqrt( 1 - v.mag2() );
    G4double temp2 = p.dot(v)/v.mag2()*(1/temp1-1);
    G4ThreeVector result( p.x() + v.x()*(temp2-tcm/temp1),
                          p.y() + v.y()*(temp2-tcm/temp1),
                          p.z() + v.z()*(temp2-tcm/temp1) );
    return result;
  }
 
 void G4PreEquilibrium::bertcem( const G4int nuclideTypeNumber, const G4int mediumNumber )
  {
    // this routine was called by cascade and by cemgeo
    //
    aNucl = GetA_FCOMON(nuclideTypeNumber,mediumNumber);
    zNucl = GetZ_FCOMON(nuclideTypeNumber,mediumNumber);
    SetT0( 3.4995 );
    //
    // neutron density = proton density = nucleon density
    //
    const G4double a[10] = { 0.95, 0.8, 0.5, 0.2, 0.1, 0.05, 0.01, 0.0, 0.0, 0.0 };
    const G4double r0n = 1.07;
    const G4double rmax = 10;
    const G4double bnr = 0.545;
    //
    G4double rn = r0n*pow(aNucl,oneThird);
    G4int i;
    //
    const G4int N = 7;
    for( i = 0; i < N; ++i )rsm.push_back( rn+bnr*(log((1-a[i])/a[i])) );
    rsm.push_back( rmax );
    //
    G4double hin[10], fi[10];
    fi[0] = fintfis( rsm[0], 0.0, bnr, rn );
    hin[0] = fints2( rsm[0], 0.0, bnr, rn );
    G4double sumhin = hin[0];
    for( i = 1; i < N; ++i )
    {
      fi[i] = fintfis( rsm[i], rsm[i-1], bnr, rn );
      hin[i] = fints2( rsm[i], rsm[i-1], bnr, rn );
      sumhin += hin[i];
    }
    //
    G4double rhon0 = aNucl/(4*pi*sumhin);
    G4double fi0 = (4*pi*zNucl*0.00144*rhon0)/aNucl;
    rhon.push_back( 3*hin[0]*rhon0/rsm[0]/rsm[0]/rsm[0] );
    af.push_back( fi[0]*fi0*3/rsm[0]/rsm[0]/rsm[0] );
    for( i = 1; i < N; ++i )
    {
      rhon.push_back( 3*hin[i]*rhon0/(rsm[i]*rsm[i]*rsm[i]-rsm[i-1]*rsm[i-1]*rsm[i-1]) );
      tfp.push_back( 0.1985*pow( rhon[i]*zNucl/aNucl, twoThirds ) );
      tfn.push_back( 0.1985*pow( rhon[i]*(aNucl-zNucl)/aNucl, twoThirds ) );
      af.push_back( 3*fi0*fi[i]/(rsm[i]*rsm[i]*rsm[i]-rsm[i-1]*rsm[i-1]*rsm[i-1]) );
      rbig.push_back( rsm[i]/rsm[N] );
    }
    rbig.push_back( rmax/rsm[N-1] );
    af.push_back( zNucl*0.00144/rmax );
    for( i = 0; i < N; ++i )
    {
      rhop.push_back( rhon[i]*zNucl/aNucl );
      rhon[i] *= (aNucl-zNucl)/aNucl;
    }
    //
    af.insert( af.begin(), 10, 0.0 );
    //
    G4ThreeVector amnucl;
    if( !cascem( amnucl ) )
    {
      SetLXYZ( amnucl );
      SetN0( exitons.protons + exitons.neutrons + exitons.hols );
      SetH0( exitons.hols );
      SetPZ0( exitons.protons );
      SetP0( exitons.protons + exitons.neutrons );
    }
    return;
  }
 
 G4double G4PreEquilibrium::fintfis( const G4double flim1, const G4double flim2,
                                     const G4double bs, const G4double rs )
  {
    const G4int N = 7;
    const G4double w[8] = { 0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
                            0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363 };
    const G4double fiks[8] = {  0.9602898565,  0.7966664774,  0.5255324099,  0.1834346425,
                               -0.1834346425, -0.5255324099, -0.7966664774, -0.9602898565 };
    //
    G4double temp = 0;
    for( G4int i = 0; i < 8; ++i )
    {
      G4double fy = ((flim1-flim2)*fiks[i]+flim1+flim2)/2;
      G4double tmp1 = 0.0;
      G4double tmp2 = 0.0;
      G4double tmp3 = 0.0;
      for( G4int j = 0; j < 8; ++j )
      {
        fy = ((flim1-flim2)*fiks[j]+flim1+flim2)/2;
        tmp1 += w[j]*(fy*fy)/(1+exp((fy-rs)/bs));
        tmp2 += w[j]*rsm[N]/(1+exp((rsm[N]-rs)/bs));
        tmp3 += w[j]*fy/(1+exp((fy-rs)/bs));
      }
      G4double f1 = tmp1*(flim1-flim2)/2;
      G4double f2 = tmp2*(flim1-flim2)/2;
      G4double f3 = tmp3*(flim1-flim2)/2;
      temp += w[i] * (f1/fy + f2 - f3)*fy*fy;
    }
    return temp*(flim1-flim2)/2;
  }
 
 G4double G4PreEquilibrium::fints2( G4double flim1, G4double flim2,
                                   G4double bs, G4double rs )
  {
    const G4double w[8] = { 0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
                            0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363 };
    const G4double fiks[8] = {  0.9602898565,  0.7966664774,  0.5255324099,  0.1834346425,
                               -0.1834346425, -0.5255324099, -0.7966664774, -0.9602898565 };
    G4double temp = 0;
    for( G4int i = 0; i < 8; ++i )
    {
      G4double fy = ((flim1-flim2)*fiks[i]+flim1+flim2)/2;
      temp += w[i] * (fy*fy)/(1+exp((fy-rs)/bs));
    }
    return temp*(flim1-flim2)/2;
  }
 
 G4double G4PreEquilibrium::arfaf( G4double *ami )
  {
    // renormalisation of emission probabilities
    //
    const G4double AMF = 0.0;
    //
    G4double sfix = 30;
    G4double smx = 0;
    for( G4int k = 0; k < 7; ++k )
    {
      G4double q8;
      if( rj[k] <= 0 )q8 = 0;
      else
      {
        G4double tmp;
        k == 6 ? tmp = AMF : tmp = ami[k];
        q8 = 2*sqrt(tmp*afj[k]*rj[k]);
      }
      smx = max( smx, q8 );
    }
    G4double result;
    smx <= sfix ? result = 0.0 : result = smx-sfix;
    return result;
  }
 
 G4double G4PreEquilibrium::fam( const G4double a, const G4double z, const G4double e )
  {
    //
    // third set of Iljinov, Mebel et al. parameters: without collective effects
    // Truran, Cameron & Hils shell corrections, without a dependence on ga
    //
    const G4double al = 0.072;
    const G4double be = 0.257;
    const G4double ga = 0.059;
    // shell is an inline function defined in the header
    return (al+be/pow(a,oneThird))*(1+(1-exp(-ga*e))*shell(a,z)/e);
  }
 
 G4double G4PreEquilibrium::gameqf( const G4int j, const G4double cc, const G4double per,
                                    const G4double am, const G4double radncl )
  {
    //
    // a possible bug in the original cem.f
    // idel set to 0 in inpcem, and nowhere else
    // idel checked in precof, bypassing setting amf
    //
    const G4double AMF = 0.0;
    //
    //                     n  p  d  t  He-3 He-4
    const G4int gam[6] = { 1, 1, 3, 3, 3,   2 }; // gammab/2
    
    G4double result;
    G4double q1, q2;
    if( j == 6 )
    {
      q1 = 2*sqrt( AMF * afj[6] * rj[6] );
      q2 = 33.4/(pi*radncl*radncl);
      result = q2/(AMF*afj[6])*((q1-1)*exp(q1-per)+exp(-per));
    }
    else
    {
      G4double alfa, beta, q3, q4, q5;
      if( j == 0 )
      {
        alfa = 0.76 + 2.2/pow(afj[0],oneThird);
        beta = (2.12/pow(afj[0],twoThirds)-0.05)/alfa;
      }
      else
      {
        alfa = 1 + cc;
        beta = 0.0;
      }
      q1 = am*afj[j];
      q2 = q1*rj[j];
      q3 = gam[j]*pow(afj[j],twoThirds)*alfa/(q1*q1);
      q4 = (2*beta*q1-3)/2 + q2;
      q5 = (2*beta*q1-3)*(sqrt(q2)-0.5) + 2*q2;
      result = q3*(q4*exp(-per)+q5*exp(2*sqrt(q2)-per));
    }
    return result;
  }
 
 G4double G4PreEquilibrium::poten( const G4int i, const IPSTRUCT ip ) const
  {
    // calculation of particle potential in i'th nuclear zone
    //
    const G4int N = 7;
    const G4double EPS = 0.007;
    const G4double VPI = 0.025;

    G4double result;
    if( i >= N+2 )
      result = 0.0;
    else if( i == N+1 )
      result = ip.proton*af[N];
    else if( ip.ip2 != 0 )
      result = ip.proton*af[i-1];
    else if( ip.nuclearZone != 0 )
      result = ip.proton*af[i-1] + tfp[i-1]*ip.proton + (1-ip.proton)*tfn[i-1] + ip.nuclearZone*EPS;
    else
      result = ip.proton*af[i-1] - VPI*(ip.ip1-1);
    return result;
  }

 G4double G4PreEquilibrium::cinema( const G4ThreeVector pstar, const G4ThreeVector v,
                                    G4ThreeVector &p, G4double &ct, G4double &st,
                                    G4double &cfi, G4double &sfi, const G4double cm ) const
  {
    // kinematic code
    //
    G4double spv = pstar.dot(v);
    G4double v2 = v.mag2();
    G4double temp1 = sqrt(1-v2);
    G4double temp2 = spv*(1/temp1-1)/v2;
    G4double pmstar = pstar.mag();
    G4double tstar = sqrt(pmstar*pmstar+cm*cm) - cm;
    p.setX( pstar.x() + v.x()*(temp2 + (tstar+cm)/temp1) );
    p.setY( pstar.y() + v.y()*(temp2 + (tstar+cm)/temp1) );
    p.setZ( pstar.z() + v.z()*(temp2 + (tstar+cm)/temp1) );
    G4double pm = p.mag();
    if( p.z()*p.z() >= pm*pm )
    {
      ct = 1;
      st = 0;
      cfi = 1;
      sfi = 0;
    }
    else
    {
      ct = p.z()/pm;
      st = sqrt(1-ct*ct);
      cfi = p.x()/(pm*st);
      sfi = p.y()/(pm*st);
    }
    return sqrt(pm*pm+cm*cm)-cm;
  }
 
 G4ThreeVector G4PreEquilibrium::rotation( const G4ThreeVector a,
                                           const G4ThreeVector b,
                                           const G4ThreeVector ps ) const
  {
    G4double sp = a.dot(b);
    //
    G4double amod = a.mag();
    G4double alpha1 = sp/amod;
    G4double alpha2 = sqrt( b.mag2() - alpha1*alpha1 );
    //
    G4ThreeVector an = a.cross(b);
    //
    G4ThreeVector pr;
    pr.setX( ps.x()*b.x()/alpha2 + (ps.z()-alpha1*ps.x()/alpha2)*a.x()/amod + (ps.y()*an.x()) );
    pr.setY( ps.x()*b.y()/alpha2 + (ps.z()-alpha1*ps.x()/alpha2)*a.y()/amod + (ps.y()*an.y()) );
    pr.setZ( ps.x()*b.z()/alpha2 + (ps.z()-alpha1*ps.x()/alpha2)*a.z()/amod + (ps.y()*an.z()) );
    pr *= 1/alpha2/amod;
    return pr;
  }
 
 G4double G4PreEquilibrium::crossSectionInterp( const G4double x, const G4int lq )
  {
    const G4double sigma[28][30] =
    { {17613.000,  330.000,  154.000,   96.000,   70.000,
          51.000,   38.400,   30.000,   23.600,   22.400,
          22.200,   22.600,   23.400,   24.700,   29.500,
          40.500,   48.500,   47.400,   47.000,   46.700,
          46.000,   45.000,   43.000,   41.200,   40.800,
          40.500,   39.800,   39.000,   39.000,   38.500},
      {17613.000,  330.000,  154.000,   96.000,   70.000,
          51.000,   38.400,   30.000,   23.600,   22.400,
          22.200,   22.600,   22.800,   22.900,   25.000,
          25.000,   25.000,   25.000,   25.000,   22.000,
          19.500,   17.500,   15.000,   13.700,   12.000,
          11.000,    9.800,    8.800,    8.500,    6.500},
      {20357.000,  950.000,  480.000,  300.000,  200.000,
         160.000,  108.000,   74.000,   50.000,   41.000,
          36.500,   34.000,   32.500,   32.000,   34.200,
          36.100,   37.800,   38.400,   39.000,   40.000,
          40.000,   40.500,   41.100,   42.300,   42.300,
          42.000,   41.100,   39.600,   39.500,   39.200},
      {20357.000,  950.000,  480.000,  300.000,  200.000,
         160.000,  108.000,   74.000,   50.000,   41.000,
          36.500,   34.000,   32.500,   31.200,   30.800,
          25.100,   19.100,   18.000,   16.500,   16.000,
          15.200,   14.000,   12.100,   10.800,   10.700,
          11.000,    9.600,    8.000,    7.000,    6.200},
      {    6.000,    6.000,    6.500,    7.000,    8.500,
          10.500,   16.500,   25.300,   57.500,   68.500,
          64.500,   52.000,   40.500,   25.700,   29.000,
          32.100,   45.600,   38.000,   44.300,   54.000,
          58.000,   45.700,   35.300,   35.000,   34.300,
          34.000,   32.400,   30.400,   26.500,   25.000},
      {    2.000,    2.000,    2.250,    2.300,    2.400,
           2.500,    4.500,    7.700,   20.000,   25.200,
          23.900,   21.000,   16.000,   10.500,   11.200,
          14.000,   20.300,   16.000,   19.300,   26.000,
          26.500,   18.700,   11.700,   10.500,    9.600,
           9.500,    6.800,    5.900,    5.000,    4.000},
      {    4.000,    4.000,    4.250,    4.700,    6.100,
           8.000,   12.000,   17.600,   37.500,   43.300,
          40.600,   31.000,   24.500,   13.100,   11.000,
           9.500,    8.300,    5.000,    5.500,    6.800,
           7.000,    3.500,    2.000,    1.900,    1.800,
           1.600,    0.220,    0.150,    0.048,    0.009},
      {    1.900,    2.300,    3.500,    5.500,    9.000,
          14.000,   28.000,   60.000,  163.000,  195.000,
         185.000,  145.000,  113.000,   45.000,   25.200,
          21.600,   15.600,   15.200,   19.500,   22.800,
          24.500,   27.600,   36.700,   41.000,   39.000,
          32.300,   28.900,   27.700,   24.900,   23.500},
      {    1.900,    2.300,    3.500,    5.500,    9.000,
          14.000,   28.000,   60.000,  163.000,  195.000,
         184.900,  144.800,  112.800,   44.400,   23.200,
          18.600,   10.800,    7.700,    9.000,   10.200,
          11.300,   13.500,   16.900,   19.000,   16.900,
          12.600,    5.700,    5.600,    4.900,    4.000},
      {   10.000,   12.000,   14.000,   16.000,   17.000,
          20.000,   30.000,   38.000,   42.000,   38.500,
          32.000,   24.000,   18.000,    1.000,    0.000,
           0.000,    0.000,    0.000,    0.000,    0.000,
           0.000,    0.000,    0.000,    0.000,    0.000,
           0.000,    0.000,    0.000,    0.000,    0.000},
      {    0.000,    0.000,    0.000,    0.050,    0.100,
           0.300,    0.600,    1.200,    2.200,    3.200,
           3.400,    3.600,    3.700,    3.800,    3.900,
           3.900,    3.900,    4.000,    4.000,    4.000,
           4.000,    4.000,    3.900,    3.800,    3.500,
           3.200,    2.900,    2.700,    2.400,    2.100},
      {    0.000,    0.000,    0.000,    0.400,    0.800,
           1.400,    2.300,    4.400,    8.000,   10.800,
          15.000,   16.100,   16.600,   17.000,   17.300,
          17.500,   17.600,   17.700,   17.800,   17.900,
          17.800,   17.500,   16.800,   16.000,   13.700,
          12.600,   11.600,   10.800,    9.200,    8.000},
      {    0.000,    0.000,    0.000,    0.300,    0.600,
           1.000,    1.600,    2.700,    5.000,    6.300,
           6.900,    7.200,    7.400,    7.500,    7.600,
           7.700,    7.800,    8.000,    8.100,    8.100,
           8.100,    8.000,    7.900,    7.700,    7.000,
           6.600,    6.200,    5.800,    4.600,    3.500},
      {    0.000,    0.000,    0.000,    0.000,    0.100,
           0.200,    0.500,    0.800,    1.800,    2.100,
           2.400,    2.600,    2.800,    3.000,    3.200,
           3.300,    3.400,    3.650,    3.800,    3.850,
           3.900,    4.000,    4.000,    3.900,    3.500,
           3.200,    2.700,    2.100,    1.000,    0.600},
      {    0.000,    0.000,    0.000,    0.000,    0.200,
           0.600,    1.200,    2.100,    3.300,    4.100,
           4.900,    6.000,    7.700,    8.900,   10.000,
          10.200,   10.100,    9.800,    9.800,   10.100,
          10.300,    8.500,    6.900,    5.500,    3.900,
           3.300,    3.000,    2.700,    2.200,    2.100},
      {    0.000,    0.000,    0.000,    0.000,    0.100,
           0.300,    0.500,    0.700,    0.900,    1.100,
           1.300,    1.500,    1.800,    2.200,    2.600,
           2.600,    2.400,    2.400,    3.000,    3.400,
           3.500,    3.100,    2.700,    2.400,    2.100,
           1.900,    1.800,    1.700,    1.500,    1.300},
      {    0.000,    0.000,    0.000,    0.200,    0.700,
           1.400,    3.000,    4.900,    5.100,    5.000,
           4.400,    3.600,    3.100,    2.700,    2.400,
           2.200,    2.000,    1.700,    1.600,    1.400,
           1.300,    1.100,    0.900,    0.800,    0.600,
           0.600,    0.600,    0.500,    0.400,    0.300},
      {    0.000,    0.000,    0.000,    0.200,    0.500,
           0.800,    1.300,    2.100,    3.400,    4.700,
           6.400,    6.800,    6.200,    5.200,    7.100,
           8.000,    8.900,    9.800,   10.300,   10.500,
          10.600,    9.900,    7.700,    5.500,    3.400,
           2.800,    2.700,    2.600,    2.200,    1.900},
      {    0.000,    0.100,    0.500,    1.200,    3.000,
           4.200,    5.000,    5.200,    6.100,    7.900,
           8.900,    9.600,    9.900,   10.000,   10.100,
          10.100,   10.000,    9.800,    9.500,    8.800,
           8.000,    6.800,    5.900,    5.200,    4.200,
           2.900,    2.700,    2.500,    2.200,    1.900},
      {    0.000,    0.011,    0.026,    0.060,    0.120,
           0.180,    0.270,    0.280,    0.240,    0.186,
           0.130,    0.076,    0.052,    0.040,    0.031,
           0.032,    0.040,    0.044,    0.039,    0.033,
           0.028,    0.025,    0.022,    0.019,    0.015,
           0.012,    0.007,    0.004,    0.002,    0.000},
      {    0.000,    0.070,    0.107,    0.143,    0.183,
           0.217,    0.240,    0.214,    0.177,    0.141,
           0.117,    0.089,    0.081,    0.080,    0.083,
           0.093,    0.103,    0.095,    0.061,    0.050,
           0.049,    0.052,    0.056,    0.050,    0.025,
           0.017,    0.008,    0.004,    0.002,    0.000},
      {    5.000,    3.000,    1.125,    0.900,    0.700,
           0.580,    0.500,    0.415,    0.350,    0.280,
           0.260,    0.270,    0.280,    0.295,    0.310,
           0.298,    0.250,    0.200,    0.145,    0.115,
           0.095,    0.065,    0.040,    0.028,    0.018,
           0.010,    0.005,    0.000,    0.000,    0.000},
      {    0.000,    0.011,    0.031,    0.055,    0.067,
           0.078,    0.084,    0.088,    0.088,    0.087,
           0.085,    0.083,    0.080,    0.078,    0.074,
           0.072,    0.069,    0.066,    0.064,    0.060,
           0.057,    0.052,    0.038,    0.031,    0.027,
           0.024,    0.023,    0.022,    0.022,    0.022},
      {    0.000,    0.002,    0.004,    0.006,    0.007,
           0.008,    0.008,    0.009,    0.011,    0.012,
           0.013,    0.014,    0.015,    0.016,    0.016,
           0.016,    0.016,    0.016,    0.016,    0.015,
           0.015,    0.014,    0.012,    0.010,    0.009,
           0.008,    0.007,    0.007,    0.007,    0.007},
      {    0.000,    0.002,    0.004,    0.006,    0.007,
           0.008,    0.008,    0.009,    0.011,    0.012,
           0.013,    0.014,    0.015,    0.016,    0.016,
           0.016,    0.016,    0.016,    0.016,    0.015,
           0.015,    0.014,    0.012,    0.010,    0.009,
           0.008,    0.007,    0.007,    0.007,    0.007},
      {    0.000,    0.010,    0.027,    0.047,    0.057,
           0.065,    0.069,    0.071,    0.069,    0.065,
           0.061,    0.057,    0.052,    0.048,    0.043,
           0.039,    0.035,    0.030,    0.027,    0.023,
           0.019,    0.015,    0.008,    0.004,    0.003,
           0.002,    0.002,    0.001,    0.001,    0.001},
      {    0.000,    0.015,    0.039,    0.067,    0.081,
           0.094,    0.100,    0.106,    0.110,    0.111,
           0.111,    0.111,    0.110,    0.110,    0.106,
           0.104,    0.101,    0.098,    0.096,    0.090,
           0.087,    0.080,    0.062,    0.051,    0.045,
           0.040,    0.037,    0.036,    0.036,    0.036},
      {    0.000,    0.081,    0.133,    0.203,    0.303,
           0.397,    0.510,    0.494,    0.417,    0.327,
           0.247,    0.180,    0.172,    0.187,    0.208,
           0.231,    0.253,    0.250,    0.211,    0.194,
           0.187,    0.187,    0.184,    0.173,    0.138,
           0.119,    0.095,    0.070,    0.049,    0.037} };
    const G4double argus[6][30] = 
    { {0.000,    0.010,    0.020,    0.030,    0.040,
       0.050,    0.070,    0.100,    0.150,    0.200,
       0.250,    0.300,    0.350,    0.400,    0.500,
       0.650,    0.850,    0.950,    1.100,    1.300,
       1.500,    2.000,    3.000,    4.000,    5.000,
       7.000,   10.000,   16.000,   22.000,   30.000},
      {0.000,    0.010,    0.020,    0.030,    0.040,
       0.050,    0.075,    0.100,    0.150,    0.175,
       0.200,    0.225,    0.250,    0.350,    0.450,
       0.500,    0.600,    0.700,    0.800,    0.850,
       0.900,    1.000,    1.200,    1.300,    1.400,
       1.600,    3.000,    4.000,   10.000,   20.000},
      {0.200,    0.250,    0.300,    0.350,    0.400,
       0.450,    0.500,    0.550,    0.600,    0.650,
       0.700,    0.750,    0.800,    0.850,    0.900,
       0.950,    1.000,    1.100,    1.200,    1.300,
       1.400,    1.600,    1.800,    2.000,    2.400,
       2.600,    2.800,    3.000,    3.500,    4.000},
      {0.150,    0.175,    0.200,    0.225,    0.250,
       0.275,    0.300,    0.325,    0.350,    0.375,
       0.400,    0.450,    0.500,    0.550,    0.600,
       0.650,    0.700,    0.750,    0.800,    0.850,
       0.900,    0.950,    1.000,    1.050,    1.150,
       1.250,    1.500,    2.000,    3.000,    4.000},
      {0.010,    0.020,    0.040,    0.050,    0.060,
       0.070,    0.080,    0.090,    0.100,    0.120,
       0.150,    0.170,    0.190,    0.210,    0.250,
       0.270,    0.300,    0.320,    0.350,    0.370,
       0.400,    0.440,    0.500,    0.600,    0.700,
       0.800,    0.900,    1.000,    2.000,    3.000},
      {0.400,    0.450,    0.500,    0.550,    0.575,
       0.600,    0.625,    0.650,    0.700,    0.750,
       0.800,    0.850,    0.900,    0.950,    1.000,
       1.050,    1.100,    1.150,    1.200,    1.250,
       1.300,    1.500,    2.000,    2.500,    3.000,
       3.500,    4.000,    4.500,    5.000,    5.500} };
    G4int i;
    if( lq < 3 )i = 0;
    else
    {
      if( lq < 9 )i = 1;
      else
      {
        if( lq < 18 )i = 2;
        else
        {
          if( lq < 20 )i = 3;
          else
          {
            if( lq < 21 )i = 4;
            else
            {
              if( lq < 26 )i = 5;
              else         i = 3;
            }
          }
        }
      }
    }
    G4int lpha = 0;
    while( x > argus[lpha][i] )++lpha;
    if( x == argus[lpha][i] )return sigma[lpha][lq];
    if( lpha == 0 )return 0.0;
    G4double phi1, psi1, phi2, psi2, phi3, psi3;
    if( lpha >= 28 )
    {
      phi1 = sigma[27][lq];
      psi1 = argus[27][i];
      phi2 = sigma[28][lq];
      psi2 = argus[28][i];
      phi3 = sigma[29][lq];
      psi3 = argus[29][i];
    }
    else
    {
      phi1 = sigma[lpha-1][lq];
      psi1 = argus[lpha-1][i];
      phi2 = sigma[lpha][lq];
      psi2 = argus[lpha][i];
      phi3 = sigma[lpha+1][lq];
      psi3 = argus[lpha+1][i];
    }
    G4double delta = (psi2-psi3)*psi1*psi1 + (psi3-psi1)*psi2*psi2 + (psi1-psi2)*psi3*psi3;
    G4double deltaa = phi1*(psi2-psi3) + phi2*(psi3-psi1) + phi3*(psi1-psi2);
    G4double deltab = (phi2-phi3)*psi1*psi1 + (phi3-phi1)*psi2*psi2 + (phi1-phi2)*psi3*psi3;
    G4double deltac = (psi2*phi3-psi3*phi2)*psi1*psi1 + (psi3*phi1-psi1*phi3)*psi2*psi2 +
                     (psi1*phi2-psi2*phi1)*psi3*psi3;
    G4double a = deltaa/delta;
    G4double b = deltab/delta;
    G4double c = deltac/delta;
    return a*x*x + b*x + c;
  }
 
 G4double G4PreEquilibrium::sigmat( const G4int i, const G4int ms, const G4int mq,
                                   const G4int ksi, const G4int iks, const G4double t )
  {
    // choose cross section type and calculate cross section value for given energy
    //
    const G4int icst[28] = {
        210,   211,   220,   221,   120,   121,   122,
        110,   111,   123,   214,   215,   224,   225,
        114,   115,   126,   124,   125, 10111, 10112,
      10113, 10115, 10114, 10116, 10118, 10117, 10110 };
    const G4int nsicst[21] = {
        112, 113, 116, 117, 127, 130, 131,
        132, 133, 134, 135, 136, 137, 212,
        213, 216, 217, 222, 223, 226, 227 };    
    
    G4int ics = 10000*i + 1000*ms + 100*mq + 10*ksi + iks;
    G4int js = 0;
    do
    {
      if( ics == icst[js] )return crossSectionInterp( t, js );
    } while ( js++ < 28 );
    G4int nsjs = 0;
    while( ics != nsicst[nsjs] )++nsjs;
    G4double result = 0;
    switch (nsjs)
    {
     case 2:
       result = crossSectionInterp( t, 9 );
       break;
     case 4:
       result = crossSectionInterp( t, 14 ) + crossSectionInterp( t, 15 );
       break;
     case 5:
       result = crossSectionInterp( t, 17 ) + crossSectionInterp( t, 18 ) +
                crossSectionInterp( t, 16 );
       break;
     case 6:
       result = (crossSectionInterp( t, 7 ) + crossSectionInterp( t, 4 ))/2;
       break;
     case 7:
       result = (crossSectionInterp( t, 8 ) + crossSectionInterp( t, 5 ) -
                 crossSectionInterp( t, 6 ))/2;
       break;
     case 8:
       result = crossSectionInterp( t, 6 );
       break;
     case 9:
       result = crossSectionInterp( t, 9 );
       break;
     case 10:
       result = (crossSectionInterp( t, 14 ) + crossSectionInterp( t, 17 ))/2;
       break;
     case 11:
       result = (crossSectionInterp( t, 15 ) + crossSectionInterp( t, 18 ))/2;
       break;
     case 12:
       result = crossSectionInterp( t, 16 )/2;
       break;
     case 13:
       result = (crossSectionInterp( t, 14 ) + crossSectionInterp( t, 15 ) +
                 crossSectionInterp( t, 17 ) + crossSectionInterp( t, 18 ) +
                 crossSectionInterp( t, 16 ))/2;
       break;
     case 17:
       result = crossSectionInterp( t, 10 ) + crossSectionInterp( t, 11 );
       break;
     case 20:
       result = crossSectionInterp( t, 12 );
       break;
     case 21:
       result = 2*crossSectionInterp( t, 12 ) + crossSectionInterp( t, 13 );
       break;
    }
    return result;
  }
 
 void G4PreEquilibrium::slqek( G4int &i, G4int &ms, G4int &mq, G4int &ksi, G4int &me,
                               const G4int lin, const G4int msin, const G4int mqin,
                               const G4int mein, const G4int ln, const G4int msn,
                               const G4int mqn, const G4int men )
  {
    // form cross-sections type
    //
    ms = msin + msn;
    i = lin + ln;    // this ln is not the data member ln
    mq = mqin + mqn;
    me = mein + men;
    if( ms != 0 )return; // strange particle
    if( i > 0 )
      ksi = 1;
    else if( mq > 1 )
      me != 1 ? ksi = 1 : ksi = 2;
    else if( me == 2 )
      ksi = 1;
    else if( me == -1 )
      ksi = 1;
    else if( me != 0 )
      mein ==  1 ? ksi = 2 : ksi = 3;
    else
      mein == -1 ? ksi = 2 : ksi = 3;
    return;
  }
 
 G4double G4PreEquilibrium::costa( const G4int j, const G4double t )
  {
    // cosinus calculation
    //
    const G4double ankj[29][4][4] = 
    { { { 2.7404, -9.6998,  10.400,  2.3882 },
        {-7.5137,  44.096, -74.379,  46.038 },
        { 7.5479, -39.274,  64.835, -41.609 },
        {-1.8369,  8.6911, -13.060,  7.1880 } },
      { {-30.853,  106.24, -129.39,  54.339 },
        { 19.465, -68.102,  96.358, -56.827 },
        {-3.4831,  12.341, -18.592,  12.024 },
        { 0.18941, -0.67880,  1.0665, -0.72910 } },
      { { 0.10258, -1.0542,  11.389, -16.638 },
        {-0.49607,  11.800, -90.857,  164.76 },
        { 1.5437, -33.769,  251.92, -450.71 },
        {-1.2021,  25.336, -186.58,  332.54 } },
      { { 0.15789,  2.9671, -5.5251,  6.8925 },
        {-7.0218, -205.34,  569.51, -898.58 },
        { 134.96,  4872.2, -14674.,  23924. },
        {-821.16, -32586.,  100980., -165530. } },
      { { 0.31531, -7.4981,  43.295, -76.360 },
        {-6.5373,  193.07, -1018.1,  1742.6 },
        { 46.864, -1303.0,  6729.1, -11075. },
        {-95.192,  2637.3, -12857.,  20294. } },
      { {-17.953,  109.72, -239.54,  228.26 },
        { 91.968, -519.63,  1126.6, -1074.0 },
        {-132.70,  741.12, -1600.0,  1524.9 },
        { 58.598, -318.74,  677.51, -640.11 } },
      { { 0.42169,  147.05, -653.35,  915.07 },
        {-3.5198, -260.19,  1225.0, -1748.1 },
        { 3.6373,  155.92, -752.01,  1079.6 },
        {-0.78041, -30.563,  147.95, -212.50 } },
      { {-0.38288,  3.7587, -6.5144,  6.7740 },
        { 103.81, -272.82,  477.59, -512.22 },
        {-1788.2,  4305.2, -7931.4,  9347.1 },
        { 7147.5, -3339.5, -4139.2, -4436.4 } },
      { { 0.24991,  32.028, -118.82,  150.99 },
        {-2.6994, -460.45,  1895.9, -2519.0 },
        { 16.268,  2138.4, -9126.2,  12431. },
        {-29.654, -3182.3,  13944., -19342. } },
      { { 3.9025, -91.126,  323.73, -400.48 },
        {-20.619,  491.70, -1715.5,  2114.3 },
        { 33.004, -766.84,  2700.3, -3352.5 },
        {-16.367,  373.94, -1320.2,  1642.3 } },
      { { 19.402, -224.46,  747.33, -935.70 },
        {-44.180,  471.94, -1485.6,  1805.5 },
        { 31.567, -301.76,  907.63, -1077.3 },
        {-6.8648,  60.476, -175.20,  203.81 } },
      { { 0.40693, -4.1404,  14.044, -17.265 },
        {-3.6799,  59.610, -162.69,  188.73 },
        { 14.556, -175.50,  458.39, -533.90 },
        {-12.621,  149.64, -381.18,  451.41 } },
      { {-0.47554,  2.2641, -12.528,  24.647 },
        { 5.1620, -9.9236,  55.623, -104.62 },
        {-8.1117,  19.315, -84.255,  139.08 },
        { 3.5187, -9.1783,  34.950, -51.243 } },
      { { 0.48173,  5.7726, -13.745,  27.125 },
        {-4.4804, -38.582,  111.59, -243.05 },
        { 16.306,  110.46, -330.45,  722.70 },
        {-15.968, -80.140,  246.16, -607.53 } },
      { {-5.1646, -6.0776,  78.989, -107.05 },
        { 21.871,  56.915, -401.59,  512.15 },
        {-27.993, -94.670,  569.28, -696.21 },
        { 11.587,  45.998, -245.6,  284.52 } },
      { {-53.067,  576.12, -1543.8,  164550. },
        { 147.50, -1638.,  4592.3, -4994.9 },
        {-134.36,  1578.0, -4446.3,  4902.2 },
        { 40.253, -488.60,  1400.1, -1560.6 } },
      { { 0.14988,  2.8753, -5.3078,  6.2233 },
        {-5.9558, -162.03,  430.79, -625.48 },
        { 128.75,  3140.2, -7918.9,  10983. },
        {-851.61, -18780.,  44607., -58790. } },
      { { 0.53689, -13.216,  81.011, -142.85 },
        {-10.550,  296.29, -1695.7,  2893.5 },
        { 69.621, -1924.5,  10620., -17468. },
        {-138.65,  3928.1, -20293.,  32058. } },
      { { 0.65288,  0.38977,  0.84078,  0.18893 },
        {-4.3964,  34.309, -73.692,  84.308 },
        { 14.889, -143.80,  312.27, -350.14 },
        {-15.658,  171.60, -372.12,  412.99 } },
      { { 0.085591,  5.0390, -13.782,  14.661 },
        { 0.054284, -9.2324,  36.397, -42.962 },
        {-0.051111,  4.6003, -20.534,  27.731 },
        { 0.0074514, -0.62529,  2.9159, -4.1101 } },
      { { 0.071622,  3.0960, -11.125,  18.130 },
        { 0.092581, -3.2186,  20.273, -33.245 },
        {-0.051531,  0.89886, -7.5084,  13.188 },
        { 0.0058258, -0.0017288,  0.70224, -1.4854 } },
      { { 0.082300,  0.15854,  3.7716, -4.0562 },
        { 0.010802, -0.33688,  1.1727, -0.67476 },
        {-0.0021798,  0.052166, -0.25816,  0.32048 },
        { 0.000065764, -0.0014711,  0.0078209, -0.010580 } },
      { { 0.11138,  0.60396,  3.0174, -4.4190 },
        {-0.017709,  0.23015, -1.8187,  3.4518 },
        { 0.0020977, -0.025458,  0.21626, -0.40692 },
        {-0.000054799,  0.00059111, -0.0055552,  0.010647 } },
      { { 0.17288,  7.1080, -17.961,  16.403 },
        {-0.14504, -13.032,  41.781, -40.799 },
        { 0.045390,  8.3515, -30.260,  32.882 },
        {-0.0047961, -1.4095,  5.3505, -6.0946 } },
      { { 0.037596,  1.4331, -3.1350,  6.4864 },
        { 0.23827,  1.8253,  1.7648, -16.735 },
        {-0.15410, -1.5201, -1.5692,  17.185 },
        { 0.025037,  0.30588,  0.32520, -3.5277 } },
      { { 0.12489,  1.3573,  0.82338, -1.4595 },
        {-0.051577, -0.35778, -1.1690,  1.8078 },
        { 0.0074864,  0.032888,  0.23744, -0.39802 },
        {-0.00029880, -0.00075117, -0.011402,  0.019505 } },
      { { 0.18470,  1.9269, -3.2979,  3.6843 },
        {-0.073932,  0.27213,  1.0600, -2.3354 },
        { 0.018907, -0.056473, -0.16487,  0.38426 },
        {-0.00092984,  0.0025506,  0.0073052, -0.017220 } },
      { {-1.0306,  32.849, -75.052,  60.255 },
        { 7.9586, -125.72,  256.04, -165.47 },
        {-14.797,  165.90, -279.91,  113.33 },
        { 8.2309, -67.871,  85.762,  5.9727 } },
      { {-237.22,  968.90, -1621.9,  1363.7 },
        { 658.00, -2694.1,  4548.0, -3846.0 },
        {-606.53,  2498.3, -4249.8,  3613.6 },
        { 186.04, -769.33,  1316.6, -1124.2 } } };
    G4double ank[4][4];
    for( G4int k = 0; k < 4; ++k )
    {
      for( G4int n = 0; n < 4; ++n )ank[n][k] = ankj[n][k][j];
    }
    G4double s1 = 0;
    G4double s2 = 0;
    G4double r = G4UniformRand();
    for( G4int n = 0; n < 4; ++n )
    {
      for( G4int k = 0; k < 4; ++k )
      {
        s1 += ank[n][k]*pow(t,k)*pow(r,n);
        s2 += ank[n][k]*pow(t,k);
      }
    }
    G4double cta = 2*sqrt(r)*(s1+(1-s2)*pow(r,4)) - 1;
    if( cta < -1 )cta = -1;
    else if( cta > 1 )cta = 1;
    return cta;
  }
 
 G4double G4PreEquilibrium::cosex( const G4int i, const G4double t, const G4double cm )
  {
    // cosinus calculation for charge exchange scattering
    //
    G4double result;
    if( i != 0 )
    {
      if( t <= 0.51 )result = costa( 13, t );
      else
        t <= 1 ? result = costa( 14, t ) : result = costa( 15, t );
    }
    else
    {
      if( t <= 0.08 )result = costa( 16, t );
      else
      {
        if( t <= 0.3 )result = costa( 17, t );
        else
        {
          if( t <= 1 )result = costa( 9, t );
          else
          {
            if( t <= 2.4 )result = costa( 10, t );
            else
            {
              G4double tmax = sqrt(t*(t+2*cm));
              result = 1+(2*log(1+G4UniformRand()*(exp(-7.5*tmax)-1)))/(7.5*tmax);
            }
          }
        }
      }
    }
    return result;
  }
 
 G4double G4PreEquilibrium::cosel( const G4int i, const G4int mq, const G4int ksi,
                                  const G4double t, const G4double cm )
  {
    // cosinus calculation for elastic scattering
    //
    G4double result;
    if( i != 0 )
    {
      t <= 0.45 ? result = costa( 11, t ) : result = costa( 12, t );
    }
    else
    {
      if( mq >= 2 )
      {
        if( ksi >= 2 )
        {
          if( t <= 0.97 )result = costa( 2, t );
          else
          {
            if( t <= 2.8 )result = (1+costa( 0, t ))/2;
            else
            {
              if( t <= 10 )result = (3-costa( 1, t ))/4;
              else
              {
                G4double tmax = sqrt(t*(t+2*cm));
                result = 1+(2*log(1+G4UniformRand()*(exp(-8.7*tmax)-1)))/(8.7*tmax);
              }
            }
          }
        }
        else
        {
          if( t <= 0.46 )result = 1-2*G4UniformRand();
          else
          {
            if( t <= 2.8 )result = (1+costa( 0, t ))/2;
            else
            {
              if( t <= 10 )result = (3-costa( 1, t ))/4;
              else
              {
                G4double tmax = sqrt(t*(t+2*cm));
                result = 1+(2*log(1+G4UniformRand()*(exp(-8.7*tmax)-1)))/(8.7*tmax);
              }
            }
          }
        }
      }
      else
      {
        if( ksi < 2 )
        {
          if( t <= 0.08 )result = costa( 3, t );
          else
          {
            if( t <= 0.3 )result = costa( 4, t );
            else
            {
              if( t <= 1 )result = costa( 5, t );
              else
              {
                if( t <= 2.4 )result = costa( 6, t );
                else
                {
                  G4double tmax = sqrt(t*(t+2*cm));
                  result = 1+(2*log(1+G4UniformRand()*(exp(-7.5*tmax)-1)))/(7.5*tmax);
                }
              }
            }
          }      
        }
        else if( ksi == 2 )
        {
          if( t <= 0.08 )result = costa( 7, t );
          else
          {
            if( t <= 0.3 )result = costa( 8, t );
            else
            {
              if( t <= 1 )result = costa( 9, t );
              else
              {
                if( t <= 2.4 )result = costa( 10, t );
                else
                {
                  G4double tmax = sqrt(t*(t+2*cm));
                  result = 1+(2*log(1+G4UniformRand()*(exp(-7.5*tmax)-1)))/(7.5*tmax);
                }
              }
            }
          }
        }
        else if( ksi > 2 )
        {
          if( G4UniformRand() <= 0.5 )
          {
            if( t <= 0.08 )result = costa( 3, t );
            else
            {
              if( t <= 0.3 )result = costa( 4, t );
              else
              {
                if( t <= 1 )result = costa( 5, t );
                else
                {
                  if( t <= 2.4 )result = costa( 6, t );
                  else
                  {
                    G4double tmax = sqrt(t*(t+2*cm));
                    result = 1+(2*log(1+G4UniformRand()*(exp(-7.5*tmax)-1)))/(7.5*tmax);
                  }
                }
              }
            }
          }
          else
          {
            if( t <= 0.08 )result = costa( 7, t );
            else
            {
              if( t <= 0.3 )result = costa( 8, t );
              else
              {
                if( t <= 1 )result = costa( 9, t );
                else
                {
                  if( t <= 2.4 )result = costa( 10, t );
                  else
                  {
                    G4double tmax = sqrt(t*(t+2*cm));
                    result = 1+(2*log(1+G4UniformRand()*(exp(-7.5*tmax)-1)))/(7.5*tmax);
                  }
                }
              }
            }
          }
        }
      }
    }
    return result;
  }

 // end of file
