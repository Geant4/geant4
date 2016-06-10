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
// $Id: G4RPGReaction.cc 94406 2015-11-13 14:52:40Z gcosmo $
//

#include <iostream>

#include "G4RPGReaction.hh"
#include "G4Log.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

G4bool G4RPGReaction::
ReactionStage(const G4HadProjectile* /*originalIncident*/,
              G4ReactionProduct& /*modifiedOriginal*/,
              G4bool& /*incidentHasChanged*/,
              const G4DynamicParticle* /*originalTarget*/,
              G4ReactionProduct& /*targetParticle*/,
              G4bool& /*targetHasChanged*/,
              const G4Nucleus& /*targetNucleus*/,
              G4ReactionProduct& /*currentParticle*/,
              G4FastVector<G4ReactionProduct,256>& /*vec*/,
              G4int& /*vecLen*/,
              G4bool /*leadFlag*/,
              G4ReactionProduct& /*leadingStrangeParticle*/)
{
  G4cout << " G4RPGReactionStage must be overridden in a derived class " 
         << G4endl;
  return false;
}


void G4RPGReaction::
AddBlackTrackParticles(const G4double epnb,            // GeV
                       const G4int npnb,
                       const G4double edta,            // GeV
                       const G4int ndta,
                       const G4ReactionProduct& modifiedOriginal,
                       G4int PinNucleus,
                       G4int NinNucleus,
                       const G4Nucleus& targetNucleus,
                       G4FastVector<G4ReactionProduct,256>& vec,
                       G4int& vecLen)
{
  // derived from original FORTRAN code in GENXPT and TWOCLU by H. Fesefeldt
  //
  // npnb is number of proton/neutron black track particles
  // ndta is the number of deuterons, tritons, and alphas produced
  // epnb is the kinetic energy available for proton/neutron black track particles
  // edta is the kinetic energy available for deuteron/triton/alpha particles

  G4ParticleDefinition* aProton = G4Proton::Proton();
  G4ParticleDefinition* aNeutron = G4Neutron::Neutron();
  G4ParticleDefinition* aDeuteron = G4Deuteron::Deuteron();
  G4ParticleDefinition* aTriton = G4Triton::Triton();
  G4ParticleDefinition* anAlpha = G4Alpha::Alpha();
    
  const G4double ekOriginal = modifiedOriginal.GetKineticEnergy()/MeV;
  const G4double atomicWeight = targetNucleus.GetA_asInt();
  const G4double atomicNumber = targetNucleus.GetZ_asInt();
    
  const G4double ika1 = 3.6;
  const G4double ika2 = 35.56;
  const G4double ika3 = 6.45;
    
  G4int i;
  G4double pp;
  G4double kinetic = 0;
  G4double kinCreated = 0;
  //  G4double cfa = 0.025*((atomicWeight-1.0)/120.0) * std::exp(-(atomicWeight-1.0)/120.0);
  G4double remainingE = 0;

  // First add protons and neutrons to final state
  if (npnb > 0) {
    //    G4double backwardKinetic = 0.0;
    G4int local_npnb = npnb;
    // DHW: does not conserve energy  for (i = 0; i < npnb; ++i) if (G4UniformRand() < sprob) local_npnb--;
    local_npnb = std::min(PinNucleus + NinNucleus , local_npnb);
    G4double local_epnb = epnb;
    if (ndta == 0) local_epnb += edta;   // Retrieve unused kinetic energy
    //    G4double ekin = local_epnb/std::max(1,local_npnb);

    remainingE = local_epnb;
    for (i = 0; i < local_npnb; ++i)
    {
      G4ReactionProduct* p1 = new G4ReactionProduct();
      //      if( backwardKinetic > local_epnb ) {
      //        delete p1;
      //        break;    
      //      }

      //      G4double ran = G4UniformRand();
      //      G4double kinetic = -ekin*std::log(ran) - cfa*(1.0+0.5*normal());
      //      if( kinetic < 0.0 )kinetic = -0.010*std::log(ran);
      //      backwardKinetic += kinetic;
      //      if( backwardKinetic > local_epnb )
      // kinetic = std::max( kineticMinimum, local_epnb-(backwardKinetic-kinetic) );

      if (G4UniformRand() > (1.0-atomicNumber/atomicWeight)) {

	// Boil off a proton if there are any left, otherwise a neutron

        if (PinNucleus > 0) {
          p1->SetDefinition( aProton );
          PinNucleus--;
	} else {
          p1->SetDefinition( aNeutron );
          NinNucleus--;
	  //        } else {
	  //          delete p1;
	  //          break;     // no nucleons left in nucleus
        }
      } else {

	// Boil off a neutron if there are any left, otherwise a proton

        if (NinNucleus > 0) {
          p1->SetDefinition( aNeutron );
          NinNucleus--;
	} else {
          p1->SetDefinition( aProton );
          PinNucleus--;
	  //        } else {
	  //          delete p1;
	  //          break;     // no nucleons left in nucleus
        }
      }

      if (i < local_npnb - 1) {
        kinetic = remainingE*G4UniformRand();
        remainingE -= kinetic;
      } else {
        kinetic = remainingE;
      }

      vec.SetElement( vecLen, p1 );
      G4double cost = G4UniformRand() * 2.0 - 1.0;
      G4double sint = std::sqrt(std::fabs(1.0-cost*cost));
      G4double phi = twopi * G4UniformRand();
      vec[vecLen]->SetNewlyAdded( true );
      vec[vecLen]->SetKineticEnergy( kinetic*GeV );
      kinCreated+=kinetic;
      pp = vec[vecLen]->GetTotalMomentum();
      vec[vecLen]->SetMomentum(pp*sint*std::sin(phi),
                               pp*sint*std::cos(phi),
                               pp*cost );
      vecLen++;
    }

    if (NinNucleus > 0) {
      if( (atomicWeight >= 10.0) && (ekOriginal <= 2.0*GeV) )
      {
        G4double ekw = ekOriginal/GeV;
        G4int ika, kk = 0;
        if( ekw > 1.0 )ekw *= ekw;
        ekw = std::max( 0.1, ekw );
        ika = G4int(ika1*G4Exp((atomicNumber*atomicNumber/
                                           atomicWeight-ika2)/ika3)/ekw);
        if( ika > 0 )
        {
          for( i=(vecLen-1); i>=0; --i )
          {
            if( (vec[i]->GetDefinition() == aProton) && vec[i]->GetNewlyAdded() )
            {
              vec[i]->SetDefinitionAndUpdateE( aNeutron );  // modified 22-Oct-97
              PinNucleus++;
              NinNucleus--;
              if( ++kk > ika )break;
            }
          }
        }
      }
    } // if (NinNucleus >0)
  } // if (npnb > 0)

  //  Next try to add deuterons, tritons and alphas to final state

  G4double ran = 0;
  if (ndta > 0) {
    //    G4double backwardKinetic = 0.0;
    G4int local_ndta=ndta;
    // DHW: does not conserve energy  for (i = 0; i < ndta; ++i) if (G4UniformRand() < sprob) local_ndta--;
    G4double local_edta = edta;
    if (npnb == 0) local_edta += epnb;  // Retrieve unused kinetic energy
    //    G4double ekin = local_edta/std::max(1,local_ndta);

    remainingE = local_edta;
    for (i = 0; i < local_ndta; ++i) {
      G4ReactionProduct* p2 = new G4ReactionProduct();
      //        if( backwardKinetic > local_edta ) {
      //          delete p2;
      //          break;
      //        }

      //	G4double ran = G4UniformRand();
      //	G4double kinetic = -ekin*std::log(ran)-cfa*(1.+0.5*normal());
      //	if( kinetic < 0.0 )kinetic = kineticFactor*std::log(ran);
      //        backwardKinetic += kinetic;
      //        if( backwardKinetic > local_edta )kinetic = local_edta-(backwardKinetic-kinetic);
      //        if( kinetic < 0.0 )kinetic = kineticMinimum;

      ran = G4UniformRand();
      if (ran < 0.60) {
        if (PinNucleus > 0 && NinNucleus > 0) {
          p2->SetDefinition( aDeuteron );
          PinNucleus--;
          NinNucleus--;
        } else if (NinNucleus > 0) {
          p2->SetDefinition( aNeutron );
          NinNucleus--;
        } else if (PinNucleus > 0) {
          p2->SetDefinition( aProton );
          PinNucleus--;
        } else {
          delete p2;
          break;
        }
      } else if (ran < 0.90) {
        if (PinNucleus > 0 && NinNucleus > 1) {
          p2->SetDefinition( aTriton );
          PinNucleus--;
          NinNucleus -= 2;
        } else if (PinNucleus > 0 && NinNucleus > 0) {
          p2->SetDefinition( aDeuteron );
          PinNucleus--;
          NinNucleus--;
        } else if (NinNucleus > 0) {
          p2->SetDefinition( aNeutron );
          NinNucleus--;
        } else if (PinNucleus > 0) {
          p2->SetDefinition( aProton );
          PinNucleus--;
        } else {
          delete p2;
          break;
        }
      } else {
        if (PinNucleus > 1 && NinNucleus > 1) {
          p2->SetDefinition( anAlpha );
          PinNucleus -= 2;
          NinNucleus -= 2;
	} else if (PinNucleus > 0 && NinNucleus > 1) {
          p2->SetDefinition( aTriton );
          PinNucleus--;
          NinNucleus -= 2;
        } else if (PinNucleus > 0 && NinNucleus > 0) {
          p2->SetDefinition( aDeuteron );
          PinNucleus--;
          NinNucleus--;
        } else if (NinNucleus > 0) {
          p2->SetDefinition( aNeutron );
          NinNucleus--;
        } else if (PinNucleus > 0) {
          p2->SetDefinition( aProton );
          PinNucleus--;
        } else {
          delete p2;
          break;
        }
      }

      if (i < local_ndta - 1) {
        kinetic = remainingE*G4UniformRand();
        remainingE -= kinetic;
      } else {
        kinetic = remainingE;
      }

      vec.SetElement( vecLen, p2 );
      G4double cost = 2.0*G4UniformRand() - 1.0;
      G4double sint = std::sqrt(std::max(0.0,(1.0-cost*cost)));
      G4double phi = twopi*G4UniformRand();
      vec[vecLen]->SetNewlyAdded( true );
      vec[vecLen]->SetKineticEnergy( kinetic*GeV );
      kinCreated+=kinetic;

      pp = vec[vecLen]->GetTotalMomentum();
      vec[vecLen]->SetMomentum( pp*sint*std::sin(phi),
                                pp*sint*std::cos(phi),
                                pp*cost );
      vecLen++;
    }
  } // if (ndta > 0)
}

 
G4double 
G4RPGReaction::GenerateNBodyEvent(const G4double totalEnergy,       // MeV
                                  const G4bool constantCrossSection,
                                  G4FastVector<G4ReactionProduct,256>& vec,
                                  G4int &vecLen)
{
  G4int i;
  const G4double expxu =  82.;           // upper bound for arg. of exp
  const G4double expxl = -expxu;         // lower bound for arg. of exp

  if (vecLen < 2) {
    G4cerr << "*** Error in G4RPGReaction::GenerateNBodyEvent" << G4endl;
    G4cerr << "    number of particles < 2" << G4endl;
    G4cerr << "totalEnergy = " << totalEnergy << "MeV, vecLen = " << vecLen << G4endl;
    return -1.0;
  }

  G4double mass[18];    // mass of each particle
  G4double energy[18];  // total energy of each particle
  G4double pcm[3][18];           // pcm is an array with 3 rows and vecLen columns
    
  G4double totalMass = 0.0;
  G4double extraMass = 0;
  G4double sm[18];
    
  for (i=0; i<vecLen; ++i) {
    mass[i] = vec[i]->GetMass()/GeV;
    if(vec[i]->GetSide() == -2) extraMass+=vec[i]->GetMass()/GeV;
    vec[i]->SetMomentum( 0.0, 0.0, 0.0 );
    pcm[0][i] = 0.0;      // x-momentum of i-th particle
    pcm[1][i] = 0.0;      // y-momentum of i-th particle
    pcm[2][i] = 0.0;      // z-momentum of i-th particle
    energy[i] = mass[i];  // total energy of i-th particle
    totalMass += mass[i];
    sm[i] = totalMass;
  }

  G4double totalE = totalEnergy/GeV;
  if (totalMass > totalE) {
    //G4cerr << "*** Error in G4RPGReaction::GenerateNBodyEvent" << G4endl;
    //G4cerr << "    total mass (" << totalMass*GeV << "MeV) > total energy ("
    //     << totalEnergy << "MeV)" << G4endl;
    totalE = totalMass;
    return -1.0;
  }

  G4double kineticEnergy = totalE - totalMass;
  G4double emm[18];
  emm[0] = mass[0];
  emm[vecLen-1] = totalE;

  if (vecLen > 2) {          // the random numbers are sorted
    G4double ran[18];
    for( i=0; i<vecLen; ++i )ran[i] = G4UniformRand();
    for (i=0; i<vecLen-2; ++i) {
      for (G4int j=vecLen-2; j>i; --j) {
        if (ran[i] > ran[j]) {
          G4double temp = ran[i];
          ran[i] = ran[j];
          ran[j] = temp;
        }
      }
    }
    for( i=1; i<vecLen-1; ++i )emm[i] = ran[i-1]*kineticEnergy + sm[i];
  }

  // Weight is the sum of logarithms of terms instead of the product of terms

  G4bool lzero = true;    
  G4double wtmax = 0.0;
  if (constantCrossSection) {
    G4double emmax = kineticEnergy + mass[0];
    G4double emmin = 0.0;
      for( i=1; i<vecLen; ++i )
      {
        emmin += mass[i-1];
        emmax += mass[i];
        G4double wtfc = 0.0;
        if( emmax*emmax > 0.0 )
        {
          G4double arg = emmax*emmax
            + (emmin*emmin-mass[i]*mass[i])*(emmin*emmin-mass[i]*mass[i])/(emmax*emmax)
            - 2.0*(emmin*emmin+mass[i]*mass[i]);
          if( arg > 0.0 )wtfc = 0.5*std::sqrt( arg );
        }
        if( wtfc == 0.0 )
        {
          lzero = false;
          break;
        }
        wtmax += G4Log( wtfc );
      }
      if( lzero )
        wtmax = -wtmax;
      else
        wtmax = expxu;
  } else {
      //   ffq(n) = pi*(2*pi)^(n-2)/(n-2)!
      const G4double ffq[18] = { 0., 3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
                                 256.3704, 268.4705, 240.9780, 189.2637,
                                 132.1308,  83.0202,  47.4210,  24.8295,
                                 12.0006,   5.3858,   2.2560,   0.8859 };
      wtmax = G4Log( std::pow( kineticEnergy, vecLen-2 ) * ffq[vecLen-1] / totalE );
  }

  // Calculate momenta for secondaries 

  lzero = true;
  G4double pd[50];

    for( i=0; i<vecLen-1; ++i )
    {
      pd[i] = 0.0;
      if( emm[i+1]*emm[i+1] > 0.0 )
      {
        G4double arg = emm[i+1]*emm[i+1]
          + (emm[i]*emm[i]-mass[i+1]*mass[i+1])*(emm[i]*emm[i]-mass[i+1]*mass[i+1])
            /(emm[i+1]*emm[i+1])
          - 2.0*(emm[i]*emm[i]+mass[i+1]*mass[i+1]);
        if( arg > 0.0 )pd[i] = 0.5*std::sqrt( arg );
      }
      if( pd[i] <= 0.0 )    //  changed from  ==  on 02 April 98
        lzero = false;
      else
        wtmax += G4Log( pd[i] );
    }
    G4double weight = 0.0;           // weight is returned by GenerateNBodyEvent
    if( lzero )weight = G4Exp( std::max(std::min(wtmax,expxu),expxl) );
    
    G4double bang, cb, sb, s0, s1, s2, c, esys, a, b, gama, beta;
    G4double ss;
    pcm[0][0] = 0.0;
    pcm[1][0] = pd[0];
    pcm[2][0] = 0.0;
    for( i=1; i<vecLen; ++i )
    {
      pcm[0][i] = 0.0;
      pcm[1][i] = -pd[i-1];
      pcm[2][i] = 0.0;
      bang = twopi*G4UniformRand();
      cb = std::cos(bang);
      sb = std::sin(bang);
      c = 2.0*G4UniformRand() - 1.0;
      ss = std::sqrt( std::fabs( 1.0-c*c ) );
      if( i < vecLen-1 )
      {
        esys = std::sqrt(pd[i]*pd[i] + emm[i]*emm[i]);
        beta = pd[i]/esys;
        gama = esys/emm[i];
        for( G4int j=0; j<=i; ++j )
        {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = std::sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*ss;                           //  rotation
          pcm[1][j] = s0*ss + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
          pcm[1][j] = gama*(pcm[1][j] + beta*energy[j]);
        }
      }
      else
      {
        for( G4int j=0; j<=i; ++j )
        {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = std::sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*s;                           //  rotation
          pcm[1][j] = s0*ss + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
        }
      }
    }

  for (i=0; i<vecLen; ++i) {
    vec[i]->SetMomentum( pcm[0][i]*GeV, pcm[1][i]*GeV, pcm[2][i]*GeV );
    vec[i]->SetTotalEnergy( energy[i]*GeV );
  }

  return weight;
}

 
G4double 
G4RPGReaction::GenerateNBodyEventT(const G4double totalEnergy,
                                   const G4bool constantCrossSection,
                                   std::vector<G4ReactionProduct*>& tempList)
{
  G4int i;
  const G4double expxu =  82.;           // upper bound for arg. of exp
  const G4double expxl = -expxu;         // lower bound for arg. of exp
  G4int listLen = tempList.size();

  if (listLen < 2) {
    G4cerr << "*** Error in G4RPGReaction::GenerateNBodyEvent" << G4endl;
    G4cerr << "    number of particles < 2" << G4endl;
    G4cerr << "totalEnergy = " << totalEnergy << "MeV, listLen = " << listLen << G4endl;
    return -1.0;
  }

  G4double mass[18];    // mass of each particle
  G4double energy[18];  // total energy of each particle
  G4double pcm[3][18];           // pcm is an array with 3 rows and listLen columns
    
  G4double totalMass = 0.0;
  G4double extraMass = 0;
  G4double sm[18];
    
  for (i=0; i<listLen; ++i) {
    mass[i] = tempList[i]->GetMass()/GeV;
    if(tempList[i]->GetSide() == -2) extraMass+=tempList[i]->GetMass()/GeV;
    tempList[i]->SetMomentum( 0.0, 0.0, 0.0 );
    pcm[0][i] = 0.0;      // x-momentum of i-th particle
    pcm[1][i] = 0.0;      // y-momentum of i-th particle
    pcm[2][i] = 0.0;      // z-momentum of i-th particle
    energy[i] = mass[i];  // total energy of i-th particle
    totalMass += mass[i];
    sm[i] = totalMass;
  }

  G4double totalE = totalEnergy/GeV;
  if (totalMass > totalE) {
    totalE = totalMass;
    return -1.0;
  }

  G4double kineticEnergy = totalE - totalMass;
  G4double emm[18];
  emm[0] = mass[0];
  emm[listLen-1] = totalE;

  if (listLen > 2) {          // the random numbers are sorted
    G4double ran[18];
    for( i=0; i<listLen; ++i )ran[i] = G4UniformRand();
    for (i=0; i<listLen-2; ++i) {
      for (G4int j=listLen-2; j>i; --j) {
        if (ran[i] > ran[j]) {
          G4double temp = ran[i];
          ran[i] = ran[j];
          ran[j] = temp;
        }
      }
    }
    for( i=1; i<listLen-1; ++i )emm[i] = ran[i-1]*kineticEnergy + sm[i];
  }

  // Weight is the sum of logarithms of terms instead of the product of terms

  G4bool lzero = true;    
  G4double wtmax = 0.0;
  if (constantCrossSection) {
    G4double emmax = kineticEnergy + mass[0];
    G4double emmin = 0.0;
      for( i=1; i<listLen; ++i )
      {
        emmin += mass[i-1];
        emmax += mass[i];
        G4double wtfc = 0.0;
        if( emmax*emmax > 0.0 )
        {
          G4double arg = emmax*emmax
            + (emmin*emmin-mass[i]*mass[i])*(emmin*emmin-mass[i]*mass[i])/(emmax*emmax)
            - 2.0*(emmin*emmin+mass[i]*mass[i]);
          if( arg > 0.0 )wtfc = 0.5*std::sqrt( arg );
        }
        if( wtfc == 0.0 )
        {
          lzero = false;
          break;
        }
        wtmax += G4Log( wtfc );
      }
      if( lzero )
        wtmax = -wtmax;
      else
        wtmax = expxu;
  } else {
      //   ffq(n) = pi*(2*pi)^(n-2)/(n-2)!
      const G4double ffq[18] = { 0., 3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
                                 256.3704, 268.4705, 240.9780, 189.2637,
                                 132.1308,  83.0202,  47.4210,  24.8295,
                                 12.0006,   5.3858,   2.2560,   0.8859 };
      wtmax = G4Log( std::pow( kineticEnergy, listLen-2 ) * ffq[listLen-1] / totalE );
  }

  // Calculate momenta for secondaries 

  lzero = true;
  G4double pd[50];

    for( i=0; i<listLen-1; ++i )
    {
      pd[i] = 0.0;
      if( emm[i+1]*emm[i+1] > 0.0 )
      {
        G4double arg = emm[i+1]*emm[i+1]
          + (emm[i]*emm[i]-mass[i+1]*mass[i+1])*(emm[i]*emm[i]-mass[i+1]*mass[i+1])
            /(emm[i+1]*emm[i+1])
          - 2.0*(emm[i]*emm[i]+mass[i+1]*mass[i+1]);
        if( arg > 0.0 )pd[i] = 0.5*std::sqrt( arg );
      }
      if( pd[i] <= 0.0 )    //  changed from  ==  on 02 April 98
        lzero = false;
      else
        wtmax += G4Log( pd[i] );
    }
    G4double weight = 0.0;           // weight is returned by GenerateNBodyEvent
    if( lzero )weight = G4Exp( std::max(std::min(wtmax,expxu),expxl) );
    
    G4double bang, cb, sb, s0, s1, s2, c, esys, a, b, gama, beta;
    G4double ss;
    pcm[0][0] = 0.0;
    pcm[1][0] = pd[0];
    pcm[2][0] = 0.0;
    for( i=1; i<listLen; ++i )
    {
      pcm[0][i] = 0.0;
      pcm[1][i] = -pd[i-1];
      pcm[2][i] = 0.0;
      bang = twopi*G4UniformRand();
      cb = std::cos(bang);
      sb = std::sin(bang);
      c = 2.0*G4UniformRand() - 1.0;
      ss = std::sqrt( std::fabs( 1.0-c*c ) );
      if( i < listLen-1 )
      {
        esys = std::sqrt(pd[i]*pd[i] + emm[i]*emm[i]);
        beta = pd[i]/esys;
        gama = esys/emm[i];
        for( G4int j=0; j<=i; ++j )
        {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = std::sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*ss;                           //  rotation
          pcm[1][j] = s0*ss + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
          pcm[1][j] = gama*(pcm[1][j] + beta*energy[j]);
        }
      }
      else
      {
        for( G4int j=0; j<=i; ++j )
        {
          s0 = pcm[0][j];
          s1 = pcm[1][j];
          s2 = pcm[2][j];
          energy[j] = std::sqrt( s0*s0 + s1*s1 + s2*s2 + mass[j]*mass[j] );
          a = s0*c - s1*ss;                           //  rotation
          pcm[1][j] = s0*ss + s1*c;
          b = pcm[2][j];
          pcm[0][j] = a*cb - b*sb;
          pcm[2][j] = a*sb + b*cb;
        }
      }
    }

  for (i=0; i<listLen; ++i) {
    tempList[i]->SetMomentum(pcm[0][i]*GeV, pcm[1][i]*GeV, pcm[2][i]*GeV);
    tempList[i]->SetTotalEnergy(energy[i]*GeV);
  }

  return weight;
}


G4double G4RPGReaction::normal()
{
  G4double ran = -6.0;
  for( G4int i=0; i<12; ++i )ran += G4UniformRand();
  return ran;
}

 
void G4RPGReaction::Defs1(const G4ReactionProduct& modifiedOriginal,
                          G4ReactionProduct& currentParticle,
                          G4ReactionProduct& targetParticle,
                          G4FastVector<G4ReactionProduct,256>& vec,
                          G4int& vecLen)
{
  // Rotate final state particle momenta by initial particle direction

  const G4double pjx = modifiedOriginal.GetMomentum().x()/MeV;
  const G4double pjy = modifiedOriginal.GetMomentum().y()/MeV;
  const G4double pjz = modifiedOriginal.GetMomentum().z()/MeV;
  const G4double p = modifiedOriginal.GetMomentum().mag()/MeV;

  if (pjx*pjx+pjy*pjy > 0.0) {
    G4double cost, sint, ph, cosp, sinp, pix, piy, piz;
    cost = pjz/p;
    sint = std::sqrt(std::abs((1.0-cost)*(1.0+cost)));
    if( pjy < 0.0 )
      ph = 3*halfpi;
    else
      ph = halfpi;
    if( std::abs( pjx ) > 0.001*MeV )ph = std::atan2(pjy,pjx);
    cosp = std::cos(ph);
    sinp = std::sin(ph);
    pix = currentParticle.GetMomentum().x()/MeV;
    piy = currentParticle.GetMomentum().y()/MeV;
    piz = currentParticle.GetMomentum().z()/MeV;
    currentParticle.SetMomentum((cost*cosp*pix - sinp*piy + sint*cosp*piz)*MeV,
                                (cost*sinp*pix + cosp*piy + sint*sinp*piz)*MeV,
                                (-sint*pix + cost*piz)*MeV);
    pix = targetParticle.GetMomentum().x()/MeV;
    piy = targetParticle.GetMomentum().y()/MeV;
    piz = targetParticle.GetMomentum().z()/MeV;
    targetParticle.SetMomentum((cost*cosp*pix - sinp*piy + sint*cosp*piz)*MeV,
                               (cost*sinp*pix + cosp*piy + sint*sinp*piz)*MeV,
                               (-sint*pix + cost*piz)*MeV);

    for (G4int i=0; i<vecLen; ++i) {
      pix = vec[i]->GetMomentum().x()/MeV;
      piy = vec[i]->GetMomentum().y()/MeV;
      piz = vec[i]->GetMomentum().z()/MeV;
      vec[i]->SetMomentum((cost*cosp*pix - sinp*piy + sint*cosp*piz)*MeV,
                          (cost*sinp*pix + cosp*piy + sint*sinp*piz)*MeV,
                          (-sint*pix + cost*piz)*MeV);
    }

  } else {
    if (pjz < 0.0) {
      currentParticle.SetMomentum( -currentParticle.GetMomentum().z() );
      targetParticle.SetMomentum( -targetParticle.GetMomentum().z() );
      for (G4int i=0; i<vecLen; ++i) vec[i]->SetMomentum( -vec[i]->GetMomentum().z() );
    }
  }
}


 void G4RPGReaction::Rotate(
  const G4double numberofFinalStateNucleons,
  const G4ThreeVector &temp,
  const G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
  const G4HadProjectile *originalIncident, // original incident particle
  const G4Nucleus &targetNucleus,
  G4ReactionProduct &currentParticle,
  G4ReactionProduct &targetParticle,
  G4FastVector<G4ReactionProduct,256> &vec,
  G4int &vecLen )
  {
    // derived from original FORTRAN code in GENXPT and TWOCLU by H. Fesefeldt
    //
    //   Rotate in direction of z-axis, this does disturb in some way our
    //    inclusive distributions, but it is necessary for momentum conservation
    //
    const G4double atomicWeight = targetNucleus.GetA_asInt();
    const G4double logWeight = G4Log(atomicWeight);
    
    G4ParticleDefinition *aPiMinus = G4PionMinus::PionMinus();
    G4ParticleDefinition *aPiPlus = G4PionPlus::PionPlus();
    G4ParticleDefinition *aPiZero = G4PionZero::PionZero();
    
    G4int i;
    G4ThreeVector pseudoParticle[4];
    for( i=0; i<4; ++i )pseudoParticle[i].set(0,0,0);
    pseudoParticle[0] = currentParticle.GetMomentum()
                        + targetParticle.GetMomentum();
    for( i=0; i<vecLen; ++i )
      pseudoParticle[0] = pseudoParticle[0] + (vec[i]->GetMomentum());
    //
    //  Some smearing in transverse direction from Fermi motion
    //
    G4double pp, pp1;
    G4double alekw, p;
    G4double r1, r2, a1, ran1, ran2, xxh, exh, pxTemp, pyTemp, pzTemp;
    
    r1 = twopi*G4UniformRand();
    r2 = G4UniformRand();
    a1 = std::sqrt(-2.0*G4Log(r2));
    ran1 = a1*std::sin(r1)*0.020*numberofFinalStateNucleons*GeV;
    ran2 = a1*std::cos(r1)*0.020*numberofFinalStateNucleons*GeV;
    G4ThreeVector fermir(ran1, ran2, 0);

    pseudoParticle[0] = pseudoParticle[0]+fermir; // all particles + fermir
    pseudoParticle[2] = temp; // original in cms system
    pseudoParticle[3] = pseudoParticle[0];
    
    pseudoParticle[1] = pseudoParticle[2].cross(pseudoParticle[3]);
    G4double rotation = 2.*pi*G4UniformRand();
    pseudoParticle[1] = pseudoParticle[1].rotate(rotation, pseudoParticle[3]);
    pseudoParticle[2] = pseudoParticle[3].cross(pseudoParticle[1]);    
    for(G4int ii=1; ii<=3; ii++)
    { 
      p = pseudoParticle[ii].mag();
      if( p == 0.0 )
        pseudoParticle[ii]= G4ThreeVector( 0.0, 0.0, 0.0 );
      else
        pseudoParticle[ii]= pseudoParticle[ii] * (1./p);
    }
    
    pxTemp = pseudoParticle[1].dot(currentParticle.GetMomentum());
    pyTemp = pseudoParticle[2].dot(currentParticle.GetMomentum());
    pzTemp = pseudoParticle[3].dot(currentParticle.GetMomentum());
    currentParticle.SetMomentum( pxTemp, pyTemp, pzTemp );
    
    pxTemp = pseudoParticle[1].dot(targetParticle.GetMomentum());
    pyTemp = pseudoParticle[2].dot(targetParticle.GetMomentum());
    pzTemp = pseudoParticle[3].dot(targetParticle.GetMomentum());
    targetParticle.SetMomentum( pxTemp, pyTemp, pzTemp );
    
    for( i=0; i<vecLen; ++i )
    {
      pxTemp = pseudoParticle[1].dot(vec[i]->GetMomentum());
      pyTemp = pseudoParticle[2].dot(vec[i]->GetMomentum());
      pzTemp = pseudoParticle[3].dot(vec[i]->GetMomentum());
      vec[i]->SetMomentum( pxTemp, pyTemp, pzTemp );
    }
    //
    //  Rotate in direction of primary particle, subtract binding energies
    //   and make some further corrections if required
    //
    Defs1( modifiedOriginal, currentParticle, targetParticle, vec, vecLen );
    G4double ekin;
    G4double dekin = 0.0;
    G4double ek1 = 0.0;
    G4int npions = 0;
    if( atomicWeight >= 1.5 )            // self-absorption in heavy molecules
    {
      // corrections for single particle spectra (shower particles)
      //
      const G4double alem[] = { 1.40, 2.30, 2.70, 3.00, 3.40, 4.60, 7.00 };
      const G4double val0[] = { 0.00, 0.40, 0.48, 0.51, 0.54, 0.60, 0.65 };
      alekw = G4Log( originalIncident->GetKineticEnergy()/GeV );
      exh = 1.0;
      if( alekw > alem[0] )   //   get energy bin
      {
        exh = val0[6];
        for( G4int j=1; j<7; ++j )
        {
          if( alekw < alem[j] ) // use linear interpolation/extrapolation
          {
            G4double rcnve = (val0[j] - val0[j-1]) / (alem[j] - alem[j-1]);
            exh = rcnve * alekw + val0[j-1] - rcnve * alem[j-1];
            break;
          }
        }
        exh = 1.0 - exh;
      }
      const G4double cfa = 0.025*((atomicWeight-1.)/120.)*G4Exp(-(atomicWeight-1.)/120.);
      ekin = currentParticle.GetKineticEnergy()/GeV - cfa*(1+normal()/2.0);
      ekin = std::max( 1.0e-6, ekin );
      xxh = 1.0;
      if( (modifiedOriginal.GetDefinition() == aPiPlus ||
           modifiedOriginal.GetDefinition() == aPiMinus) &&
           currentParticle.GetDefinition() == aPiZero &&
           G4UniformRand() <= logWeight) xxh = exh;
      dekin += ekin*(1.0-xxh);
      ekin *= xxh;
      if (currentParticle.GetDefinition()->GetParticleSubType() == "pi") {
        ++npions;
        ek1 += ekin;
      }
      currentParticle.SetKineticEnergy( ekin*GeV );
      pp = currentParticle.GetTotalMomentum()/MeV;
      pp1 = currentParticle.GetMomentum().mag()/MeV;
      if( pp1 < 0.001*MeV )
      {
        G4double costheta = 2.*G4UniformRand() - 1.;
        G4double sintheta = std::sqrt(1. - costheta*costheta);
        G4double phi = twopi*G4UniformRand();
        currentParticle.SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                                     pp*sintheta*std::sin(phi)*MeV,
                                     pp*costheta*MeV ) ;
      }
      else
        currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
      ekin = targetParticle.GetKineticEnergy()/GeV - cfa*(1+normal()/2.0);
      ekin = std::max( 1.0e-6, ekin );
      xxh = 1.0;
      if( (modifiedOriginal.GetDefinition() == aPiPlus ||
           modifiedOriginal.GetDefinition() == aPiMinus) &&
           targetParticle.GetDefinition() == aPiZero &&
           G4UniformRand() < logWeight) xxh = exh;
      dekin += ekin*(1.0-xxh);
      ekin *= xxh;
      if (targetParticle.GetDefinition()->GetParticleSubType() == "pi") {
        ++npions;
        ek1 += ekin;
      }
      targetParticle.SetKineticEnergy( ekin*GeV );
      pp = targetParticle.GetTotalMomentum()/MeV;
      pp1 = targetParticle.GetMomentum().mag()/MeV;
      if( pp1 < 0.001*MeV )
      {
        G4double costheta = 2.*G4UniformRand() - 1.;
        G4double sintheta = std::sqrt(1. - costheta*costheta);
        G4double phi = twopi*G4UniformRand();
        targetParticle.SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                                    pp*sintheta*std::sin(phi)*MeV,
                                    pp*costheta*MeV ) ;
      }
      else
        targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
      for( i=0; i<vecLen; ++i )
      {
        ekin = vec[i]->GetKineticEnergy()/GeV - cfa*(1+normal()/2.0);
        ekin = std::max( 1.0e-6, ekin );
        xxh = 1.0;
        if( (modifiedOriginal.GetDefinition() == aPiPlus ||
             modifiedOriginal.GetDefinition() == aPiMinus) &&
             vec[i]->GetDefinition() == aPiZero &&
             G4UniformRand() < logWeight) xxh = exh;
        dekin += ekin*(1.0-xxh);
        ekin *= xxh;
        if (vec[i]->GetDefinition()->GetParticleSubType() == "pi") {
          ++npions;
          ek1 += ekin;
        }
        vec[i]->SetKineticEnergy( ekin*GeV );
        pp = vec[i]->GetTotalMomentum()/MeV;
        pp1 = vec[i]->GetMomentum().mag()/MeV;
        if( pp1 < 0.001*MeV )
        {
          G4double costheta = 2.*G4UniformRand() - 1.;
          G4double sintheta = std::sqrt(1. - costheta*costheta);
          G4double phi = twopi*G4UniformRand();
          vec[i]->SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                               pp*sintheta*std::sin(phi)*MeV,
                               pp*costheta*MeV ) ;
        }
        else
          vec[i]->SetMomentum( vec[i]->GetMomentum() * (pp/pp1) );
      }
    }
    if( (ek1 != 0.0) && (npions > 0) )
    {
      dekin = 1.0 + dekin/ek1;
      //
      //  first do the incident particle
      //
      if (currentParticle.GetDefinition()->GetParticleSubType() == "pi") 
      {
        currentParticle.SetKineticEnergy(
         std::max( 0.001*MeV, dekin*currentParticle.GetKineticEnergy() ) );
        pp = currentParticle.GetTotalMomentum()/MeV;
        pp1 = currentParticle.GetMomentum().mag()/MeV;
        if( pp1 < 0.001 )
        {
          G4double costheta = 2.*G4UniformRand() - 1.;
          G4double sintheta = std::sqrt(1. - costheta*costheta);
          G4double phi = twopi*G4UniformRand();
          currentParticle.SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                                       pp*sintheta*std::sin(phi)*MeV,
                                       pp*costheta*MeV ) ;
        } else {
          currentParticle.SetMomentum( currentParticle.GetMomentum() * (pp/pp1) );
	}
      }

      if (targetParticle.GetDefinition()->GetParticleSubType() == "pi") 
      {
        targetParticle.SetKineticEnergy(
         std::max( 0.001*MeV, dekin*targetParticle.GetKineticEnergy() ) );
        pp = targetParticle.GetTotalMomentum()/MeV;
        pp1 = targetParticle.GetMomentum().mag()/MeV;
        if( pp1 < 0.001 )
        {
          G4double costheta = 2.*G4UniformRand() - 1.;
          G4double sintheta = std::sqrt(1. - costheta*costheta);
          G4double phi = twopi*G4UniformRand();
          targetParticle.SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                                      pp*sintheta*std::sin(phi)*MeV,
                                      pp*costheta*MeV ) ;
        } else {
          targetParticle.SetMomentum( targetParticle.GetMomentum() * (pp/pp1) );
	}
      }

      for( i=0; i<vecLen; ++i )
      {
        if (vec[i]->GetDefinition()->GetParticleSubType() == "pi")
        {
          vec[i]->SetKineticEnergy( std::max( 0.001*MeV, dekin*vec[i]->GetKineticEnergy() ) );
          pp = vec[i]->GetTotalMomentum()/MeV;
          pp1 = vec[i]->GetMomentum().mag()/MeV;
          if( pp1 < 0.001 )
          {
            G4double costheta = 2.*G4UniformRand() - 1.;
            G4double sintheta = std::sqrt(1. - costheta*costheta);
            G4double phi = twopi*G4UniformRand();
            vec[i]->SetMomentum( pp*sintheta*std::cos(phi)*MeV,
                                 pp*sintheta*std::sin(phi)*MeV,
                                 pp*costheta*MeV ) ;
          } else {
            vec[i]->SetMomentum( vec[i]->GetMomentum() * (pp/pp1) );
	  }
        }
      } // for i
    } // if (ek1 != 0)
  }
 
  std::pair<G4int, G4int> G4RPGReaction::GetFinalStateNucleons(
   const G4DynamicParticle* originalTarget,
   const G4FastVector<G4ReactionProduct,256>& vec, 
   const G4int& vecLen)
  {
    // Get number of protons and neutrons removed from the target nucleus
 
    G4int protonsRemoved = 0;
    G4int neutronsRemoved = 0;
    if (originalTarget->GetDefinition()->GetParticleName() == "proton")
      protonsRemoved++;
    else
      neutronsRemoved++;
 
    G4String secName;
    for (G4int i = 0; i < vecLen; i++) {
      secName = vec[i]->GetDefinition()->GetParticleName();
      if (secName == "proton") {
        protonsRemoved++;
      } else if (secName == "neutron") {
        neutronsRemoved++;
      } else if (secName == "anti_proton") {
        protonsRemoved--;
      } else if (secName == "anti_neutron") {
        neutronsRemoved--;
      }
    }

    return std::pair<G4int, G4int>(protonsRemoved, neutronsRemoved);
  }


 G4ThreeVector G4RPGReaction::Isotropic(const G4double& pp)
  {
    G4double costheta = 2.*G4UniformRand() - 1.;
    G4double sintheta = std::sqrt(1. - costheta*costheta);
    G4double phi = twopi*G4UniformRand();
    return G4ThreeVector(pp*sintheta*std::cos(phi),
                         pp*sintheta*std::sin(phi),
                         pp*costheta);
  }


 void G4RPGReaction::MomentumCheck(
   const G4ReactionProduct &modifiedOriginal,
   G4ReactionProduct &currentParticle,
   G4ReactionProduct &targetParticle,
   G4FastVector<G4ReactionProduct,256> &vec,
   G4int &vecLen )
  {
    const G4double pOriginal = modifiedOriginal.GetTotalMomentum()/MeV;
    G4double testMomentum = currentParticle.GetMomentum().mag()/MeV;
    G4double pMass;
    if( testMomentum >= pOriginal )
    {
      pMass = currentParticle.GetMass()/MeV;
      currentParticle.SetTotalEnergy(
       std::sqrt( pMass*pMass + pOriginal*pOriginal )*MeV );
      currentParticle.SetMomentum(
       currentParticle.GetMomentum() * (pOriginal/testMomentum) );
    }
    testMomentum = targetParticle.GetMomentum().mag()/MeV;
    if( testMomentum >= pOriginal )
    {
      pMass = targetParticle.GetMass()/MeV;
      targetParticle.SetTotalEnergy(
       std::sqrt( pMass*pMass + pOriginal*pOriginal )*MeV );
      targetParticle.SetMomentum(
       targetParticle.GetMomentum() * (pOriginal/testMomentum) );
    }
    for( G4int i=0; i<vecLen; ++i )
    {
      testMomentum = vec[i]->GetMomentum().mag()/MeV;
      if( testMomentum >= pOriginal )
      {
        pMass = vec[i]->GetMass()/MeV;
        vec[i]->SetTotalEnergy(
         std::sqrt( pMass*pMass + pOriginal*pOriginal )*MeV );
        vec[i]->SetMomentum( vec[i]->GetMomentum() * (pOriginal/testMomentum) );
      }
    }
  }

 void G4RPGReaction::NuclearReaction(
   G4FastVector<G4ReactionProduct,4> &vec,
   G4int &vecLen,
   const G4HadProjectile *originalIncident,
   const G4Nucleus &targetNucleus,
   const G4double theAtomicMass,
   const G4double *mass )
  {
    // derived from original FORTRAN code NUCREC by H. Fesefeldt (12-Feb-1987)
    //
    // Nuclear reaction kinematics at low energies
    //
    G4ParticleDefinition *aGamma = G4Gamma::Gamma();
    G4ParticleDefinition *aProton = G4Proton::Proton();
    G4ParticleDefinition *aNeutron = G4Neutron::Neutron();
    G4ParticleDefinition *aDeuteron = G4Deuteron::Deuteron();
    G4ParticleDefinition *aTriton = G4Triton::Triton();
    G4ParticleDefinition *anAlpha = G4Alpha::Alpha();
    
    const G4double aProtonMass = aProton->GetPDGMass()/MeV;
    const G4double aNeutronMass = aNeutron->GetPDGMass()/MeV;
    const G4double aDeuteronMass = aDeuteron->GetPDGMass()/MeV;
    const G4double aTritonMass = aTriton->GetPDGMass()/MeV;
    const G4double anAlphaMass = anAlpha->GetPDGMass()/MeV;

    G4ReactionProduct currentParticle;
    currentParticle = *originalIncident;
    //
    // Set beam particle, take kinetic energy of current particle as the
    // fundamental quantity.  Due to the difficult kinematic, all masses have to
    // be assigned the best measured values
    //
    G4double p = currentParticle.GetTotalMomentum();
    G4double pp = currentParticle.GetMomentum().mag();
    if( pp <= 0.001*MeV )
    {
      G4double phinve = twopi*G4UniformRand();
      G4double rthnve = std::acos( std::max( -1.0, std::min( 1.0, -1.0 + 2.0*G4UniformRand() ) ) );
      currentParticle.SetMomentum( p*std::sin(rthnve)*std::cos(phinve),
                                   p*std::sin(rthnve)*std::sin(phinve),
                                   p*std::cos(rthnve) );
    }
    else
      currentParticle.SetMomentum( currentParticle.GetMomentum() * (p/pp) );
    //
    // calculate Q-value of reactions
    //
    G4double currentKinetic = currentParticle.GetKineticEnergy()/MeV;
    G4double currentMass = currentParticle.GetDefinition()->GetPDGMass()/MeV;
    G4double qv = currentKinetic + theAtomicMass + currentMass;
    
    G4double qval[9];
    qval[0] = qv - mass[0];
    qval[1] = qv - mass[1] - aNeutronMass;
    qval[2] = qv - mass[2] - aProtonMass;
    qval[3] = qv - mass[3] - aDeuteronMass;
    qval[4] = qv - mass[4] - aTritonMass;
    qval[5] = qv - mass[5] - anAlphaMass;
    qval[6] = qv - mass[6] - aNeutronMass - aNeutronMass;
    qval[7] = qv - mass[7] - aNeutronMass - aProtonMass;
    qval[8] = qv - mass[8] - aProtonMass  - aProtonMass;
    
    if( currentParticle.GetDefinition() == aNeutron )
    {
      const G4double A = targetNucleus.GetA_asInt();    // atomic weight
      if( G4UniformRand() > ((A-1.0)/230.0)*((A-1.0)/230.0) )
        qval[0] = 0.0;
      if( G4UniformRand() >= currentKinetic/7.9254*A )
        qval[2] = qval[3] = qval[4] = qval[5] = qval[8] = 0.0;
    }
    else
      qval[0] = 0.0;
    
    G4int i;
    qv = 0.0;
    for( i=0; i<9; ++i )
    {
      if( mass[i] < 500.0*MeV )qval[i] = 0.0;
      if( qval[i] < 0.0 )qval[i] = 0.0;
      qv += qval[i];
    }
    G4double qv1 = 0.0;
    G4double ran = G4UniformRand();
    G4int index;
    for( index=0; index<9; ++index )
    {
      if( qval[index] > 0.0 )
      {
        qv1 += qval[index]/qv;
        if( ran <= qv1 )break;
      }
    }
    if( index == 9 )  // loop continued to the end
    {
      throw G4HadronicException(__FILE__, __LINE__,
           "G4RPGReaction::NuclearReaction: inelastic reaction kinematically not possible");
    }
    G4double ke = currentParticle.GetKineticEnergy()/GeV;
    G4int nt = 2;
    if( (index>=6) || (G4UniformRand()<std::min(0.5,ke*10.0)) )nt = 3;
    
    G4ReactionProduct **v = new G4ReactionProduct * [3];
    v[0] =  new G4ReactionProduct;
    v[1] =  new G4ReactionProduct;
    v[2] =  new G4ReactionProduct;
    
    v[0]->SetMass( mass[index]*MeV );
    switch( index )
    {
     case 0:
       v[1]->SetDefinition( aGamma );
       v[2]->SetDefinition( aGamma );
       break;
     case 1:
       v[1]->SetDefinition( aNeutron );
       v[2]->SetDefinition( aGamma );
       break;
     case 2:
       v[1]->SetDefinition( aProton );
       v[2]->SetDefinition( aGamma );
       break;
     case 3:
       v[1]->SetDefinition( aDeuteron );
       v[2]->SetDefinition( aGamma );
       break;
     case 4:
       v[1]->SetDefinition( aTriton );
       v[2]->SetDefinition( aGamma );
       break;
     case 5:
       v[1]->SetDefinition( anAlpha );
       v[2]->SetDefinition( aGamma );
       break;
     case 6:
       v[1]->SetDefinition( aNeutron );
       v[2]->SetDefinition( aNeutron );
       break;
     case 7:
       v[1]->SetDefinition( aNeutron );
       v[2]->SetDefinition( aProton );
       break;
     case 8:
       v[1]->SetDefinition( aProton );
       v[2]->SetDefinition( aProton );
       break;
    }
    //
    // calculate centre of mass energy
    //
    G4ReactionProduct pseudo1;
    pseudo1.SetMass( theAtomicMass*MeV );
    pseudo1.SetTotalEnergy( theAtomicMass*MeV );
    G4ReactionProduct pseudo2 = currentParticle + pseudo1;
    pseudo2.SetMomentum( pseudo2.GetMomentum() * (-1.0) );
    //
    // use phase space routine in centre of mass system
    //
    G4FastVector<G4ReactionProduct,256> tempV;
    tempV.Initialize( nt );
    G4int tempLen = 0;
    tempV.SetElement( tempLen++, v[0] );
    tempV.SetElement( tempLen++, v[1] );
    if( nt == 3 )tempV.SetElement( tempLen++, v[2] );
    G4bool constantCrossSection = true;
    GenerateNBodyEvent( pseudo2.GetMass()/MeV, constantCrossSection, tempV, tempLen );
    v[0]->Lorentz( *v[0], pseudo2 );
    v[1]->Lorentz( *v[1], pseudo2 );
    if( nt == 3 )v[2]->Lorentz( *v[2], pseudo2 );
    
    G4bool particleIsDefined = false;
    if( v[0]->GetMass()/MeV - aProtonMass < 0.1 )
    {
      v[0]->SetDefinition( aProton );
      particleIsDefined = true;
    }
    else if( v[0]->GetMass()/MeV - aNeutronMass < 0.1 )
    {
      v[0]->SetDefinition( aNeutron );
      particleIsDefined = true;
    }
    else if( v[0]->GetMass()/MeV - aDeuteronMass < 0.1 )
    {
      v[0]->SetDefinition( aDeuteron );
      particleIsDefined = true;
    }
    else if( v[0]->GetMass()/MeV - aTritonMass < 0.1 )
    {
      v[0]->SetDefinition( aTriton );
      particleIsDefined = true;
    }
    else if( v[0]->GetMass()/MeV - anAlphaMass < 0.1 )
    {
      v[0]->SetDefinition( anAlpha );
      particleIsDefined = true;
    }
    currentParticle.SetKineticEnergy(
     std::max( 0.001, currentParticle.GetKineticEnergy()/MeV ) );
    p = currentParticle.GetTotalMomentum();
    pp = currentParticle.GetMomentum().mag();
    if( pp <= 0.001*MeV )
    {
      G4double phinve = twopi*G4UniformRand();
      G4double rthnve = std::acos( std::max( -1.0, std::min( 1.0, -1.0 + 2.0*G4UniformRand() ) ) );
      currentParticle.SetMomentum( p*std::sin(rthnve)*std::cos(phinve),
                                   p*std::sin(rthnve)*std::sin(phinve),
                                   p*std::cos(rthnve) );
    }
    else
      currentParticle.SetMomentum( currentParticle.GetMomentum() * (p/pp) );
    
    if( particleIsDefined )
    {
      v[0]->SetKineticEnergy(
       std::max( 0.001, 0.5*G4UniformRand()*v[0]->GetKineticEnergy()/MeV ) );
      p = v[0]->GetTotalMomentum();
      pp = v[0]->GetMomentum().mag();
      if( pp <= 0.001*MeV )
      {
        G4double phinve = twopi*G4UniformRand();
        G4double rthnve = std::acos( std::max(-1.0,std::min(1.0,-1.0+2.0*G4UniformRand())) );
        v[0]->SetMomentum( p*std::sin(rthnve)*std::cos(phinve),
                          p*std::sin(rthnve)*std::sin(phinve),
                          p*std::cos(rthnve) );
      }
      else
        v[0]->SetMomentum( v[0]->GetMomentum() * (p/pp) );
    }
    if( (v[1]->GetDefinition() == aDeuteron) ||
        (v[1]->GetDefinition() == aTriton)   ||
        (v[1]->GetDefinition() == anAlpha) ) 
      v[1]->SetKineticEnergy(
       std::max( 0.001, 0.5*G4UniformRand()*v[1]->GetKineticEnergy()/MeV ) );
    else
      v[1]->SetKineticEnergy( std::max( 0.001, v[1]->GetKineticEnergy()/MeV ) );
    
    p = v[1]->GetTotalMomentum();
    pp = v[1]->GetMomentum().mag();
    if( pp <= 0.001*MeV )
    {
      G4double phinve = twopi*G4UniformRand();
      G4double rthnve = std::acos( std::max(-1.0,std::min(1.0,-1.0+2.0*G4UniformRand())) );
      v[1]->SetMomentum( p*std::sin(rthnve)*std::cos(phinve),
                        p*std::sin(rthnve)*std::sin(phinve),
                        p*std::cos(rthnve) );
    }
    else
      v[1]->SetMomentum( v[1]->GetMomentum() * (p/pp) );
    
    if( nt == 3 )
    {
      if( (v[2]->GetDefinition() == aDeuteron) ||
          (v[2]->GetDefinition() == aTriton)   ||
          (v[2]->GetDefinition() == anAlpha) ) 
        v[2]->SetKineticEnergy(
         std::max( 0.001, 0.5*G4UniformRand()*v[2]->GetKineticEnergy()/MeV ) );
      else
        v[2]->SetKineticEnergy( std::max( 0.001, v[2]->GetKineticEnergy()/MeV ) );
      
      p = v[2]->GetTotalMomentum();
      pp = v[2]->GetMomentum().mag();
      if( pp <= 0.001*MeV )
      {
        G4double phinve = twopi*G4UniformRand();
        G4double rthnve = std::acos( std::max(-1.0,std::min(1.0,-1.0+2.0*G4UniformRand())) );
        v[2]->SetMomentum( p*std::sin(rthnve)*std::cos(phinve),
                          p*std::sin(rthnve)*std::sin(phinve),
                          p*std::cos(rthnve) );
      }
      else
        v[2]->SetMomentum( v[2]->GetMomentum() * (p/pp) );
    }
    G4int del;
    for(del=0; del<vecLen; del++) delete vec[del];
    vecLen = 0;
    if( particleIsDefined )
    {
      vec.SetElement( vecLen++, v[0] );
    }
    else
    {
      delete v[0];
    }
    vec.SetElement( vecLen++, v[1] );
    if( nt == 3 )
    {
      vec.SetElement( vecLen++, v[2] );
    }
    else
    {
      delete v[2];
    }
    delete [] v;
    return;
  }
 
 /* end of file */
 
