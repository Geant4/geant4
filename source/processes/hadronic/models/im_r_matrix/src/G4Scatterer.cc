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
// $Id: G4Scatterer.cc,v 1.16 2010-03-12 15:45:18 gunter Exp $ //
//

#include <vector>

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Scatterer.hh"
#include "G4KineticTrack.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"

#include "G4CollisionNN.hh"
#include "G4CollisionPN.hh"
#include "G4CollisionMesonBaryon.hh"

#include "G4CollisionInitialState.hh"
#include "G4HadTmpUtil.hh"
#include "G4Pair.hh"
#include "G4AutoLock.hh"

//Mutex for control of shared resource
namespace  {
    G4Mutex collisions_mutex = G4MUTEX_INITIALIZER;
    G4bool setupDone = false;
}

// Declare the categories of collisions the Scatterer can handle
typedef GROUP2(G4CollisionNN, G4CollisionMesonBaryon) theChannels;

G4CollisionVector G4Scatterer::collisions;

//----------------------------------------------------------------------------

G4Scatterer::G4Scatterer()
{
  G4AutoLock l(&collisions_mutex);
  if ( ! setupDone )
  {
      Register aR;
      G4ForEach<theChannels>::Apply(&aR, &collisions);
      setupDone = true;
  }
}

//----------------------------------------------------------------------------

G4Scatterer::~G4Scatterer()
{
  G4AutoLock l(&collisions_mutex);
  std::for_each(collisions.begin(), collisions.end(), G4Delete());
  collisions.clear();
}

//----------------------------------------------------------------------------

G4double G4Scatterer::GetTimeToInteraction(const G4KineticTrack& trk1,
					   const G4KineticTrack& trk2) const
{
  G4double time = DBL_MAX;
    G4double distance_fast;
  G4LorentzVector mom1 = trk1.GetTrackingMomentum();
//  G4cout << "zcomp=" << std::abs(mom1.vect().unit().z() -1 ) << G4endl;
  G4double collisionTime;

  if ( std::abs(mom1.vect().unit().z() -1 ) < 1e-6 )
  {
     G4ThreeVector position = trk2.GetPosition() - trk1.GetPosition();
     G4double deltaz=position.z();
     G4double velocity = mom1.z()/mom1.e() * c_light;

     collisionTime=deltaz/velocity;
     distance_fast=position.x()*position.x() + position.y()*position.y();
  } else {

    //  The nucleons of the nucleus are FROZEN, ie. do not move..

    G4ThreeVector position = trk2.GetPosition() - trk1.GetPosition();

    G4ThreeVector velocity = mom1.vect()/mom1.e() * c_light;  // mom1.boostVector() will exit on slightly negative mass
    collisionTime = (position * velocity) / velocity.mag2();    // can't divide by /c_light;
    position -= velocity * collisionTime;
    distance_fast=position.mag2();

//    if ( collisionTime>0 ) G4cout << " dis1/2 square" << dis1 <<" "<< dis2 << G4endl;
//     collisionTime = GetTimeToClosestApproach(trk1,trk2);
  }
     if (collisionTime > 0)
	{
	   static const G4double maxCrossSection = 500*millibarn;
	   if(0.7*pi*distance_fast>maxCrossSection) return time;


           G4LorentzVector mom2(0,0,0,trk2.Get4Momentum().mag());

// 	   G4ThreeVector momLab = mom1.vect();// frozen Nucleus - mom2.vect();
// 	   G4ThreeVector posLab = trk1.GetPosition() - trk2.GetPosition();
// 	   G4double disLab=posLab * posLab - (posLab*momLab) * (posLab*momLab) /(momLab.mag2());

	   G4LorentzRotation toCMSFrame((-1)*(mom1 + mom2).boostVector());
	   mom1 = toCMSFrame * mom1;
	   mom2 = toCMSFrame * mom2;

	   G4LorentzVector coordinate1(trk1.GetPosition(), 100.);
	   G4LorentzVector coordinate2(trk2.GetPosition(), 100.);
	   G4ThreeVector pos = ((toCMSFrame * coordinate1).vect() -
				(toCMSFrame * coordinate2).vect());

	   G4ThreeVector mom = mom1.vect() - mom2.vect();

	  // Calculate the impact parameter

	   G4double distance = pos * pos - (pos*mom) * (pos*mom) / (mom.mag2());

//     G4cout << " disDiff " << distance-disLab << " " << disLab
//            << " " << std::abs(distance-disLab)/distance << G4endl
//	    << " mom/Lab " << mom << " " << momLab << G4endl
//	    << " pos/Lab " << pos << " " << posLab
//	    << G4endl;

	   if(pi*distance>maxCrossSection) return time;

	   // charged particles special
	   static const G4double maxChargedCrossSection = 200*millibarn;
	   if(std::abs(trk1.GetDefinition()->GetPDGCharge())>0.1 &&
	      std::abs(trk2.GetDefinition()->GetPDGCharge())>0.1 &&
	      pi*distance>maxChargedCrossSection) return time;

           G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
	   // neutrons special  pn is largest cross-section, but above 1.91 GeV is less than 200 mb
	   if(( trk1.GetDefinition() == G4Neutron::Neutron() ||
	        trk2.GetDefinition() == G4Neutron::Neutron() ) &&
		     sqrtS>1.91*GeV && pi*distance>maxChargedCrossSection) return time;

/*
 * 	  if(distance <= sqr(1.14*fermi))
 * 	  {
 * 	    time = collisionTime;
 *
 * *
 *  * 	     G4cout << "Scatter distance/time: " << std::sqrt(distance)/fermi <<
 *  * 	         " / "<< time/ns << G4endl;
 *  * 	      G4ThreeVector pos1=trk1.GetPosition();
 *  * 	      G4ThreeVector pos2=trk2.GetPosition();
 *  * 	      G4LorentzVector xmom1 = trk1.Get4Momentum();
 *  * 	      G4LorentzVector xmom2 = trk2.Get4Momentum();
 *  * 	      G4cout << "position1: " <<  pos1.x() << " " << pos1.y() << " "
 *  * 	      		<< pos1.z();
 *  * 	      pos1+=(collisionTime*c_light/xmom1.e())*xmom1.vect();
 *  * 	      G4cout << " straight line trprt: "
 *  * 	      		<<  pos1.x() << " " << pos1.y() << " "
 *  * 			<<  pos1.z()  << G4endl;
 *  * 	      G4cout << "position2: " <<  pos2.x() << " " << pos2.y() << " "
 *  * 	      		<< pos2.z()  << G4endl;
 *  * 	      G4cout << "straight line distance 2 fixed:" << (pos1-pos2).mag()/fermi << G4endl;
 *  * 	      pos2+= (collisionTime*c_light/xmom2.e())*xmom2.vect();
 *  * 	      G4cout<< " straight line trprt: "
 *  * 	      		<<  pos2.x() << " " << pos2.y() << " "
 *  * 			<<  pos2.z() << G4endl;
 *  * 	      G4cout << "straight line distance :" << (pos1-pos2).mag()/fermi << G4endl;
 *  *
 * 	  }
 *
 * 	  if(1)
 * 	    return time;
 */

	   if ((trk1.GetActualMass()+trk2.GetActualMass()) > sqrtS) return time;



	  const G4VCollision* collision = FindCollision(trk1,trk2);
	  G4double totalCrossSection;
	  // The cross section is interpreted geometrically as an area
	  // Two particles are assumed to collide if their distance is < (totalCrossSection/pi)

	  if (collision != 0)
	    {
	      totalCrossSection = collision->CrossSection(trk1,trk2);
	      if ( totalCrossSection > 0 )
	        {
/*		    G4cout << " totalCrossection = "<< totalCrossSection << ", trk1/2, s, e-m: "
 *		           << trk1.GetDefinition()->GetParticleName()
 *			   << " / "
 *		           << trk2.GetDefinition()->GetParticleName()
 *			   << ", "
 *			   << (trk1.Get4Momentum()+trk2.Get4Momentum()).mag()
 *			   << ", "
 *			   << (trk1.Get4Momentum()+trk2.Get4Momentum()).mag()-
 *			       trk1.Get4Momentum().mag() - trk2.Get4Momentum().mag()
 *			   << G4endl;
 */
		 if (distance <= totalCrossSection / pi)
		   {
		     time = collisionTime;
		   }
	        } else
		{

 		 // For debugging...
 //		    G4cout << " totalCrossection = 0, trk1/2, s, e-m: "
 //		           << trk1.GetDefinition()->GetParticleName()
 //			   << " / "
 //		           << trk2.GetDefinition()->GetParticleName()
 //			   << ", "
 //			   << (trk1.Get4Momentum()+trk2.Get4Momentum()).mag()
 //			   << ", "
 //			   << (trk1.Get4Momentum()+trk2.Get4Momentum()).mag()-
 //			       trk1.Get4Momentum().mag() - trk2.Get4Momentum().mag()
 //			   << G4endl;

		}
/*
 * 	      if(distance <= sqr(5.*fermi))
 * 	       {
 * 		  G4cout << " distance,xsect, std::sqrt(xsect/pi) : " << std::sqrt(distance)/fermi
 * 			 << " " << totalCrossSection/sqr(fermi)
 * 			 << " " << std::sqrt(totalCrossSection / pi)/fermi << G4endl;
 * 	       }
 */

	    }
	  else
	    {
	      time = DBL_MAX;
//	      /*
	      // For debugging
//hpw	      	      G4cout << "G4Scatterer - collision not found: "
//hpw		     << trk1.GetDefinition()->GetParticleName()
//hpw		     << " - "
//hpw		     << trk2.GetDefinition()->GetParticleName()
//hpw		     << G4endl;
	      // End of debugging
//	      */
	    }
	}

      else
	{
	      /*
	  // For debugging
	  G4cout << "G4Scatterer - negative collisionTime"
		 << ": collisionTime = " << collisionTime
		 << ", position = " << position
		 << ", velocity = " << velocity
		 << G4endl;
	  // End of debugging
	      */
	}

  return time;
}

//----------------------------------------------------------------------------

G4KineticTrackVector* G4Scatterer::Scatter(const G4KineticTrack& trk1,
					      const G4KineticTrack& trk2) const
{
//   G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();
   G4LorentzVector pInitial=trk1.Get4Momentum() + trk2.Get4Momentum();
   G4double energyBalance = pInitial.t();
   G4double pxBalance = pInitial.vect().x();
   G4double pyBalance = pInitial.vect().y();
   G4double pzBalance = pInitial.vect().z();
   G4int chargeBalance = G4lrint(trk1.GetDefinition()->GetPDGCharge()
                       + trk2.GetDefinition()->GetPDGCharge());
   G4int baryonBalance = trk1.GetDefinition()->GetBaryonNumber()
                       + trk2.GetDefinition()->GetBaryonNumber();

   const G4VCollision* collision = FindCollision(trk1,trk2);
   if (collision != 0)
   {
     G4double aCrossSection = collision->CrossSection(trk1,trk2);
     if (aCrossSection > 0.0)
     {


  	#ifdef debug_G4Scatterer
	G4cout << "be4 FinalState 1(p,e,m): "
        << trk1.Get4Momentum() << " "
        << trk1.Get4Momentum().mag()
	<< ", 2: "
        << trk2.Get4Momentum()<< " "
        << trk2.Get4Momentum().mag() << " "
        << G4endl;
	#endif


       G4KineticTrackVector* products = collision->FinalState(trk1,trk2);
       if(!products || products->size() == 0) return products;

  	#ifdef debug_G4Scatterer
       G4cout << "size of FS: "<<products->size()<<G4endl;
	#endif

       G4KineticTrack *final= products->operator[](0);


  	#ifdef debug_G4Scatterer
        G4cout << "    FinalState 1: "
		<< final->Get4Momentum()<< " "
		<< final->Get4Momentum().mag() ;
	#endif

        if(products->size() == 1) return products;
	final=products->operator[](1);
  	#ifdef debug_G4Scatterer
	G4cout << ", 2: "
		<< final->Get4Momentum() << " "
        	<< final->Get4Momentum().mag() << " " << G4endl;
	#endif

       final= products->operator[](0);
       G4LorentzVector pFinal=final->Get4Momentum();
       if(products->size()==2)
       {
         final=products->operator[](1);
         pFinal +=final->Get4Momentum();
       }

       #ifdef debug_G4Scatterer
       if ( (pInitial-pFinal).mag() > 0.1*MeV )
       {
          G4cout << "G4Scatterer: momentum imbalance, pInitial= " <<pInitial << " pFinal= " <<pFinal<< G4endl;
       }
       G4cout << "Scatterer costh= " << trk1.Get4Momentum().vect().unit() *(products->operator[](0))->Get4Momentum().vect().unit()<< G4endl;
       #endif

       for(size_t hpw=0; hpw<products->size(); hpw++)
       {
         energyBalance-=products->operator[](hpw)->Get4Momentum().t();
         pxBalance-=products->operator[](hpw)->Get4Momentum().vect().x();
         pyBalance-=products->operator[](hpw)->Get4Momentum().vect().y();
         pzBalance-=products->operator[](hpw)->Get4Momentum().vect().z();
	 chargeBalance-=G4lrint(products->operator[](hpw)->GetDefinition()->GetPDGCharge());
         baryonBalance-=products->operator[](hpw)->GetDefinition()->GetBaryonNumber();
       }
       if(getenv("ScattererEnergyBalanceCheck"))
         std::cout << "DEBUGGING energy balance A: "
	           <<energyBalance<<" "
	           <<pxBalance<<" "
	           <<pyBalance<<" "
	           <<pzBalance<<" "
		   <<chargeBalance<<" "
		   <<baryonBalance<<" "
		   <<G4endl;
       if(chargeBalance !=0 )
       {
         G4cout << "track 1"<<trk1.GetDefinition()->GetParticleName()<<G4endl;
         G4cout << "track 2"<<trk2.GetDefinition()->GetParticleName()<<G4endl;
         for(size_t hpw=0; hpw<products->size(); hpw++)
         {
            G4cout << products->operator[](hpw)->GetDefinition()->GetParticleName()<<G4endl;
         }
         G4Exception("G4Scatterer", "im_r_matrix001", FatalException,
             "Problem in ChargeBalance");
       }
       return products;
     }
   }

   return NULL;
}

//----------------------------------------------------------------------------

const G4VCollision* G4Scatterer::FindCollision(const G4KineticTrack& trk1,
					 const G4KineticTrack& trk2) const
{
  G4VCollision* collisionInCharge = 0;

  size_t i;
  for (i=0; i<collisions.size(); i++)
    {
      G4VCollision* component = collisions[i];
      if (component->IsInCharge(trk1,trk2))
	{
	  collisionInCharge = component;
	  break;
	}
    }
//    if(collisionInCharge)
//    {
//      G4cout << "found collision : "
//         << collisionInCharge->GetName()<< " "
// 	<< "for "
// 	<< trk1.GetDefinition()->GetParticleName()<<" + "
// 	<< trk2.GetDefinition()->GetParticleName()<<" "
// 	<< G4endl;;
//    }
  return collisionInCharge;
}

//----------------------------------------------------------------------------

G4double G4Scatterer::GetCrossSection(const G4KineticTrack& trk1,
				      const G4KineticTrack& trk2) const
{
   const G4VCollision* collision = FindCollision(trk1,trk2);
   G4double aCrossSection = 0;
   if (collision != 0)
   {
     aCrossSection = collision->CrossSection(trk1,trk2);
   }
   return aCrossSection;
}

//----------------------------------------------------------------------------

const std::vector<G4CollisionInitialState *> & G4Scatterer::
GetCollisions(G4KineticTrack * aProjectile,
              std::vector<G4KineticTrack *> & someCandidates,
	      G4double aCurrentTime)
{
  theCollisions.clear();
  std::vector<G4KineticTrack *>::iterator j=someCandidates.begin();
  for(; j != someCandidates.end(); ++j)
  {
    G4double collisionTime = GetTimeToInteraction(*aProjectile, **j);
    if(collisionTime == DBL_MAX)  // no collision
    {
      continue;
    }
    G4KineticTrackVector aTarget;
    aTarget.push_back(*j);
    theCollisions.push_back(
      new G4CollisionInitialState(collisionTime+aCurrentTime, aProjectile, aTarget, this) );
//      G4cerr <<" !!!!!! debug collisions "<<collisionTime<<" "<<pkt->GetDefinition()->GetParticleName()<<G4endl;
   }
   return theCollisions;
}


G4KineticTrackVector * G4Scatterer::
GetFinalState(G4KineticTrack * aProjectile,
	      std::vector<G4KineticTrack *> & theTargets)
{
    G4KineticTrack target_reloc(*(theTargets[0]));
    return Scatter(*aProjectile, target_reloc);
}
//----------------------------------------------------------------------------
