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
// $Id: G4Scatterer.cc,v 1.13.2.1 2004/03/24 13:18:45 hpw Exp $ //
//

#include "globals.hh"
#include <vector>
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


// Declare the categories of collisions the Scatterer can handle
typedef GROUP2(G4CollisionNN, G4CollisionMesonBaryon) theChannels;


G4Scatterer::G4Scatterer()
{
  Register aR;
  G4ForEach<theChannels>::Apply(&aR, &collisions);
}


G4Scatterer::~G4Scatterer()
{
  std::for_each(collisions.begin(), collisions.end(), G4Delete());
  collisions.clear();
}


G4double G4Scatterer::GetTimeToInteraction(const G4KineticTrack& trk1, 
					   const G4KineticTrack& trk2)
{
  G4double time = DBL_MAX;
  
  G4double sqrtS = (trk1.Get4Momentum() + trk2.Get4Momentum()).mag();

  // Check whether there is enough energy for elastic scattering 
  // (to put the particles on to mass shell
  
//  if (trk1.GetDefinition()->GetPDGMass() + trk2.GetDefinition()->GetPDGMass() < sqrtS)
    if (trk1.GetActualMass() + trk2.GetActualMass() < sqrtS)
    {
      G4LorentzVector mom1 = trk1.GetTrackingMomentum();
      //  The nucleons of the nucleus are FROZEN, ie. do not move..
      
      G4ThreeVector position = trk1.GetPosition() - trk2.GetPosition();    

      // G4ThreeVector velocity = (mom1.boostVector() - mom2.boostVector()) * c_light;
      if ( mom1.mag2() < -1.*eV )
      {
        G4cout << "G4Scatterer::GetTimeToInteraction(): negative m2:" << mom1.mag2() << G4endl;
      } 
      G4ThreeVector velocity = mom1.vect()/mom1.e() * c_light; 
      G4double collisionTime = - (position * velocity) / (velocity * velocity);    // can't divide by /c_light;
      
     if (collisionTime > 0)
	{ 
           G4LorentzVector mom2(0,0,0,trk2.Get4Momentum().mag());
	   G4LorentzRotation toCMSFrame((-1)*(mom1 + mom2).boostVector());
	   mom1 = toCMSFrame * mom1;
	   mom2 = toCMSFrame * mom2;

	   G4LorentzVector coordinate1(trk1.GetPosition(), 100.);
	   G4LorentzVector coordinate2(trk2.GetPosition(), 100.);
	   G4ThreeVector pos = ((toCMSFrame * coordinate1).vect() - 
				(toCMSFrame * coordinate2).vect());

	   G4ThreeVector mom = mom1.vect() - mom2.vect();

	  // Calculate the impact parameter

	   G4double distance = pos * pos - (pos*mom) * (pos*mom) / (mom*mom);
	   // global optimization
	   static const G4double maxCrossSection = 500*millibarn;
	   if(pi*distance>maxCrossSection) return time;
	   
	   // charged particles special
	   static const G4double maxChargedCrossSection = 200*millibarn;
	   if(abs(trk1.GetDefinition()->GetPDGCharge())>0.1 && 
	      abs(trk2.GetDefinition()->GetPDGCharge())>0.1 &&
	      pi*distance>maxChargedCrossSection) return time;
	      
	   // neutrons special   
	   if(( trk1.GetDefinition() == G4Neutron::Neutron() ||
	        trk1.GetDefinition() == G4Neutron::Neutron() ) &&
		sqrtS>1.91*GeV && pi*distance>maxChargedCrossSection) return time;

/*
 * 	  if(distance <= sqr(1.14*fermi))
 * 	  {
 * 	    time = collisionTime;
 * 	  
 * *
 *  * 	     G4cout << "Scatter distance/time: " << sqrt(distance)/fermi <<
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
	    
	  
	  G4VCollision* collision = FindCollision(trk1,trk2);
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
 * 		  G4cout << " distance,xsect, sqrt(xsect/pi) : " << sqrt(distance)/fermi
 * 			 << " " << totalCrossSection/sqr(fermi)
 * 			 << " " << sqrt(totalCrossSection / pi)/fermi << G4endl;
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
    }
  else
    {
	      /*
      // For debugging
      G4cout << "G4Scatterer - Infinite time to interaction"
	     << ": sqrtS " << sqrtS
	     << ", mass1 = " << trk1.GetDefinition()->GetPDGMass()
	     << ", mass2 = " << trk2.GetDefinition()->GetPDGMass()
	     << G4endl;
	throw G4HadronicException(__FILE__, __LINE__, "G4Scatterer TimeToInteraction is INF");
      // End of debugging
	      */
    }

  return time;
}

G4KineticTrackVector* G4Scatterer::Scatter(const G4KineticTrack& trk1, 
					      const G4KineticTrack& trk2)  
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
   
   G4VCollision* collision = FindCollision(trk1,trk2);
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
	 G4Exception("We have the problem");
       }
       return products;
     } 
   } 
   
   return NULL;
}


G4VCollision* G4Scatterer::FindCollision(const G4KineticTrack& trk1, 
					 const G4KineticTrack& trk2)
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
  
G4double G4Scatterer::GetCrossSection(const G4KineticTrack& trk1, 
				      const G4KineticTrack& trk2)  
{
   G4VCollision* collision = FindCollision(trk1,trk2);
   G4double aCrossSection = 0;
   if (collision != 0)
   {
     aCrossSection = collision->CrossSection(trk1,trk2);
   }
   return aCrossSection;
}


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
