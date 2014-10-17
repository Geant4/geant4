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
// $Id: G4CollisionComposite.cc,v 1.9 2010-03-12 15:45:18 gunter Exp $ //

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4CollisionComposite.hh"
#include "G4VCollision.hh"
#include "G4CollisionVector.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4VCrossSectionSource.hh"
#include "G4HadTmpUtil.hh"
#include "G4AutoLock.hh"

const G4int G4CollisionComposite::nPoints = 32;

const G4double G4CollisionComposite::theT[nPoints] =
{.01, .03, .05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10., 15, 20, 50, 100};

G4CollisionComposite::G4CollisionComposite()
{ 
  G4MUTEXINIT( bufferMutex );
}


G4CollisionComposite::~G4CollisionComposite()
{
  G4MUTEXDESTROY(bufferMutex);
  std::for_each(components.begin(), components.end(), G4Delete());
}


G4double G4CollisionComposite::CrossSection(const G4KineticTrack& trk1, 
					    const G4KineticTrack& trk2) const
{
  G4double crossSect = 0.;
  const G4VCrossSectionSource* xSource = GetCrossSectionSource();
  if (xSource != 0)
  // There is a total cross section for this Collision
  {
    crossSect = xSource->CrossSection(trk1,trk2);
  }
  else
  {
    G4AutoLock l(&bufferMutex);
    // waiting for mutable to enable buffering.
    const_cast<G4CollisionComposite *>(this)->BufferCrossSection(trk1.GetDefinition(), trk2.GetDefinition());
//    G4cerr << "Buffer filled, reying with sqrts = "<< (trk1.Get4Momentum()+trk2.Get4Momentum()).mag() <<G4endl;
    crossSect = BufferedCrossSection(trk1,trk2);
  }
  return crossSect;
}


G4KineticTrackVector* G4CollisionComposite::FinalState(const G4KineticTrack& trk1, 
							  const G4KineticTrack& trk2) const
{
  std::vector<G4double> cxCache;
  G4double partialCxSum = 0.0;

  size_t i;
  for (i=0; i<components.size(); i++) 
  {
    G4double partialCx;
//    cout << "comp" << i << " " << components[i]()->GetName();
    if (components[i]->IsInCharge(trk1,trk2)) 
    {
      partialCx = components[i]->CrossSection(trk1,trk2);
    } 
    else 
    {
      partialCx = 0.0;
    }
//    cout << "   cx=" << partialCx << endl;
    partialCxSum += partialCx;
    cxCache.push_back(partialCx);
  }

  G4double random = G4UniformRand()*partialCxSum;
  G4double running = 0;
  for (i=0; i<cxCache.size(); i++) 
  {
    running += cxCache[i];
    if (running > random) 
    {
      return components[i]->FinalState(trk1, trk2);
    }
  }
//  G4cerr <<"in charge = "<<IsInCharge(trk1, trk2)<<G4endl;
//  G4cerr <<"Cross-section = "<<CrossSection(trk1, trk2)/millibarn<<" "<<running<<" "<<cxCache.size()<<G4endl;
//  G4cerr <<"Names = "<<trk1.GetDefinition()->GetParticleName()<<", "<<trk2.GetDefinition()->GetParticleName()<<G4endl;
//  throw G4HadronicException(__FILE__, __LINE__, "G4CollisionComposite: no final state found!");
  return NULL;
}


G4bool G4CollisionComposite::IsInCharge(const G4KineticTrack& trk1, 
					const G4KineticTrack& trk2) const
{
  G4bool isInCharge = false;

  // The composite is in charge if any of its components is in charge

  const G4CollisionVector* comps = GetComponents();
  if (comps)
    {
      G4CollisionVector::const_iterator iter;
      for (iter = comps->begin(); iter != comps->end(); ++iter)
	{
	 if ( ((*iter))->IsInCharge(trk1,trk2) ) isInCharge = true;
	}
    }

  return isInCharge;
}

void G4CollisionComposite::
BufferCrossSection(const G4ParticleDefinition * aP, const G4ParticleDefinition * bP)
{
   // check if already buffered
   size_t i;
   for(i=0; i<theBuffer.size(); i++)
   {
     if(theBuffer[i].InCharge(aP, bP)) return;
   }
//   G4cerr << "Buffering for "<<aP->GetParticleName()<<" "<<bP->GetParticleName()<<G4endl;
   
   // buffer the new one.
   G4CrossSectionBuffer aNewBuff(aP, bP);
   size_t maxE=nPoints;
   for(size_t tt=0; tt<maxE; tt++)
   {
     G4double aT = theT[tt]*GeV;
     G4double crossSect = 0;
     // The total cross-section is summed over all the component channels
     
     //A.R. 28-Sep-2012 Fix reproducibility problem
     //                 Assign the kinetic energy to the lightest of the
     //                 two particles, instead to the first one always.
     G4double atime = 0;
     G4double btime = 0;
     G4ThreeVector aPosition(0,0,0);
     G4ThreeVector bPosition(0,0,0);
     G4double aM = aP->GetPDGMass();
     G4double bM = bP->GetPDGMass();
     G4double aE = aM;
     G4double bE = bM;
     G4ThreeVector aMom(0,0,0);
     G4ThreeVector bMom(0,0,0);
     if ( aM <= bM ) {
      aE += aT;
      aMom = G4ThreeVector(0,0,std::sqrt(aE*aE-aM*aM));
     } else {
      bE += aT;
      bMom = G4ThreeVector(0,0,std::sqrt(bE*bE-bM*bM));
     }
     G4LorentzVector a4Momentum(aE, aMom);
     G4LorentzVector b4Momentum(bE, bMom);
     G4KineticTrack a(aP, atime, aPosition, a4Momentum);
     G4KineticTrack b(bP, btime, bPosition, b4Momentum);
     
     for (i=0; i<components.size(); i++)
     {
       if(components[i]->IsInCharge(a,b))
       {
	 crossSect += components[i]->CrossSection(a,b);
       }
     }
     G4double sqrts = (a4Momentum+b4Momentum).mag();
     aNewBuff.push_back(sqrts, crossSect);
   }
   theBuffer.push_back(aNewBuff);
//   theBuffer.back().Print();
}


G4double G4CollisionComposite::
BufferedCrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
{
   for(size_t i=0; i<theBuffer.size(); i++)
   {
     if(theBuffer[i].InCharge(trk1.GetDefinition(), trk2.GetDefinition())) 
     {
       return theBuffer[i].CrossSection(trk1, trk2);
     }
   }
   throw G4HadronicException(__FILE__, __LINE__, "G4CollisionComposite::BufferedCrossSection - Blitz !!");
   return 0;
}

