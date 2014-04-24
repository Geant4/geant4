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
#ifndef G4GeneralNNCollision_h
#define G4GeneralNNCollision_h

#include "G4CollisionComposite.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

class G4GeneralNNCollision : public G4CollisionComposite
{
  public:

  G4bool 
  IsInCharge(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const
  {
    G4bool result = false;
    const G4ParticleDefinition * aD = trk1.GetDefinition();
    const G4ParticleDefinition * bD = trk2.GetDefinition();
    if(  (aD==G4Proton::Proton() || aD == G4Neutron::Neutron())
       &&(bD==G4Proton::Proton() || bD == G4Neutron::Neutron()) ) result = true;
    return result;
  }
  
  protected:

  template <int dm, int d0, int dp, int dpp, class channelType> 
  struct MakeNNToNDelta {
  static void Make(G4CollisionComposite * aC)
  {
    typedef INT4<channelType, NeutronPC, NeutronPC, NeutronPC, d0>  theC1;
    typedef INT4<channelType, NeutronPC, NeutronPC, ProtonPC,  dm>  theC2;
    typedef INT4<channelType, NeutronPC, ProtonPC,  ProtonPC,  d0>  theC3;
    typedef INT4<channelType, NeutronPC, ProtonPC,  NeutronPC, dp>  theC4;
    typedef INT4<channelType, ProtonPC,  ProtonPC,  NeutronPC, dpp> theC5;
    typedef INT4<channelType, ProtonPC,  ProtonPC,  ProtonPC,  dp>  theC6;
    typedef GROUP6(theC1, theC2, theC3, theC4, theC5, theC6) theChannels;
    G4CollisionComposite::Resolve aR;
    G4ForEach<theChannels>::Apply(&aR, aC); 
  }};

  template <int Np, int Nn, class channelType> 
  struct MakeNNToNNStar{
  static void Make(G4CollisionComposite * aC)
  {
    typedef INT4<channelType, NeutronPC, NeutronPC, NeutronPC, Nn>  theC1;
    typedef INT4<channelType, ProtonPC,  ProtonPC,  ProtonPC,  Np>  theC2;
    typedef INT4<channelType, NeutronPC, ProtonPC,  NeutronPC, Np>  theC3;
    typedef INT4<channelType, NeutronPC, ProtonPC,  ProtonPC,  Nn>  theC4;
    typedef GROUP4(theC1, theC2, theC3, theC4) theChannels;
    G4CollisionComposite::Resolve aR;
    G4ForEach<theChannels>::Apply(&aR, aC); 
  }};
  
  template <class channelType, int Np, int Nn> 
  struct MakeNNStarToNN{
  static void Make(G4CollisionComposite * aC)
  {
    typedef INT4<channelType, Nn, NeutronPC, NeutronPC, NeutronPC>  theC1;
    typedef INT4<channelType, Np, ProtonPC,  ProtonPC,  ProtonPC>   theC2;
    typedef INT4<channelType, Np, NeutronPC, NeutronPC, ProtonPC>   theC3;
    typedef INT4<channelType, Nn, ProtonPC, NeutronPC, ProtonPC>    theC4;
    typedef GROUP4(theC1, theC2, theC3, theC4) theChannels;
    G4CollisionComposite::Resolve aR;
    G4ForEach<theChannels>::Apply(&aR, aC); 
  }};
  
  template <int Np, class channelType, int Nn> 
  struct MakeNNToDeltaNstar{
  static void Make(G4CollisionComposite * aC)
  {
    typedef INT4<channelType, NeutronPC, NeutronPC, D1232::D0,  Nn>  theC1;
    typedef INT4<channelType, NeutronPC, NeutronPC, D1232::Dm,  Np>  theC2;
    typedef INT4<channelType, ProtonPC,  ProtonPC,  D1232::Dp,  Np>  theC3;
    typedef INT4<channelType, ProtonPC,  ProtonPC,  D1232::Dpp, Nn>  theC4;
    typedef INT4<channelType, NeutronPC, ProtonPC,  D1232::D0,  Np>  theC5;
    typedef INT4<channelType, NeutronPC, ProtonPC,  D1232::Dp,  Nn>  theC6;
    typedef GROUP6(theC1, theC2, theC3, theC4, theC5, theC6) theChannels;
    G4CollisionComposite::Resolve aR;
    G4ForEach<theChannels>::Apply(&aR, aC); 
  }};
  
  template <int dm, int d0, int dp, int dpp, class channelType> 
  struct MakeNNToDeltaDelta{
  static void Make(G4CollisionComposite * aC)
  {
    typedef INT4<channelType, NeutronPC, NeutronPC, DeltamPC, dp>  theC1;
    typedef INT4<channelType, NeutronPC, NeutronPC, Delta0PC, d0>  theC2;
    typedef INT4<channelType, NeutronPC, NeutronPC, DeltapPC, dm>  theC3;
    typedef INT4<channelType, NeutronPC, ProtonPC, DeltapPC, d0>  theC4;
    typedef INT4<channelType, NeutronPC, ProtonPC, Delta0PC, dp>  theC5;
    typedef INT4<channelType, NeutronPC, ProtonPC, DeltamPC, dpp>  theC6;
    typedef INT4<channelType, NeutronPC, ProtonPC, DeltappPC, dm>  theC7;
    typedef INT4<channelType, ProtonPC, ProtonPC, Delta0PC, dpp>  theC8;
    typedef INT4<channelType, ProtonPC, ProtonPC, DeltapPC, dp>  theC9;
    typedef INT4<channelType, ProtonPC, ProtonPC, DeltappPC, d0>  theC10;
    typedef GROUP10(theC1, theC2, theC3, theC4, theC5, theC6, theC7, theC8, theC9, theC10) theChannels;
    G4CollisionComposite::Resolve aR;
    G4ForEach<theChannels>::Apply(&aR, aC); 
  }};
};

#endif
