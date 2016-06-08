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
// $Id: G4CollisionMesonBaryonToResonance.hh,v 1.7 2002/12/12 19:17:37 gunter Exp $ //

#ifndef G4CollisionMesonBaryonToResonance_h
#define G4CollisionMesonBaryonToResonance_h

#include "globals.hh"
#include "G4CollisionComposite.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include "g4std/vector"
#include "G4XpipNTotal.hh"
#include "G4XpimNTotal.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4PionMinus.hh"


class G4CollisionMesonBaryonToResonance : public G4CollisionComposite
{

public:

  G4CollisionMesonBaryonToResonance();
  virtual ~G4CollisionMesonBaryonToResonance(){};
  virtual G4String GetName() const { return "mN -> baryon Resonance Collision"; }  
  virtual G4double CrossSection(const G4KineticTrack& trk1, 
				const G4KineticTrack& trk2) const
  {
    G4double result=0;
//     if(trk1.GetDefinition() == G4PionPlus::PionPlus() && trk2.GetDefinition()==G4Proton::Proton())
//     {
//       result=thepipp.CrossSection(trk1, trk2);
//     }
//     else if(trk2.GetDefinition() == G4PionPlus::PionPlus() && trk1.GetDefinition()==G4Proton::Proton())
//     {
//       result=thepipp.CrossSection(trk2, trk1);
//     }
//     else if(trk1.GetDefinition() == G4PionMinus::PionMinus() && trk2.GetDefinition()==G4Proton::Proton())
//     {
//       result=thepimp.CrossSection(trk1, trk2);
//     }
//     else if(trk2.GetDefinition() == G4PionMinus::PionMinus() && trk1.GetDefinition()==G4Proton::Proton())
//     {
//       result=thepimp.CrossSection(trk2, trk1);
//     }
//     else 
    {
      result = G4CollisionComposite::CrossSection(trk1, trk2);
    }
    return result;
  }

protected:
  
  virtual const G4std::vector<G4String>& GetListOfColliders(G4int whichOne) const
  {
    G4Exception("Tried to call G4CollisionMesonBaryonToResonance::GetListOfColliders. Please find out why!");
    G4std::vector<G4String> * aList = new G4std::vector<G4String>;
    return *aList;
  } 
private:
  G4XpipNTotal thepipp;
  G4XpimNTotal thepimp;
};

#endif
