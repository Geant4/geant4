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
// $Id: G4CollisionMesonBaryonToResonance.hh,v 1.3 2006-06-29 20:32:52 gunter Exp $ //

#ifndef G4CollisionMesonBaryonToResonance_h
#define G4CollisionMesonBaryonToResonance_h

#include "globals.hh"
#include "G4CollisionComposite.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include <vector>
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

private:
  G4CollisionMesonBaryonToResonance(const G4CollisionMesonBaryonToResonance &);
  G4CollisionMesonBaryonToResonance & operator= (const G4CollisionMesonBaryonToResonance &);


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
  
  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    throw G4HadronicException(__FILE__, __LINE__, "Tried to call G4CollisionMesonBaryonToResonance::GetListOfColliders. Please find out why!");
    std::vector<G4String> * aList = new std::vector<G4String>;
    return *aList;
  } 
private:
  G4XpipNTotal thepipp;
  G4XpimNTotal thepimp;
};

#endif
