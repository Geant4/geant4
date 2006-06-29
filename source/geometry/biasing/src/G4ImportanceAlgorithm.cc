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
// $Id: G4ImportanceAlgorithm.cc,v 1.14 2006-06-29 18:17:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceAlgorithm.cc
//
// ----------------------------------------------------------------------

#include "G4Types.hh"
#include <sstream>
#include "Randomize.hh"

#include "G4ImportanceAlgorithm.hh"

G4ImportanceAlgorithm::G4ImportanceAlgorithm(): fWorned(false)
{}

G4ImportanceAlgorithm::~G4ImportanceAlgorithm()
{
  if(fWorned) {
    G4cout << G4endl;
    Warning("~G4ImportanceAlgorithm: ipre_over_ipost ! in [0.25, 4] seen");
    G4cout << G4endl;
  }
}

G4Nsplit_Weight
G4ImportanceAlgorithm::Calculate(G4double ipre,
				 G4double ipost,
                                 G4double init_w) const
{
  G4Nsplit_Weight nw = {0,0};
  if (ipost>0.){
    if (!(ipre>0.)){
      Error("G4ImportanceAlgorithm::Calculate() - ipre==0.");
    }
    G4double ipre_over_ipost = ipre/ipost;
    if ((ipre_over_ipost<0.25 || ipre_over_ipost> 4) && !fWorned) {
      std::ostringstream os;
      os << "Calculate: ipre_over_ipost ! in [0.25, 4]: ipre_over_ipost = "
	 << ipre_over_ipost << '\0' << G4endl;
      Warning(os.str());
      fWorned = true;
      if (ipre_over_ipost<=0) {
	Error("G4ImportanceAlgorithm::Calculate() - ipre_over_ipost<=0");
      }
    }
    if (init_w<=0.) {
      Error("G4ImportanceAlgorithm::Calculate() - iniitweight<= 0. found");
    }

    // default geometrical splitting 
    // in integer mode 
    // for ipre_over_ipost <= 1
    G4double inv = 1./ipre_over_ipost;
    nw.fN = static_cast<G4int>(inv);
    nw.fW = init_w * ipre_over_ipost;
    
    // geometrical splitting for double mode
    if (ipre_over_ipost<1) {
      if ( static_cast<G4double>(nw.fN) != inv) {
	// double mode
	// probability p for splitting into n+1 tracks
	G4double p = inv - nw.fN;
	// get a random number out of [0,1)
	G4double r = G4UniformRand();
	if (r<p) {
	  nw.fN++;
	} 
      }  
    }
    // ipre_over_ipost > 1
    //  russian roulett
    else if (ipre_over_ipost>1) {
      // probabiity for killing track
      G4double p = 1-inv;
      // get a random number out of [0,1)
      G4double r = G4UniformRand();
      if (r<p) {
	// kill track
	nw.fN = 0;
	nw.fW = 0;
      }
      else {
	nw.fN = 1;     
      }
    }
  }
  return nw;
}

void G4ImportanceAlgorithm::Error(const G4String &m) const
{
  G4cerr << "ERROR - G4ImportanceAlgorithm: " << m << G4endl;
  G4Exception("G4ImportanceAlgorithm::Error()",
              "FatalException", FatalException, m);
}

void G4ImportanceAlgorithm::Warning(const G4String &m) const
{
  G4Exception("G4ImportanceAlgorithm::Warning()",
              "Notification", JustWarning, m);
}
