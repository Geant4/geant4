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
// $Id: G4ImportanceAlgorithm.cc,v 1.9 2002-11-04 10:43:07 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceAlgorithm.cc
//
// ----------------------------------------------------------------------

#include "g4std/strstream"
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
      G4Exception("Error: G4ImportanceAlgorithm::Calculate: ipre==0.");
    }
    G4double ipre_over_ipost = ipre/ipost;
    if ((ipre_over_ipost<0.25 || ipre_over_ipost> 4) && !fWorned) {
      G4std::ostrstream os;
      os << "Calculate: ipre_over_ipost ! in [0.25, 4]: ipre_over_ipost = "
	 << ipre_over_ipost << '\0' << G4endl;
      Warning(os.str());
      fWorned = true;
      if (ipre_over_ipost<=0) {
	Error("Calculate: ipre_over_ipost<=0");
      }
    }
    if (init_w<=0.) {
      Error("Calculate:  iniitweight<= 0. found");
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
  G4cout << "ERROR - G4ImportanceAlgorithm::" << m << G4endl;
  G4Exception("Program aborted.");
}

void G4ImportanceAlgorithm::Warning(const G4String &m) const
{
  G4cout << "WARNING - G4ImportanceAlgorithm::" << m << G4endl;
}
