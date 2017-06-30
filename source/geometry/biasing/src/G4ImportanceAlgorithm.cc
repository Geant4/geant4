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
// $Id: G4ImportanceAlgorithm.cc 102994 2017-03-07 16:31:28Z gcosmo $
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
#include "G4Threading.hh"

#include "G4ImportanceAlgorithm.hh"

#ifdef G4MULTITHREADED
G4Mutex G4ImportanceAlgorithm::ImportanceMutex = G4MUTEX_INITIALIZER;
#endif

G4ImportanceAlgorithm::G4ImportanceAlgorithm(): fWorned(false)
{
}

G4ImportanceAlgorithm::~G4ImportanceAlgorithm()
{
}

G4Nsplit_Weight
G4ImportanceAlgorithm::Calculate(G4double ipre,
				 G4double ipost,
                                 G4double init_w) const
{

#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ImportanceAlgorithm::ImportanceMutex);
#endif

  G4Nsplit_Weight nw = {0,0};
  if (ipost>0.){
    if (!(ipre>0.)){
      Error("Calculate() - ipre==0.");
    }
    G4double ipre_over_ipost = ipre/ipost;
    if ((ipre_over_ipost<0.25 || ipre_over_ipost> 4) && !fWorned) {
      std::ostringstream os;
      os << "Calculate() - ipre_over_ipost ! in [0.25, 4]." << G4endl
         << "ipre_over_ipost = " << ipre_over_ipost << ".";
      Warning(os.str());
      fWorned = true;
      if (ipre_over_ipost<=0) {
	Error("Calculate() - ipre_over_ipost<=0.");
      }
    }
    if (init_w<=0.) {
      Error("Calculate() - iniitweight<= 0. found!");
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
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ImportanceAlgorithm::ImportanceMutex);
#endif
  return nw;
}

void G4ImportanceAlgorithm::Error(const G4String &msg) const
{
  G4Exception("G4ImportanceAlgorithm::Error()",
              "GeomBias0002", FatalException, msg);
}

void G4ImportanceAlgorithm::Warning(const G4String &msg) const
{
  G4Exception("G4ImportanceAlgorithm::Warning()",
              "GeomBias1001", JustWarning, msg);
}
