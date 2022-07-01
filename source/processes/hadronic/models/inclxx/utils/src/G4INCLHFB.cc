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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/*
 * HFB.cc
 *
 *  \date Oct 25, 2017
 * \author Jose Luis Rodriguez Sanchez
 */

#include "G4INCLHFB.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"
#include "G4Threading.hh"
#include <algorithm>
#include <istream>

namespace G4INCL {

  namespace {
   #ifdef INCLXX_IN_GEANT4_MODE
         G4ThreadLocal G4double radiusP[TableZSize][TableASize];
         G4ThreadLocal G4double radiusN[TableZSize][TableASize];
         G4ThreadLocal G4double diffusenessP[TableZSize][TableASize];
         G4ThreadLocal G4double diffusenessN[TableZSize][TableASize];
   #else
         G4double radiusP[TableZSize][TableASize];
         G4double radiusN[TableZSize][TableASize];
         G4double diffusenessP[TableZSize][TableASize];
         G4double diffusenessN[TableZSize][TableASize];
   #endif

   void cleanTable(){
    for(G4int i=0;i<TableZSize;++i)
     for(G4int j=0;j<TableASize;++j){
      radiusP[i][j]=-1.;
      radiusN[i][j]=-1.;
      diffusenessP[i][j]=-1.;
      diffusenessN[i][j]=-1.;
     }
   }
  }

  namespace HFB {

#ifdef INCLXX_IN_GEANT4_MODE
      void initialize() {
#else
      void initialize(const std::string &path) {
#endif

      // Clear the existing tables, if any
      cleanTable();

#ifdef INCLXX_IN_GEANT4_MODE
       if(!G4FindDataDir("G4INCLDATA")) {
        G4ExceptionDescription ed;
        ed << " Data missing: set environment variable G4INCLDATA\n"
           << " to point to the directory containing data files needed\n"
           << " by the INCL++ model" << G4endl;
           G4Exception("G4INCLDataFile::readData()","table_radius_hfb.dat",
                FatalException, ed);
      }
      G4String dataPath0(G4FindDataDir("G4INCLDATA"));
      G4String dataPath(dataPath0 + "/table_radius_hfb.dat");
#else
      // File name
      std::string dataPath(path + "/table_radius_hfb.dat");
      INCL_DEBUG("Reading radius and diffuseness parameters from file " << dataPath << '\n');
#endif

      // Open the file stream
      std::ifstream hfbTableIn(dataPath.c_str());
      if(!hfbTableIn.good()) {
        std::cerr << "Cannot open " << dataPath << " data file." << '\n';
        std::abort();
        return;
      }

      // read the file
      G4int z, a, nbnuclei=0;
      G4double rp, rn, dp, dn;
      while(hfbTableIn.good()) { /* Loop checking, 22.01.2018, J.L. Rodriguez */
        hfbTableIn >> z >> a >> rp >> rn >> dp >> dn;
        radiusP[z][a] = rp;
        radiusN[z][a] = rn;
        diffusenessP[z][a] = dp;
        diffusenessN[z][a] = dn;
        nbnuclei++;
      }
      hfbTableIn.close();
      INCL_DEBUG("Read " << nbnuclei << " nuclei" << '\n');

    }

    G4double getRadiusParameterHFB(const ParticleType t, const G4int A, const G4int Z){
        // HFB calculations
      G4double r0=0.;
        if(t==Neutron)
         if(radiusN[Z][A]>0.)r0=radiusN[Z][A];
        if(t==Proton)
         if(radiusP[Z][A]>0.)r0=radiusP[Z][A];
      return r0;
    }

    G4double getSurfaceDiffusenessHFB(const ParticleType t, const G4int A, const G4int Z){
        // HFB calculations
      G4double a=0.;
        if(t==Neutron)
          if(diffusenessN[Z][A]>0.)a=diffusenessN[Z][A];
        if(t==Proton)
          if(diffusenessP[Z][A]>0.)a=diffusenessP[Z][A];
      return a;
    }
  }
}
