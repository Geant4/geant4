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

/** \file G4INCLNuclearMassTable.cc
 * \brief Functions that encapsulate a mass table
 *
 * \date 22nd October 2013
 * \author Davide Mancusi
 */

#ifndef INCLXX_IN_GEANT4_MODE

#include "G4INCLNuclearMassTable.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLGlobals.hh"
#include <algorithm>
#include <istream>

namespace G4INCL {

  namespace {

    G4ThreadLocal G4double **theTable = NULL;
    G4ThreadLocal G4int AMax = 0;
    G4ThreadLocal G4int *ZMaxArray = NULL;
    G4ThreadLocal G4double protonMass = 0.;
    G4ThreadLocal G4double neutronMass = 0.;

    const G4double amu = 931.494061; // atomic mass unit in MeV/c^2
    const G4double eMass = 0.5109988; // electron mass in MeV/c^2

    G4double getWeizsaeckerMass(const G4int A, const G4int Z) {
      const G4int Npairing = (A-Z)%2;                  // pairing
      const G4int Zpairing = Z%2;
      const G4double fA = (G4double) A;
      const G4double fZ = (G4double) Z;
      G4double binding =
        - 15.67*fA                          // nuclear volume
        + 17.23*Math::pow23(fA)                // surface energy
        + 93.15*((fA/2.-fZ)*(fA/2.-fZ))/fA       // asymmetry
        + 0.6984523*fZ*fZ*Math::powMinus13(fA);      // coulomb
      if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(fA);  // pairing

      return fZ*::G4INCL::ParticleTable::getRealMass(Proton)+((G4double)(A-Z))
        *::G4INCL::ParticleTable::getRealMass(Neutron)+binding;
    }

    void setMass(const G4int A, const G4int Z, const G4double mass) {
      theTable[A][Z] = mass;
    }

    class MassRecord {
      public:
        MassRecord() :
          A(0),
          Z(0),
          excess(0.)
      {}

        MassRecord(const G4int a, const G4int z, const G4double e) :
          A(a),
          Z(z),
          excess(e)
      {}

        friend std::istream &operator>>(std::istream &in, MassRecord &record);

        G4int A;
        G4int Z;
        G4double excess;
    };

    std::istream &operator>>(std::istream &in, MassRecord &record) {
      return (in >> record.A >> record.Z >> record.excess);
    }

    G4bool compareA(const MassRecord &lhs, const MassRecord &rhs) {
      return (lhs.A < rhs.A);
    }

  }

  namespace NuclearMassTable {

    void initialize(const std::string &path, const G4double pMass, const G4double nMass) {
      protonMass = pMass;
      neutronMass = nMass;

      // Clear the existing tables, if any
      deleteTable();

      // File name
      std::string fileName(path + "/walletlifetime.dat");
      INCL_DEBUG("Reading real nuclear masses from file " << fileName << '\n');

      // Open the file stream
      std::ifstream massTableIn(fileName.c_str());
      if(!massTableIn.good()) {
        std::cerr << "Cannot open " << fileName << " data file." << '\n';
        std::abort();
        return;
      }

      // read the file
      std::vector<MassRecord> records;
      MassRecord record;
      while(massTableIn.good()) { /* Loop checking, 10.07.2015, D.Mancusi */
        massTableIn >> record;
        records.push_back(record);
      }
      massTableIn.close();
      INCL_DEBUG("Read " << records.size() << " nuclear masses" << '\n');

      // determine the max A
      AMax = std::max_element(records.begin(), records.end(), compareA)->A;
      INCL_DEBUG("Max A in nuclear-mass table = " << AMax << '\n');
      ZMaxArray = new G4int[AMax+1];
      std::fill(ZMaxArray, ZMaxArray+AMax+1, 0);
      theTable = new G4double*[AMax+1];
      std::fill(theTable, theTable+AMax+1, static_cast<G4double*>(NULL));

      // determine the max A per Z
      for(std::vector<MassRecord>::const_iterator i=records.begin(), e=records.end(); i!=e; ++i) {
        ZMaxArray[i->A] = std::max(ZMaxArray[i->A], i->Z);
      }

      // allocate the arrays
      for(G4int A=1; A<=AMax; ++A) {
        theTable[A] = new G4double[ZMaxArray[A]+1];
        std::fill(theTable[A], theTable[A]+ZMaxArray[A]+1, -1.);
      }

      // fill the actual masses
      for(std::vector<MassRecord>::const_iterator i=records.begin(), e=records.end(); i!=e; ++i) {
        setMass(i->A, i->Z, i->A*amu + i->excess - i->Z*eMass);
      }
    }

    G4double getMass(const G4int A, const G4int Z) {
      if(A>AMax || Z>ZMaxArray[A]) {
        INCL_DEBUG("Real mass unavailable for isotope A=" << A << ", Z=" << Z
                   << ", using Weizsaecker's formula"
                   << '\n');
        return getWeizsaeckerMass(A,Z);
      }

      const G4double mass = theTable[A][Z];
      if(mass<0.) {
        INCL_DEBUG("Real mass unavailable for isotope A=" << A << ", Z=" << Z
                   << ", using Weizsaecker's formula"
                   << '\n');
        return getWeizsaeckerMass(A,Z);
      } else
        return mass;
    }

    void deleteTable() {
      delete[] ZMaxArray;
      ZMaxArray = NULL;
      for(G4int A=1; A<=AMax; ++A)
        delete[] theTable[A];
      delete[] theTable;
      theTable = NULL;
    }
  }

}

#endif // INCLXX_IN_GEANT4_MODE
