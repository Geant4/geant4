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
 * G4INCLRanecu3.cc
 *
 *  \date 7 juin 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLRanecu3.hh"

namespace G4INCL {

  Ranecu3::Ranecu3() :
    iseed1(666),
    iseed2(777),
    iseed3(1234),
    i1(0), i2(0), i3(0), iz(0),
    uscale(1.0/2.147483563e9),
    m1(2147483563), m2(2147483399), m3(2147482739),
    a1(40014), a2(40692), a3(45742),
    q1(m1/a1), q2(m2/a2), q3(m3/a3),
    r1(m1%a1), r2(m2%a2), r3(m3%a3)
  {
  }

  Ranecu3::Ranecu3(const Random::SeedVector &sv) :
    i1(0), i2(0), i3(0), iz(0),
    uscale(1.0/2.147483563e9),
    m1(2147483563), m2(2147483399), m3(2147482739),
    a1(53668), a2(52774), a3(46947),
    q1(m1/a1), q2(m2/a2), q3(m3/a3),
    r1(m1%a1), r2(m2%a2), r3(m3%a3)
  {
    setSeeds(sv);
  }

  Ranecu3::~Ranecu3() {}

  G4double Ranecu3::flat()
  {
    i1=iseed1/q1;
    iseed1=a1*(iseed1-i1*q1)-i1*r1;
    if(iseed1 < 0) iseed1 = iseed1 + m1;

    i2=iseed2/q2;
    iseed2=a2*(iseed2-i2*q2)-i2*r2;
    if(iseed2 < 0) iseed2 = iseed2 + m2;

    i3=iseed3/q3;
    iseed3=a3*(iseed3-i3*q3)-i3*r3;
    if(iseed3 < 0) iseed3 = iseed3 + m3;

    iz = iseed1 - iseed2 + iseed3;
    if(iz < 1) iz = iz + 2147483562;

    return iz*uscale;
  }

}
