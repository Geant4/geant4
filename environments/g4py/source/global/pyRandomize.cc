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
// $Id: pyRandomize.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyRandomize.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "Randomize.hh"
#include <vector>

using namespace boost::python;
using namespace CLHEP;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyRandomize {

// problem with BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS for static method
// setTheSeed
void f1_setTheSeed(long seed)
{
  HepRandom::setTheSeed(seed);
}

void f2_setTheSeed(long seed, int lux)
{
  HepRandom::setTheSeed(seed, lux);
}

// setTheSeeds
void f1_setTheSeeds(const list& seedList)
{
  // check size...
  int idx=0;
  while(1) {
    long val= extract<long>(seedList[idx]);
    if(val==0) break;
    idx++;
  }
  int nsize= idx+1;

  long* seedArray= new long[nsize]; // should not be deleted!!
                                    // (this is a problem with CLHEP.)
  for (int i=0; i< nsize; i++) {
    seedArray[i]= extract<long>(seedList[i]);
  }

  HepRandom::setTheSeeds(seedArray);
}

void f2_setTheSeeds(const list& seedList, int aux)
{
  // check size...
  int idx=0;
  while(1) {
    long val= extract<long>(seedList[idx]);
    if(val==0) break;
    idx++;
  }
  int nsize= idx+1;

  long* seedArray= new long[nsize];
  for (int i=0; i< nsize; i++) {
    seedArray[i]= extract<long>(seedList[i]);
  }

  HepRandom::setTheSeeds(seedArray, aux);
}

// getTheSeeds
list f_getTheSeeds()
{
  list seedList;
  const long* seeds= HepRandom::getTheSeeds();
  int idx=0;
  while(1) {
    if( seeds[idx]==0) break;
    seedList.append(seeds[idx]);
    idx++;
  }
  return seedList;
}

// getTheTableSeeds
list f_getTheTableSeeds(int index)
{
  long seedPair[2];
  HepRandom::getTheTableSeeds(seedPair, index);

  list seedList;
  seedList.append(seedPair[0]);
  seedList.append(seedPair[1]);

  return seedList;
}


// saveEngineStatus
void f1_saveEngineStatus()
{
  HepRandom::saveEngineStatus();
}

void f2_saveEngineStatus(const char* filename)
{
  HepRandom::saveEngineStatus(filename);
}

// restoreEngineStatus
void f1_restoreEngineStatus()
{
  HepRandom::restoreEngineStatus();
}

void f2_restoreEngineStatus(const char* filename)
{
  HepRandom::restoreEngineStatus(filename);
}

// RandBit::shootBit
int f1_RandBit_shootBit()
{
  return RandBit::shootBit();
}

// RandGaussQ::shoot
double f1_RandGaussQ_shoot()
{
  return RandGaussQ::shoot();
}

double f2_RandGaussQ_shoot(double mean, double stdDev)
{
  return RandGaussQ::shoot(mean, stdDev);
}


// G4UniformRand
double f_G4UniformRand()
{
  return G4UniformRand();
}

}

using namespace pyRandomize;

// ====================================================================
// module definition
// ====================================================================
void export_Randomize()
{
  class_<HepRandom>("HepRandom", "generate random number")
    // ---
    .def(init<long>())
    .def(init<HepRandomEngine&>())
    .def(init<HepRandomEngine*>())
    // ---
    .def("setTheSeed",       f1_setTheSeed)
    .def("setTheSeed",       f2_setTheSeed)
    .staticmethod("setTheSeed")
    .def("getTheSeed",           &HepRandom::getTheSeed)
    .staticmethod("getTheSeed")
    .def("setTheSeeds",      f1_setTheSeeds)
    .def("setTheSeeds",      f2_setTheSeeds)
    .staticmethod("setTheSeeds")
    .def("getTheSeeds",      f_getTheSeeds)
    .staticmethod("getTheSeeds")
    .def("getTheTableSeeds", f_getTheTableSeeds)
    .staticmethod("getTheTableSeeds")
    // ---
    .def("getTheGenerator",     &HepRandom::getTheGenerator,
         return_value_policy<reference_existing_object>())
    .staticmethod("getTheGenerator")
    .def("setTheEngine",        &HepRandom::setTheEngine)
    .staticmethod("setTheEngine")
    .def("getTheEngine",        &HepRandom::getTheEngine,
         return_value_policy<reference_existing_object>())
    .staticmethod("getTheEngine")
    .def("saveEngineStatus",    f1_saveEngineStatus)
    .def("saveEngineStatus",    f2_saveEngineStatus)
    .staticmethod("saveEngineStatus")
    .def("restoreEngineStatus", f1_restoreEngineStatus)
    .def("restoreEngineStatus", f2_restoreEngineStatus)
    .staticmethod("restoreEngineStatus")
    .def("showEngineStatus",    &HepRandom::showEngineStatus)
    .staticmethod("showEngineStatus")
    .def("createInstance",      &HepRandom::createInstance)
    .staticmethod("createInstance")
    ;

  // ---
  class_<RandBit, boost::noncopyable>
    ("RandBit", "generate bit random number", no_init)
    .def("shootBit", f1_RandBit_shootBit)
    .staticmethod("shootBit")
    ;

  // ---
  class_<G4RandGauss, boost::noncopyable>
    ("G4RandGauss", "generate gaussian random number", no_init)
    .def("shoot", f1_RandGaussQ_shoot)
    .def("shoot", f2_RandGaussQ_shoot)
    .staticmethod("shoot")
    ;

  // ---
  def("G4UniformRand", f_G4UniformRand);

}

