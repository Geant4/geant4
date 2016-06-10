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
// $Id: pyPhysicsLists.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyPhysicsLists.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include <vector>
#include <algorithm>
#include "FTF_BIC.hh"
#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "FTFP_BERT_TRV.hh"
#include "FTFP_INCLXX.hh"
#include "FTFP_INCLXX_HP.hh"
#include "LBE.hh"
#include "QBBC.hh"
#include "QGS_BIC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_FTFP_BERT.hh"
#include "QGSP_INCLXX.hh"
#include "QGSP_INCLXX_HP.hh"
#include "Shielding.hh"

// macro for adding physics lists
#define ADD_PHYSICS_LIST(plname) \
  class_<plname, plname*, bases<G4VUserPhysicsList>, boost::noncopyable> \
    (#plname, #plname " physics list") \
    ; \
  AddPhysicsList(#plname);

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyPhysicsLists {

  static std::vector<std::string> plList;

void AddPhysicsList(const G4String& plname) {
  plList.push_back(plname);
}

void ListPhysicsList() {
  for (G4int i=0; i< plList.size(); i++) {
    G4cout << plList[i] << G4endl;
  }
}

}

using namespace pyPhysicsLists;

// ====================================================================
// module definition
// ====================================================================
void export_PhysicsLists()
{
  def("ListPhysicsList",   ListPhysicsList);

  ADD_PHYSICS_LIST(FTF_BIC);
  ADD_PHYSICS_LIST(FTFP_BERT);
  ADD_PHYSICS_LIST(FTFP_BERT_HP);
  ADD_PHYSICS_LIST(FTFP_BERT_TRV);
  ADD_PHYSICS_LIST(FTFP_INCLXX);
  ADD_PHYSICS_LIST(FTFP_INCLXX_HP);
  ADD_PHYSICS_LIST(LBE);
  ADD_PHYSICS_LIST(QBBC);
  ADD_PHYSICS_LIST(QGS_BIC);
  ADD_PHYSICS_LIST(QGSP_BERT);
  ADD_PHYSICS_LIST(QGSP_BERT_HP);
  ADD_PHYSICS_LIST(QGSP_BIC);
  ADD_PHYSICS_LIST(QGSP_BIC_HP);
  ADD_PHYSICS_LIST(QGSP_FTFP_BERT);
  ADD_PHYSICS_LIST(QGSP_INCLXX);
  ADD_PHYSICS_LIST(QGSP_INCLXX_HP);
  ADD_PHYSICS_LIST(Shielding);

  // sort PL vector
  std::sort(plList.begin(), plList.end());
}
