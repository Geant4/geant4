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
// $Id: pyPhysicsLists.cc,v 1.9 2009/11/20 03:36:51 kmura Exp $
// $Name: geant4-09-03 $
// ====================================================================
//   pyPhysicsLists.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include <vector>
#include <algorithm>

#if G4VERSION_NUMBER >= 820

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronInelasticQLHEP.hh"

#include "FTFC.hh"
#include "FTFP.hh"
#include "LBE.hh"
#include "LHEP.hh"
#include "LHEP_BERT.hh"
#include "LHEP_BERT_HP.hh"
#include "LHEP_EMV.hh"
#include "LHEP_PRECO_HP.hh"
#include "QBBC.hh"
#include "QGSC.hh"
#include "QGSC_EFLOW.hh"
#include "QGSC_EMV.hh"
#include "QGSP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_EMV.hh"
#include "QGSP_EMX.hh"
#include "QGSP_QEL.hh"

#if G4VERSION_NUMBER <= 830
#include "LHEP_BIC.hh"
#include "LHEP_BIC_HP.hh"
#include "LHEP_HP.hh"
#include "LHEP_LEAD.hh"
#include "LHEP_LEAD_HP.hh"
#include "LHEP_PRECO.hh"
#include "QGSC_LEAD.hh"
#include "QGSC_LEAD_HP.hh"
#include "QGSP_HP.hh"
#endif

#endif

#if G4VERSION_NUMBER >= 830
#include "FTFP_EMV.hh"
#include "QGSP_BERT_EMV.hh"
#include "QGSP_BERT_NQE.hh"
#include "QGSP_BERT_TRV.hh"
#include "QGSP_EMV_NQE.hh"
#include "QGSP_NQE.hh"
#endif

#if G4VERSION_NUMBER >= 910
#include "FTFP_BERT.hh"
#include "FTF_BIC.hh"
#include "QGSC_BERT.hh"
#include "QGS_BIC.hh"
#include "QGSP_DIF.hh"
#include "QGSP_BERT_DIF.hh"
#endif

#if G4VERSION_NUMBER >= 930
#include "CHIPS.hh"
#include "FTFP_BERT_EMV.hh"
#include "FTFP_BERT_EMX.hh"
#include "FTFP_BERT_TRV.hh"
#include "QGSC_CHIPS.hh"
#include "QGSC_QGSC.hh"
#include "QGSP_BERT_EMX.hh"
#include "QGSP_BERT_NOLEP.hh"
#include "QGSP_BIC_EMY.hh"
#include "QGSP_FTFP_BERT.hh"
#include "QGSP_INCL_ABLA.hh"
#endif

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

};

using namespace pyPhysicsLists;

// ====================================================================
// module definition
// ====================================================================
void export_PhysicsLists()
{
  def("ListPhysicsList",   ListPhysicsList);

#if G4VERSION_NUMBER >= 820
  ADD_PHYSICS_LIST(FTFC);
  ADD_PHYSICS_LIST(FTFP);
  ADD_PHYSICS_LIST(LBE);
  ADD_PHYSICS_LIST(LHEP);
  ADD_PHYSICS_LIST(LHEP_BERT);
  ADD_PHYSICS_LIST(LHEP_BERT_HP);
  ADD_PHYSICS_LIST(LHEP_EMV);
  ADD_PHYSICS_LIST(LHEP_PRECO_HP);
  ADD_PHYSICS_LIST(QBBC);
  ADD_PHYSICS_LIST(QGSC);
  ADD_PHYSICS_LIST(QGSC_EFLOW);
  ADD_PHYSICS_LIST(QGSC_EMV);
  ADD_PHYSICS_LIST(QGSP);
  ADD_PHYSICS_LIST(QGSP_BERT);
  ADD_PHYSICS_LIST(QGSP_BERT_HP);
  ADD_PHYSICS_LIST(QGSP_BIC);
  ADD_PHYSICS_LIST(QGSP_BIC_HP);
  ADD_PHYSICS_LIST(QGSP_EMV);
  ADD_PHYSICS_LIST(QGSP_EMX);
  ADD_PHYSICS_LIST(QGSP_QEL);

#if G4VERSION_NUMBER <= 830
  ADD_PHYSICS_LIST(LHEP_BIC);
  ADD_PHYSICS_LIST(LHEP_BIC_HP);
  ADD_PHYSICS_LIST(LHEP_HP);
  ADD_PHYSICS_LIST(LHEP_LEAD);
  ADD_PHYSICS_LIST(LHEP_LEAD_HP);
  ADD_PHYSICS_LIST(LHEP_PRECO);
  ADD_PHYSICS_LIST(QGSC_LEAD);
  ADD_PHYSICS_LIST(QGSC_LEAD_HP);
  ADD_PHYSICS_LIST(QGSP_HP);
#endif

#endif

#if G4VERSION_NUMBER >= 830
  ADD_PHYSICS_LIST(FTFP_EMV);
  ADD_PHYSICS_LIST(QGSP_BERT_EMV);
  ADD_PHYSICS_LIST(QGSP_BERT_NQE);
  ADD_PHYSICS_LIST(QGSP_BERT_TRV);
  ADD_PHYSICS_LIST(QGSP_EMV_NQE);
  ADD_PHYSICS_LIST(QGSP_NQE);
#endif

#if G4VERSION_NUMBER >= 910
  ADD_PHYSICS_LIST(FTFP_BERT);
  ADD_PHYSICS_LIST(FTF_BIC);
  ADD_PHYSICS_LIST(QGSC_BERT);
  ADD_PHYSICS_LIST(QGS_BIC);
  ADD_PHYSICS_LIST(QGSP_DIF);
  ADD_PHYSICS_LIST(QGSP_BERT_DIF);
#endif

#if G4VERSION_NUMBER >= 930
  ADD_PHYSICS_LIST(CHIPS);
  ADD_PHYSICS_LIST(FTFP_BERT_EMV);
  ADD_PHYSICS_LIST(FTFP_BERT_EMX);
  ADD_PHYSICS_LIST(FTFP_BERT_TRV);
  ADD_PHYSICS_LIST(QGSC_CHIPS);
  ADD_PHYSICS_LIST(QGSC_QGSC);
  ADD_PHYSICS_LIST(QGSP_BERT_EMX);
  ADD_PHYSICS_LIST(QGSP_BERT_NOLEP);
  ADD_PHYSICS_LIST(QGSP_BIC_EMY);
  ADD_PHYSICS_LIST(QGSP_FTFP_BERT);
  ADD_PHYSICS_LIST(QGSP_INCL_ABLA);
#endif

  // sort PL vector
  std::sort(plList.begin(), plList.end());

}

