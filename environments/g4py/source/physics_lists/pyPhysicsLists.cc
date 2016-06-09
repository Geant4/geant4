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
// $Id: pyPhysicsLists.cc,v 1.1 2006/11/20 08:53:05 kmura Exp $
// $Name: geant4-08-02 $
// ====================================================================
//   pyPhysicsLists.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "pyG4Version.hh"

#if G4VERSION_NUMBER >= 820

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronInelasticQLHEP.hh"

#include "FTFC.hh"
#include "FTFP.hh"
#include "LBE.hh"
#include "LHEP.hh"
#include "LHEP_BERT.hh"
#include "LHEP_BERT_HP.hh"
#include "LHEP_BIC.hh"
#include "LHEP_BIC_HP.hh"
#include "LHEP_EMV.hh"
#include "LHEP_HP.hh"
#include "LHEP_LEAD.hh"
#include "LHEP_LEAD_HP.hh"
#include "LHEP_PRECO.hh"
#include "LHEP_PRECO_HP.hh"
#include "QBBC.hh"
#include "QGSC.hh"
#include "QGSC_LEAD.hh"
#include "QGSC_LEAD_HP.hh"
#include "QGSP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"
#include "QGSP_EMV.hh"
#include "QGSP_EMX.hh"
#include "QGSP_HP.hh"

#endif

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_PhysicsLists()
{
#if G4VERSION_NUMBER >= 820

  class_<FTFC, FTFC*, bases<G4VUserPhysicsList> >
    ("FTFC", "FTFC physics list")
    ;

  class_<FTFP, FTFP*, bases<G4VUserPhysicsList> >
    ("FTFP", "FTFP physics list")
    ;

  class_<LBE, LBE*, bases<G4VUserPhysicsList> >
    ("LBE", "LBE physics list")
    ;

  class_<LHEP, LHEP*, bases<G4VUserPhysicsList> >
    ("LHEP", "LHEP physics list")
    ;

  class_<LHEP_BERT, LHEP_BERT*, bases<G4VUserPhysicsList> >
    ("LHEP_BERT", "LHEP_BERT physics list")
    ;

  class_<LHEP_BERT_HP, LHEP_BERT_HP*, bases<G4VUserPhysicsList> >
    ("LHEP_BERT_HP", "LHEP_BERT_HP physics list")
    ;

  class_<LHEP_BIC, LHEP_BIC*, bases<G4VUserPhysicsList> >
    ("LHEP_BIC", "LHEP_BIC physics list")
    ;

  class_<LHEP_BIC_HP, LHEP_BIC_HP*, bases<G4VUserPhysicsList> >
    ("LHEP_BIC_HP", "LHEP_BIC_HP physics list")
    ;

  class_<LHEP_EMV, LHEP_EMV*, bases<G4VUserPhysicsList> >
    ("LHEP_EMV", "LHEP_EMV physics list")
    ;

  class_<LHEP_HP, LHEP_HP*, bases<G4VUserPhysicsList> >
    ("LHEP_HP", "LHEP_HP physics list")
    ;

  class_<LHEP_LEAD, LHEP_LEAD*, bases<G4VUserPhysicsList> >
    ("LHEP_LEAD", "LHEP_LEAD physics list")
    ;

  class_<LHEP_LEAD_HP, LHEP_LEAD_HP*, bases<G4VUserPhysicsList> >
    ("LHEP_LEAD_HP", "LHEP_LEAD_HP physics list")
    ;

  class_<LHEP_PRECO, LHEP_PRECO*, bases<G4VUserPhysicsList> >
    ("LHEP_PRECO", "LHEP_PRECO physics list")
    ;

  class_<LHEP_PRECO_HP, LHEP_PRECO_HP*, bases<G4VUserPhysicsList> >
    ("LHEP_PRECO_HP", "LHEP_PRECO_HP physics list")
    ;

  //class_<QBBC, QBBC*, bases<G4VUserPhysicsList> >
  //("QBBC", "QBBC physics list")
  //  ;

  class_<QGSC, QGSC*, bases<G4VUserPhysicsList> >
    ("QGSC", "QGSC physics list")
    ;

  class_<QGSC_LEAD, QGSC_LEAD*, bases<G4VUserPhysicsList> >
    ("QGSC_LEAD", "QGSC_LEAD physics list")
    ;

  class_<QGSC_LEAD_HP, QGSC_LEAD_HP*, bases<G4VUserPhysicsList> >
    ("QGSC_LEAD_HP", "QGSC_LEAD_HP physics list")
    ;

  class_<QGSP, QGSP*, bases<G4VUserPhysicsList> >
    ("QGSP", "QGSP physics list")
    ;

  class_<QGSP_BERT, QGSP_BERT*, bases<G4VUserPhysicsList> >
    ("QGSP_BERT", "QGSP_BERT physics list")
    ;

  class_<QGSP_BERT_HP, QGSP_BERT_HP*, bases<G4VUserPhysicsList> >
    ("QGSP_BERT_HP", "QGSP_BERT_HP physics list")
    ;

  class_<QGSP_BIC, QGSP_BIC*, bases<G4VUserPhysicsList> >
    ("QGSP_BIC", "QGSP_BIC physics list")
    ;

  class_<QGSP_EMV, QGSP_EMV*, bases<G4VUserPhysicsList> >
    ("QGSP_EMV", "QGSP_EMV physics list")
    ;

  class_<QGSP_EMX, QGSP_EMX*, bases<G4VUserPhysicsList> >
    ("QGSP_EMX", "QGSP_EMX physics list")
    ;

  class_<QGSP_HP, QGSP_HP*, bases<G4VUserPhysicsList> >
    ("QGSP_HP", "QGSP_HP physics list")
    ;

#endif

}

