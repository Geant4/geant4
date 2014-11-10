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
// $Id: RegisterPhysLists.icc 66704 2013-01-10 18:20:17Z rhatcher $
//
//---------------------------------------------------------------------------
//
// ClassName:   RegisterPhysLists
//
// Author: R. Hatcher 2014-10-15
//
// Modified:
//
//----------------------------------------------------------------------------

#include "G4PhysListStamper.hh"

#include "FTFP_BERT.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTFP_BERT);

#include "FTFP_BERT_HP.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTFP_BERT_HP);

#include "FTFP_BERT_TRV.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTFP_BERT_TRV);

#include "FTFP_INCLXX.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTFP_INCLXX);

#include "FTFP_INCLXX_HP.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTFP_INCLXX_HP);

#include "FTF_BIC.hh"
G4_DECLARE_PHYSLIST_FACTORY(FTF_BIC);

#include "LBE.hh"
G4_DECLARE_PHYSLIST_FACTORY(LBE);

#include "QBBC.hh"
G4_DECLARE_PHYSLIST_FACTORY(QBBC);

#include "QGSP_BERT.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_BERT);

#include "QGSP_BERT_HP.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_BERT_HP);

#include "QGSP_BIC.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_BIC);

#include "QGSP_BIC_HP.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_BIC_HP);

#include "QGSP_FTFP_BERT.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_FTFP_BERT);

#include "QGS_BIC.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGS_BIC);

#include "QGSP_INCLXX.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_INCLXX);

#include "QGSP_INCLXX_HP.hh"
G4_DECLARE_PHYSLIST_FACTORY(QGSP_INCLXX_HP);

#include "Shielding.hh"
G4_DECLARE_PHYSLIST_FACTORY(Shielding);

/***********************************************************************
// rhatcher 2014-11-07 
// investigating why this is problematic for  mac/clang & linux/icc
// so for now don't support ShieldingLEND and ShieldingM

// some extra hoops because the physlist factory expects to be able
// to construct a list using no args (or just a verbosity)
// but 
//   "ShieldingLEND" is Shielding(verbose,"LEND");
//   "ShieldingM"    is Shielding(verbose,"HP","M");
template<class T>
class TShieldingLEND : public Shielding
{
public:
  explicit TShieldingLEND(G4int verbose = 1) : TShielding(verbose,"LEND") { ; }
  virtual ~TShieldingLEND() { ; }
};
typedef TShieldingLEND<G4VModularPhysicsList> ShieldingLEND;
G4_DECLARE_PHYSLIST_FACTORY(ShieldingLEND);

template<class T>
class TShieldingM : public Shielding
{
public:
  explicit TShieldingM(G4int verbose = 1) : TShielding(verbose,"HP","M") { ; }
  virtual ~TShieldingM() { ; }
};
typedef TShieldingM<G4VModularPhysicsList> ShieldingM;
G4_DECLARE_PHYSLIST_FACTORY(ShieldingM);

***********************************************************************/

#include "NuBeam.hh"
G4_DECLARE_PHYSLIST_FACTORY(NuBeam);

/***********************************************************************
// rhatcher 2014-11-07 
// investigating why this is problematic for  mac/clang ( & linux/icc ? )

#include "G4GenericPhysicsList.hh"
G4_DECLARE_PHYSLIST_FACTORY(G4GenericPhysicsList);

***********************************************************************/

