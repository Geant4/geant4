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
// $Id: A01HodoscopeHit.cc,v 1.7 2005/06/07 10:50:02 perl Exp $
// --------------------------------------------------------------
//
#include "A01HodoscopeHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Allocator<A01HodoscopeHit> A01HodoscopeHitAllocator;

A01HodoscopeHit::A01HodoscopeHit(G4int i,G4double t)
{
  id = i;
  time = t;
  pLogV = 0;
}

A01HodoscopeHit::~A01HodoscopeHit()
{;}

A01HodoscopeHit::A01HodoscopeHit(const A01HodoscopeHit &right)
    : G4VHit() {
  id = right.id;
  time = right.time;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
}

const A01HodoscopeHit& A01HodoscopeHit::operator=(const A01HodoscopeHit &right)
{
  id = right.id;
  time = right.time;
  pos = right.pos;
  rot = right.rot;
  pLogV = right.pLogV;
  return *this;
}

int A01HodoscopeHit::operator==(const A01HodoscopeHit &/*right*/) const
{
  return 0;
}

void A01HodoscopeHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Transform3D trans(rot.inverse(),pos);
    G4VisAttributes attribs;
    const G4VisAttributes* pVA = pLogV->GetVisAttributes();
    if(pVA) attribs = *pVA;
    G4Colour colour(0.,1.,1.);
    attribs.SetColour(colour);
    attribs.SetForceSolid(true);
    pVVisManager->Draw(*pLogV,attribs,trans);
  }
}

const std::map<G4String,G4AttDef>* A01HodoscopeHit::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("A01HodoscopeHit",isNew);
  if (isNew) {
    G4String HitType("HitType");
    (*store)[HitType] = G4AttDef(HitType,"Hit Type","Bookkeeping","","G4String");

    G4String ID("ID");
    (*store)[ID] = G4AttDef(ID,"ID","Bookkeeping","","G4int");

    G4String Column("Column");
    (*store)[Column] = G4AttDef(Column,"Column ID","Bookkeeping","","G4int");

    G4String Row("Row");
    (*store)[Row] = G4AttDef(Row,"Row ID","Bookkeeping","","G4int");

    //G4String Time("Time");
    //(*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

    G4String Energy("Energy");
    (*store)[Energy] = G4AttDef(Energy,"Energy Deposited","Physics","G4BestUnit","G4double");

    G4String Pos("Pos");
    (*store)[Pos] = G4AttDef(Pos, "Position",
		      "Physics","G4BestUnit","G4ThreeVector");

    G4String LVol("LVol");
    (*store)[LVol] = G4AttDef(LVol,"Logical Volume","Bookkeeping","","G4String");
  }
  return store;
}

std::vector<G4AttValue>* A01HodoscopeHit::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back(G4AttValue("HitType","HodoscopeHit",""));

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(id),""));

  values->push_back
    (G4AttValue("Column"," ",""));

  values->push_back
    (G4AttValue("Row"," ",""));

  //values->push_back
  //  (G4AttValue("Time",G4BestUnit(time,"Time"),""));

  G4double noEnergy = 0.*MeV;
  values->push_back
    (G4AttValue("Energy",G4BestUnit(noEnergy,"Energy"),""));

  values->push_back
    (G4AttValue("Pos",G4BestUnit(pos,"Length"),""));

  if (pLogV)
    values->push_back
      (G4AttValue("LVol",pLogV->GetName(),""));
  else
    values->push_back
      (G4AttValue("LVol"," ",""));

  return values;
}

void A01HodoscopeHit::Print()
{
    G4cout << "  Hodoscope[" << id << "] " << time/ns << " (nsec)" << G4endl;
}



