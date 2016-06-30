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
// $Id: G4MultiFunctionalDetector.cc 94771 2015-12-09 09:44:05Z gcosmo $
//
// G4MultiFunctionalDetector
//
// 2010-07-23 T.Aso Call PS if the step length or energy deposit is not zero.
//       
//
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"

G4MultiFunctionalDetector::G4MultiFunctionalDetector(G4String name)
:G4VSensitiveDetector(name)
{;}

G4MultiFunctionalDetector::~G4MultiFunctionalDetector()
{
  for(auto pr : primitives) delete pr;
  primitives.clear();
}

G4bool G4MultiFunctionalDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* aTH)
{
  if(aStep->GetStepLength()>0. || aStep->GetTotalEnergyDeposit()>0.){
    for(auto pr : primitives)
    { pr->HitPrimitive(aStep,aTH); }
  }
  return true;
}

G4bool G4MultiFunctionalDetector::RegisterPrimitive(G4VPrimitiveScorer* aPS)
{
   for(auto pr : primitives)
   {
     if(pr==aPS)
     {
       G4ExceptionDescription ED;
       ED << "Primitive <" << aPS->GetName() << "> is already defined in <" << SensitiveDetectorName
              << ">." << G4endl << "Method RegisterPrimitive() is ignored." << G4endl;
       G4Exception("G4MultiFunctionalDetector::RegisterPrimitive","Det0101",
                   JustWarning,ED);
       return false;
     }
   }
   primitives.push_back(aPS);
   aPS->SetMultiFunctionalDetector(this);
   collectionName.insert(aPS->GetName());
   if(G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName,false))
   {
    // This G4MultiFunctionalDetector has already been registered to G4SDManager.
    // Make sure this new primitive is registered as well.
    G4SDManager::GetSDMpointer()->AddNewCollection(SensitiveDetectorName,aPS->GetName());
   }
   return true;
}

G4bool G4MultiFunctionalDetector::RemovePrimitive(G4VPrimitiveScorer* aPS)
{
   auto pr = std::find(primitives.begin(),primitives.end(),aPS);
   if(pr!=primitives.end())
   { 
     primitives.erase(pr);
     aPS->SetMultiFunctionalDetector(nullptr);
     return true;
   }
   G4cerr << "Primitive <" << aPS->GetName() << "> is not defined in <" << SensitiveDetectorName
          << ">." << G4endl << "Method RemovePrimitive() is ignored." << G4endl;
   return false;
}   

void G4MultiFunctionalDetector::Initialize(G4HCofThisEvent* HC)
{
   for(auto pr : primitives)
   { pr->Initialize(HC); }
}

void G4MultiFunctionalDetector::EndOfEvent(G4HCofThisEvent* HC)
{
   for(auto pr : primitives)
   { pr->EndOfEvent(HC); }
}

void G4MultiFunctionalDetector::clear()
{
   for(auto pr : primitives)
   { pr->clear(); }
}

void G4MultiFunctionalDetector::DrawAll()
{
   for(auto pr : primitives)
   { pr->DrawAll(); }
}

void G4MultiFunctionalDetector::PrintAll()
{
   for(auto pr : primitives)
   { pr->PrintAll(); }
}

