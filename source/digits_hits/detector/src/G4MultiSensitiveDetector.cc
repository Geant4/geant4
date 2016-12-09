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
// G4MultiSensitiveDetector

#include "G4MultiSensitiveDetector.hh"
#include "G4SDManager.hh"
#include <sstream>

//#define MSDDEBUG
#ifdef MSDDEBUG
#define DBG( msg ) G4cout<<msg<<G4endl
#else
#define DBG( msg )
#endif
#define VDBG( vl , msg ) if ( vl<=verboseLevel ) G4cout<<msg<<G4endl

G4MultiSensitiveDetector::G4MultiSensitiveDetector(G4String name)
	: G4VSensitiveDetector(name)
{
#ifdef MSDDEBUG
	verboseLevel = 3;
#endif
	VDBG(1,"Creating G4MultiSenstiveDetector with name: "<<name);
}

G4MultiSensitiveDetector::~G4MultiSensitiveDetector()
{
	VDBG(2,GetName()<<" : Destructing G4MultiSensitiveDetector");
	ClearSDs();
}

G4MultiSensitiveDetector::G4MultiSensitiveDetector(const G4MultiSensitiveDetector& rhs)
: G4VSensitiveDetector(rhs) ,
  fSensitiveDetectors(rhs.fSensitiveDetectors)
{
	VDBG(3,GetName()<<" : Copy constructor called.");
}

G4MultiSensitiveDetector&
G4MultiSensitiveDetector::operator=(const G4MultiSensitiveDetector& rhs)
{
	if ( this != &rhs ) {
		//G4VSensitiveDetector::operator=(static_cast<const G4VSensitiveDetector&>(rhs));
		G4VSensitiveDetector::operator=(static_cast<const G4VSensitiveDetector&>(rhs));
		fSensitiveDetectors = rhs.fSensitiveDetectors;
	}
	return *this;
}

void
G4MultiSensitiveDetector::Initialize(G4HCofThisEvent* )
{
    //SDManager is resposnsible for calling this since the granular SDs
  // are also registered
	//for ( auto sd : fSensitiveDetectors ) sd->Initialize(hcte);
}

void
G4MultiSensitiveDetector::EndOfEvent(G4HCofThisEvent* )
{
  //SDManager is resposnsible for calling this since the granular SDs
// are also registered
	//for ( auto sd : fSensitiveDetectors ) sd->EndOfEvent(hcte);
}

void
G4MultiSensitiveDetector::clear()
{
	for ( auto sd : fSensitiveDetectors ) sd->clear();
}

void
G4MultiSensitiveDetector::DrawAll()
{
	for ( auto sd : fSensitiveDetectors ) sd->DrawAll();
}

void
G4MultiSensitiveDetector::PrintAll()
{
	for ( auto sd : fSensitiveDetectors ) sd->PrintAll();
}

G4bool
G4MultiSensitiveDetector::ProcessHits(G4Step*aStep,G4TouchableHistory*)
{
	VDBG(2,GetName()<<" : Called processHits: "<<aStep<<" with Edep: "<<aStep->GetTotalEnergyDeposit());
	G4bool result = true;
	for (auto sd : fSensitiveDetectors )
		result &= sd->Hit(aStep);
	return result;
}

G4int G4MultiSensitiveDetector::GetCollectionID(G4int)
{
	G4ExceptionDescription msg;
	msg << GetName()<<" : This method cannot be called for an instance of type G4MultiSensitiveDetector."
		<< " First retrieve a contained G4VSensitiveDetector with. i.e. GetSD and then "
		<< " call this method.";
	G4Exception("G4MultiSensitiveDetector::GetCollectionID","Det0011",FatalException,msg);
	return -1;
}

//This method requires all contained SD to be clonable
G4VSensitiveDetector* G4MultiSensitiveDetector::Clone() const
{
	VDBG(2,GetName()<<"Cloning an instance of G4MultiSensitiveDetector");
	G4MultiSensitiveDetector* newInst = new G4MultiSensitiveDetector(this->GetName());
	for ( auto sd : fSensitiveDetectors )
		newInst->AddSD( sd->Clone() );
	return newInst;
}
