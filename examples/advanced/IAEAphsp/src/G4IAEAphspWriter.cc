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
// Author: M.A. Cortes-Giraldo, Universidad de Sevilla
//
// History changelog prior creation of this example:
// - 17/10/2009: inheritance removed (not needed)
// - 17/10/2009: version 1.0
// - 07/10/2024: version 2.0 for Geant4 example (MT compliant)
// - 18/10/2025: version 3.0
//   - Creation of IAEASourceIdRegistry for thread-safe source_id assignation
//


#include "G4IAEAphspWriter.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"

#include "iaea_phsp.h"
#include "IAEASourceIdRegistry.hh"
#include "G4IAEAphspWriterStack.hh"

#include <map>
#include <sstream>
#include <vector>

//==============================================================================

G4IAEAphspWriter::G4IAEAphspWriter(const G4String filename)
{
  fFileName = filename;
  fZphspVec = new std::vector<G4double>;
  fConstVariables = new std::map<G4int, G4double>;
  fOrigHistories = new std::vector<G4int>;
  fIAEASourcesOpen = 0;
  G4cout << "G4IAEAphspWriter object constructed." << G4endl;
}


G4IAEAphspWriter::~G4IAEAphspWriter()
{
  if (fZphspVec) delete fZphspVec;
  if (fConstVariables) delete fConstVariables;
  if (fOrigHistories) delete fOrigHistories;
}


//==============================================================================

void G4IAEAphspWriter::AddZphsp(const G4double zphsp)
{
  fZphspVec->push_back(zphsp);
  G4cout << "G4IAEAphspWriter: Registered phase-space plane at z = "
	 << zphsp/cm << " cm." << G4endl;

  // fOrigHistories must incorporate a new element (equal size as fZphspVec)
  fOrigHistories->push_back(0);
}


//==============================================================================

void G4IAEAphspWriter::SetDataFromIAEAStack(const G4IAEAphspWriterStack* stack)
{
  if (!stack) {
    G4ExceptionDescription msg;
    msg <<  "No G4IAEAphspWriterStack has been constructed!" << G4endl;
    G4Exception("G4IAEAphspWriter::SetDataFromIAEAStack()",
		"IAEAphspWriter001", FatalException, msg );
    return;
  }
  else {
    fFileName = stack->GetFileName();

    if (stack->GetZphspVec()->size() > 0) {
      (*fZphspVec) = *(stack->GetZphspVec()); // copy objects, not pointers

      G4cout << "G4IAEAphspWriter::fFileName = " << fFileName << G4endl;
      G4cout << "G4IAEAphspWriter::fZphspVec->size() = "
	     << fZphspVec->size() << G4endl;

      // fOrigHistories must have as many element as fZphspVec
      fOrigHistories->assign(fZphspVec->size(), 0);
    }
    else {
      G4ExceptionDescription msg;
      msg << "No phsp plane z-coordinate has been defined!" << G4endl;
      G4Exception("G4IAEAphspWriter::SetDataFromIAEAStack()",
		  "IAEAphspWriter002", FatalErrorInArgument, msg );
      return;
    }
  }
}


//==============================================================================

void G4IAEAphspWriter::SetConstVariable(const G4int idx, const G4double val)
{
  if (idx < 0 || idx > 6) {
    G4cout << "No constant variable applies for index " << idx
	   << ". Doing nothing." << G4endl;
    return;
  }

  fConstVariables->insert( std::pair<G4int, G4double>(idx, val) );
  switch (idx) {
  case 0:
    G4cout << "Variable 'x' set to constant value " << val << " cm" << G4endl;
    break;
  case 1:
    G4cout << "Variable 'y' set to constant value " << val << " cm" << G4endl;
    break;
  case 2:
    G4cout << "Variable 'z' set to constant value " << val << " cm" << G4endl;
    break;
  case 3:
    G4cout << "Variable 'u' set to constant value " << val << G4endl;
    break;
  case 4:
    G4cout << "Variable 'v' set to constant value " << val << G4endl;
    break;
  case 5:
    G4cout << "Variable 'w' set to constant value " << val << G4endl;
    break;
  case 6:
    G4cout << "Variable 'wt' set to constant value " << val << G4endl;
  }
}


//==============================================================================

void G4IAEAphspWriter::WriteIAEAParticle(const size_t idx, const G4int incHist,
					 const G4int pdg, const G4double kinE,
					 const G4double wt,
					 const G4ThreeVector pos,
					 const G4ThreeVector momDir)
{
  IAEA_I32 partType;
  switch(pdg) {
  case 22:
    partType = 1;  // gamma
    break;
  case 11:
    partType = 2;  // electron
    break;
  case -11:
    partType = 3;  // positron
    break;
  case 2112:
    partType = 4;  // neutron
    break;
  case 2212:
    partType = 5;  // proton
    break;
  default:
    G4ExceptionDescription msg;
    msg << "PDG " << pdg
	<< " is not supported by IAEAphsp format and will not be recorded."
	<< G4endl;
    G4Exception("G4IAEAphspWriter::WriteIAEAParticle()",
		"IAEAphspWriter003", JustWarning, msg);
    return;
  }

  IAEA_I32 sourceID = static_cast<IAEA_I32>(idx+fIAEASourcesOpen);

  IAEA_I32 nStat = static_cast<IAEA_I32>(incHist);
  IAEA_Float energy = static_cast<IAEA_Float>(kinE/MeV);
  IAEA_Float weight = static_cast<IAEA_Float>(wt);
  IAEA_Float x = static_cast<IAEA_Float>( pos.x()/cm );
  IAEA_Float y = static_cast<IAEA_Float>( pos.y()/cm );
  IAEA_Float z = static_cast<IAEA_Float>( pos.z()/cm );
  IAEA_Float u = static_cast<IAEA_Float>( momDir.x() );
  IAEA_Float v = static_cast<IAEA_Float>( momDir.y() );
  IAEA_Float w = static_cast<IAEA_Float>( momDir.z() );

  // Extra variables
  IAEA_Float extraFloat = -1; // no extra floats stored
  IAEA_I32 extraInt = nStat;

  // And finally store the particle following the IAEA routines
  iaea_write_particle(&sourceID, &nStat, &partType,
		      &energy, &weight,
		      &x, &y, &z, &u, &v, &w, &extraFloat, &extraInt);
}



//==============================================================================

void G4IAEAphspWriter::WriteIAEAParticle(const G4Step* aStep,
					 const G4int zStopIdx)
{
  IAEA_I32 sourceID =
    static_cast<IAEA_I32>(zStopIdx+fIAEASourcesOpen); // beware
  G4double zStop = (*fZphspVec)[zStopIdx];

  // The particle type and kinetic energy
  // -------------------------------------
  const G4Track* aTrack = aStep->GetTrack();
  G4int PDGCode = aTrack->GetDefinition()->GetPDGEncoding();
  IAEA_I32 partType;
  G4double postE = aStep->GetPostStepPoint()->GetKineticEnergy();
  G4double preE = aStep->GetPreStepPoint()->GetKineticEnergy();
  IAEA_Float kinEnergyMeV;
  G4ThreeVector postR = aStep->GetPostStepPoint()->GetPosition();
  G4ThreeVector preR = aStep->GetPreStepPoint()->GetPosition();
  G4double postZ = postR.z();
  G4double preZ = preR.z();

  switch(PDGCode) {
  case 22:
    partType = 1;  // gamma
    kinEnergyMeV = static_cast<IAEA_Float>(preE/MeV);
    break;
  case 11:
    partType = 2;  // electron
    kinEnergyMeV =
      static_cast<IAEA_Float>( (preE+
				(postE-preE)*(zStop-preZ)/(postZ-preZ))/MeV );
    break;
  case -11:
    partType = 3;  // positron
    kinEnergyMeV =
      static_cast<IAEA_Float>( (preE+
				(postE-preE)*(zStop-preZ)/(postZ-preZ))/MeV );
    break;
  case 2112:
    partType = 4;  // neutron
    kinEnergyMeV = static_cast<IAEA_Float>(preE/MeV);
    break;
  case 2212:
    partType = 5;  // proton
    kinEnergyMeV =
      static_cast<IAEA_Float>( (preE+
				(postE-preE)*(zStop-preZ)/(postZ-preZ))/MeV );
    break;
  default:
    G4String pname = aTrack->GetDefinition()->GetParticleName();
    G4String errmsg = "'" + pname + "' is not supported by IAEAphsp format"
      + " and will not be recorded.";
    G4Exception("G4IAEAphspWriter::WriteIAEAParticle()",
		"IAEAphspWriter004", JustWarning, errmsg.c_str() );
    return;
  }

  // Track weight
  IAEA_Float wt = static_cast<IAEA_Float>(aTrack->GetWeight());

  // Position
  G4double postX = postR.x();
  G4double preX = preR.x();
  G4double postY = postR.y();
  G4double preY = preR.y();
  IAEA_Float x =
    static_cast<IAEA_Float>( (preX+
			      (postX-preX)*(zStop-preZ)/(postZ-preZ))/cm );
  IAEA_Float y =
    static_cast<IAEA_Float>( (preY+
			      (postY-preY)*(zStop-preZ)/(postZ-preZ))/cm );
  IAEA_Float z =
    static_cast<IAEA_Float>( zStop/cm );

  // Momentum direction
  G4ThreeVector momDir = aStep->GetPreStepPoint()->GetMomentumDirection();
  IAEA_Float u = static_cast<IAEA_Float>(momDir.x());
  IAEA_Float v = static_cast<IAEA_Float>(momDir.y());
  IAEA_Float w = static_cast<IAEA_Float>(momDir.z());

  // Extra variables
  IAEA_Float extraFloat = -1; // no extra floats stored
  // IAEA_I32 extraInt = static_cast<IAEA_I32>((*fIncrNumberVector)[zStopIdx]);
  // IAEA_I32 nStat = extraInt;
  IAEA_I32 extraInt = -1;
  IAEA_I32 nStat = 0;  //MACG Admitted values, this MUST change

  // And finally store the particle following the IAEA routines
  iaea_write_particle(&sourceID, &nStat, &partType,
		      &kinEnergyMeV, &wt,
		      &x, &y, &z, &u, &v, &w, &extraFloat, &extraInt);
}



//==============================================================================

void G4IAEAphspWriter::OpenIAEAphspOutFiles(const G4Run* aRun)
{
  // Open all the files intended to store
  // the phase spaces following the IAEA format.

  const IAEA_I32 accessWrite = 2;  // 2 = Writing mode in IAEA routines

  size_t nZphsps = fZphspVec->size();
  for (size_t ii = 0; ii < nZphsps; ii++) {
    // Set the source ID and file name in a unique way
    std::stringstream sstr;
    sstr << ((*fZphspVec)[ii]/cm);
    G4String zphsp(sstr.str());
    G4String fullName = fFileName + "_" + zphsp + "cm";

    // This part is only added when running several runs
    // during the simulation.
    G4int runID = aRun->GetRunID();
    if (runID > 0) {
      std::stringstream sstr2;
      sstr2 << runID;
      G4String runIDStr(sstr2.str());
      fullName += "_";
      fullName += runIDStr;
    }

    
    // Create the file to store the IAEA phase space

    // Reserve a global ID and request it explicitly
    G4int reserved = IAEASourceIdRegistry::Instance().ReserveNextLowest();
    if (reserved < 0) {
      G4ExceptionDescription ed;
      ed << "No free IAEA source IDs available for writer" << G4endl;
      G4Exception("G4IAEAphspWriter::OpenIAEAphspOutFiles",
		  "IAEAphspWriter005",FatalException, ed);
    }
    IAEA_I32 sourceWrite = static_cast<IAEA_I32>(reserved);
    char* filename = const_cast<char*>(fullName.data());
    IAEA_I32 result = 0;
    
    iaea_new_source( &sourceWrite, filename, &accessWrite,
		     &result, fullName.size()+1 );

    if (result < 0 || sourceWrite < 0) {
      IAEASourceIdRegistry::Instance().Release(reserved);
      G4ExceptionDescription ed;
      ed << "IAEAphsp output file opening operation failed!" << G4endl;
      G4Exception("G4IAEAphspWriter::OpenIAEAphspOutFiles()",
		  "IAEAphspWriter006", FatalException, ed);
    }
    
    // This difference tells about the number of IAEA files already open.
    fIAEASourcesOpen = sourceWrite - ii;
    G4cout << "G4IAEAphspWriter::OpenOutputIAEAphspFiles() ==> "
	   << "\"" << fullName << "\"   IAEAphsp id = " << sourceWrite << "."
	   << G4endl;


    // Set the global information and options.

    // Set constant variables
    std::map<G4int, G4double>::iterator itmap;
    for (itmap = fConstVariables->begin();
	 itmap != fConstVariables->end(); itmap++) {
      IAEA_I32 varIdx = static_cast<IAEA_I32>( (*itmap).first );
      IAEA_Float varValue = static_cast<IAEA_I32>( (*itmap).second );
      iaea_set_constant_variable(&sourceWrite, &varIdx, &varValue);
    }
    //MACG
    // Set constant Z
    // IAEA_I32 varIdx = 2;  // 0=x, 1=y, 2=z, 3=u, 4=v, 5=w, 6=wt
    // IAEA_Float varValue = static_cast<IAEA_Float>((*fZphspVec)[ii]/cm);
    // iaea_set_constant_variable(&sourceWrite, &varIdx, &varValue);

    // Extra variables
    IAEA_I32 extraFloats = 0;
    IAEA_I32 extraInts = 1;
    iaea_set_extra_numbers(&sourceWrite, &extraFloats, &extraInts);

    // Extra variables types
    IAEA_I32 longIdx = 0;
    IAEA_I32 longType = 1; // incremental history number
    iaea_set_type_extralong_variable(&sourceWrite, &longIdx, &longType);
  }
}


//==============================================================================

void G4IAEAphspWriter::CloseIAEAphspOutFiles()
{
  // Close the IAEA files
  G4int nZphsps = fZphspVec->size();
  for (G4int ii = 0; ii < nZphsps; ii++) {
    const IAEA_I32 sourceID = static_cast<IAEA_I32>(ii+fIAEASourcesOpen);
    IAEA_I64 nEvts = static_cast<IAEA_I64>( fOrigHistories->at(ii) );
    iaea_set_total_original_particles(&sourceID, &nEvts);

    IAEA_I32 result = 0;
    iaea_print_header(&sourceID, &result);
    if (result < 0) {
      G4Exception("G4IAEAphspWriter::EndOfRunAction()",
		  "IAEAphspWriter007", JustWarning,
		  "IAEA phsp source not found");
    }

    iaea_destroy_source(&sourceID, &result);
    if (result > 0) {
      IAEASourceIdRegistry::Instance().Release(static_cast<G4int>(sourceID));
      G4cout << "Phase-space file at z_phsp = " << (*fZphspVec)[ii]/cm
	     << " cm (IAEA source id #" << sourceID << ") closed successfully!"
	     << G4endl << G4endl;
    }
    else {
      G4Exception("G4IAEAphspWriter::EndOfRunAction()",
		  "IAEAphspWriter008", JustWarning,
		  "IAEA file not closed properly");
    }
  }
}
