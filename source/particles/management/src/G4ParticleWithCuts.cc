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
//
// $Id: G4ParticleWithCuts.cc,v 1.19 2003-03-10 08:43:53 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//  History: 
//   first implementation, based on object model of Hisaya Kurashige,
//                                                  21 Oct 1996
//   calculation of Range Table is based on implementeation for Muon 
//                                         by L.Urban, 10 May 1996
//   modify CalcEnergyCuts                 09 Nov. 1998, L.Urban
//   added  RestoreCuts  H.Kurashige 09 Mar. 2001
//   modify for material-V03-02-02 (STL migration)  H.Kurashige 19 Sep. 2001
//   introduced material dependent range cuts   08 Oct. 2001
//   restructuring for Cuts per Region  by Hisaya    11 MAr.2003 
// ----------------------------------------------------------------
// Class Description
//  Dummy to be removed in future 
// ------------------------------------------------------------
#include "globals.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"

#include "G4ios.hh"

#include "g4std/strstream"

G4ParticleWithCuts::G4ParticleWithCuts(
		const G4String&  aName,  
                G4double         mass,     
                G4double         width,
                G4double         charge,   
                G4int            iSpin,
                G4int            iParity,
                G4int            iConjugation,
                G4int            iIsospin,   
                G4int            iIsospinZ, 
                G4int            gParity,
                const G4String&  pType,
                G4int            lepton,
                G4int            baryon,
                G4int            encoding,
                G4bool           stable,
                G4double         lifetime,
                G4DecayTable     *decaytable,
		G4bool           shortlived)
	: G4ParticleDefinition(aName, mass, width, charge, iSpin, iParity,
                               iConjugation, iIsospin, iIsospinZ, gParity,
                               pType, lepton, baryon, encoding, stable,
                               lifetime, decaytable, shortlived)
{
}

G4ParticleWithCuts::~G4ParticleWithCuts()
{ 
}
