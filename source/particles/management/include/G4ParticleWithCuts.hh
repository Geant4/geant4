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
// $Id: G4ParticleWithCuts.hh,v 1.16 2003-03-10 08:43:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//       first implementation, based on object model of Hisaya Kurashige, 
//                                                      21 Oct 1996
//       calculation of Range Table is based on implementeation for Muon 
//                                           by L.Urban, 10 May 1996
//       added  RestoreCuts  H.Kurashige 09 Mar. 2001
//       introduced material dependent range cuts   08 Sep. 2001
//       
//       restructuring for Cuts per Region  by Hisaya    07 Oct.2002 
//       restructuring for Cuts per Region  by Hisaya    11 MAr.2003 
// ----------------------------------------------------------------
// Class Description
//  Dummy to be removed in future 
//

#ifndef G4ParticleWithCuts_h
#define G4ParticleWithCuts_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"

class G4PhysicsLogVector;
class G4ProductionCutsTable;

class G4ParticleWithCuts : public G4ParticleDefinition
{
  public:
     G4ParticleWithCuts(const G4String&  aName,  
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
		G4bool           resonance = false);
      virtual ~G4ParticleWithCuts();
   
  //--------------for SetCuts-------------------------------------------

};

#endif
