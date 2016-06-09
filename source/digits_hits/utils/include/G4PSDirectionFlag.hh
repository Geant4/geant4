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
// $Id: G4PSDirectionFlag.hh,v 1.1 2005/11/17 22:53:38 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
//---------------------------------------------------------------
//
// G4PSDirectionFlag.hh
//
// Class Description:
//   This is an enumerator to define the direction of
// the particle's Surface Current or Flux for PrimitiveScorer
//
//
//---------------------------------------------------------------

#ifndef G4PSDirectionFlag_h
#define G4PSDirectionFlag_h 1

//////////////////
enum G4PSFluxFlag
//////////////////
{ 
    fFlux_InOut, 
    // For both direction In / Out.
    fFlux_In, 
    // IN : Direction which comes into the geometry
    fFlux_Out 
};

/////////////////////
enum G4PSCurrentFlag
/////////////////////
{ 
    fCurrent_InOut, 
    // For both direction In / Out.
    fCurrent_In, 
    // IN : Direction which comes into the geometry
    fCurrent_Out 
    // OUT : Direction which goes out fromthe geometry
};

#endif


