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
// $Id: G4ProcessType.hh,v 1.4 2001-07-11 10:08:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4ProcessType.hh
//
// Class Description:
//   This is an enumerator to define process type
//
//
//---------------------------------------------------------------

#ifndef G4ProcessType_h
#define G4ProcessType_h 1

enum G4ProcessType
{
  fNotDefined,
  fTransportation,
  fElectromagnetic,
  fOptical,             
  fHadronic,
  fPhotolepton_hadron,
  fDecay,
  fGeneral,
  fParameterisation,
  fUserDefined
};
#endif


