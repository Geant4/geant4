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
//
// Modified version of the G4Nucleus class by T. Lampen, HIP, June 2000
 // original G4Nucleusby H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 19-Nov-1996
 // last modified: 27-Mar-1997
 // J.P.Wellisch: 23-Apr-97: minor simplifications
 // modified by J.L.Chuma 24-Jul-97  to set the total momentum in Cinema and
 //                                  EvaporationEffects
 // modified by J.L.Chuma 21-Oct-97  put std::abs() around the totalE^2-mass^2
 //                                  in calculation of total momentum in
 //                                  Cinema and EvaporationEffects
 // Chr. Volcker, 10-Nov-1997: new methods and class variables.
 // HPW added utilities for low energy neutron transport. (12.04.1998)
 // M.G. Pia, 2 Oct 1998: modified GetFermiMomentum to avoid memory leaks
 
#include "G4LayeredNucleus.hh"
#include "Randomize.hh"

G4ThreeVector G4LayeredNucleus::GetMomentum()
{
  return  momentumVector;
}


void G4LayeredNucleus::SetMomentum( const G4ThreeVector& mom )
  {
    momentumVector = mom;
  }



/* end of file */

