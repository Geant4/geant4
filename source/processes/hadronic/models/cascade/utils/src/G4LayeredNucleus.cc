// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Modified version of the G4Nucleus class by T. Lampen, HIP, June 2000
 // original G4Nucleusby H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 19-Nov-1996
 // last modified: 27-Mar-1997
 // J.P.Wellisch: 23-Apr-97: minor simplifications
 // modified by J.L.Chuma 24-Jul-97  to set the total momentum in Cinema and
 //                                  EvaporationEffects
 // modified by J.L.Chuma 21-Oct-97  put abs() around the totalE^2-mass^2
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


void G4LayeredNucleus::SetMomentum( const G4ThreeVector mom )
  {
    momentumVector = mom;
  }



/* end of file */

