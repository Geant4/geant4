// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IrregularXTRdEdx.cc,v 1.1 2001-02-27 15:23:48 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//



#include "G4IrregularXTRdEdx.hh"
#include "Randomize.hh"

#include "G4Gamma.hh"

////////////////////////////////////////////////////////////////////////////
//
// Constructor, destructor

G4IrregularXTRdEdx::G4IrregularXTRdEdx(G4Envelope *anEnvelope, 
                                                  G4double a, G4double b) :
  G4VXTRdEdx(anEnvelope,a,b)
{
  G4cout<<"Irregular X-ray TR dE/dx model is called"<<G4endl ;

  // Build energy and angular integral spectra of X-ray TR photons inside
  // a radiator

    BuildTable() ;
}

///////////////////////////////////////////////////////////////////////////

G4IrregularXTRdEdx::~G4IrregularXTRdEdx()
{
  ;
}




///////////////////////////////////////////////////////////////////////////
//
// Very rough approximation for radiator interference factor for the case of
// fully irregular radiator. The plate and gas gap thicknesses are distributed 
// according to exponent. The mean values of the plate and gas gap thicknesses 
// are supposed to be much more than XTR formation zones but much less than 
// mean absorption length of XTR photons in coresponding material.

G4double 
G4IrregularXTRdEdx::GetStackFactor( G4double energy, 
                                         G4double gamma, G4double varAngle )
{
  G4double result ;

  G4complex R  = OneInterfaceXTRdEdx(energy,gamma,varAngle) ;
  R           *= G4double(fPlateNumber) ;
  result       = 2.0*G4std::real(R) ;

  return      result ;
}


//
//
////////////////////////////////////////////////////////////////////////////








