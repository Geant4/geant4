// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VXTRdEdx.hh,v 1.2 2001-02-27 09:34:02 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// base class for 'fast' parametrisation model describing X-ray transition
// created in some G4Envelope. Anglur distribuiton is very rough !!! (see DoIt
// method
// 
// History:
// 26.02.01 V. Grichine first version 
// 26.02.01 V. Grichine, DoIt was transformed from virtual
//


#ifndef G4VXrayTRadModel_h
#define G4VXrayTRadModel_h 1


#include "globals.hh"
#include "templates.hh"
#include "g4std/complex"

#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Gamma.hh"

#include "G4VXrayTRmodel.hh"


class G4VXTRdEdx : public G4VXrayTRmodel
{
public:

   G4VXTRdEdx (G4LogicalVolume *anEnvelope,G4double,G4double);
   virtual  ~G4VXTRdEdx ();

  // Pure virtual functions from base class
 
  void DoIt(const G4FastTrack&, G4FastStep&)  ;

  // Pure virtuals must be implemented in inherited particular TR radiators

  virtual  G4double GetStackFactor( G4double energy, G4double gamma,
                                                     G4double varAngle ) = 0  ;

protected:

  void BuildTable() ;
  void BuildEnergyTable() ;
  void BuildAngleTable() ;

  G4complex OneInterfaceXTRdEdx( G4double energy, 
                                G4double gamma,
                                G4double varAngle ) ;

  G4double SpectralAngleXTRdEdx(G4double varAngle) ;

  G4double SpectralXTRdEdx(G4double energy) ;

  G4double AngleSpectralXTRdEdx(G4double energy) ;

  G4double G4VXTRdEdx::AngleXTRdEdx(G4double varAngle) ;
};

#endif

