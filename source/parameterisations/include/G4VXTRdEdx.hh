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
// $Id: G4VXTRdEdx.hh,v 1.5 2001-09-18 09:02:00 gcosmo Exp $
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


#ifndef G4VXTRdEdx_h
#define G4VXTRdEdx_h 1


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

  G4double AngleXTRdEdx(G4double varAngle) ;
};

#endif

