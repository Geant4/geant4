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
// $Id: G4NURBStubesector.hh,v 1.8 2006-06-29 19:05:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tubesector builder prototype
// OC 290896

// Class Description:
// Tubesector builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.

#ifndef __C_G4NURBStubesector__
#define __C_G4NURBStubesector__ 1 

#include "G4NURBS.hh"

class  G4NURBStubesector : public G4NURBS
{
  public:  
    // angle in radians

    // If PHI2 smaller (or equal) than PHI1 , it is incremented
    // by 2pi as necessary to become strictly greater.

    // Except that, you can use any value for the arguments,
    // it's the renderer or you that will have troubles. 

    G4NURBStubesector(G4double RMIN, G4double RMAX,
                      G4double DZ, G4double PHI1, G4double PHI2);
    virtual ~G4NURBStubesector();
    virtual const char* Whoami() const;

  protected:
    char *  mpwhoami;

  private:
    static t_inddCtrlPt DecideNbrCtrlPts(G4double PHI1, G4double PHI2);
};

#endif
// end of __C_G4NURBStubesector__
