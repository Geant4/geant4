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
// $Id: G4SurfaceOfRevolution.hh,v 1.4 2001-07-11 09:59:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SurfaceOfRevolution
//
// Class description:
// 
// Definition of a surface of revolution.
// (not implemented yet).

// ----------------------------------------------------------------------
#ifndef included_G4SurfaceOfRevolution
#define included_G4SurfaceOfRevolution

#include "G4Surface.hh"

class G4SurfaceOfRevolution : public G4Surface
{

public:  // with description

  G4SurfaceOfRevolution();
  virtual ~G4SurfaceOfRevolution();
    // Constructor & destructor.

private:

  G4SurfaceOfRevolution(const G4SurfaceOfRevolution &);
  G4SurfaceOfRevolution& operator=(const G4SurfaceOfRevolution &);
    // Private copy constructor and assignment operator.
};

#endif
