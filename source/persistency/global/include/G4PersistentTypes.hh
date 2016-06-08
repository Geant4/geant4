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
// $Id: G4PersistentTypes.hh,v 1.4 2001/07/11 10:02:25 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

// file description:
//   Persistent-capable typedefs for Geant4/Persistency category

// History:
// 15.06.98 Y.Morita - Created

#ifndef G4_PERSISTENT_TYPES_HH
#define G4_PERSISTENT_TYPES_HH

// Typedefs to decouple from library classes
#include <HepODBMS/odbms/HepODBMS.h>

typedef d_String G4PString;
typedef d_Double G4Pdouble;
typedef d_Float  G4Pfloat;
typedef d_Long   G4Pint;
typedef d_Long   G4Plong;

#endif /* G4_PERSISTENT_TYPES_HH */
