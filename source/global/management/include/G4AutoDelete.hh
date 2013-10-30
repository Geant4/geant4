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
// $Id$
//
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//   This function implements a simplified "garbage collection" mechanism
//   for G4-MT model. Objects are registered when created on the heap and they 
//   will be deleted at the end of
//   the program (like they would be if marked as "static").
// 
// Limitation:
//   The registered object, should not
//   contain any G4ThreadLocal data member. Note that in general,
//   if object is to be thread-private it is unnecessary to mark
//   any data-member as G4ThreadLocal.
//
// Performance issues:
//   This function uses G4ThreadLocalSingleton that on its own uses
//   locks and mutexes. Thus its use should be limited to only when 
//   really necessary.
//
// Example:
//   class G4SharedByThreads {
//     void calledByThreads() {
//          G4Something* anObject = new G4Something;
//          G4AutoDelete::Register( anObject );
//     }
//   };
//
// History:
//  28 October 2013: A. Dotti - First implementation

#ifndef G4AUTODELETE_HH
#define G4AUTODELETE_HH

#include "G4ThreadLocalSingleton.hh"
namespace G4AutoDelete {
  template<class T>
  void Register( T* inst ) {
    static G4ThreadLocalSingleton<T> container;
    container.Register(inst);
  }
}

#endif //G4AUTODELETE_HH
