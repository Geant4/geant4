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
// DUMMY CLASS
#ifndef AIDA_ITUPLEENTRY_H
#define AIDA_ITUPLEENTRY_H 1

#include <string>
#include <vector>

namespace AIDA {

/**
 * User level interface to a Tuple.
 *
 * @author The AIDA team (http://aida.freehep.org/)
 *
 */

class ITupleEntry {

public:
    /// Destructor.
    virtual ~ITupleEntry() { /* nop */; }

}; // class
}; // namespace AIDA
#endif /* ifndef AIDA_ITUPLEENTRY_H */
