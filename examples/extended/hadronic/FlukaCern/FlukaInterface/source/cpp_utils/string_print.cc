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
#ifdef G4_USE_FLUKA


#include "string_print.h"

#include <stdlib.h>
#include <stdarg.h>
#include <cstring>
#include <memory>


namespace cpp_utils {

  // ***************************************************************************
  // String format like printf() for strings.
  // ***************************************************************************
  std::string sformat(const char* fmt, ...) {
    std::unique_ptr<char[]> formatted;
    int n = (int)strlen(fmt) * 2;	// Reserve two times as much as the length of the fmt_str
    while (true) {
      // Wrap the plain char array into the unique_ptr
      formatted.reset(new char[n]);
      strcpy(&formatted[0], fmt);
      va_list ap;
      va_start(ap, fmt);
      int final_n = vsnprintf(&formatted[0], n, fmt, ap);
      va_end(ap);
      if (final_n < 0 || final_n >= n) {
        n += std::abs(final_n - n + 1);
      }
      else {
        break;
      }
    }
    return std::string(formatted.get());
  }

} // namespace cpp_utils


#endif // G4_USE_FLUKA
