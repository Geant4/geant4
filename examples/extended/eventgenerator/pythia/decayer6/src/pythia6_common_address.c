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
// $Id: pythia6_common_address.c 100687 2016-10-31 11:20:33Z gcosmo $
//
// According to pythia6_common_address.c provided in Root
// Pythia6 distribution:
// http://root.cern.ch/
// see http://root.cern.ch/root/License.html
// ----------------------------------------------------------------------------

#include <string.h>
// declaration of PYTHIA6 common clocks
#ifndef WIN32
# define pycomp pycomp_
# define py1ent py1ent_
# define pyjets pyjets_
# define pydat1 pydat1_
# define pydat3 pydat3_
# define type_of_call
#else
# define pycomp PYCOMP
# define py1ent PY1ENT
# define pyjets PYJETS
# define pydat1 PYDAT1
# define pydat3 PYDAT3
# define type_of_call _stdcall
#endif

extern int pyjets[2+5*4000+2*2*5*4000];
extern int pydat1[200+2*200+200+2*200];
extern int pydat3[3*500+2*8000+2*8000+5*8000];  /* KNDCAY=8000 */

void* pythia6_common_address(const char* name)
{
   if      (!strcmp(name,"PYJETS")) return pyjets;
   else if (!strcmp(name,"PYDAT1")) return pydat1;
   else if (!strcmp(name,"PYDAT3")) return pydat3;
   return 0;
}
