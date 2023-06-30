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
// Gathers all dependencies to FLUKA INCLUDE FILES,
// as needed by the G4 <-> FLUKA interface.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef FLUKA_COMMON_DEPENDENCIES_HH
#define FLUKA_COMMON_DEPENDENCIES_HH


#include "dblprc.h"
#include "dimpar.h"
#include "iounit.h"
#include "beamcm.h"
#include "blnkcm.h"
#include "caslim.h"
#include "cmelds.h"
#include "cmphnu.h"
#include "ctitle.h"
#include "currpt.h"
#include "evaflg.h"
#include "evapix.h"
#include "fheavy.h"
#include "flkmat.h"
#include "genflg.h"
#include "genstk.h"
#include "genthr.h"
#include "isotop.h"
#include "ncsfta.h"
#include "ndnicm.h"
#include "nucdat.h"
#include "nucflg.h"
#include "nucgeo.h"
#include "nucpot.h"
#include "nuinfo.h"
#include "paprop.h"
#include "parevt.h"
#include "part2.h"
#include "phnccm.h"
#include "resnuc.h"
#include "sgtbcm.h"
#include "sumcou.h"
#include "thrscm.h"
#include "usryld.h"


#endif
#endif // G4_USE_FLUKA
