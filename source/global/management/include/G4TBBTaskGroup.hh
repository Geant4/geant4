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
//
// This file wraps a TBB task_group into a G4TaskGroup
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (May 28st 2020)
// ---------------------------------------------------------------

#pragma once

#ifndef G4GMAKE
#  include "G4GlobalConfig.hh"
#endif

#if defined(GEANT4_USE_TBB)
#  include "PTL/TBBTaskGroup.hh"
#else
#  include "PTL/TaskGroup.hh"
#endif

#if defined(GEANT4_USE_TBB)
template <typename Tp, typename Arg = Tp>
using G4TBBTaskGroup = PTL::TBBTaskGroup<Tp, Arg>;
#else
template <typename Tp, typename Arg = Tp>
using G4TBBTaskGroup = PTL::TaskGroup<Tp, Arg>;
#endif
