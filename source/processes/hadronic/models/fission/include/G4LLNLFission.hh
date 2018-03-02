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
// This software was developed by Lawrence Livermore National Laboratory.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 3. The name of the author may not be used to endorse or promote products
//   derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Copyright (c) 2006 The Regents of the University of California.
// All rights reserved.
// UCRL-CODE-224807
//
// $Id: G4LLNLFission.hh 67966 2013-03-13 09:38:38Z gcosmo $
//

#include "G4Types.hh"

//
// This file is a copy of Fission.hh, made for use with Geant4.
//

  extern void genspfissevt_(G4int *isotope, G4double *time);
/*
 * This function is called to trigger a spontaneous fission.
 * Multiple neutrons and photons are generated and stored
 * in a stack along with their energies, directions and 
 * emission times.
 * The arguments of this function are
 *      isotope:        94239 for Pu-239 for instance
 *      time:           the time of the spontaneous fission
 */

  extern void genfissevt_(G4int *isotope, G4double *time, G4double *nubar, G4double *eng);
/*
 * This function is called to trigger a neutron-induced fission.
 * Multiple neutrons and photons are generated and stored
 * in a stack along with their energies, directions and 
 * emission times. In addition to the arguments above, this
 * function needs
 *      nubar:          user-specified average number of neutrons emitted 
 *                      per fission (e.g. as tabulated in the cross-section 
 *                      libraries used by the particle transport code)
 *      eng:            energy of the neutron inducing fission
 */

  extern G4int getnnu_();
/*
 * This function returns the number of neutrons emitted by the
 * fission, -1 if there is no neutron data for that isotope in 
 * the fission library.
 */

  extern G4int getpnu_();
/*
 * This function returns the number of photons emitted by the
 * fission, -1 if there is no photon data for that isotope in 
 * the fission library.
 */

  extern G4double getneng_(G4int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * its energy, -1 if index isout of range.
 */

  extern G4double getnvel_(G4int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * the amplitude of its velocity, -1 if index is out of range.
 */

  extern G4double getndircosu_(G4int *index);
  extern G4double getndircosv_(G4int *index);
  extern G4double getndircosw_(G4int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * the direction cosines of its velocity vector on the x, y and z 
 * axes.
 */

  extern G4double getpeng_(G4int *index);
/*
 * Given the index of the emitted photon, this function returns
 * its energy, -1 if index is out of range.
 */

  extern G4double getpvel_(G4int *index);
/*
 * Given the index of the emitted photon, this function returns
 * the amplitude of its velocity, -1 if index is out of range.
 */

  extern G4double getpdircosu_(G4int *index);
  extern G4double getpdircosv_(G4int *index);
  extern G4double getpdircosw_(G4int *index);
/*
 * Given the index of the emitted photon, this function returns
 * the direction cosines of its velocity.
 */

  extern G4double getnage_(G4int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * its age, -1 if index is out of range.
 * This age will be different from the time specified
 * in generateFissionEvent and generateSpontaneousFissionEvent
 * for non-prompt neutrons, i.e. delayed neutrons. 
 */

  extern G4double getpage_(G4int *index);
/*
 * Given the index of the emitted photon, this function returns
 * its age, -1 of index is out of range.
 *  This age will be different from the time specified
 * in generateFissionEvent and generateSpontaneousFissionEvent
 * for photons that are emitted by beta-decay of the fission
 * fragments.
 */

  extern void setdelay_(G4int *delay);
/*
 * This function is called to enable delayed neutrons and photons
 * Input
 *      delay:
 *              0 (default) for strictly prompt neutrons and 
 *                photons
 *              1 (n/a) for prompt neutrons, prompt and delayed 
 *                photons
 *              2 (n/a) for prompt and delayed neutrons, prompt 
 *                photons
 *              3 (n/a) for prompt and delayed neutrons, prompt 
 *                and delayed photons
 */

  extern void setcorrel_(G4int *correlation);
/*
 * This function is called to set the type of neutron photon correlation
 * Input
 *      correlation:
 *              0 (default) for no correlation between neutrons and
 *                photons
 *              1 (n/a) for number correlation between neutrons and 
 *                photons
 *              2 (n/a) for number and energy correlation between 
 *                neutrons and photons
 */

  extern void setnudist_(G4int *nudist);
/*
 * This function is called to set the data to be sampled for the neutron
 * number distributions in induced fissions
 * Input
 *      nudist:
 *               0 to use the fit to the Zucker and Holden tabulated 
 *                 P(nu) distributions as a function of energy for 
 *                 U235, U238 and Pu239. Terrell for other isotopes.
 *               1 to use fits to the Zucker and Holden tabulated 
 *                 P(nu) distribution as a function of energy for 
 *                 U238 and Pu239, and a fit to the Zucker and Holden 
 *                 data as well as the Gwin, Spencer and Ingle data 
 *                 (at thermal energies) as a function of energy for 
 *                 U235. Terrell for other isotopes.
 *               2 (default) to use the fit to the Zucker and Holden 
 *                 tabulated P(nu) distributions as a function of nubar. 
 *                 The U238 fit is used for the U232, U234, U236 and 
 *                 U238 isotopes, the U235 fit for U233 and U235, the 
 *                 Pu239 fit for Pu239 and Pu241. Terrell for other 
 *                 isotopes.
 */


  extern void setcf252_(G4int *ndist, G4int *neng);
/*
 * This function is called to set the data to be sampled for the 
 * (a) Cf252 spontaneous fission number distribution, and 
 * (b) Cf252 spontaneous fission neutron energy spectrum
 * Input
 *      ndist:
 *              0 (default) to sample the number of neutrons from the 
 *                tabulated data measured by Spencer
 *              1 to sample the number of neutrons from Boldeman's data
 *      neng:
 *              0 to sample the spontaneous fission neutron energy from 
 *                Mannhart corrected Maxwellian spectrum
 *              1 to sample the spontaneous fission neutron energy from 
 *                Madland-Nix theoretical spectrum
 *              2 to sample the spontaneous fission neutron energy from 
 *                the Froehner Watt spectrum
 */

  extern void setrngf_(G4float (*funcptr) (void));
/*
 * This function sets the random number generator to the user-defined
 * one specified in the argument. If either setrngf_ or setrngd_ are
 * not specified, the default system call srand48 will be called.
 * Input
 *      funcptr:
 *               a random number generator function that returns a
 *               variable of type G4float
 */

  extern void setrngd_(G4double (*funcptr) (void));
/*
 * This function sets the random number generator to the user-defined
 * one specified in the argument. If either setrngf_ or setrngd_ are
 * not specified, the default system call srand48 will be called.
 * Input
 *      funcptr:
 *               a random number generator function that returns a
 *               variable of type G4double
 */
// }
