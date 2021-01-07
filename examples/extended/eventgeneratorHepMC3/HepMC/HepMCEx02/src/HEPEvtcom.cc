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
/// \file eventgenerator/HepMC/HepMCEx02/src/HEPEvtcom.cc
/// \brief Implementation of the HEPEvtcom class
//
// ======================================================================
//      PARAMETER (NMXHEP=4000)
//      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
//     &        JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
// ======================================================================
///**********************************************************/
///*           D E S C R I P T I O N :                      */
///*--------------------------------------------------------*/
///* NEVHEP          - event number (or some special meaning*/
///*                    (see documentation for details)     */
///* NHEP            - actual number of entries in current  */
///*                    event.                              */
///* ISTHEP[IHEP]    - status code for IHEP'th entry - see  */
///*                    documentation for details           */
///* IDHEP [IHEP]    - IHEP'th particle identifier according*/
///*                    to PDG.                             */
///* JMOHEP[IHEP][0] - pointer to position of 1st mother    */
///* JMOHEP[IHEP][1] - pointer to position of 2nd mother    */
///* JDAHEP[IHEP][0] - pointer to position of 1st daughter  */
///* JDAHEP[IHEP][1] - pointer to position of 2nd daughter  */
///* PHEP  [IHEP][0] - X momentum [Gev/c]                   */
///* PHEP  [IHEP][1] - Y momentum [Gev/c]                   */
///* PHEP  [IHEP][2] - Z momentum [Gev/c]                   */
///* PHEP  [IHEP][3] - Energy [Gev]                         */
///* PHEP  [IHEP][4] - Mass[Gev/c^2]                        */
///* VHEP  [IHEP][0] - X vertex [mm]                        */
///* VHEP  [IHEP][1] - Y vertex [mm]                        */
///* VHEP  [IHEP][2] - Z vertex [mm]                        */
///* VHEP  [IHEP][3] - production time [mm/c]               */
///*========================================================*/
//
// This interface to HEPEVT common block treats the block as
// an array of bytes --- the precision and number of entries
// is determined "on the fly" by the wrapper and used to decode
// each entry.
//
// HEPEVT_EntriesAllocation is the maximum size of the HEPEVT common block
// that can be interfaced. It is NOT the actual size of the HEPEVT common
// used in each individual application. The actual size can be changed on
// the fly using HepMC::HEPEVT_Wrapper::set_max_number_entries().
// Thus HEPEVT_EntriesAllocation should typically be set
// to the maximum possible number of entries --- 10000 is a good choice
// (and is the number used by ATLAS versions of Pythia).

#include <ctype.h>

enum {HEPEVT_EntriesAllocation=4000};

const unsigned int hepevt_bytes_allocation =
sizeof(long int) * ( 2 + 4 * HEPEVT_EntriesAllocation )
  + sizeof(double) * ( 9 * HEPEVT_EntriesAllocation );

extern "C" struct hepevt{
  char data[hepevt_bytes_allocation];
};

hepevt hepevt_;
