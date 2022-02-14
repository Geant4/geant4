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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4IsotopeList
//
// Author  V. Ivantchenko, 22 October 2020
//
 
#ifndef G4IsotopeList_h
#define G4IsotopeList_h 1

static const G4int amin[] = {
   0,
   1,   3,   6,   9,  10,  12,  14,  16,  19,  20, //1-10
  23,  24,  27,  27,  31,  32,  35,  36,  39,  40, //11-20
  45,  46,  50,  50,  55,  54,  59,  58,  63,  64, //21-30
  69,  70,  75,  74,  79,  78,  85,  84,  89,  90, //31-40
  93,  92,  98,  96, 103, 102, 107, 106, 113, 112, //41-50
 121, 120, 127, 124, 133, 130, 138, 136, 141, 142, //51-60
 145, 144, 151, 152, 158, 156, 165, 162, 169, 168, //61-70
 175, 174, 180, 180, 185, 184, 191, 190, 197, 196, //71-80
 203, 204, 209, 209, 210, 222, 223, 226, 227, 232, //81-90
 231, 233, 237, 238};

static const G4int amax[] = {
   0,
   3,   4,   7,   9,  11,  14,  15,  18,  19,  22,  //1-10
  23,  26,  27,  30,  31,  36,  37,  40,  41,  48,  //11-20
  45,  50,  51,  54,  55,  58,  59,  64,  65,  70,  //21-30
  71,  76,  75,  82,  81,  86,  87,  90,  89,  96,  //31-40
  94, 100,  98, 104, 103, 110, 109, 116, 115, 124,  //41-50
 123, 130, 129, 136, 137, 138, 139, 142, 141, 150,  //51-60
 145, 154, 153, 160, 159, 164, 165, 170, 169, 176,  //61-70
 176, 180, 181, 186, 187, 192, 193, 198, 197, 204,  //71-80
 205, 208, 209, 209, 210, 222, 223, 226, 227, 232,  //81-90
 231, 238, 237, 244};

static const G4double aeff[] = {
  0.,
 1.00794, 4.00264, 6.94003, 9.01218,  10.811, 12.0107, 14.0068, 15.9994, 18.9984,   20.18,  //1-10
 22.9898,  24.305, 26.9815, 28.0854, 30.9738, 32.0661, 35.4526, 39.9477, 39.0983,  40.078,  //11-20
 44.9559, 47.8667, 50.9415, 51.9961,  54.938, 55.8451, 58.9332, 58.6933, 63.5456, 65.3955,  //21-30
 69.7231, 72.6128, 74.9216, 78.9594, 79.9035, 83.7993, 85.4677, 87.6166, 88.9058, 91.2236,  //31-40
 92.9064, 95.9313, 97.9072, 101.065, 102.906, 106.415, 107.868, 112.411, 114.818,  118.71,  //41-50
  121.76, 127.603, 126.904, 131.292, 132.905, 137.327, 138.905, 140.115, 140.908, 144.236,  //51-60
 144.913, 150.366, 151.964, 157.252, 158.925, 162.497,  164.93, 167.256, 168.934, 173.038,  //61-70
 174.967, 178.485, 180.948, 183.842, 186.207, 190.225, 192.216, 195.078, 196.967, 200.599,  //71-80
 204.383, 207.217,  208.98, 208.982, 209.987, 222.018,  223.02, 226.025, 227.028, 232.038,  //81-90
 231.036, 238.029, 237.048, 244.064};

#endif
