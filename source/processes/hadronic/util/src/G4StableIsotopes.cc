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
//
#include "G4StableIsotopes.hh"

const G4int G4StableIsotopes::protonCount[92] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92};
const G4String G4StableIsotopes::elementName[92] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "", "Nb", "Mo", "", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "", "", "", "", "", "", "Th", "", "U"};
const G4int G4StableIsotopes::nIsotopes[92] = {2, 2, 2, 1, 2, 2, 2, 3, 1, 3, 1, 3, 1, 3, 1, 4, 2, 3, 3, 6, 1, 5, 2, 4, 1, 4, 1, 5, 2, 5, 2, 5, 1, 6, 2, 6, 2, 4, 1, 5, 1, 7, 0, 7, 1, 6, 2, 8, 2, 10, 2, 8, 1, 9, 1, 7, 2, 4, 1, 7, 0, 7, 2, 7, 1, 7, 1, 6, 1, 7, 2, 6, 2, 5, 2, 7, 2, 6, 1, 7, 2, 4, 1, 0, 0, 0, 0, 0, 0, 1, 0, 3};
const G4int G4StableIsotopes::start[92] = {0, 2, 4, 6, 7, 9, 11, 13, 16, 17, 20, 21, 24, 25, 28, 29, 33, 35, 38, 41, 47, 48, 53, 55, 59, 60, 64, 65, 70, 72, 77, 79, 84, 85, 91, 93, 99, 101, 105, 106, 111, 112, 119, 119, 126, 127, 133, 135, 143, 145, 155, 157, 165, 166, 175, 176, 183, 185, 189, 190, 197, 197, 204, 206, 213, 214, 221, 222, 228, 229, 236, 238, 244, 246, 251, 253, 260, 262, 268, 269, 276, 278, 282, 283, 283, 283, 283, 283, 283, 283, 284, 284};
const G4int G4StableIsotopes::nucleonCount[287] = {1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 36, 35, 37, 36, 38, 40, 39, 40, 41, 40, 42, 43, 44, 46, 48, 45, 46, 47, 48, 49, 50, 50, 51, 50, 52, 53, 54, 55, 54, 56, 57, 58, 59, 58, 60, 61, 62, 64, 63, 65, 64, 66, 67, 68, 70, 69, 71, 70, 72, 73, 74, 76, 75, 74, 76, 77, 78, 80, 82, 79, 81, 78, 80, 82, 83, 84, 86, 85, 87, 84, 86, 87, 88, 89, 90, 91, 92, 94, 96, 93, 92, 94, 95, 96, 97, 98, 100, 96, 98, 99, 100, 101, 102, 104, 103, 102, 104, 105, 106, 108, 110, 107, 109, 106, 108, 110, 111, 112, 113, 114, 116, 113, 115, 112, 114, 115, 116, 117, 118, 119, 120, 122, 124, 121, 123, 120, 122, 123, 124, 125, 126, 128, 130, 127, 124, 126, 128, 129, 130, 131, 132, 134, 136, 133, 130, 132, 134, 135, 136, 137, 138, 138, 139, 136, 138, 140, 142, 141, 142, 143, 144, 145, 146, 148, 150, 144, 147, 148, 149, 150, 152, 154, 151, 153, 152, 154, 155, 156, 157, 158, 160, 159, 156, 158, 160, 161, 162, 163, 164, 165, 162, 164, 166, 167, 168, 170, 169, 168, 170, 171, 172, 173, 174, 176, 175, 176, 174, 176, 177, 178, 179, 180, 180, 181, 180, 182, 183, 184, 186, 185, 187, 184, 186, 187, 188, 189, 190, 192, 191, 193, 190, 192, 194, 195, 196, 198, 197, 196, 198, 199, 200, 201, 202, 204, 203, 205, 204, 206, 207, 208, 209, 232, 234, 235, 238};
const G4double G4StableIsotopes::abundance[287] = {99.985, 0.015, 0.000137, 99.9999, 7.5, 92.5, 100, 19.9, 80.1, 98.9, 1.1, 99.634, 0.366, 99.762, 0.038, 0.2, 100, 90.48, 0.27, 9.25, 100, 78.99, 10, 11.01, 100, 92.23, 4.67, 3.1, 100, 95.02, 0.75, 4.21, 0.02, 75.77, 24.23, 0.337, 0.063, 99.6, 93.2581, 0.0117, 6.7302, 96.941, 0.647, 0.135, 2.086, 0.004, 0.187, 100, 8, 7.2, 73.8, 5.5, 5.4, 0.25, 99.75, 4.345, 83.789, 9.501, 2.365, 100, 5.8, 91.72, 2.2, 0.28, 100, 68.277, 26.223, 1.14, 3.634, 0.926, 69.17, 30.83, 48.6, 27.9, 4.1, 18.8, 0.6, 60.108, 39.892, 21.23, 27.66, 7.73, 36.94, 7.44, 100, 0.89, 9.36, 7.63, 23.78, 49.61, 8.73, 50.69, 49.31, 0.35, 2.25, 11.6, 11.5, 57, 17.3, 72.165, 27.835, 0.56, 9.86, 7, 82.58, 100, 51.45, 11.22, 17.15, 17.38, 2.8, 100, 14.84, 9.25, 15.92, 16.68, 9.55, 24.13, 9.63, 5.52, 1.88, 12.7, 12.6, 17, 31.6, 18.7, 100, 1.02, 11.14, 22.33, 27.33, 26.46, 11.72, 51.839, 48.161, 1.25, 0.89, 12.49, 12.8, 24.13, 12.22, 28.73, 7.49, 4.3, 95.7, 0.97, 0.65, 0.36, 14.53, 7.68, 24.22, 8.58, 32.59, 4.63, 5.79, 57.36, 42.64, 0.096, 2.6, 0.908, 4.816, 7.14, 18.95, 31.69, 33.8, 100, 0.1, 0.09, 1.91, 26.4, 4.1, 21.2, 26.9, 10.4, 8.9, 100, 0.106, 0.101, 2.417, 6.592, 7.854, 11.23, 71.7, 0.0902, 99.9098, 0.19, 0.25, 88.48, 11.08, 100, 27.13, 12.18, 23.8, 8.3, 17.19, 5.76, 5.64, 3.1, 15, 11.3, 13.8, 7.4, 26.7, 22.7, 47.8, 52.2, 0.2, 2.18, 14.8, 20.47, 15.65, 24.84, 21.86, 100, 0.06, 0.1, 2.34, 18.9, 25.5, 24.9, 28.2, 100, 0.14, 1.61, 33.6, 22.95, 26.8, 14.9, 100, 0.13, 3.05, 14.3, 21.9, 16.12, 31.8, 12.7, 97.41, 2.59, 0.162, 5.206, 18.606, 27.297, 13.629, 35.1, 0.012, 99.988, 0.13, 26.3, 14.3, 30.67, 28.6, 37.4, 62.6, 0.02, 1.58, 1.6, 13.3, 16.1, 26.4, 41, 37.3, 62.7, 0.01, 0.79, 32.9, 33.8, 25.3, 7.2, 100, 0.15, 9.97, 16.87, 23.1, 13.18, 29.86, 6.87, 29.524, 70.476, 1.4, 24.1, 22.1, 52.4, 100, 100, 0.005, 0.72, 99.275}; 
