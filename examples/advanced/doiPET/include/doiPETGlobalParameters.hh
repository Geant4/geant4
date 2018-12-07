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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan

//The header file consists of specification if the scanner geometry

//Number of crystals (scintillators) in the DOI (x), tangential (y) and axial (z) direction. 
#define numberOfCrystal_DOI 4
#define numberOfCrystal_tangential 16
#define numberOfCrystal_axial 16

//size of the crystals in the DOI (x), tangential (y) and axial (z) direction. 
#define sizeOfCrystal_DOI 7.5*mm
#define sizeOfCrystal_tangential 2.8*mm
#define sizeOfCrystal_axial 2.8*mm

//Inter crystal gap between crystals. 
#define crystalGap_DOI 0.0*mm
#define crystalGap_tangential 0.08*mm
#define crystalGap_axial 0.08*mm



//Define Aluminum thickness cover for the detector
#define AluminumCoverThickness 0.3 * mm

//Number of PET detectors per ring
#define numberOfDetector_perRing 40

//Number of rings in the PET system 
#define numberOfRings 4

//Radius of the PET system. Note that the radius is defined between to opposing sensitive detectors (It does not include detector block cover).
#define scannerRadius 330 * mm

//Gap between two adjacent rings. There are three gaps in four ring PET system
#define ringGap 10.4 * mm