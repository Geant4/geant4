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
// File:        G4MaterialPropertiesIndex.hh
// Description: Indices and Names for G4MaterialProperties
// Created:     29-06-2017
// Author:      Soon Yung Jun
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertiesIndex_h
#define G4MaterialPropertiesIndex_h 1

#include <vector>
#include "G4String.hh"

enum G4MaterialPropertyIndex   {
  kNullPropertyIndex = -1,     // the number of G4MaterialPropertyIndex
  kRINDEX,                     // index of refraction                  
  kREFLECTIVITY,               // reflectivity         
  kREALRINDEX,                 // real part of the refractive index
  kIMAGINARYRINDEX,            // imaginary part of the refractive index
  kEFFICIENCY,                 // efficiency 
  kTRANSMITTANCE,              // transmittance of a dielectric surface
  kSPECULARLOBECONSTANT,       // reflection probability about the normal of a micro facet. 
  kSPECULARSPIKECONSTANT,      // reflection probability about the average surface normal
  kBACKSCATTERCONSTANT,        // for the case of several reflections within a deep groove
  kGROUPVEL,                   // group velocity
  kMIEHG,                      // Mie scattering length
  kRAYLEIGH,                   // Rayleigh scattering attenuation length
  kWLSCOMPONENT,               // the relative emission spectrum of the material as a function of the photon's momentum
  kWLSABSLENGTH,               // the absorption length of the material as a function of the photon's momentum
  kABSLENGTH,                  // the absorption length
  kFASTCOMPONENT,              // fast component of scintillation
  kSLOWCOMPONENT,              // slow component of scintillation
  kPROTONSCINTILLATIONYIELD,   // scintillation light yield by protons  
  kDEUTERONSCINTILLATIONYIELD, // scintillation light yield by deuterons
  kTRITONSCINTILLATIONYIELD,   // scintillation light yield by tritons
  kALPHASCINTILLATIONYIELD,    // scintillation light yield by alphas
  kIONSCINTILLATIONYIELD,      // scintillation light yield by ions
  kELECTRONSCINTILLATIONYIELD, // scintillation light yield by electrons
  kNumberOfPropertyIndex       // the number of G4MaterialPropertyIndex
} ;

enum G4MaterialConstPropertyIndex 
{ 
  kNullConstPropertyIndex = -1, // the number of G4MaterialPropertyIndex
  kSURFACEROUGHNESS,            // surface microroughness      
  kISOTHERMAL_COMPRESSIBILITY,  // isothermal compressibility
  kRS_SCALE_FACTOR,             // Rayleigh scattering scale factor
  kWLSMEANNUMBERPHOTONS,        // WLS mean number of photons
  kWLSTIMECONSTANT,             // any time delay which may occur between absorption and re-emission of the photon
  kMIEHG_FORWARD,               // forward angle of Mie scattering based on Henyey-Greenstein phase function
  kMIEHG_BACKWARD,              // backward angle of Mie scattering based on Henyey-Greenstein phase function
  kMIEHG_FORWARD_RATIO,	        // ratio of the MIEHG forward scattering 
  kSCINTILLATIONYIELD,	        // scintillation light yield
  kRESOLUTIONSCALE,	        // resolution scale
  kFASTTIMECONSTANT,	        // fast scintillation time constant 
  kFASTSCINTILLATIONRISETIME,   // rise time of fast scintillation 
  kSLOWTIMECONSTANT,	        // slow scintillation time constant
  kSLOWSCINTILLATIONRISETIME,   // rise time of slow scintillation 
  kYIELDRATIO,                  // yield ratio
  kFERMIPOT,                    // the Fermi potential (in neV)
  kDIFFUSION,                   // diffusion
  kSPINFLIP,		        // spin flip
  kLOSS,		        // loss
  kLOSSCS,		        // loss cross-section
  kABSCS,		        // 1/v energy dependent absorption cross section
  kSCATCS,                      // incoherent elastic scattering cross-section
  kMR_NBTHETA,                  // number of theta bins of microroughness (MR)
  kMR_NBE,                      // number of energy bins 
  kMR_RRMS,                     // RMS of roughness
  kMR_CORRLEN,                  // correlation length
  kMR_THETAMIN,                 // minimum value of theta
  kMR_THETAMAX,                 // maximum value of theta
  kMR_EMIN,                     // mininum value of energy
  kMR_EMAX,                     // maximum value of energy
  kMR_ANGNOTHETA,               // number of theta angles in the look-up table
  kMR_ANGNOPHI,                 // number of phi angles in the look-up table
  kMR_ANGCUT,                   // angular cut
  kNumberOfConstPropertyIndex   // the number of G4MaterialConstPropertyIndex
};

#endif /* G4MaterialPropertiesIndex_h */
