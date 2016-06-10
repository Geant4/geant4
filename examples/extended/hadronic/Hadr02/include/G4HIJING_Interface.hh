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
#ifndef G4HIJING_Interface_hh
#define G4HIJING_Interface_hh

//
// MODULE:           G4HIJING_Model.hh
//
// Version:          1.B
// Date:       10/09/2013
// Author:     Khaled Abdel-Waged 
// Institute:  Umm Al-Qura University
// Country:    SAUDI ARABIA       
//

// First: 1-COMMON BLOCK FOR OPTIONS AND PARAMETERS
//-----
// HIPARNT
// input parameters 
// HIPR1, HIPR2 for event options
// HINT1, HINT2 of current event
//---------------------------------
struct cchijinghiparnt
{

float hipr1[100];
G4int ihpr2[50];
float hint1[100];
G4int ihnt2[50];

};

//****************************************************************************
// Second: 5-COMMON BLOCKS FOR EVENT INFORMATION
//----------------------------------------------------------------------------

//
// HIMAIN1-> Global information of the events are defined
//
struct cchijinghimain1
{
G4int natt;
float eatt;
G4int jatt, nt, np,n0,n01,n10,n11;
};
//-----------
// HIMAIN2->information of produced stable and undecayed particles
//-------
struct cchijinghimain2
{
G4int katt[4][130000];
float patt[4][130000];
};
//--------
// HIJJET1-> information about produced partons which 
//           are connected with the valence quarks, diquarks,...
//--------
struct cchijinghijjet1
{
G4int npj[300], kfpj[500][300];
float pjpx[500][300],pjpy[500][300], pjpz[500][300],pjpe[500][300];
float pjpm[500][300];
G4int   ntj[300], kftj[500][300];
float pjtx[500][300], pjty[500][300], pjtz[500][300];
float pjte[500][300], pjtm[500][300];
};
//--------
// HIJJET1-> information about produced partons which 
//           will form string systems without being connected with
//           valence quarks, diquarks,...
// ------
struct cchijinghijjet2
{
G4int nsg, njsg[900], iasg[3][900], k1sg[100][900];
G4int k2sg[100][900];
float pxsg[100][900], pysg[100][900], pzsg[100][900];
float pesg[100][900], pmsg[100][900];
};

//------
// HISTRNG
// contain information about the projectile and target nucleons
//-----
struct cchijinghistrng
{
G4int nfp[15][300];
float pp[15][300];
G4int nft[15][300];
float pt[15][300];
};
//****************************************************************************
// third: 2-COMMON BLOCKS which contain specific information
//----------------------------------------------------------------------------

struct cchijinghijjet4
{
G4int ndr, iadr[2][900], kfdr[900];
float pdr[5][900];
};

struct cchijinghijcrdn
{
float yp[300][3], yt[300][3];
};
//---------------------------------------------------------
// fourth: 5-Other common blocks
// --------------------------------------------------------


struct cchijingbveg1
{
float xl, xu, acc;
G4int ndim, ncall, itmx, nprn;
};




struct cchijingseedvax
{
G4int num1;
};



struct cchijingranseed
{
float nseed;
};



struct cchijinghijdat
{
float hidat0[10][10],hidat[10];
};


struct cchijinghipyint
{
G4int mint4, mint5;
float atco[20][200], atxs[200+1];
};

// hijing

extern "C"
{
// initialize HIJING for specified event type, 
// collision frame and energy 

//extern void hijset_ (float*,
////                     const char*, const char*,const char*,
//                     G4int*, G4int*, G4int*, G4int*);

extern void hijset_ (float*);

// to generate a complete event as specified by sybroutine HIJSET

//extern void hijing_ (const char*,
//                    float*, float*);

extern void hijing_ (float*,float*);

extern float ulmass_ (G4int*);



// reset all relevant common blocks and variables and initialize HIJING 
// for each event

extern void hijini_ ();

// calculate cross sections for minijet production, cross section of
// the triggered processes, elastic, inelastic, total cross section..
extern void hijcrs_ ();

// 
// initialize program for generating hard scattering 
//as specified by  parameters and options

extern void jetini_ (G4int*, G4int*, G4int*);
//

// re-initiate PYTHIA for the triggered hard processe
// or simulate one hard scattering among the multiple jet production 
//per NN-collision

extern void hijhrd_ (G4int*, G4int*, G4int*, G4int*, G4int*);

// 
//generate soft interaction for each binary NN-collision

extern void hijsft_ (G4int*, G4int*, G4int*);

// rearrange gluon jets in a string system according to their rapidities

extern void hijsrt_ (G4int*, G4int*);

// perform jet quenching by allowing final state interaction of produced jet
// inside excited strings

extern void quench_ (G4int*, G4int*);

// 
// arrange produced partons together with the valence quarks and diquarks

extern void hijfrg_ (G4int*, G4int*, G4int*);

// perform soft radiation according to the Lund dipole approx.

extern void attrad_ (G4int*);

// generate flavor codes of the valence quark (diquark) 
//inside a given nucleon (hadron).

extern void attflv_ (G4int*, G4int*, G4int*);

// perform elastic scattering and possible elastic NN cascading

extern void hijcsc_ (G4int*, G4int*);

// three parameter Wood-Sax distribution

extern void hijwds_ (G4int*, G4int*, float*);

// gives profile function of 2 colliding nuclei at a given impact parameter

extern float profile_ (float*);

// transform the produced particles from c.m to lab frame

extern void hiboost_ ();

//----------------------------------------
// the default values of the parametrs and options to initialize
// the event record common blocks

extern void g4hijingblockdata_ ();
// ----------------------
// random generator

extern void rlu_ (G4int*);

//-----------------------------------------
extern struct cchijinghiparnt hiparnt_;

extern struct cchijinghimain1 himain1_;
extern struct cchijinghimain2 himain2_;

extern struct cchijinghijjet1 hijjet1_;
extern struct cchijinghijjet2 hijjet2_;

extern struct cchijinghistrng histrng_;

extern struct cchijinghijjet4 hijjet4_;

extern struct cchijinghijcrdn hijcrdn_;

extern struct cchijingbveg1 bveg1_;

extern struct cchijingseedvax seedvax_;

extern struct cchijingranseed ranseed_;

extern struct cchijinghijdat hijdat_;

extern struct cchijinghipyint hipyint_;

}

#endif

