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
// Code developed by:
// S.Guatelli, M. Large and A. Malaroda, University of Wollongong
//
//
#include "ICRP110PhantomMaterial_Female.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

ICRP110PhantomMaterial_Female::ICRP110PhantomMaterial_Female(): 
  fLung(nullptr), fTeeth(nullptr), fBone(nullptr), fHumeri_upper(nullptr), fHumeri_lower(nullptr),
  fArm_lower(nullptr), fHand(nullptr), fClavicle(nullptr), fCranium(nullptr), fFemora_upper(nullptr),
  fFemora_lower(nullptr), fLeg_lower(nullptr), fFoot(nullptr), fMandible(nullptr),
  fPelvis(nullptr), fRibs(nullptr), fScapulae(nullptr), fSpine_cervical(nullptr),
  fSpine_lumbar(nullptr), fSpine_thoratic(nullptr), fSacrum(nullptr),
  fSternum(nullptr), fHf_upper(nullptr), fHf_lower(nullptr), fMed_lowerleg(nullptr),
  fMed_lowerarm(nullptr), fCartilage(nullptr), fSkin(nullptr), fBlood(nullptr),
  fMuscle(nullptr), fLiver(nullptr), fPancreas(nullptr), fBrain(nullptr), fHeart(nullptr), fEye(nullptr),
  fKidney(nullptr), fStomach(nullptr), fIntestine_sml(nullptr), fIntestine_lrg(nullptr),
  fSpleen(nullptr), fThyroid(nullptr), fBladder(nullptr), fOvaries_testes(nullptr), fAdrenals(nullptr),
  fOesophagus(nullptr), fMisc(nullptr), fUterus_prostate(nullptr), fLymph(nullptr),
  fBreast_glandular(nullptr), fBreast_adipose(nullptr), fGastro_content(nullptr),
  fUrine(nullptr)
{;}

ICRP110PhantomMaterial_Female::~ICRP110PhantomMaterial_Female()
{;}

void ICRP110PhantomMaterial_Female::DefineMaterials()
{
  // Define required materials

  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4double d;  // density
 
  // General elements
 
  A = 1.01*g/mole;
  auto elH = new G4Element ("Hydrogen","H",Z = 1.,A);

  A = 12.011*g/mole;
  auto elC = new G4Element("Carbon","C",Z = 6.,A);  

  A = 14.01*g/mole;
  auto elN = new G4Element("Nitrogen","N",Z = 7.,A);

  A = 16.00*g/mole;
  auto elO = new G4Element("Oxygen","O",Z = 8.,A);

  A = 22.99*g/mole;
  auto elNa = new G4Element("Sodium","Na",Z = 11.,A);

  A = 24.305*g/mole;
  auto elMg = new G4Element("Magnesium","Mg",Z = 12.,A);

  A = 30.974*g/mole;
  auto elP = new G4Element("Phosphorus","P",Z = 15.,A);
 
  A = 32.064*g/mole;
  auto elS = new G4Element("Sulfur","S",Z = 16.,A);
 
  A = 35.453*g/mole;
  auto elCl = new G4Element("Chlorine","Cl",Z = 17.,A);
 
  A = 39.098*g/mole;
  auto elK = new G4Element("Potassium","K",Z = 19.,A);

  A = 40.08*g/mole;
  auto elCa = new G4Element("Calcium","Ca",Z = 20.,A);

  A = 55.85*g/mole;
  auto elFe  = new G4Element("Iron","Fe",Z = 26.,A);
 
  A = 126.90447 *g/mole;
  auto elI = new G4Element("Iodine","I", Z = 53.,A);

  //Added tissues for phantom including their tissue ID
  //Teeth -> Tissue ID 1
  d = 2.750 *g/cm3;
  fTeeth = new G4Material("teeth",d,7);
  fTeeth -> AddElement(elH,0.022);
  fTeeth -> AddElement(elC,0.095);
  fTeeth -> AddElement(elN,0.029);
  fTeeth -> AddElement(elO,0.421);
  fTeeth -> AddElement(elMg,0.007);
  fTeeth -> AddElement(elP,0.137);
  fTeeth -> AddElement(elCa,0.289);
  
  //Mineral Bone -> Tissue ID 2
  d = 1.920 *g/cm3;
  fBone = new G4Material("bone",d,9);
  fBone -> AddElement(elH,0.036);
  fBone -> AddElement(elC,0.159);
  fBone -> AddElement(elN,0.042);
  fBone -> AddElement(elO,0.448);
  fBone -> AddElement(elNa,0.003);
  fBone -> AddElement(elMg,0.002);
  fBone -> AddElement(elP,0.094);
  fBone -> AddElement(elS,0.003);
  fBone -> AddElement(elCa,0.213);
  
  //Humeri, upper half, spongiosa -> ID 3
  d = 1.185 *g/cm3;
  fHumeri_upper = new G4Material("humeri_upper",d,11);
  fHumeri_upper -> AddElement(elH,0.087) ;
  fHumeri_upper -> AddElement(elC,0.366);
  fHumeri_upper -> AddElement(elN,0.025);
  fHumeri_upper -> AddElement(elO,0.422);
  fHumeri_upper -> AddElement(elNa,0.002);
  fHumeri_upper -> AddElement(elMg,0.001);
  fHumeri_upper -> AddElement(elP,0.030);
  fHumeri_upper -> AddElement(elS,0.003);
  fHumeri_upper -> AddElement(elCl,0.001);
  fHumeri_upper -> AddElement(elK, 0.001);
  fHumeri_upper -> AddElement(elCa,0.062);  
  
  //Humeri, lower half, spongiosa -> ID 4
  d = 1.117 *g/cm3;
  fHumeri_lower = new G4Material("humeri_lower",d,9);
  fHumeri_lower -> AddElement(elH,0.096);
  fHumeri_lower -> AddElement(elC,0.473);
  fHumeri_lower -> AddElement(elN,0.017);
  fHumeri_lower -> AddElement(elO,0.341);
  fHumeri_lower -> AddElement(elNa,0.002);
  fHumeri_lower -> AddElement(elP,0.022);
  fHumeri_lower -> AddElement(elS,0.002);
  fHumeri_lower -> AddElement(elCl,0.001);
  fHumeri_lower -> AddElement(elCa,0.046); 
  
  //Lower arm bones, spongiosa -> ID 5
  d = 1.117 *g/cm3;
  fArm_lower = new G4Material("arm_lower",d,9);
  fArm_lower -> AddElement(elH,0.096);
  fArm_lower -> AddElement(elC,0.473);
  fArm_lower -> AddElement(elN,0.017);
  fArm_lower -> AddElement(elO,0.341);
  fArm_lower -> AddElement(elNa,0.002);
  fArm_lower -> AddElement(elP,0.022);
  fArm_lower -> AddElement(elS,0.002);
  fArm_lower -> AddElement(elCl,0.001);
  fArm_lower -> AddElement(elCa,0.046); 
  
  //Hand Bones, Spongiosa ->ID 6
  d = 1.117 *g/cm3;
  fHand = new G4Material("hand",d,9);
  fHand -> AddElement(elH,0.096);
  fHand -> AddElement(elC,0.473);
  fHand -> AddElement(elN,0.017);
  fHand -> AddElement(elO,0.341);
  fHand -> AddElement(elNa,0.002);
  fHand -> AddElement(elP,0.022);
  fHand -> AddElement(elS,0.002);
  fHand -> AddElement(elCl,0.001);
  fHand -> AddElement(elCa,0.046); 
  
  //Clavicles, spongiosa -> ID 7
  d = 1.191 *g/cm3;
  fClavicle = new G4Material("clavicle",d,11);
  fClavicle -> AddElement(elH,0.087);
  fClavicle -> AddElement(elC,0.361);
  fClavicle -> AddElement(elN,0.025);
  fClavicle -> AddElement(elO,0.424);
  fClavicle -> AddElement(elNa,0.002);
  fClavicle -> AddElement(elMg,0.001);
  fClavicle -> AddElement(elP,0.031);
  fClavicle -> AddElement(elS,0.003);
  fClavicle -> AddElement(elCl,0.001);
  fClavicle -> AddElement(elK,0.001);
  fClavicle -> AddElement(elCa,0.064);
  
  // Cranium, spongiosa -> ID 8
  d = 1.245 *g/cm3;
  fCranium = new G4Material("cranium",d,11);
  fCranium -> AddElement(elH,0.081);
  fCranium -> AddElement(elC,0.317);
  fCranium -> AddElement(elN,0.028);
  fCranium -> AddElement(elO,0.451);
  fCranium -> AddElement(elNa,0.002);
  fCranium -> AddElement(elMg, 0.001);
  fCranium -> AddElement(elP,0.037);
  fCranium -> AddElement(elS,0.003);
  fCranium -> AddElement(elCl,0.001);
  fCranium -> AddElement(elK,0.001);
  fCranium -> AddElement(elCa,0.078);
  
  //Femora, upper half, spongiosa -> ID 9
  d = 1.046 *g/cm3;
  fFemora_upper = new G4Material("femora_upper",d,10);
  fFemora_upper -> AddElement(elH,0.104);
  fFemora_upper -> AddElement(elC,0.496);
  fFemora_upper -> AddElement(elN,0.018);
  fFemora_upper -> AddElement(elO,0.349);
  fFemora_upper -> AddElement(elNa,0.001);
  fFemora_upper -> AddElement(elP,0.009);
  fFemora_upper -> AddElement(elS,0.002);
  fFemora_upper -> AddElement(elCl,0.001);
  fFemora_upper -> AddElement(elK,0.001);
  fFemora_upper -> AddElement(elCa,0.019);
  
  //Femora, lower half, spongiosa -> ID 10
  d = 1.117 *g/cm3;
  fFemora_lower = new G4Material("femora_lower",d,9);
  fFemora_lower -> AddElement(elH,0.096);
  fFemora_lower -> AddElement(elC,0.473);
  fFemora_lower -> AddElement(elN,0.017);
  fFemora_lower -> AddElement(elO,0.341);
  fFemora_lower -> AddElement(elNa,0.002);
  fFemora_lower -> AddElement(elP,0.022);
  fFemora_lower -> AddElement(elS,0.002);
  fFemora_lower -> AddElement(elCl,0.001);
  fFemora_lower -> AddElement(elCa,0.046); 
  
  //Lower leg bones, spongiosa -> ID 11
  d = 1.117 *g/cm3;
  fLeg_lower = new G4Material("leg_lower",d,9);
  fLeg_lower -> AddElement(elH,0.096);
  fLeg_lower -> AddElement(elC,0.473);
  fLeg_lower -> AddElement(elN,0.017);
  fLeg_lower -> AddElement(elO,0.341);
  fLeg_lower -> AddElement(elNa,0.002);
  fLeg_lower -> AddElement(elP,0.022);
  fLeg_lower -> AddElement(elS,0.002);
  fLeg_lower -> AddElement(elCl,0.001);
  fLeg_lower -> AddElement(elCa,0.046); 
  
  //Foot bones, spongiosa ->ID 12
  d = 1.117 *g/cm3;
  fFoot = new G4Material("foot",d,9);
  fFoot -> AddElement(elH,0.096);
  fFoot -> AddElement(elC,0.473);
  fFoot -> AddElement(elN,0.017);
  fFoot -> AddElement(elO,0.341);
  fFoot -> AddElement(elNa,0.002);
  fFoot -> AddElement(elP,0.022);
  fFoot -> AddElement(elS,0.002);
  fFoot -> AddElement(elCl,0.001);
  fFoot -> AddElement(elCa,0.046); 
  
  //Mandible, spongiosa -> ID 13
  d = 1.189 *g/cm3;
  fMandible = new G4Material("mandible",d,11);
  fMandible -> AddElement(elH,0.087);
  fMandible -> AddElement(elC,0.357);
  fMandible -> AddElement(elN,0.026);
  fMandible -> AddElement(elO,0.429);
  fMandible -> AddElement(elNa,0.002);
  fMandible -> AddElement(elMg,0.001);
  fMandible -> AddElement(elP,0.030);
  fMandible -> AddElement(elS,0.003);
  fMandible -> AddElement(elCl,0.001);
  fMandible -> AddElement(elK,0.001);
  fMandible -> AddElement(elCa,0.063);
  
  //Pelvis, Spongiosa -> ID 14
  d = 1.109 *g/cm3;
  fPelvis = new G4Material("pelvis",d,10);
  fPelvis -> AddElement(elH,0.096);
  fPelvis -> AddElement(elC,0.406);
  fPelvis -> AddElement(elN,0.025);
  fPelvis -> AddElement(elO,0.412);
  fPelvis -> AddElement(elNa,0.001);
  fPelvis -> AddElement(elP,0.018);
  fPelvis -> AddElement(elS,0.002);
  fPelvis -> AddElement(elCl,0.001);
  fPelvis -> AddElement(elK,0.001);
  fPelvis -> AddElement(elCa,0.038);
  
  //Ribs, spongiosa -> ID 15
  d = 1.092 *g/cm3;
  fRibs = new G4Material("ribs",d,11);
  fRibs -> AddElement(elH,0.097);
  fRibs -> AddElement(elC,0.381);
  fRibs -> AddElement(elN,0.028);
  fRibs -> AddElement(elO,0.445);
  fRibs -> AddElement(elNa,0.001);
  fRibs -> AddElement(elP,0.014);
  fRibs -> AddElement(elS,0.002);
  fRibs -> AddElement(elCl,0.002);
  fRibs -> AddElement(elK,0.001);
  fRibs -> AddElement(elCa,0.028);
  fRibs -> AddElement(elFe,0.001);
  
  //Scapulae, spongiosa -> ID 16
  d = 1.128 *g/cm3;
  fScapulae = new G4Material("scapulae",d,10);
  fScapulae -> AddElement(elH,0.094);
  fScapulae -> AddElement(elC,0.406);
  fScapulae -> AddElement(elN,0.024);
  fScapulae -> AddElement(elO,0.404);
  fScapulae -> AddElement(elNa,0.001);
  fScapulae -> AddElement(elP,0.022);
  fScapulae -> AddElement(elS,0.002);
  fScapulae -> AddElement(elCl,0.001);
  fScapulae -> AddElement(elK,0.001);
  fScapulae -> AddElement(elCa,0.045);
  
  //Cervical Spine, spongiosa -> ID 17
  d = 1.135 *g/cm3;
  fSpine_cervical = new G4Material("spine_cervical",d,10);
  fSpine_cervical -> AddElement(elH,0.092);
  fSpine_cervical -> AddElement(elC,0.351);
  fSpine_cervical -> AddElement(elN,0.029);
  fSpine_cervical -> AddElement(elO,0.458);
  fSpine_cervical -> AddElement(elNa,0.001);
  fSpine_cervical -> AddElement(elP,0.021);
  fSpine_cervical -> AddElement(elS,0.002);
  fSpine_cervical -> AddElement(elCl,0.002);
  fSpine_cervical -> AddElement(elK,0.002);
  fSpine_cervical -> AddElement(elCa,0.043);
  
  //Thoratic Spine, spongiosa -> ID 18
  d = 1.084 *g/cm3;
  fSpine_thoratic = new G4Material("spine_thoratic",d,11);
  fSpine_thoratic -> AddElement(elH,0.098);
  fSpine_thoratic -> AddElement(elC,0.386);
  fSpine_thoratic -> AddElement(elN,0.028);
  fSpine_thoratic -> AddElement(elO,0.442);
  fSpine_thoratic -> AddElement(elNa,0.001);
  fSpine_thoratic -> AddElement(elP,0.013);
  fSpine_thoratic -> AddElement(elS,0.002);
  fSpine_thoratic -> AddElement(elCl,0.002);
  fSpine_thoratic -> AddElement(elK,0.001);
  fSpine_thoratic -> AddElement(elCa,0.026);
  fSpine_thoratic -> AddElement(elFe,0.001);
  
  //Lumbar Spine, spongiosa -> ID 19
  d = 1.171 *g/cm3;
  fSpine_lumbar = new G4Material("spine_lumbar",d,11);
  fSpine_lumbar -> AddElement(elH,0.088);
  fSpine_lumbar -> AddElement(elC,0.329);
  fSpine_lumbar -> AddElement(elN,0.030);
  fSpine_lumbar -> AddElement(elO,0.466);
  fSpine_lumbar -> AddElement(elNa,0.001);
  fSpine_lumbar -> AddElement(elMg,0.001);
  fSpine_lumbar -> AddElement(elP,0.026);
  fSpine_lumbar -> AddElement(elS,0.003);
  fSpine_lumbar -> AddElement(elCl,0.001);
  fSpine_lumbar -> AddElement(elK,0.001);
  fSpine_lumbar -> AddElement(elCa,0.054);
  
  //Sacrum, spongiosa -> ID 20
  d = 1.052 *g/cm3;
  fSacrum = new G4Material("sacrum",d,11);
  fSacrum -> AddElement(elH,0.102);
  fSacrum -> AddElement(elC,0.410);
  fSacrum -> AddElement(elN,0.027);
  fSacrum -> AddElement(elO,0.433);
  fSacrum -> AddElement(elNa,0.001);
  fSacrum -> AddElement(elP,0.007);
  fSacrum -> AddElement(elS,0.002);
  fSacrum -> AddElement(elCl,0.002);
  fSacrum -> AddElement(elK,0.001);
  fSacrum -> AddElement(elCa,0.014);
  fSacrum -> AddElement(elFe,0.001);
  
  //Sternum, spongiosa -> ID 21
  d = 1.076 *g/cm3;
  fSternum = new G4Material("sternum",d,11);
  fSternum -> AddElement(elH,0.099);
  fSternum -> AddElement(elC,0.392);
  fSternum -> AddElement(elN,0.028);
  fSternum -> AddElement(elO,0.439);
  fSternum -> AddElement(elNa,0.001);
  fSternum -> AddElement(elP,0.012);
  fSternum -> AddElement(elS,0.002);
  fSternum -> AddElement(elCl,0.002);
  fSternum -> AddElement(elK,0.001);
  fSternum -> AddElement(elCa,0.023);
  fSternum -> AddElement(elFe,0.001);
  
  //Humeri and femora, upper halves, medullary cavity -> ID 22
  d = 0.980 *g/cm3;
  fHf_upper = new G4Material("hf_upper",d,7);
  fHf_upper -> AddElement(elH,0.115);
  fHf_upper -> AddElement(elC,0.637);
  fHf_upper -> AddElement(elN,0.007);
  fHf_upper -> AddElement(elO,0.238);
  fHf_upper -> AddElement(elNa,0.001);
  fHf_upper -> AddElement(elS,0.001);
  fHf_upper -> AddElement(elCl,0.001);
  
  //Humeri and femora, lower halves, medullary cavity -> ID 23
  d = 0.980 *g/cm3;
  fHf_lower = new G4Material("hf_lower",d,7);
  fHf_lower -> AddElement(elH,0.115);
  fHf_lower -> AddElement(elC,0.637);
  fHf_lower -> AddElement(elN,0.007);
  fHf_lower -> AddElement(elO,0.238);
  fHf_lower -> AddElement(elNa,0.001);
  fHf_lower -> AddElement(elS,0.001);
  fHf_lower -> AddElement(elCl,0.001);
  
  //Lower arm bones, medullary cavity -> ID 24
  d = 0.980 *g/cm3;
  fMed_lowerarm = new G4Material("med_lowerarm",d,7);
  fMed_lowerarm -> AddElement(elH,0.115);
  fMed_lowerarm -> AddElement(elC,0.637);
  fMed_lowerarm -> AddElement(elN,0.007);
  fMed_lowerarm -> AddElement(elO,0.238);
  fMed_lowerarm -> AddElement(elNa,0.001);
  fMed_lowerarm -> AddElement(elS,0.001);
  fMed_lowerarm -> AddElement(elCl,0.001);
  
  //Lower leg bones, medullary cavity -> ID 25
  d = 0.980 *g/cm3;
  fMed_lowerleg = new G4Material("med_lowerleg",d,7);
  fMed_lowerleg -> AddElement(elH,0.115);
  fMed_lowerleg -> AddElement(elC,0.637);
  fMed_lowerleg -> AddElement(elN,0.007);
  fMed_lowerleg -> AddElement(elO,0.238);
  fMed_lowerleg -> AddElement(elNa,0.001);
  fMed_lowerleg -> AddElement(elS,0.001);
  fMed_lowerleg -> AddElement(elCl,0.001);
  
  //Cartilage -> ID 26
  d = 1.100 *g/cm3;
  fCartilage = new G4Material("cartilage",d,8);
  fCartilage -> AddElement(elH,0.096);
  fCartilage -> AddElement(elC,0.099);
  fCartilage -> AddElement(elN,0.022);
  fCartilage -> AddElement(elO,0.744);
  fCartilage -> AddElement(elNa,0.005);
  fCartilage -> AddElement(elP,0.022);
  fCartilage -> AddElement(elS,0.009);
  fCartilage -> AddElement(elCl,0.003);
  
  //Skin -> Id 27
  d = 1.090 *g/cm3;
  fSkin = new G4Material("skin",d,9);
  fSkin -> AddElement(elH,0.100);
  fSkin -> AddElement(elC,0.199);
  fSkin -> AddElement(elN,0.042);
  fSkin -> AddElement(elO,0.650);
  fSkin -> AddElement(elNa,0.002);
  fSkin -> AddElement(elP,0.001);
  fSkin -> AddElement(elS,0.002);
  fSkin -> AddElement(elCl,0.003);
  fSkin -> AddElement(elK,0.001);
  
  //Blood -> ID 28
  d = 1.060 *g/cm3;
  fBlood = new G4Material("blood",d,10);
  fBlood -> AddElement(elH,0.102);
  fBlood -> AddElement(elC,0.110);
  fBlood -> AddElement(elN,0.033);
  fBlood -> AddElement(elO,0.745);
  fBlood -> AddElement(elNa,0.001);
  fBlood -> AddElement(elP,0.001);
  fBlood -> AddElement(elS,0.002);
  fBlood -> AddElement(elCl,0.003);
  fBlood -> AddElement(elK,0.002);
  fBlood -> AddElement(elFe,0.001);
  
  //Muscular Tissue -> ID 29
  d = 1.050 *g/cm3;
  fMuscle = new G4Material("muscle",d,9);
  fMuscle -> AddElement(elH,0.102);
  fMuscle -> AddElement(elC,0.142);
  fMuscle -> AddElement(elN,0.034);
  fMuscle -> AddElement(elO,0.711);
  fMuscle -> AddElement(elNa,0.001);
  fMuscle -> AddElement(elP,0.002);
  fMuscle -> AddElement(elS,0.003);
  fMuscle -> AddElement(elCl,0.001);
  fMuscle -> AddElement(elK,0.004);
  
  //Liver -> ID 30
  d = 1.050 *g/cm3;
  fLiver = new G4Material("liver",d,9);
  fLiver -> AddElement(elH,0.102);
  fLiver -> AddElement(elC,0.131);
  fLiver -> AddElement(elN,0.031);
  fLiver -> AddElement(elO,0.724);
  fLiver -> AddElement(elNa,0.002);
  fLiver -> AddElement(elP,0.002);
  fLiver -> AddElement(elS,0.003);
  fLiver -> AddElement(elCl,0.002);
  fLiver -> AddElement(elK,0.003);
  
  //Pancreas ->ID 31
  d = 1.050 *g/cm3;
  fPancreas = new G4Material("pancreas",d,9);
  fPancreas -> AddElement(elH,0.105);
  fPancreas -> AddElement(elC,0.157);
  fPancreas -> AddElement(elN,0.024);
  fPancreas -> AddElement(elO,0.705);
  fPancreas -> AddElement(elNa,0.002);
  fPancreas -> AddElement(elP,0.002);
  fPancreas -> AddElement(elS,0.001);
  fPancreas -> AddElement(elCl,0.002);
  fPancreas -> AddElement(elK,0.002);
  
  //Brain -> ID 32
  d = 1.050 *g/cm3;
  fBrain = new G4Material("brain",d,9);
  fBrain -> AddElement(elH,0.107);
  fBrain -> AddElement(elC,0.144);
  fBrain -> AddElement(elN,0.022);
  fBrain -> AddElement(elO,0.713);
  fBrain -> AddElement(elNa,0.002);
  fBrain -> AddElement(elP,0.004);
  fBrain -> AddElement(elS,0.002);
  fBrain -> AddElement(elCl,0.003);
  fBrain -> AddElement(elK,0.003);
  
  //Heart -> ID 33
  d = 1.050 *g/cm3;
  fHeart = new G4Material("heart",d,9);
  fHeart -> AddElement(elH,0.104);
  fHeart -> AddElement(elC,0.138);
  fHeart -> AddElement(elN,0.029);
  fHeart -> AddElement(elO,0.719);
  fHeart -> AddElement(elNa,0.001);
  fHeart -> AddElement(elP,0.002);
  fHeart -> AddElement(elS,0.002);
  fHeart -> AddElement(elCl,0.002);
  fHeart -> AddElement(elK,0.003);
  
  //Eye ->ID 34
  d = 1.050 *g/cm3;
  fEye = new G4Material("eye",d,8);
  fEye -> AddElement(elH,0.097);
  fEye -> AddElement(elC,0.183);
  fEye -> AddElement(elN,0.054);
  fEye -> AddElement(elO,0.660);
  fEye -> AddElement(elCa,0.001);
  fEye -> AddElement(elP,0.001);
  fEye -> AddElement(elS,0.003);
  fEye -> AddElement(elCl,0.001);
  
  //Kidneys -> ID 35
  d = 1.050 *g/cm3;
  fKidney = new G4Material("kidney",d,10);
  fKidney -> AddElement(elH,0.103);
  fKidney -> AddElement(elC,0.125);
  fKidney -> AddElement(elN,0.031);
  fKidney -> AddElement(elO,0.730);
  fKidney -> AddElement(elNa,0.002);
  fKidney -> AddElement(elP,0.002);
  fKidney -> AddElement(elS,0.002);
  fKidney -> AddElement(elCl,0.002);
  fKidney -> AddElement(elK,0.002);
  fKidney -> AddElement(elCa,0.001);
  
  //Stomach ->ID 36
  d = 1.040 *g/cm3;
  fStomach = new G4Material("stomach",d,9);
  fStomach -> AddElement(elH,0.105);
  fStomach -> AddElement(elC,0.114);
  fStomach -> AddElement(elN,0.025);
  fStomach -> AddElement(elO,0.750);
  fStomach -> AddElement(elNa,0.001);
  fStomach -> AddElement(elP,0.001);
  fStomach -> AddElement(elS,0.001);
  fStomach -> AddElement(elCl,0.002);
  fStomach -> AddElement(elK,0.001);
  
  //Small intestine ->ID 37
  d = 1.040 *g/cm3;
  fIntestine_sml = new G4Material("intestine_sml",d,9);
  fIntestine_sml -> AddElement(elH,0.105);
  fIntestine_sml -> AddElement(elC,0.114);
  fIntestine_sml -> AddElement(elN,0.025);
  fIntestine_sml -> AddElement(elO,0.750);
  fIntestine_sml -> AddElement(elNa,0.001);
  fIntestine_sml -> AddElement(elP,0.001);
  fIntestine_sml -> AddElement(elS,0.001);
  fIntestine_sml -> AddElement(elCl,0.002);
  fIntestine_sml -> AddElement(elK,0.001);
  
  //Large intestine ->ID 38
  d = 1.040 *g/cm3;
  fIntestine_lrg = new G4Material("intestine_lrg",d,9);
  fIntestine_lrg -> AddElement(elH,0.105);
  fIntestine_lrg -> AddElement(elC,0.114);
  fIntestine_lrg -> AddElement(elN,0.025);
  fIntestine_lrg -> AddElement(elO,0.750);
  fIntestine_lrg -> AddElement(elNa,0.001);
  fIntestine_lrg -> AddElement(elP,0.001);
  fIntestine_lrg -> AddElement(elS,0.001);
  fIntestine_lrg -> AddElement(elCl,0.002);
  fIntestine_lrg -> AddElement(elK,0.001);
  
  //Spleen -> ID 39
  d = 1.040 *g/cm3;
  fSpleen = new G4Material("spleen",d,9);
  fSpleen -> AddElement(elH,0.103);
  fSpleen -> AddElement(elC,0.112);
  fSpleen -> AddElement(elN,0.032);
  fSpleen -> AddElement(elO,0.743);
  fSpleen -> AddElement(elNa,0.001);
  fSpleen -> AddElement(elP,0.002);
  fSpleen -> AddElement(elS,0.002);
  fSpleen -> AddElement(elCl,0.002);
  fSpleen -> AddElement(elK,0.003);

  //Thyroid -> ID 40
  d = 1.040 *g/cm3;
  fThyroid = new G4Material("thyroid",d,10);
  fThyroid -> AddElement(elH,0.104);
  fThyroid -> AddElement(elC,0.118);
  fThyroid -> AddElement(elN,0.025);
  fThyroid -> AddElement(elO,0.745);
  fThyroid -> AddElement(elNa,0.002);
  fThyroid -> AddElement(elP,0.001);
  fThyroid -> AddElement(elS,0.001);
  fThyroid -> AddElement(elCl,0.002);
  fThyroid -> AddElement(elK,0.001);
  fThyroid -> AddElement(elI,0.001);
  
  //Urinary Bladder -> ID 41
  d = 1.040 *g/cm3;
  fBladder = new G4Material("bladder",d,9);
  fBladder -> AddElement(elH,0.105);
  fBladder -> AddElement(elC,0.096);
  fBladder -> AddElement(elN,0.026);
  fBladder -> AddElement(elO,0.761);
  fBladder -> AddElement(elNa,0.002);
  fBladder -> AddElement(elP,0.002);
  fBladder -> AddElement(elS,0.002);
  fBladder -> AddElement(elCl,0.003);
  fBladder -> AddElement(elK,0.003);
  
  //  Ovaries (Defined as ovaries_testes to be called upon by ColourMap.dat 
  //for visualisation purposes) -> ID 42
  d = 1.040 *g/cm3;
  fOvaries_testes = new G4Material("ovaries_testes",d,9);
  fOvaries_testes -> AddElement(elH,0.105);
  fOvaries_testes -> AddElement(elC,0.094);
  fOvaries_testes -> AddElement(elN,0.025);
  fOvaries_testes -> AddElement(elO,0.766);
  fOvaries_testes -> AddElement(elNa,0.002);
  fOvaries_testes -> AddElement(elP,0.002);
  fOvaries_testes -> AddElement(elS,0.002);
  fOvaries_testes -> AddElement(elCl,0.002);
  fOvaries_testes -> AddElement(elK,0.002);
  
  //Adrenals -> ID 43
  d = 1.030 *g/cm3;
  fAdrenals = new G4Material("adrenals",d,9);
  fAdrenals -> AddElement(elH,0.104);
  fAdrenals -> AddElement(elC,0.228);
  fAdrenals -> AddElement(elN,0.028);
  fAdrenals -> AddElement(elO,0.630);
  fAdrenals -> AddElement(elNa,0.001);
  fAdrenals -> AddElement(elP,0.002);
  fAdrenals -> AddElement(elS,0.003);
  fAdrenals -> AddElement(elCl,0.002);
  fAdrenals -> AddElement(elK,0.002);
  
  //Oesophagus -> ID 44
  d = 1.030 *g/cm3;
  fOesophagus = new G4Material("oesophagus",d,9);
  fOesophagus -> AddElement(elH,0.104);
  fOesophagus -> AddElement(elC,0.222);
  fOesophagus -> AddElement(elN,0.028);
  fOesophagus -> AddElement(elO,0.636);
  fOesophagus -> AddElement(elNa,0.001);
  fOesophagus -> AddElement(elP,0.002);
  fOesophagus -> AddElement(elS,0.003);
  fOesophagus -> AddElement(elCl,0.002);
  fOesophagus -> AddElement(elK,0.002);
  
  //Miscillaneous (Gallbladder, Trachea, Thymus, Tonsils, Ureters, ...) -> ID 45
  d = 1.030 *g/cm3;
  fMisc = new G4Material("misc",d,9);
  fMisc -> AddElement(elH,0.105);
  fMisc -> AddElement(elC,0.235);
  fMisc -> AddElement(elN,0.028);
  fMisc -> AddElement(elO,0.622);
  fMisc -> AddElement(elNa,0.001);
  fMisc -> AddElement(elP,0.002);
  fMisc -> AddElement(elS,0.003);
  fMisc -> AddElement(elCl,0.002);
  fMisc -> AddElement(elK,0.002);
  
  //Uterus (Defined as fUterus_prostate to be called upon by ColourMap.dat 
  // for visualisation purposes) -> ID 46
  d = 1.030 *g/cm3;
  fUterus_prostate = new G4Material("uterus_prostate",d,9);
  fUterus_prostate -> AddElement(elH,0.105);
  fUterus_prostate -> AddElement(elC,0.286);
  fUterus_prostate -> AddElement(elN,0.025);
  fUterus_prostate -> AddElement(elO,0.576);
  fUterus_prostate -> AddElement(elNa,0.001);
  fUterus_prostate -> AddElement(elP,0.002);
  fUterus_prostate -> AddElement(elS,0.002);
  fUterus_prostate -> AddElement(elCl,0.001);
  fUterus_prostate -> AddElement(elK,0.002); 
  
  //Lymph -> ID 47
  d = 1.030 *g/cm3;
  fLymph = new G4Material("lymph",d,7);
  fLymph -> AddElement(elH,0.108);
  fLymph -> AddElement(elC,0.042);
  fLymph -> AddElement(elN,0.011);
  fLymph -> AddElement(elO,0.831);
  fLymph -> AddElement(elNa,0.003);
  fLymph -> AddElement(elS,0.001);
  fLymph -> AddElement(elCl,0.004);

  //Breast (Mammary Gland)-> ID 48
  d = 1.020 *g/cm3;
  fBreast_glandular = new G4Material("breast_glandular",d,4);
  fBreast_glandular -> AddElement(elH,0.114);
  fBreast_glandular -> AddElement(elC,0.461);
  fBreast_glandular -> AddElement(elN,0.005);
  fBreast_glandular -> AddElement(elO,0.420);

  //Adipose tissue (breast) -> ID 49
  d = 0.950 *g/cm3;
  fBreast_adipose = new G4Material("breast_adipose",d,7);
  fBreast_adipose -> AddElement(elH,0.114);
  fBreast_adipose -> AddElement(elC,0.589);
  fBreast_adipose -> AddElement(elN,0.007);
  fBreast_adipose -> AddElement(elO,0.287);
  fBreast_adipose -> AddElement(elNa,0.001);
  fBreast_adipose -> AddElement(elS,0.001);
  fBreast_adipose -> AddElement(elCl,0.001);
  
  //Lung Tissue (Compressed Lung) -> ID 50
  d = 0.385 *g/cm3;
  fLung = new G4Material("lung",d,9);
  fLung -> AddElement(elH,0.103);
  fLung -> AddElement(elC,0.107);
  fLung -> AddElement(elN,0.032);
  fLung -> AddElement(elO,0.746);
  fLung -> AddElement(elNa,0.002);
  fLung -> AddElement(elP,0.002);
  fLung -> AddElement(elS,0.003);
  fLung -> AddElement(elCl,0.003);
  fLung -> AddElement(elK,0.002);
  
  //Contents of gastro-intestinal tract -> ID 51
  d = 1.040 *g/cm3;
  fGastro_content = new G4Material("gastro_content",d,10);
  fGastro_content -> AddElement(elH,0.100);
  fGastro_content -> AddElement(elC,0.222);
  fGastro_content -> AddElement(elN,0.022);
  fGastro_content -> AddElement(elO,0.644);
  fGastro_content -> AddElement(elNa,0.001);
  fGastro_content -> AddElement(elP,0.002);
  fGastro_content -> AddElement(elS,0.003);
  fGastro_content -> AddElement(elCl,0.001);
  fGastro_content -> AddElement(elK,0.004);
  fGastro_content -> AddElement(elCa,0.001);
  
  //Urine -> ID 52
  d = 1.040 *g/cm3;
  fUrine = new G4Material("urine",d,7);
  fUrine -> AddElement(elH,0.107);
  fUrine -> AddElement(elC,0.003);
  fUrine -> AddElement(elN,0.010);
  fUrine -> AddElement(elO,0.873);
  fUrine -> AddElement(elNa,0.004);
  fUrine -> AddElement(elP,0.001);
  fUrine -> AddElement(elK,0.002);

}

G4Material* ICRP110PhantomMaterial_Female::GetMaterial(G4String material)
{
  // Returns a material 
  G4Material* pttoMaterial = G4Material::GetMaterial(material); 
  if (!pttoMaterial) G4cout << "WARNING: material '" << material << "' is not defined!" << G4endl;
  return pttoMaterial; 
}
