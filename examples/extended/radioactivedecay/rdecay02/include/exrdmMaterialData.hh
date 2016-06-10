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
// $Id: exrdmMaterialData.hh 68007 2013-03-13 11:28:03Z gcosmo $
//
/// \file radioactivedecay/rdecay02/include/exrdmMaterialData.hh
/// \brief Definition of the exrdmMaterialData class
//

#include "G4Material.hh"

const G4String exrdmMaterial::fELU[110] =
{
  " H","He","Li","Be"," B"," C"," N"," O"," F","Ne",
  "Na","Mg","Al","Si"," P"," S","Cl","Ar"," K","Ca",
  "Sc","Ti"," V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr","Rb","Sr"," Y","Zr",
  "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
  "Sb","Te"," I","Xe","Cs","Ba","La","Ce","Pr","Nd",
  "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
  "Lu","Hf","Ta"," W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
  "Pa"," U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
  "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","UN"
};
const G4String exrdmMaterial::fELL[110] =
{
  " h","he","li","be"," b"," c"," n"," o"," f","ne",
  "na","mg","al","si"," p"," s","cl","ar"," k","ca",
  "sc","ti"," v","cr","mn","fe","co","ni","cu","zn",
  "ga","ge","as","se","br","kr","rb","sr"," y","zr",
  "nb","mo","tc","ru","rh","pd","ag","cd","in","sn",
  "sb","te"," i","xe","cs","ba","la","ce","pr","nd",
  "pm","sm","eu","gd","tb","dy","ho","er","tm","yb",
  "lu","hf","ta"," w","re","os","ir","pt","au","hg",
  "tl","pb","bi","po","at","rn","fr","ra","ac","th",
  "pa"," u","np","pu","am","cm","bk","cf","es","fm",
  "md","no","lr","rf","db","sg","bh","hs","mt","un"
};
const G4String exrdmMaterial::fEUU[110] =
{
  " H","HE","LI","BE"," B"," C"," N"," O"," F","NE",
  "NA","MG","AL","SI"," P"," S","CL","AR"," K","CA",
  "SC","TI"," V","CR","MN","FE","CO","NI","CU","ZN",
  "GA","GE","AS","SE","BR","KR","RB","SR"," Y","ZR",
  "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN",
  "SB","TE"," I","XE","CS","BA","LA","CE","PR","ND",
  "PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
  "LU","HF","TA"," W","RE","OS","IR","PT","AU","HG",
  "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
  "PA"," U","NP","PU","AM","CM","BK","CF","ES","FM",
  "MD","NO","LR","RF","DB","SG","BH","HS","MT","UN"
};

const G4double exrdmMaterial::fA[110] =
{
  1.00794,4.002602,6.941,9.012182,10.811,12.0107,14.00674,15.9994,18.9984,
  20.1797,22.98977,24.3050,26.9815,28.0855,30.973761,32.066,35.4527,39.948,
  39.0983,40.078,44.95591,47.867,50.9415,51.9961,54.938049,55.845,58.9332,
  58.6934,63.546,65.39,69.723,72.61,74.9216,78.96,79.904,83.8,85.4678,87.62,
  88.90585,91.224,92.90638,95.94,97.9072,101.07,102.9055,106.42,107.8682,
  112.411,114.818,118.71,121.76,127.6,126.90447,131.29,132.90545,137.327,
  138.9055,140.116,140.90765,144.24,144.9127,150.36,151.964,157.25,158.92534,
  162.5,164.93032,167.26,168.93421,173.04,174.967,178.49,180.9479,183.84,
  186.207,190.23,192.217,195.078,196.96655,200.59,
  204.3833,207.2,208.98038,208.9824,209.9871,222.0176,223.0197,226.0254,
  227.0277,232.038,231.03588,238.0289,237.0482,244.0642,243.0614,247.0703,
  247.0703,251.0796,252.083,257.0951,258.0984,259.1011,262.11,263.1125,262.1144,
  266.1219,264.1247,269.1341,268.1388,272.1463
};
