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
// $Id: G4GoudsmitSaundersonTable.cc 94933 2015-12-18 09:22:52Z gcosmo $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class implementation file
//
// File name:     G4GoudsmitSaundersonTable
//
// Author:        Mihaly Novak / (Omrane Kadri)
//
// Creation date: 20.02.2009
//
// Class description:
//   Class to handle multiple scattering angular distributions precomputed by
//   using Kawrakow-Bielajew Goudsmit-Saunderson MSC model based on the screened
//   Rutherford DCS for elastic scattering of electrons/positrons [1,2]. This
//   class is used by G4GoudsmitSaundersonMscModel to sample the angular
//   deflection of electrons/positrons after travelling a given path.
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
// 26.08.2009 O.Kadri: avoiding unuseful calculations and optimizing the root
//                     finding parameter error's within SampleTheta method
// 08.02.2010 O.Kadri: reduce delared variables; reduce error of finding root
//                     in secant method
// 26.03.2010 O.Kadri: minimum of used arrays in computation within the dichotomie
//                     finding method the error was the lowest value of uvalues
// 12.05.2010 O.Kadri: changing of sqrt((b-a)*(b-a)) with fabs(b-a)
// 18.05.2015 M. Novak This class has been completely replaced (only the original
//            class name was kept; class description was also inserted):
//            A new version of Kawrakow-Bielajew Goudsmit-Saunderson MSC model
//            based on the screened Rutherford DCS for elastic scattering of
//            electrons/positrons has been introduced[1,2]. The corresponding MSC
//            angular distributions over a 2D parameter grid have been recomputed
//            and the CDFs are now stored in a variable transformed (smooth) form
//            together with the corresponding rational interpolation parameters.
//            The new version is several times faster, more robust and accurate
//            compared to the earlier version (G4GoudsmitSaundersonMscModel class
//            that use these data has been also completely replaced)
//
// References:
//   [1] A.F.Bielajew, NIMB, 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//
// -----------------------------------------------------------------------------

#include "G4GoudsmitSaundersonTable.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>


const G4double G4GoudsmitSaundersonTable::fgLambdaValues[] ={
    1.000000000000e+00, 1.165914401180e+00, 1.359356390879e+00, 1.584893192461e+00,
    1.847849797422e+00, 2.154434690032e+00, 2.511886431510e+00, 2.928644564625e+00,
    3.414548873834e+00, 3.981071705535e+00, 4.641588833613e+00, 5.411695265465e+00,
    6.309573444802e+00, 7.356422544596e+00, 8.576958985909e+00, 1.000000000000e+01,
    1.165914401180e+01, 1.359356390879e+01, 1.584893192461e+01, 1.847849797422e+01,
    2.154434690032e+01, 2.511886431510e+01, 2.928644564625e+01, 3.414548873834e+01,
    3.981071705535e+01, 4.641588833613e+01, 5.411695265465e+01, 6.309573444802e+01,
    7.356422544596e+01, 8.576958985909e+01, 1.000000000000e+02, 1.165914401180e+02,
    1.359356390879e+02, 1.584893192461e+02, 1.847849797422e+02, 2.154434690032e+02,
    2.511886431510e+02, 2.928644564625e+02, 3.414548873834e+02, 3.981071705535e+02,
    4.641588833613e+02, 5.411695265465e+02, 6.309573444802e+02, 7.356422544596e+02,
    8.576958985909e+02, 1.000000000000e+03, 1.165914401180e+03, 1.359356390879e+03,
    1.584893192461e+03, 1.847849797422e+03, 2.154434690032e+03, 2.511886431510e+03,
    2.928644564625e+03, 3.414548873834e+03, 3.981071705535e+03, 4.641588833613e+03,
    5.411695265465e+03, 6.309573444802e+03, 7.356422544596e+03, 8.576958985909e+03,
    1.000000000000e+04, 1.165914401180e+04, 1.359356390879e+04, 1.584893192461e+04,
    1.847849797422e+04, 2.154434690032e+04, 2.511886431510e+04, 2.928644564625e+04,
    3.414548873834e+04, 3.981071705535e+04, 4.641588833613e+04, 5.411695265465e+04,
    6.309573444802e+04, 7.356422544596e+04, 8.576958985909e+04, 1.000000000000e+05
};

const G4double G4GoudsmitSaundersonTable::fgLamG1Values[]={
       0.0010, 0.0509, 0.1008, 0.1507, 0.2006, 0.2505, 0.3004, 0.3503, 0.4002,
       0.4501, 0.5000, 0.5499, 0.5998, 0.6497, 0.6996, 0.7495, 0.7994, 0.8493,
       0.8992, 0.9491, 0.9990
    };
const G4double G4GoudsmitSaundersonTable::fgLamG1ValuesII[]={
        0.999, 1.332, 1.665, 1.998, 2.331, 2.664, 2.997, 3.33, 3.663,
        3.996, 4.329, 4.662, 4.995, 5.328, 5.661, 5.994, 6.327,
        6.66, 6.993, 7.326, 7.659, 7.992
    };

const G4double G4GoudsmitSaundersonTable::fgUValues[]={
    0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
    0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59,
    0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69,
    0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
    0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
    0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
    1.00
  };

const G4double G4GoudsmitSaundersonTable::fgG1Values[]={
    1.00000000000000e-10, 1.15478198468946e-10, 1.33352143216332e-10, 1.53992652605949e-10, 1.77827941003892e-10,
    2.05352502645715e-10, 2.37137370566166e-10, 2.73841963426436e-10, 3.16227766016838e-10, 3.65174127254838e-10,
    4.21696503428582e-10, 4.86967525165863e-10, 5.62341325190349e-10, 6.49381631576211e-10, 7.49894209332456e-10,
    8.65964323360065e-10, 1.00000000000000e-09, 1.15478198468946e-09, 1.33352143216332e-09, 1.53992652605949e-09,
    1.77827941003892e-09, 2.05352502645715e-09, 2.37137370566166e-09, 2.73841963426436e-09, 3.16227766016838e-09,
    3.65174127254838e-09, 4.21696503428582e-09, 4.86967525165863e-09, 5.62341325190349e-09, 6.49381631576211e-09,
    7.49894209332456e-09, 8.65964323360065e-09, 1.00000000000000e-08, 1.15478198468946e-08, 1.33352143216332e-08,
    1.53992652605949e-08, 1.77827941003892e-08, 2.05352502645715e-08, 2.37137370566166e-08, 2.73841963426436e-08,
    3.16227766016838e-08, 3.65174127254838e-08, 4.21696503428582e-08, 4.86967525165863e-08, 5.62341325190349e-08,
    6.49381631576211e-08, 7.49894209332456e-08, 8.65964323360065e-08, 1.00000000000000e-07, 1.15478198468946e-07,
    1.33352143216332e-07, 1.53992652605949e-07, 1.77827941003892e-07, 2.05352502645715e-07, 2.37137370566166e-07,
    2.73841963426436e-07, 3.16227766016838e-07, 3.65174127254838e-07, 4.21696503428582e-07, 4.86967525165863e-07,
    5.62341325190349e-07, 6.49381631576211e-07, 7.49894209332456e-07, 8.65964323360065e-07, 1.00000000000000e-06,
    1.15478198468946e-06, 1.33352143216332e-06, 1.53992652605949e-06, 1.77827941003892e-06, 2.05352502645715e-06,
    2.37137370566166e-06, 2.73841963426436e-06, 3.16227766016838e-06, 3.65174127254838e-06, 4.21696503428582e-06,
    4.86967525165863e-06, 5.62341325190349e-06, 6.49381631576211e-06, 7.49894209332456e-06, 8.65964323360065e-06,
    1.00000000000000e-05, 1.15478198468946e-05, 1.33352143216332e-05, 1.53992652605949e-05, 1.77827941003892e-05,
    2.05352502645715e-05, 2.37137370566166e-05, 2.73841963426436e-05, 3.16227766016838e-05, 3.65174127254838e-05,
    4.21696503428582e-05, 4.86967525165863e-05, 5.62341325190349e-05, 6.49381631576211e-05, 7.49894209332456e-05,
    8.65964323360065e-05, 1.00000000000000e-04, 1.15478198468946e-04, 1.33352143216332e-04, 1.53992652605949e-04,
    1.77827941003892e-04, 2.05352502645715e-04, 2.37137370566166e-04, 2.73841963426436e-04, 3.16227766016838e-04,
    3.65174127254838e-04, 4.21696503428582e-04, 4.86967525165863e-04, 5.62341325190349e-04, 6.49381631576211e-04,
    7.49894209332456e-04, 8.65964323360065e-04, 1.00000000000000e-03, 1.15478198468946e-03, 1.33352143216332e-03,
    1.53992652605949e-03, 1.77827941003892e-03, 2.05352502645715e-03, 2.37137370566166e-03, 2.73841963426436e-03,
    3.16227766016838e-03, 3.65174127254838e-03, 4.21696503428582e-03, 4.86967525165863e-03, 5.62341325190349e-03,
    6.49381631576211e-03, 7.49894209332456e-03, 8.65964323360065e-03, 1.00000000000000e-02, 1.15478198468946e-02,
    1.33352143216332e-02, 1.53992652605949e-02, 1.77827941003892e-02, 2.05352502645715e-02, 2.37137370566166e-02,
    2.73841963426436e-02, 3.16227766016838e-02, 3.65174127254838e-02, 4.21696503428582e-02, 4.86967525165863e-02,
    5.62341325190349e-02, 6.49381631576211e-02, 7.49894209332456e-02, 8.65964323360065e-02, 1.00000000000000e-01,
    1.15478198468946e-01, 1.33352143216332e-01, 1.53992652605949e-01, 1.77827941003892e-01, 2.05352502645715e-01,
    2.37137370566166e-01, 2.73841963426436e-01, 3.16227766016838e-01, 3.65174127254838e-01, 4.21696503428582e-01,
    4.86967525165863e-01, 5.62341325190349e-01, 6.49381631576211e-01, 7.49894209332456e-01, 8.65964323360065e-01
};

const G4double G4GoudsmitSaundersonTable::fgScreeningParam[]={
    1.92484052135703e-12, 2.23565438533804e-12, 2.59674772669370e-12, 3.01627006395395e-12, 3.50369464193945e-12,
    4.07003394459337e-12, 4.72809038421179e-12, 5.49274792465712e-12, 6.38131134142232e-12, 7.41390092242269e-12,
    8.61391169587491e-12, 1.00085477656079e-11, 1.16294440746525e-11, 1.35133899458129e-11, 1.57031711107699e-11,
    1.82485496926619e-11, 2.12074048158664e-11, 2.46470602564846e-11, 2.86458299060786e-11, 3.32948169025090e-11,
    3.87000082054842e-11, 4.49847133009623e-11, 5.22924037716594e-11, 6.07900198618435e-11, 7.06718211166503e-11,
    8.21638709501378e-11, 9.55292598967889e-11, 1.11074189683990e-10, 1.29155060543825e-10, 1.50186727847036e-10,
    1.74652121757794e-10, 2.03113455838333e-10, 2.36225288153009e-10, 2.74749742338471e-10, 3.19574247380381e-10,
    3.71732214707185e-10, 4.32427141127727e-10, 5.03060707798248e-10, 5.85265540790081e-10, 6.80943410264647e-10,
    7.92309775465429e-10, 9.21945734889531e-10, 1.07285861883013e-09, 1.24855266934879e-09, 1.45311149575389e-09,
    1.69129427781576e-09, 1.96864802125653e-09, 2.29163855873534e-09, 2.66780344424806e-09, 3.10593042087591e-09,
    3.61626576440636e-09, 4.21075753406168e-09, 4.90333961465076e-09, 5.71026343332307e-09, 6.65048540388811e-09,
    7.74611952188878e-09, 9.02296613895651e-09, 1.05111298261663e-08, 1.22457414410270e-08, 1.42678020976703e-08,
    1.66251697709241e-08, 1.93737128201322e-08, 2.25786588894194e-08, 2.63161725354421e-08, 3.06752006784906e-08,
    3.57596317176864e-08, 4.16908220722131e-08, 4.86105532157964e-08, 5.66844932060409e-08, 6.61062495628185e-08,
    7.71021154618766e-08, 8.99366289841022e-08, 1.04919086073377e-07, 1.22411172469263e-07, 1.42835908860059e-07,
    1.66688137634119e-07, 1.94546819824371e-07, 2.27089458246136e-07, 2.65109018729386e-07, 3.09533787294144e-07,
    3.61450678951755e-07, 4.22132605720001e-07, 4.93070620012029e-07, 5.76011677884208e-07, 6.73003018378312e-07,
    7.86444334741884e-07, 9.19149125868381e-07, 1.07441686808177e-06, 1.25611794581916e-06, 1.46879363370693e-06,
    1.71777384258668e-06, 2.00931584092477e-06, 2.35076775595475e-06, 2.75076136411589e-06, 3.21943951980544e-06,
    3.76872457154335e-06, 4.41263530713955e-06, 5.16766139268237e-06, 6.05320597042372e-06, 7.09210911391032e-06,
    8.31126727283062e-06, 9.74236675732094e-06, 1.14227528119565e-05, 1.33964600351938e-05, 1.57154349593153e-05,
    1.84409877007598e-05, 2.16455169438847e-05, 2.54145614063097e-05, 2.98492416879238e-05, 3.50691694441802e-05,
    4.12159166620983e-05, 4.84571570931433e-05, 5.69916154060348e-05, 6.70549883575744e-05, 7.89270374851336e-05,
    9.29400960653008e-05, 1.09489286334512e-04, 1.29044808732337e-04, 1.52166746391861e-04, 1.79522929336019e-04,
    2.11910529073029e-04, 2.50282212267455e-04, 2.95777880652745e-04, 3.49763274771614e-04, 4.13877036477045e-04,
    4.90088229201730e-04, 5.80766832133355e-04, 6.88770389861787e-04, 8.17550860329037e-04, 9.71286825640136e-04,
    1.15504770108335e-03, 1.37499852012182e-03, 1.63865645831619e-03, 1.95521372859357e-03, 2.33594617839223e-03,
    2.79473334292108e-03, 3.34872458412049e-03, 4.01919834698563e-03, 4.83267910882034e-03, 5.82240174726978e-03,
    7.03024963447929e-03, 8.50934682352796e-03, 1.03275659802088e-02, 1.25723383030196e-02, 1.53573467148336e-02,
    1.88319961908362e-02, 2.31950693494851e-02, 2.87148467953825e-02, 3.57594981860148e-02, 4.48443278665925e-02,
    5.67077407341705e-02, 7.24383639141981e-02, 9.36982315562396e-02, 1.23138322807982e-01, 1.65231234966025e-01,
    2.28105620915395e-01, 3.28136673526381e-01, 5.03725917697149e-01, 8.70438998827130e-01, 2.00702876580162e+00
};

const G4double G4GoudsmitSaundersonTable::fgSrcAValues[] = {
    1.04015860967600e+00, 1.04040150744807e+00, 1.04064741953810e+00, 1.04089640311900e+00, 1.04114851683194e+00,
    1.04140382083415e+00, 1.04166237684853e+00, 1.04192424821535e+00, 1.04218949994586e+00, 1.04245819877822e+00,
    1.04273041323563e+00, 1.04300621368681e+00, 1.04328567240902e+00, 1.04356886365371e+00, 1.04385586371482e+00,
    1.04414675100008e+00, 1.04444160610522e+00, 1.04474051189143e+00, 1.04504355356608e+00, 1.04535081876704e+00,
    1.04566239765056e+00, 1.04597838298311e+00, 1.04629887023725e+00, 1.04662395769178e+00, 1.04695374653644e+00,
    1.04728834098127e+00, 1.04762784837111e+00, 1.04797237930520e+00, 1.04832204776253e+00, 1.04867697123289e+00,
    1.04903727085420e+00, 1.04940307155639e+00, 1.04977450221209e+00, 1.05015169579471e+00, 1.05053478954418e+00,
    1.05092392514088e+00, 1.05131924888815e+00, 1.05172091190400e+00, 1.05212907032243e+00, 1.05254388550507e+00,
    1.05296552426367e+00, 1.05339415909403e+00, 1.05382996842239e+00, 1.05427313686456e+00, 1.05472385549904e+00,
    1.05518232215471e+00, 1.05564874171422e+00, 1.05612332643389e+00, 1.05660629628136e+00, 1.05709787929208e+00,
    1.05759831194585e+00, 1.05810783956480e+00, 1.05862671673423e+00, 1.05915520774789e+00, 1.05969358707925e+00,
    1.06024213988081e+00, 1.06080116251315e+00, 1.06137096310598e+00, 1.06195186215339e+00, 1.06254419314585e+00,
    1.06314830324155e+00, 1.06376455398007e+00, 1.06439332204139e+00, 1.06503500005380e+00, 1.06568999745437e+00,
    1.06635874140595e+00, 1.06704167777519e+00, 1.06773927217628e+00, 1.06845201108575e+00, 1.06918040303385e+00,
    1.06992497987890e+00, 1.07068629817131e+00, 1.07146494061474e+00, 1.07226151763259e+00, 1.07307666904867e+00,
    1.07391106589198e+00, 1.07476541233634e+00, 1.07564044778664e+00, 1.07653694912487e+00, 1.07745573313039e+00,
    1.07839765909007e+00, 1.07936363161610e+00, 1.08035460369071e+00, 1.08137157995935e+00, 1.08241562029607e+00,
    1.08348784366764e+00, 1.08458943232554e+00, 1.08572163635876e+00, 1.08688577864359e+00, 1.08808326023102e+00,
    1.08931556621724e+00, 1.09058427214773e+00, 1.09189105101191e+00, 1.09323768089211e+00, 1.09462605333836e+00,
    1.09605818254966e+00, 1.09753621545270e+00, 1.09906244278036e+00, 1.10063931126624e+00, 1.10226943708649e+00,
    1.10395562069822e+00, 1.10570086324418e+00, 1.10750838471683e+00, 1.10938164410283e+00, 1.11132436176008e+00,
    1.11334054431726e+00, 1.11543451242832e+00, 1.11761093176532e+00, 1.11987484769228e+00, 1.12223172413234e+00,
    1.12468748722310e+00, 1.12724857445230e+00, 1.12992199008194e+00, 1.13271536780716e+00, 1.13563704176072e+00,
    1.13869612717300e+00, 1.14190261223537e+00, 1.14526746300462e+00, 1.14880274353680e+00, 1.15252175386784e+00,
    1.15643918898390e+00, 1.16057132257199e+00, 1.16493622014373e+00, 1.16955398712386e+00, 1.17444705874537e+00,
    1.17964054016868e+00, 1.18516260723776e+00, 1.19104498083299e+00, 1.19732349105124e+00, 1.20403875167602e+00,
    1.21123697091984e+00, 1.21897093167774e+00, 1.22730118415524e+00, 1.23629750662007e+00, 1.24604070744941e+00,
    1.25662486545524e+00, 1.26816013838519e+00, 1.28077631555756e+00, 1.29462735591193e+00, 1.30989724673401e+00,
    1.32680765564328e+00, 1.34562805256148e+00, 1.36668928751138e+00, 1.39040208794754e+00, 1.41728269492652e+00,
    1.44798908275733e+00, 1.48337325073139e+00, 1.52455859525395e+00, 1.57305765486728e+00, 1.63095721573069e+00,
    1.70122060331290e+00, 1.78820418157127e+00, 1.89858943787341e+00, 2.04318263724065e+00, 2.24070141295433e+00,
    2.52670054341548e+00, 2.97823240973856e+00, 3.80070470423559e+00, 5.80504410098720e+00, 0.0
};

const G4double G4GoudsmitSaundersonTable::fgSrcBValues[] = {
    4.85267101555493e-02, 4.87971711674324e-02, 4.90707875373474e-02, 4.93476163203902e-02, 4.96277159833787e-02,
    4.99111464494480e-02, 5.01979691443985e-02, 5.04882470448502e-02, 5.07820447282697e-02, 5.10794284249799e-02,
    5.13804660722384e-02, 5.16852273704654e-02, 5.19937838417706e-02, 5.23062088908123e-02, 5.26225778681911e-02,
    5.29429681364105e-02, 5.32674591386191e-02, 5.35961324701870e-02, 5.39290719533364e-02, 5.42663637149114e-02,
    5.46080962674657e-02, 5.49543605938995e-02, 5.53052502357015e-02, 5.56608613851103e-02, 5.60212929813128e-02,
    5.63866468109177e-02, 5.67570276129776e-02, 5.71325431886598e-02, 5.75133045160478e-02, 5.78994258700939e-02,
    5.82910249481984e-02, 5.86882230016298e-02, 5.90911449731285e-02, 5.94999196410672e-02, 5.99146797704773e-02,
    6.03355622714283e-02, 6.07627083650703e-02, 6.11962637578703e-02, 6.16363788244832e-02, 6.20832087997576e-02,
    6.25369139805019e-02, 6.29976599373900e-02, 6.34656177379503e-02, 6.39409641809610e-02, 6.44238820432201e-02,
    6.49145603393270e-02, 6.54131945953447e-02, 6.59199871371941e-02, 6.64351473947443e-02, 6.69588922226319e-02,
    6.74914462388500e-02, 6.80330421823362e-02, 6.85839212907651e-02, 6.91443337000157e-02, 6.97145388666155e-02,
    7.02948060148960e-02, 7.08854146105257e-02, 7.14866548622191e-02, 7.20988282536816e-02, 7.27222481079186e-02,
    7.33572401862953e-02, 7.40041433248063e-02, 7.46633101103861e-02, 7.53351076002229e-02, 7.60199180872846e-02,
    7.67181399156741e-02, 7.74301883495570e-02, 7.81564964999089e-02, 7.88975163136540e-02, 7.96537196301176e-02,
    8.04255993102895e-02, 8.12136704447979e-02, 8.20184716471344e-02, 8.28405664392583e-02, 8.36805447373290e-02,
    8.45390244462388e-02, 8.54166531722937e-02, 8.63141100644354e-02, 8.72321077953886e-02, 8.81713946953474e-02,
    8.91327570520271e-02, 9.01170215924217e-02, 9.11250581632225e-02, 9.21577826286684e-02, 9.32161600065677e-02,
    9.43012078656787e-02, 9.54140000100007e-02, 9.65556704785876e-02, 9.77274178926762e-02, 9.89305101855936e-02,
    1.00166289755146e-01, 1.01436179082793e-01, 1.02741686869338e-01, 1.04084414742952e-01, 1.05466064602181e-01,
    1.06888446664513e-01, 1.08353488300132e-01, 1.09863243740694e-01, 1.11419904764858e-01, 1.13025812475804e-01,
    1.14683470301675e-01, 1.16395558367910e-01, 1.18164949411119e-01, 1.19994726428692e-01, 1.21888202285954e-01,
    1.23848941535858e-01, 1.25880784744097e-01, 1.27987875657392e-01, 1.30174691605324e-01, 1.32446077587840e-01,
    1.34807284573851e-01, 1.37264012622799e-01, 1.39822459544273e-01, 1.42489375933601e-01, 1.45272127568431e-01,
    1.48178766328292e-01, 1.51218111012368e-01, 1.54399839689173e-01, 1.57734595526106e-01, 1.61234108430976e-01,
    1.64911335308903e-01, 1.68780622319388e-01, 1.72857893239124e-01, 1.77160868934300e-01, 1.81709324071779e-01,
    1.86525388617776e-01, 1.91633903472561e-01, 1.97062841888240e-01, 2.02843811271485e-01, 2.09012653799732e-01,
    2.15610169273952e-01, 2.22682990203222e-01, 2.30284647840434e-01, 2.38476879578463e-01, 2.47331243935506e-01,
    2.56931130997193e-01, 2.67374286122045e-01, 2.78776006654390e-01, 2.91273230921782e-01, 3.05029824532737e-01,
    3.20243494424876e-01, 3.37154947792243e-01, 3.56060196110919e-01, 3.77327342746089e-01, 4.01419886817026e-01,
    4.28929703940887e-01, 4.60624750236354e-01, 4.97519791888632e-01, 5.40984294142651e-01, 5.92912498016053e-01,
    6.56002088728957e-01, 7.34232298458703e-01, 8.33731313414106e-01, 9.64463147693015e-01, 1.14381339640412e+00,
    1.40517335796248e+00, 1.82227062151045e+00, 2.59912058457152e+00, 4.62773895951686e+00, 0.0
};


G4double G4GoudsmitSaundersonTable::fgInverseQ2CDFs[fgNumLambdas*fgNumLamG1*fgNumUvalues] = {0.0};
G4double G4GoudsmitSaundersonTable::fgInterParamsA2[fgNumLambdas*fgNumLamG1*fgNumUvalues] = {0.0};
G4double G4GoudsmitSaundersonTable::fgInterParamsB2[fgNumLambdas*fgNumLamG1*fgNumUvalues] = {0.0};
G4double G4GoudsmitSaundersonTable::fgInverseQ2CDFsII[fgNumLambdas*fgNumLamG1II*fgNumUvalues] = {0.0};
G4double G4GoudsmitSaundersonTable::fgInterParamsA2II[fgNumLambdas*fgNumLamG1II*fgNumUvalues] = {0.0};
G4double G4GoudsmitSaundersonTable::fgInterParamsB2II[fgNumLambdas*fgNumLamG1II*fgNumUvalues] = {0.0};

std::vector<G4double>* G4GoudsmitSaundersonTable::fgMoliereBc  = 0;
std::vector<G4double>* G4GoudsmitSaundersonTable::fgMoliereXc2 = 0;

G4bool G4GoudsmitSaundersonTable::fgIsInitialised = FALSE;

////////////////////////////////////////////////////////////////////////////////
void G4GoudsmitSaundersonTable::Initialise(){
   if(!fgIsInitialised){
     LoadMSCData();
     LoadMSCDataII();
     fgIsInitialised = TRUE;
   }

   InitMoliereMSCParams();
}

////////////////////////////////////////////////////////////////////////////////
G4GoudsmitSaundersonTable::~G4GoudsmitSaundersonTable(){
  if(fgMoliereBc){
    delete fgMoliereBc;
    fgMoliereBc = 0;
  }
  if(fgMoliereXc2){
    delete fgMoliereXc2;
    fgMoliereXc2 = 0;
  }
  fgIsInitialised = FALSE;
}

////////////////////////////////////////////////////////////////////////////////
// sample cos(theta) i.e. angular deflection from the precomputed angular
// distributions in the real multiple scattering case
G4double G4GoudsmitSaundersonTable::SampleCosTheta(G4double lambdavalue, G4double lamG1value,
                            G4double screeningparam, G4double rndm1, G4double rndm2,
                            G4double rndm3){
  const G4double lnl0   = 0.0;              //log(fgLambdaValues[0]);
  const G4double invlnl = 6.5144172285e+00; //1./log(fgLambdaValues[i+1]/fgLambdaValues[i])
  G4double logLambdaVal = G4Log(lambdavalue);
  G4int lambdaindx = (G4int)((logLambdaVal-lnl0)*invlnl);

  const G4double dd = 1.0/0.0499; // delta
  G4int lamg1indx  = (G4int)((lamG1value-fgLamG1Values[0])*dd);

  G4double probOfLamJ = G4Log(fgLambdaValues[lambdaindx+1]/lambdavalue)*invlnl;
  if(rndm1 > probOfLamJ)
     ++lambdaindx;

  // store 1.0/(fgLamG1Values[lamg1indx+1]-fgLamG1Values[lamg1indx]) const value
  const G4double onePerLamG1Diff = dd;
  G4double probOfLamG1J = (fgLamG1Values[lamg1indx+1]-lamG1value)*onePerLamG1Diff;
  if(rndm2 > probOfLamG1J)
     ++lamg1indx;

  G4int begin = lambdaindx*(fgNumLamG1*fgNumUvalues)+lamg1indx*fgNumUvalues;

  // sample bin
  G4int indx =  (G4int)(rndm3*100); // value/delatU ahol deltaU = 0.01; find u-indx
  G4double cp1 = fgUValues[indx];
//  G4double cp2 = fgUValues[indx+1];
  indx +=begin;

  G4double u1  = fgInverseQ2CDFs[indx];
  G4double u2  = fgInverseQ2CDFs[indx+1];

  G4double dum1 = rndm3 - cp1;
  const G4double dum2  = 0.01;//cp2-cp1;  // this is constans 0.01
  const G4double dum22 = 0.0001;
  // rational interpolation of the inverse CDF
  G4double sample= u1+(1.0+fgInterParamsA2[indx]+fgInterParamsB2[indx])*dum1*dum2/(dum22+
               dum1*dum2*fgInterParamsA2[indx]+fgInterParamsB2[indx]*dum1*dum1)*(u2-u1);

  // transform back u to cos(theta) :
  // -compute parameter value a = w2*screening_parameter
  G4double t = logLambdaVal;
  G4double dum;
  if(lambdavalue >= 10.0)
    dum = 0.5*(-2.77164+t*( 2.94874-t*(0.1535754-t*0.00552888) ));
  else
    dum = 0.5*(1.347+t*(0.209364-t*(0.45525-t*(0.50142-t*0.081234))));

  G4double a = dum*(lambdavalue+4.0)*screeningparam;

  // this is the sampled cos(theta)
  return (2.0*a*sample)/(1.0-sample+a);
}

////////////////////////////////////////////////////////////////////////////////
// sample cos(theta) i.e. angular deflection from the precomputed angular
// distributions in the real multiple scattering case (For the second grid)
G4double G4GoudsmitSaundersonTable::SampleCosThetaII(G4double lambdavalue, G4double lamG1value,
                            G4double screeningparam, G4double rndm1, G4double rndm2,
                            G4double rndm3){
  const G4double lnl0   = 0.0;              //log(fgLambdaValues[0]);
  const G4double invlnl = 6.5144172285e+00; //1./log(fgLambdaValues[i+1]/fgLambdaValues[i])
  G4double logLambdaVal = G4Log(lambdavalue);
  G4int lambdaindx = (G4int)((logLambdaVal-lnl0)*invlnl);

  const G4double dd = 1.0/0.333; // delta
  G4int lamg1indx  = (G4int)((lamG1value-fgLamG1ValuesII[0])*dd);

  G4double probOfLamJ = G4Log(fgLambdaValues[lambdaindx+1]/lambdavalue)*invlnl;
  if(rndm1 > probOfLamJ)
     ++lambdaindx;

  // store 1.0/(fgLamG1Values[lamg1indx+1]-fgLamG1Values[lamg1indx]) const value
  const G4double onePerLamG1Diff = dd;
  G4double probOfLamG1J = (fgLamG1ValuesII[lamg1indx+1]-lamG1value)*onePerLamG1Diff;
  if(rndm2 > probOfLamG1J)
     ++lamg1indx;

  // protection against cases when the sampled lamG1 values are not possible due
  // to the fact that G1<=1 (G1 = lam_el/lam_tr1)
  // -- in theory it should never happen when lamg1indx=0 but check it just to be sure
  if(fgLamG1ValuesII[lamg1indx]>lambdavalue && lamg1indx>0) 
    --lamg1indx;
    

  G4int begin = lambdaindx*(fgNumLamG1II*fgNumUvalues)+lamg1indx*fgNumUvalues;

  // sample bin
  G4int indx =  (G4int)(rndm3*100); // value/delatU ahol deltaU = 0.01; find u-indx
  G4double cp1 = fgUValues[indx];
//  G4double cp2 = fgUValues[indx+1];
  indx +=begin;

  G4double u1  = fgInverseQ2CDFsII[indx];
  G4double u2  = fgInverseQ2CDFsII[indx+1];

  G4double dum1 = rndm3 - cp1;
  const G4double dum2  = 0.01;//cp2-cp1;  // this is constans 0.01
  const G4double dum22 = 0.0001;
  // rational interpolation of the inverse CDF
  G4double sample= u1+(1.0+fgInterParamsA2II[indx]+fgInterParamsB2II[indx])*dum1*dum2/(dum22+
               dum1*dum2*fgInterParamsA2II[indx]+fgInterParamsB2II[indx]*dum1*dum1)*(u2-u1);

  // transform back u to cos(theta) :
  // -compute parameter value a = w2*screening_parameter
  G4double t = logLambdaVal;
  G4double dum = 0.;
  if(lambdavalue >= 10.0)
    dum = 0.5*(-2.77164+t*( 2.94874-t*(0.1535754-t*0.00552888) ));
  else
    dum = 0.5*(1.347+t*(0.209364-t*(0.45525-t*(0.50142-t*0.081234))));

  G4double a = dum*(lambdavalue+4.0)*screeningparam;

  // this is the sampled cos(theta)
  return (2.0*a*sample)/(1.0-sample+a);
}

////////////////////////////////////////////////////////////////////////////////
// samplig multiple scattering angles cos(theta) and sin(thata)
// including no-scattering, single, "few" scattering cases as well
void G4GoudsmitSaundersonTable::Sampling(G4double lambdavalue, G4double lamG1value,
                            G4double scrPar, G4double &cost, G4double &sint){
  G4double rand0 = G4UniformRand();
  G4double expn  = G4Exp(-lambdavalue);

  //
  // no scattering case
  //
  if(rand0 < expn) {
    cost = 1.0;
    sint = 0.0;
    return;
  }

  //
  // single scattering case : sample (analitically) from the single scattering PDF
  //                          that corresponds to the modified-Wentzel DCS i.e.
  //                          Wentzel DCS that gives back the PWA first-transport mfp
  //
  if( rand0 < (1.+lambdavalue)*expn) {
    G4double rand1 = G4UniformRand();
    // sampling 1-cos(theta)
    G4double dum0  = 2.0*scrPar*rand1/(1.0 - rand1 + scrPar);
    // add protections
    if(dum0 < 0.0)       dum0 = 0.0;
    else if(dum0 > 2.0 ) dum0 = 2.0;
    // compute cos(theta) and sin(theta)
    cost = 1.0 - dum0;
    sint = std::sqrt(dum0*(2.0 - dum0));
    return;
  }

  //
  // handle this case:
  //      -lambdavalue < 1 i.e. mean #elastic events along the step is < 1 but
  //       the currently sampled case is not 0 or 1 scattering. [Our minimal
  //       lambdavalue (that we have precomputed, transformed angular distributions
  //       stored in a form of equally probabe intervalls together with rational
  //       interp. parameters) is 1.]
  //      -probability of having n elastic events follows Poisson stat. with
  //       lambdavalue parameter.
  //      -the max. probability (when lambdavalue=1) of having more than one
  //       elastic events is 0.2642411 and the prob of having 2,3,..,n elastic
  //       events decays rapidly with n. So set a max n to 10.
  //      -sampling of this cases is done in a one-by-one single elastic event way
  //       where the current #elastic event is sampled from the Poisson distr.
  //      -(other possibility would be to use lambdavalue=1 distribution because
  //        they are very close)
  if(lambdavalue < 1.0) {
     G4double prob, cumprob;
     prob = cumprob = expn;
     G4double curcost,cursint;
     // init cos(theta) and sin(theta) to the zero scattering values
     cost = 1.0;
     sint = 0.0;

     for(G4int iel = 1; iel < 10; ++iel){
        // prob of having iel scattering from Poisson
        prob    *= lambdavalue/(G4double)iel;
        cumprob += prob;
        //
        //sample cos(theta) from the singe scattering pdf
        //
        G4double rand1 = G4UniformRand();
        // sampling 1-cos(theta)
        G4double dum0  = 2.0*scrPar*rand1/(1.0 - rand1 + scrPar);
        // compute cos(theta) and sin(theta)
        curcost = 1.0 - dum0;
        cursint = dum0*(2.0 - dum0);
        //
        // if we got current deflection that is not too small
        // then update cos(theta) sin(theta)
        if(cursint > 1.0e-20){
          cursint = std::sqrt(cursint);
          G4double curphi = CLHEP::twopi*G4UniformRand();
          cost = cost*curcost - sint*cursint*std::cos(curphi);
          sint = std::sqrt(std::max(0.0, (1.0-cost)*(1.0+cost)));
        }
        //
        // check if we have done enough scattering i.e. sampling from the Poisson
        if( rand0 < cumprob)
          return;
     }
     // if reached the max iter i.e. 10
     return;
  }


  /// **** let it izotropic if we are above the grid i.e. true path length is long
  //  IT IS DONE NOW IN THE MODEL
  //if(lamG1value>fgLamG1Values[fgNumLamG1-1]) {
  //  cost = 1.-2.*G4UniformRand();
  //  sint = std::sqrt((1.-cost)*(1.+cost));
  //  return;
  //}


  //
  // multiple scattering case with lambdavalue >= 1:
  //   - use the precomputed and transformed Goudsmit-Saunderson angular
  //     distributions that are stored in a form of equally probabe intervalls
  //     together with the corresponding rational interp. parameters
  //
  // add protections
  if(lamG1value<fgLamG1Values[0]) lamG1value = fgLamG1Values[0];
  if(lamG1value>fgLamG1ValuesII[fgNumLamG1II-1]) lamG1value = fgLamG1ValuesII[fgNumLamG1II-1];

  // sample 1-cos(theta) :: depending on the grids
  G4double dum0 = 0.0;
  if(lamG1value<fgLamG1Values[fgNumLamG1-1]) {
      dum0 = SampleCosTheta(lambdavalue, lamG1value, scrPar, G4UniformRand(),
                                 G4UniformRand(), G4UniformRand());

  } else {
    dum0 = SampleCosThetaII(lambdavalue, lamG1value, scrPar, G4UniformRand(),
                               G4UniformRand(), G4UniformRand());
  }

  // add protections
  if(dum0 < 0.0)  dum0 = 0.0;
  if(dum0 > 2.0 ) dum0 = 2.0;

  // compute cos(theta) and sin(theta)
     cost = 1.0-dum0;
     sint = std::sqrt(dum0*(2.0 - dum0));
}

////////////////////////////////////////////////////////////////////////////////
G4double G4GoudsmitSaundersonTable::GetScreeningParam(G4double G1){
   if(G1<fgG1Values[0])                      return fgScreeningParam[0];
   if(G1>fgG1Values[fgNumScreeningParams-1]) return fgScreeningParam[fgNumScreeningParams-1];
   // normal case
   const G4double lng10   = -2.30258509299405e+01;   //log(fgG1Values[0]);
   const G4double invlng1 =  6.948711710452e+00;     //1./log(fgG1Values[i+1]/fgG1Values[i])
   G4double logG1 = G4Log(G1);
   G4int k = (G4int)((logG1-lng10)*invlng1);
 return G4Exp(logG1*fgSrcAValues[k])*fgSrcBValues[k]; // log-log linear
}

////////////////////////////////////////////////////////////////////////////////
void G4GoudsmitSaundersonTable::LoadMSCData(){
 G4double dum;
 char  basename[512];
 char* path = getenv("G4LEDATA");
 if (!path) {
    G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
 }

 std::string pathString(path);
 sprintf(basename,"inverseq2CDF");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		 FatalException,
		 msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInverseQ2CDFs[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }


 sprintf(basename,"inverseq2CDFRatInterpPA");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		 FatalException,
		 msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInterParamsA2[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }

 sprintf(basename,"inverseq2CDFRatInterpPB");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		FatalException,
		msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInterParamsB2[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }
}

////////////////////////////////////////////////////////////////////////////////
// Load the second grid data
////////////////////////////////////////////////////////////////////////////////
void G4GoudsmitSaundersonTable::LoadMSCDataII(){
G4double dum;
char  basename[512];
char* path = getenv("G4LEDATA");
if (!path) {
   G4Exception("G4GoudsmitSaundersonTable::LoadMSCDataII()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
   return;
}


 std::string pathString(path);
 sprintf(basename,"inverseq2CDFII");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCDataII()","em0006",
		 FatalException,
		 msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1II;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1II+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInverseQ2CDFsII[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }


 sprintf(basename,"inverseq2CDFRatInterpPAII");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCDataII()","em0006",
		 FatalException,
		 msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1II;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1II+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInterParamsA2II[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }

 sprintf(basename,"inverseq2CDFRatInterpPBII");
 for(int k = 0; k < fgNumLambdas; ++k){
   char fname[512];
   sprintf(fname,"%s/msc_GS/%s_%d",path,basename,k);
   std::ifstream infile(fname,std::ios::in);
   if(!infile.is_open()){
     char msgc[512];
     sprintf(msgc,"Cannot open file: %s .",fname);
     G4Exception("G4GoudsmitSaundersonTable::LoadMSCDataII()","em0006",
		FatalException,
		msgc);
     return;
   }

   int limk = k*fgNumUvalues*fgNumLamG1II;
   for(int i = 0; i < fgNumUvalues; ++i)
    for(int j = 0; j < fgNumLamG1II+1; ++j)
    if(j==0) infile >> dum;
    else     infile >> fgInterParamsB2II[limk+(j-1)*fgNumUvalues+i];

   infile.close();
 }
}

////////////////////////////////////////////////////////////////////////////////
// compute material dependent Moliere MSC parameters at initialisation
void G4GoudsmitSaundersonTable::InitMoliereMSCParams(){
   const G4double const1   = 7821.6;      // [cm2/g]
   const G4double const2   = 0.1569;      // [cm2 MeV2 / g]
   const G4double finstrc2 = 5.325135453E-5; // fine-structure const. square

   G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
   // get number of materials in the table
   unsigned int numMaterials = theMaterialTable->size();
   // make sure that we have long enough vectors
   if(!fgMoliereBc) {
     fgMoliereBc  = new std::vector<G4double>(numMaterials);
     fgMoliereXc2 = new std::vector<G4double>(numMaterials);
   } else {
     fgMoliereBc->resize(numMaterials);
     fgMoliereXc2->resize(numMaterials);
   }

   const G4double xims = 0.0; // use xi_MS = 0.0 value now. We will see later on.
   for(unsigned int imat = 0; imat < numMaterials; ++imat) {
     const G4Material      *theMaterial  = (*theMaterialTable)[imat];
     const G4ElementVector *theElemVect  = theMaterial->GetElementVector();

     const G4int     numelems        = theMaterial->GetNumberOfElements();
     const G4double* theFractionVect = theMaterial->GetFractionVector();
     G4double zs = 0.0;
     G4double zx = 0.0;
     G4double ze = 0.0;

     for(G4int ielem = 0; ielem < numelems; ielem++) {
       G4double izet = (*theElemVect)[ielem]->GetZ();
       G4double iwa  = (*theElemVect)[ielem]->GetN();
//       G4int ipz     = theAtomsVect[ielem];
       G4double ipz  = theFractionVect[ielem];

       G4double dum  = ipz/iwa*izet*(izet+xims);
       zs += dum;
       ze += dum*(-2.0/3.0)*G4Log(izet);
       zx += dum*G4Log(1.0+3.34*finstrc2*izet*izet);
//       wa += ipz*iwa;
     }
     G4double density = theMaterial->GetDensity()*CLHEP::cm3/CLHEP::g; // [g/cm3]

     (*fgMoliereBc)[theMaterial->GetIndex()]  = const1*density*zs*G4Exp(ze/zs)/G4Exp(zx/zs);  //[1/cm]
     (*fgMoliereXc2)[theMaterial->GetIndex()] = const2*density*zs;  // [MeV2/cm]

     // change to Geant4 internal units of 1/length and energ2/length
     (*fgMoliereBc)[theMaterial->GetIndex()]  *= 1.0/CLHEP::cm;
     (*fgMoliereXc2)[theMaterial->GetIndex()] *= CLHEP::MeV*CLHEP::MeV/CLHEP::cm;
   }

}
