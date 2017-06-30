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
// $Id: G4GoudsmitSaundersonTable.cc 103884 2017-05-03 08:04:50Z gcosmo $
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
// 28.04.2017 M. Novak: New representation of the angular distribution data with
//            significantly reduced data size.  
//
// References:
//   [1] A.F.Bielajew, NIMB, 111 (1996) 195-208
//   [2] I.Kawrakow, A.F.Bielajew, NIMB 134(1998) 325-336
//
// -----------------------------------------------------------------------------

#include "G4GoudsmitSaundersonTable.hh"

#include <fstream>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4Log.hh"
#include "G4Exp.hh"



const G4double G4GoudsmitSaundersonTable::gG1Values1[]={
  1.00000000000000e-10,	1.11813656534780e-10,	1.25022937876777e-10,	1.39792718347230e-10,	1.56307349953404e-10,
  1.74772963415515e-10,	1.95420041029080e-10,	2.18506293476381e-10,	2.44319876494558e-10,	2.73182987549824e-10,
  3.05455887410409e-10,	3.41541396814338e-10,	3.81889924358074e-10,	4.27005088362666e-10,	4.77450002887863e-10,
  5.33854306354332e-10,	5.96922020503163e-10,	6.67440337785873e-10,	7.46289446866469e-10,	8.34453518874582e-10,
  9.33032991536807e-10,	1.04325830451314e-09,	1.16650525737890e-09,	1.30431218194578e-09,	1.45839914326215e-09,
  1.63068940895331e-09,	1.82333345487608e-09,	2.03873580671887e-09,	2.27958505257621e-09,	2.54888740110574e-09,
  2.85000420413065e-09,	3.18669391203342e-09,	3.56315898561579e-09,	3.98409834996457e-09,	4.45476604503721e-09,
  4.98103680502590e-09,	5.56947938504262e-09,	6.22743855036691e-09,	6.96312675162171e-09,	7.78572663013967e-09,
  8.70550563296124e-09,	9.73394416805517e-09,	1.08838788993564e-08,	1.21696629701878e-08,	1.36073451549260e-08,
  1.52148701750309e-08,	1.70123026797217e-08,	1.90220776869612e-08,	2.12692806106778e-08,	2.37819603694417e-08,
  2.65914794847249e-08,	2.97329055385667e-08,	3.32454488767036e-08,	3.71729520204430e-08,	4.15644368959766e-08,
  4.64747167114826e-08,	5.19650801192889e-08,	5.81040562026047e-08,	6.49682698351557e-08,	7.26433980900701e-08,
  8.12252396356236e-08,	9.08209104657278e-08,	1.01550180889909e-07,	1.13546970470690e-07,	1.26961019567745e-07,
  1.41959758352533e-07,	1.58730396621904e-07,	1.77482260495109e-07,	1.98449405160165e-07,	2.21893536281100e-07,
  2.48107276530226e-07,	2.77417818017302e-07,	3.10191006204147e-07,	3.46835906278882e-07,	3.87809908985959e-07,
  4.33624439641402e-07,	4.84851341591500e-07,	5.42130013791391e-07,	6.06175391592660e-07,	6.77786870353773e-07,
  7.57858283255199e-07,	8.47389057859347e-07,	9.47496690668154e-07,	1.05943069538209e-06,	1.18458819895856e-06,
  1.32453138013506e-06,	1.48100696807959e-06,	1.65596804454467e-06,	1.85159842165288e-06,	2.07033989959036e-06,
  2.31492274443046e-06,	2.58839976650298e-06,	2.89418442466467e-06,	3.23609343207765e-06,	3.61839439528787e-06,
  4.04585908122089e-06,	4.52382297695752e-06,	5.05825188569672e-06,	5.65581639013695e-06,	6.32397511270551e-06,
  7.07106781186548e-06,	7.90641947650063e-06,	8.84045671765333e-06,	9.88483791038276e-06,	1.10525987101351e-05,
  1.23583147599179e-05,	1.38182836191416e-05,	1.54507281849087e-05,	1.72760241447962e-05,	1.93169543001281e-05,
  2.15989929341255e-05,	2.41506237743345e-05,	2.70036955180411e-05,	3.01938193582402e-05,	3.37608134719546e-05,
  3.77492000188789e-05,	4.22087608537362e-05,	4.71951588885831e-05,	5.27706328607239e-05,	5.90047741781194e-05,
  6.59753955386448e-05,	7.37695021650427e-05,	8.24843777782376e-05,	9.22287988638087e-05,	1.03124392387732e-04,
  1.15307153907997e-04,	1.28929145030718e-04,	1.44160391397874e-04,	1.61191004896813e-04,	1.80233556580283e-04,
  2.01525729915095e-04,	2.25333287476472e-04,	2.51953388117470e-04,	2.81718296017409e-04,	3.14999527904539e-04,
  3.52212490217359e-04,	3.93821664084232e-04,	4.40346402838697e-04,	4.92367414433318e-04,	5.50534009663646e-04,
  6.15572206672459e-04,	6.88293792892308e-04,	7.69606457534812e-04,	8.60525121097460e-04,	9.62184603299410e-04,
  1.07585378756374e-03,	1.20295145884293e-03,	1.34506401247076e-03,	1.50396525507698e-03,	1.68163854471420e-03,
  1.88030154654320e-03,	2.10243391306996e-03,	2.35080823443078e-03,	2.62852464503775e-03,	2.93904951853454e-03,
  3.28625873404131e-03,	3.67448605372515e-03,	4.10857721553061e-03,	4.59395041623963e-03,	5.13666393979225e-03,
  5.74349177498518e-03,	6.42200816638524e-03,	7.18068215379752e-03,	8.02898328030137e-03,	8.97749978827105e-03,
  1.00380707786680e-02,	1.12239339831779e-02,	1.25498909936409e-02,	1.40324920111190e-02,	1.56902424205829e-02,
  1.75438337696249e-02,	1.96164020342010e-02,	2.19338163950031e-02,	2.45250021288779e-02,	2.74223016455310e-02,
  3.06618781758653e-02,	3.42841671506745e-02,	3.83343809036650e-02,	4.28630729983582e-02,	4.79267692226361e-02,
  5.35886731268147e-02,	5.99194549115623e-02,	6.69981335123266e-02,	7.49130628901859e-02,	8.37630348397159e-02,
  9.36585120787880e-02,	1.04723007011361e-01,	1.17094623372576e-01,	1.30927779998507e-01,	1.46395138236142e-01,
  1.63689757050976e-01,	1.83027502731593e-01,	2.04649743268488e-01,	2.28826361037536e-01,	2.55859121391545e-01,
  2.86085439205647e-01,	3.19882590389418e-01,	3.57672420932580e-01,	3.99926612261187e-01,	4.47172568624903e-01,
  5.00000000000001e-01
};

const G4double G4GoudsmitSaundersonTable::gScrAValues1[] = {
  1.04157965108795e+00,	1.00489574216333e+00,	1.04054171766980e+00,	1.04070044848564e+00,	1.06610486501067e+00,
  1.03212411797639e+00,	1.04186181391655e+00,	1.03806695397030e+00,	1.02828633376454e+00,	1.04189017675023e+00,
  1.05394563184103e+00,	1.03804859577890e+00,	1.04287935130576e+00,	1.04089067509304e+00,	1.03680659160640e+00,
  1.04315290812772e+00,	1.04885551056127e+00,	1.04149443516180e+00,	1.04415435067166e+00,	1.04310819097932e+00,
  1.04146814384517e+00,	1.04640481749889e+00,	1.04530214331355e+00,	1.04383868894705e+00,	1.04554873719815e+00,
  1.04499781690510e+00,	1.04441460529014e+00,	1.04683911728964e+00,	1.04545554920892e+00,	1.04663018492595e+00,
  1.04703168231144e+00,	1.04683581255684e+00,	1.04653923434779e+00,	1.04829979258202e+00,	1.04687530707532e+00,
  1.04798795480120e+00,	1.04856435291023e+00,	1.04852473343216e+00,	1.04851363909937e+00,	1.04945890378783e+00,
  1.04870544718128e+00,	1.04996471037899e+00,	1.05005469174707e+00,	1.05041253331647e+00,	1.05027955019518e+00,
  1.05094931288355e+00,	1.05105026713329e+00,	1.05146143471074e+00,	1.05186311864079e+00,	1.05211792402488e+00,
  1.05229626092483e+00,	1.05273676509171e+00,	1.05306578843975e+00,	1.05339905039215e+00,	1.05376106738709e+00,
  1.05410304876713e+00,	1.05440774351329e+00,	1.05474466789715e+00,	1.05513191826058e+00,	1.05549243912963e+00,
  1.05585783258782e+00,	1.05624821688601e+00,	1.05659792043794e+00,	1.05696998024445e+00,	1.05737019747710e+00,
  1.05776152902401e+00,	1.05815838472034e+00,	1.05857004546784e+00,	1.05897169614415e+00,	1.05937160731340e+00,
  1.05980349982220e+00,	1.06022984673212e+00,	1.06066249265147e+00,	1.06110945722573e+00,	1.06153938214552e+00,
  1.06199967846382e+00,	1.06245899725372e+00,	1.06292537694282e+00,	1.06339898465598e+00,	1.06388339726801e+00,
  1.06436517542617e+00,	1.06486492927630e+00,	1.06536923137948e+00,	1.06588168246212e+00,	1.06640248548639e+00,
  1.06693282043666e+00,	1.06746902378861e+00,	1.06801714089088e+00,	1.06857352353498e+00,	1.06913938237093e+00,
  1.06971496644464e+00,	1.07030090511985e+00,	1.07089636394310e+00,	1.07150231293887e+00,	1.07211985687749e+00,
  1.07274812967738e+00,	1.07338782431134e+00,	1.07403926162142e+00,	1.07470277473663e+00,	1.07537870967102e+00,
  1.07606742595665e+00,	1.07676929731457e+00,	1.07748471236607e+00,	1.07821407538698e+00,	1.07895789130915e+00,
  1.07971626136770e+00,	1.08049014701186e+00,	1.08127968685614e+00,	1.08208546070942e+00,	1.08290798546009e+00,
  1.08374784764054e+00,	1.08460542143151e+00,	1.08548157805625e+00,	1.08637674308402e+00,	1.08729160600829e+00,
  1.08822683862055e+00,	1.08918317336784e+00,	1.09016122856967e+00,	1.09116195193296e+00,	1.09218603392848e+00,
  1.09323435159941e+00,	1.09430779436841e+00,	1.09540730659958e+00,	1.09653382728267e+00,	1.09768845084086e+00,
  1.09887221998465e+00,	1.10008628700310e+00,	1.10133186213536e+00,	1.10261018775782e+00,	1.10392262316543e+00,
  1.10527059084102e+00,	1.10665553492655e+00,	1.10807908026622e+00,	1.10954286235417e+00,	1.11104865114364e+00,
  1.11259831365381e+00,	1.11419383300353e+00,	1.11583729384280e+00,	1.11753094916587e+00,	1.11927715602890e+00,
  1.12107842962296e+00,	1.12293746859341e+00,	1.12485712546544e+00,	1.12684044667885e+00,	1.12889070084775e+00,
  1.13101136387157e+00,	1.13320615891507e+00,	1.13547908296162e+00,	1.13783440834654e+00,	1.14027673014285e+00,
  1.14281098480473e+00,	1.14544248142819e+00,	1.14817694097055e+00,	1.15102053514709e+00,	1.15397992901364e+00,
  1.15706233242302e+00,	1.16027555423021e+00,	1.16362806276602e+00,	1.16712906153756e+00,	1.17078856139551e+00,
  1.17461747454111e+00,	1.17862771404547e+00,	1.18283230939909e+00,	1.18724553705940e+00,	1.19188307189529e+00,
  1.19676215742827e+00,	1.20190180478471e+00,	1.20732302052005e+00,	1.21304907048578e+00,	1.21910578803603e+00,
  1.22552193133196e+00,	1.23232960387200e+00,	1.23956474857061e+00,	1.24726773187585e+00,	1.25548403737872e+00,
  1.26426509489660e+00,	1.27366927496387e+00,	1.28376309013730e+00,	1.29462265524641e+00,	1.30633547329637e+00,
  1.31900263594822e+00,	1.33274155598412e+00,	1.34768938750939e+00,	1.36400734614952e+00,	1.38188621839712e+00,
  1.40155346152780e+00,	1.42328245813182e+00,	1.44740473052999e+00,	1.47432628433781e+00,	1.50454981066593e+00,
  1.53870535800673e+00,	1.57759350750855e+00,	1.62224744405266e+00,	1.67402435052153e+00,	1.73474370149111e+00,
  1.80690321932344e+00,	1.89402874011488e+00,	2.00126625381261e+00,	2.13643769500363e+00,	2.31204919591699e+00,
  0.00000000000000e+00
};

const G4double G4GoudsmitSaundersonTable::gScrBValues1[] = {
  5.03312053141067e-02,	2.17158694570135e-02,	4.89530420350471e-02,	4.91296760086723e-02,	8.71886858692111e-02,
  4.06342377881512e-02,	5.05167856049234e-02,	4.64275030486448e-02,	3.73907006070029e-02,	5.04503967256328e-02,
  6.57011334316898e-02,	4.64603705597810e-02,	5.15916031665536e-02,	4.94249215340636e-02,	4.52770321697444e-02,
  5.18471229639412e-02,	5.85229132211994e-02,	5.00937731823507e-02,	5.29737893758232e-02,	5.18278719763304e-02,
  5.00902942334804e-02,	5.54744303244664e-02,	5.42303662705311e-02,	5.26308403011457e-02,	5.44942276528391e-02,
  5.38901284419719e-02,	5.32613848445252e-02,	5.59091658773638e-02,	5.43908715707870e-02,	5.56698979976760e-02,
  5.61114230117329e-02,	5.58968126339045e-02,	5.55752599197682e-02,	5.75002338073153e-02,	5.59465721884279e-02,
  5.71493691866209e-02,	5.77788824961718e-02,	5.77356464303655e-02,	5.77236166819877e-02,	5.87514247126034e-02,
  5.79355837467423e-02,	5.92972073210311e-02,	5.93951223491295e-02,	5.97837285914114e-02,	5.96399021521681e-02,
  6.03632962616234e-02,	6.04724111049718e-02,	6.09160596930517e-02,	6.13498664951689e-02,	6.16248954003208e-02,
  6.18168887323794e-02,	6.22906285773529e-02,	6.26445419525010e-02,	6.30027198820527e-02,	6.33915610812053e-02,
  6.37586513609228e-02,	6.40853272813854e-02,	6.44460823428755e-02,	6.48604273634117e-02,	6.52459401901135e-02,
  6.56363234017371e-02,	6.60531083887963e-02,	6.64261143028855e-02,	6.68225028982192e-02,	6.72485258208146e-02,
  6.76647594484705e-02,	6.80864826339858e-02,	6.85235661778902e-02,	6.89496332396662e-02,	6.93733890763592e-02,
  6.98305909274069e-02,	7.02815314979391e-02,	7.07386938681932e-02,	7.12105552923658e-02,	7.16639570586871e-02,
  7.21488807066093e-02,	7.26323200186316e-02,	7.31226973850434e-02,	7.36201694780953e-02,	7.41284824122965e-02,
  7.46334967838100e-02,	7.51567953920159e-02,	7.56843139369856e-02,	7.62197880088967e-02,	7.67634067721601e-02,
  7.73163807961576e-02,	7.78748607678250e-02,	7.84451178308181e-02,	7.90233355480754e-02,	7.96107426405302e-02,
  8.02075698638796e-02,	8.08144421214313e-02,	8.14304645410927e-02,	8.20566078378425e-02,	8.26939838829410e-02,
  8.33416671123110e-02,	8.40003369095949e-02,	8.46702885152005e-02,	8.53518281751840e-02,	8.60452736689284e-02,
  8.67509548673856e-02,	8.74692143245562e-02,	8.82004079042808e-02,	8.89449054446883e-02,	8.97031776652742e-02,
  9.04752806639226e-02,	9.12621449420644e-02,	9.20638618073932e-02,	9.28809677060177e-02,	9.37139327647979e-02,
  9.45632950069344e-02,	9.54293742753128e-02,	9.63129898767921e-02,	9.72145083896325e-02,	9.81345586346619e-02,
  9.90737481888985e-02,	1.00032741939832e-01,	1.01012086104919e-01,	1.02012653000163e-01,	1.03035053408570e-01,
  1.04080079737091e-01,	1.05148532047009e-01,	1.06241260133595e-01,	1.07359103069986e-01,	1.08503048771995e-01,
  1.09674028081727e-01,	1.10873074047149e-01,	1.12101270532853e-01,	1.13359726193764e-01,	1.14649657617391e-01,
  1.15972336002916e-01,	1.17329045635105e-01,	1.18721239438085e-01,	1.20150371810475e-01,	1.21618018926522e-01,
  1.23125841585931e-01,	1.24675602992041e-01,	1.26269153716434e-01,	1.27908495768747e-01,	1.29595719341508e-01,
  1.31333054059601e-01,	1.33122892445295e-01,	1.34967760239165e-01,	1.36870353646818e-01,	1.38833565227297e-01,
  1.40860468186044e-01,	1.42954353084258e-01,	1.45118751619748e-01,	1.47357436315975e-01,	1.49674463733404e-01,
  1.52074190179311e-01,	1.54561299160811e-01,	1.57140836162498e-01,	1.59818242883263e-01,	1.62599394679300e-01,
  1.65490646080663e-01,	1.68498878580249e-01,	1.71631553951876e-01,	1.74896781476726e-01,	1.78303381453875e-01,
  1.81860968644217e-01,	1.85580040380620e-01,	1.89472079055570e-01,	1.93549667959742e-01,	1.97826625868490e-01,
  2.02318158342527e-01,	2.07041034807091e-01,	2.12013791462807e-01,	2.17256966510187e-01,	2.22793375200011e-01,
  2.28648428943953e-01,	2.34850511284297e-01,	2.41431420028350e-01,	2.48426890453852e-01,	2.55877217174183e-01,
  2.63827998186181e-01,	2.72331028202796e-01,	2.81445378786391e-01,	2.91238712555341e-01,	3.01788891991870e-01,
  3.13185963581234e-01,	3.25534623955993e-01,	3.38957309662964e-01,	3.53598103677331e-01,	3.69627721982226e-01,
  3.87249946066427e-01,	4.06710015890883e-01,	4.28305718666080e-01,	4.52402242352685e-01,	4.79452376902252e-01,
  5.10024456502661e-01,	5.44841745996517e-01,	5.84839150981769e-01,	6.31246863293948e-01,	6.85717185564307e-01,
  7.50523056016889e-01,	8.28880627894856e-01,	9.25497179158956e-01,	1.04755297673620e+00,	1.20658126625305e+00,
  0.00000000000000e+00
};

const G4double G4GoudsmitSaundersonTable::gG1Values2[]={
  5.00000000000000e-01,	5.06918722752298e-01,	5.13933182953643e-01,	5.21044705365768e-01,	5.28254633081725e-01,
  5.35564327779545e-01,	5.42975169979400e-01,	5.50488559304339e-01,	5.58105914744617e-01,	5.65828674925688e-01,
  5.73658298379911e-01,	5.81596263822002e-01,	5.89644070428317e-01,	5.97803238119977e-01,	6.06075307849934e-01,
  6.14461841893989e-01,	6.22964424145851e-01,	6.31584660416271e-01,	6.40324178736321e-01,	6.49184629664860e-01,
  6.58167686600270e-01,	6.67275046096487e-01,	6.76508428183425e-01,	6.85869576691813e-01,	6.95360259582547e-01,
  7.04982269280583e-01,	7.14737423013460e-01,	7.24627563154504e-01,	7.34654557570783e-01,	7.44820299975873e-01,
  7.55126710287506e-01,	7.65575734990175e-01,	7.76169347502743e-01,	7.86909548551151e-01,	7.97798366546274e-01,
  8.08837857967014e-01,	8.20030107748688e-01,	8.31377229676789e-01,	8.42881366786204e-01,	8.54544691765948e-01,
  8.66369407369501e-01,	8.78357746830827e-01,	8.90511974286139e-01,	9.02834385201514e-01,	9.15327306806416e-01,
  9.27993098533219e-01,	9.40834152462814e-01,	9.53852893776382e-01,	9.67051781213414e-01,	9.80433307536078e-01,
  9.94000000000000e-01
};

const G4double G4GoudsmitSaundersonTable::gScrAValues2[] = {
  2.43260802029458e+00,	2.46262608568646e+00,	2.49389298606006e+00,	2.52648830233721e+00,	2.56049853754171e+00,
  2.59601788683344e+00,	2.63314911269562e+00,	2.67200454242540e+00,	2.71270720837972e+00,	2.75539215535386e+00,
  2.80020794442213e+00,	2.84731838855627e+00,	2.89690456283251e+00,	2.94916714123393e+00,	3.00432912375239e+00,
  3.06263903196253e+00,	3.12437466980803e+00,	3.18984756971552e+00,	3.25940827440272e+00,	3.33345264348985e+00,
  3.41242942478934e+00,	3.49684939631883e+00,	3.58729647295325e+00,	3.68444128860530e+00,	3.78905792262525e+00,
  3.90204465382849e+00,	4.02444992107417e+00,	4.15750508109666e+00,	4.30266613502606e+00,	4.46166742598720e+00,
  4.63659151671749e+00,	4.82996123729109e+00,	5.04486256966139e+00,	5.28511113787701e+00,	5.55548149760216e+00,
  5.86202872540949e+00,	6.21254879271518e+00,	6.61725305332142e+00,	7.08978285849237e+00,	7.64878291654446e+00,
  8.32042898816053e+00,	9.14266198370974e+00,	1.01726447874941e+01,	1.15007249602827e+01,	1.32786704661903e+01,
  1.57827651696577e+01,	1.95753182250328e+01,	2.60088916622360e+01,	3.94153107765581e+01,	8.75131588110543e+01,
  0.00000000000000e+00
};

const G4double G4GoudsmitSaundersonTable::gScrBValues2[] = {
  1.31174192812636e+00,	1.33876879475328e+00,	1.36692483335422e+00,	1.39628213980526e+00,	1.42691908588269e+00,
  1.45892101775522e+00,	1.49238104990491e+00,	1.52740097004838e+00,	1.56409227363037e+00,	1.60257735002838e+00,
  1.64299084710729e+00,	1.68548124620309e+00,	1.73021268642906e+00,	1.77736708555737e+00,	1.82714661536618e+00,
  1.87977660250105e+00,	1.93550894278961e+00,	1.99462613821718e+00,	2.05744609330159e+00,	2.12432784287395e+00,
  2.19567842949722e+00,	2.27196120902661e+00,	2.35370594287901e+00,	2.44152114220022e+00,	2.53610927301215e+00,
  2.63828562730580e+00,	2.74900193477337e+00,	2.86937616592981e+00,	3.00073050806567e+00,	3.14464025535566e+00,
  3.30299745868968e+00,	3.47809481239111e+00,	3.67273770943219e+00,	3.89039616290446e+00,	4.13541419825107e+00,
  4.41330381199485e+00,	4.73116626310813e+00,	5.09831013668125e+00,	5.52718261974326e+00,	6.03481658964974e+00,
  6.64516146679752e+00,	7.39300060933648e+00,	8.33088077994393e+00,	9.54216525157779e+00,	1.11676544350690e+01,
  1.34658259138613e+01,	1.69701651097459e+01,	2.29982069382437e+01,	3.60380269802614e+01,	9.32253464839576e+01,
  0.00000000000000e+00
};

bool G4GoudsmitSaundersonTable::gIsInitialised = false;

std::vector<G4double>* G4GoudsmitSaundersonTable::fgMoliereBc  = 0;
std::vector<G4double>* G4GoudsmitSaundersonTable::fgMoliereXc2 = 0;

G4GoudsmitSaundersonTable::G4GoudsmitSaundersonTable() {}

G4GoudsmitSaundersonTable::~G4GoudsmitSaundersonTable() {
  for (size_t i=0; i<fGSMSCAngularDistributions1.size(); ++i) {
    if (fGSMSCAngularDistributions1[i]) {
      delete [] fGSMSCAngularDistributions1[i]->fUValues;
      delete [] fGSMSCAngularDistributions1[i]->fParamA;
      delete [] fGSMSCAngularDistributions1[i]->fParamB;
      delete fGSMSCAngularDistributions1[i];
    }  
  }
  fGSMSCAngularDistributions2.clear();
  for (size_t i=0; i<fGSMSCAngularDistributions2.size(); ++i) {
    if (fGSMSCAngularDistributions2[i]) {
      delete [] fGSMSCAngularDistributions2[i]->fUValues;
      delete [] fGSMSCAngularDistributions2[i]->fParamA;
      delete [] fGSMSCAngularDistributions2[i]->fParamB;
      delete fGSMSCAngularDistributions2[i];
    }
  }
  fGSMSCAngularDistributions2.clear();
  gIsInitialised = false;
}

void G4GoudsmitSaundersonTable::Initialise() {
  // load precomputed angular distributions and set up several values used during the sampling
  if (!gIsInitialised) {
    G4double lLambdaMin  = G4Log(gLAMBMIN);
    G4double lLambdaMax  = G4Log(gLAMBMAX);
    fLogLambda0        = lLambdaMin;
    fLogDeltaLambda    = (lLambdaMax-lLambdaMin)/(gLAMBNUM-1.);
    fInvLogDeltaLambda = 1./fLogDeltaLambda;
    fInvDeltaQ1        = 1./((gQMAX1-gQMIN1)/(gQNUM1-1.));
    fDeltaQ2           = (gQMAX2-gQMIN2)/(gQNUM2-1.);
    fInvDeltaQ2        = 1./fDeltaQ2;
    // for the A(G1) function
    fLogG1FuncMin1      = G4Log(gG1Values1[0]);
    fInvLogDeltaG1Func1 = 1./((G4Log(gG1Values1[gNUMSCR1-1]/gG1Values1[0]))/(gNUMSCR1-1.));
    fLogG1FuncMin2      = G4Log(gG1Values2[0]);
    fInvLogDeltaG1Func2 = 1./((G4Log(gG1Values2[gNUMSCR2-1]/gG1Values2[0]))/(gNUMSCR2-1.));
    LoadMSCData();
  }
  InitMoliereMSCParams();
}


// samplig multiple scattering angles cos(theta) and sin(thata)
// including no-scattering, single, "few" scattering cases as well
// lambdaval : s/lambda_el
// qval      : s/lambda_el G1
// scra      : screening parameter
void G4GoudsmitSaundersonTable::Sampling(G4double lambdaval, G4double qval, G4double scra, G4double &cost, G4double &sint){
  G4double rand0 = G4UniformRand();
  G4double expn  = G4Exp(-lambdaval);
  //
  // no scattering case
  if (rand0<expn) {
    cost = 1.0;
    sint = 0.0;
    return;
  }
  //
  // single scattering case : sample from the single scattering PDF
  if (rand0<(1.+lambdaval)*expn) {
    G4double rand1 = G4UniformRand();
    // sampling 1-cos(theta)
    G4double dum0  = 2.0*scra*rand1/(1.0-rand1+scra);
    // add protections
    if (dum0 < 0.0) {
      dum0 = 0.0;
    }
    if (dum0>2.0) {
      dum0 = 2.0;
    }
    // compute cos(theta) and sin(theta) from the sampled 1-cos(theta)
    cost = 1.0-dum0;
    sint = std::sqrt(dum0*(2.0-dum0));
    return;
  }
  //
  // handle this case:
  //      -lambdaval < 1 i.e. mean #elastic events along the step is < 1 but
  //       the currently sampled case is not 0 or 1 scattering. [Our minimal
  //       lambdaval (that we have precomputed, transformed angular distributions
  //       stored in a form of equally probabe intervalls together with rational
  //       interp. parameters) is 1.]
  //      -probability of having n elastic events follows Poisson stat. with
  //       lambdaval parameter.
  //      -the max. probability (when lambdaval=1) of having more than one
  //       elastic events is 0.2642411 and the prob of having 2,3,..,n elastic
  //       events decays rapidly with n. So set a max n to 10.
  //      -sampling of this cases is done in a one-by-one single elastic event way
  //       where the current #elastic event is sampled from the Poisson distr.
  if (lambdaval<1.0) {
    G4double prob, cumprob;
    prob = cumprob = expn;
    G4double curcost,cursint;
    // init cos(theta) and sin(theta) to the zero scattering values
    cost = 1.0;
    sint = 0.0;
    for (G4int iel=1; iel<10; ++iel) {
      // prob of having iel scattering from Poisson
      prob    *= lambdaval/(G4double)iel;
      cumprob += prob;
      //
      //sample cos(theta) from the singe scattering pdf
      //
      G4double rand1 = G4UniformRand();
      // sampling 1-cos(theta)
      G4double dum0  = 2.0*scra*rand1/(1.0-rand1+scra);
      // compute cos(theta) and sin(theta) from the sampled 1-cos(theta)
      curcost = 1.0-dum0;        // cos(theta)
      cursint = dum0*(2.0-dum0); // sin^2(theta)
      //
      // if we got current deflection that is not too small
      // then update cos(theta) sin(theta)
      if (cursint>1.0e-20) {
        cursint         = std::sqrt(cursint);
        G4double curphi = CLHEP::twopi*G4UniformRand();
        cost            = cost*curcost-sint*cursint*std::cos(curphi);
        sint            = std::sqrt(std::max(0.0, (1.0-cost)*(1.0+cost)));
      }
      //
      // check if we have done enough scattering i.e. sampling from the Poisson
      if (rand0<cumprob) {
        return;
      }
    }
    // if reached the max iter i.e. 10
    return;
  }
  //
  // multiple scattering case with lambdavalue >= 1:
  //   - use the precomputed and transformed Goudsmit-Saunderson angular
  //     distributions to sample 1.-cos(theta)
  G4double vrand[3];
  G4Random::getTheEngine()->flatArray(3,vrand);
  //double dum0 = SampleCosTheta(lambdaval, qval, scra, r1, r2, r3);
  double dum0 = SampleCosTheta(lambdaval, qval, scra, vrand[0], vrand[1], vrand[2]);
  // add protections
  if (dum0<0.0) dum0 = 0.0;
  if (dum0>2.0) dum0 = 2.0;
  // compute cos(theta) and sin(theta) from the sampled 1-cos(theta)
  cost = 1.0-dum0;
  sint = std::sqrt(dum0*(2.0-dum0));
}


G4double G4GoudsmitSaundersonTable::SampleCosTheta(G4double lambdaval, G4double qval, G4double scra, G4double rndm1, G4double rndm2, G4double rndm3) {
  // make sure that lambda = s/lambda_el is in [gLAMBMIN,gLAMBMAX)
  // : note that lambda<gLAMBMIN=1 is already handeled before so lambda>= gLAMBMIN for sure
  if (lambdaval>=gLAMBMAX) {
    lambdaval = gLAMBMAX-1.e-8;
  }
  // if we need the second grid
  if (qval>=gQMIN2) {
    return SampleCosTheta2(lambdaval, qval, scra, rndm1, rndm2, rndm3);
  } else {
    return SampleCosTheta1(lambdaval, qval, scra, rndm1, rndm2, rndm3);
  }
}



// for Q Grid 1
G4double G4GoudsmitSaundersonTable::SampleCosTheta1(G4double lambdaval, G4double qval, G4double scra, G4double rndm1, G4double rndm2, G4double rndm3) {
  // determine lower lambda index
  // logspacing in lambdaval => L=s/lambda_el
  G4double logLambda       = G4Log(lambdaval);
  G4double pLambdaIndxPlus = (logLambda-fLogLambda0)*fInvLogDeltaLambda;
  G4int    lambdaIndx      = (G4int)(pLambdaIndxPlus);    // lower index of the lambda bin
  pLambdaIndxPlus        = pLambdaIndxPlus-lambdaIndx;    // probability of taking the higher index distribution
  if (rndm1<pLambdaIndxPlus) {
    ++lambdaIndx;
  }
  // determine lower Q index
  // linear spacing in qval => Q=s/lambda_el G1
  // assuming that qval is in [gQMIN1, gQMAX1) !
  if (qval<gQMIN1) {
    qval = gQMIN1;
  }
  G4double pQIndxPlus = (qval-gQMIN1)*fInvDeltaQ1;
  G4int    qIndx      = (int)(pQIndxPlus);        // lower index of the Q bin
  pQIndxPlus        = pQIndxPlus-qIndx;
  if (rndm2<pQIndxPlus) {
    ++qIndx;
  }
  G4int indx      = lambdaIndx*gQNUM1+qIndx;
  GSMSCAngularDtr *gsDtr = fGSMSCAngularDistributions1[indx];
  // sampling form the selected distribution
  G4double delta  = 1.0/(gsDtr->fNumData-1.);
  // determine lower cumulative bin inidex
  G4int indxl     = rndm3/delta;
  G4double  aval  = rndm3-indxl*delta;
  G4double  dum0  = delta*aval;
  G4double  dum1  = (1.0+gsDtr->fParamA[indxl]+gsDtr->fParamB[indxl])*dum0;
  G4double  dum2  = delta*delta + gsDtr->fParamA[indxl]*dum0 + gsDtr->fParamB[indxl]*aval*aval;
  G4double sample = gsDtr->fUValues[indxl] +  dum1/dum2 *(gsDtr->fUValues[indxl+1]-gsDtr->fUValues[indxl]);
  // transform back u to cos(theta) :
  // -compute parameter value a = w2*screening_parameter
  if (lambdaval>10.0) {
    dum0 = 0.5*(-2.77164+logLambda*( 2.94874-logLambda*(0.1535754-logLambda*0.00552888) ));
  } else {
    dum0 = 0.5*(1.347+logLambda*(0.209364-logLambda*(0.45525-logLambda*(0.50142-logLambda*0.081234))));
  }
  G4double para = dum0*(lambdaval+4.0)*scra;
  // this is the sampled 1-cos(theta) = (2.0*para*sample)/(1.0-sample+para)
  return (2.0*para*sample)/(1.0-sample+para);
}

// for Q Grid 2
G4double G4GoudsmitSaundersonTable::SampleCosTheta2(G4double lambdaval, G4double qval, G4double scra, G4double rndm1, G4double rndm2, G4double rndm3) {
  // determine lower lambda index
  // logspacing in lambdaval => L=s/lambda_el
  G4double logLambda       = G4Log(lambdaval);
  G4double pLambdaIndxPlus = (logLambda-fLogLambda0)*fInvLogDeltaLambda;
  G4int    lambdaIndx      = (G4int)(pLambdaIndxPlus);        // lower index of the lambda bin
  pLambdaIndxPlus        = pLambdaIndxPlus-lambdaIndx;    // probability of taking the higher index distribution
  // protect: lower almbda might be smaller than Q (G1=Q/L should be <1)
  G4double lambdaScale     = 1.0;
  G4double logQ            = G4Log(qval);
  G4double lambdaLow       = fLogLambda0+lambdaIndx*fLogDeltaLambda;
  if (lambdaLow<logQ) {
    G4double min  = (logQ-lambdaLow);
    lambdaScale = (fLogDeltaLambda-min/pLambdaIndxPlus)/(fLogDeltaLambda-min);
  }
  pLambdaIndxPlus *= lambdaScale;
  // select lower or upper lambda grid point
  if (rndm1<pLambdaIndxPlus) {
    ++lambdaIndx;
  }
  // determine lower Q index
  // linear spacing in qval => Q=s/lambda_el G1
  // assuming that qval is in [gQMIN2, gQMAX2) !
  GSMSCAngularDtr *gsDtr = nullptr;
  if (qval<gQMAX2) {
    G4double pQIndxPlus = (qval-gQMIN2)*fInvDeltaQ2;
    G4int    qIndx      = (G4int)(pQIndxPlus);        // lower index of the Q bin
    pQIndxPlus        = pQIndxPlus-qIndx;
    G4int indx          = lambdaIndx*gQNUM2+qIndx;
    if (fGSMSCAngularDistributions2[indx]) {
      pQIndxPlus     *= fGSMSCAngularDistributions2[indx]->fQScale;
      if (rndm2<pQIndxPlus) {
        ++indx;
      }
      gsDtr = fGSMSCAngularDistributions2[indx];
    }
  }
  // sample from the selected distribution: nullptr means that cos(theta) is isotropic
  if (gsDtr) {
    //
    // sampling form the selected distribution
    G4double delta  = 1.0/(gsDtr->fNumData-1.);
    // determine lower cumulative bin inidex
    G4int indxl     = rndm3/delta;
    G4double  aval  = rndm3-indxl*delta;
    G4double  dum0  = delta*aval;
    G4double  dum1  = (1.0+gsDtr->fParamA[indxl]+gsDtr->fParamB[indxl])*dum0;
    G4double  dum2  = delta*delta + gsDtr->fParamA[indxl]*dum0 + gsDtr->fParamB[indxl]*aval*aval;
    G4double sample = gsDtr->fUValues[indxl] +  dum1/dum2 *(gsDtr->fUValues[indxl+1]-gsDtr->fUValues[indxl]);
    // transform back u to cos(theta) :
    // -compute parameter value a = w2*screening_parameter
    if (lambdaval>10.0) {
      dum0 = 0.5*(-2.77164+logLambda*( 2.94874-logLambda*(0.1535754-logLambda*0.00552888) ));
    } else {
      dum0 = 0.5*(1.347+logLambda*(0.209364-logLambda*(0.45525-logLambda*(0.50142-logLambda*0.081234))));
    }
    G4double para = dum0*(lambdaval+4.0)*scra;
    // this is the sampled 1-cos(theta) = (2.0*para*sample)/(1.0-sample+para)
    // so returns with 1-cos(theta)
    return (2.0*para*sample)/(1.0-sample+para);
  } else {
    // isotropic 1-cos(theta)
    return rndm3*2.0;
  }
}



G4double G4GoudsmitSaundersonTable::GetScreeningParam(G4double G1){
   if (G1<gG1Values1[0]) {
    return gSCRMIN1;
  } else if (G1<gG1Values2[0]) {
    G4double logG1 = G4Log(G1);
    G4int k = (G4int)((logG1-fLogG1FuncMin1)*fInvLogDeltaG1Func1);
    return G4Exp(gScrAValues1[k]*logG1)*gScrBValues1[k]; // log-log linear
  } else if (G1<gG1Values2[gNUMSCR2-1]) {
    G4double logG1 = G4Log(G1);
    G4int k = (G4int)((logG1-fLogG1FuncMin2)*fInvLogDeltaG1Func2);
    return G4Exp(gScrAValues2[k]*logG1)*gScrBValues2[k]; // log-log linear
  } else {
    return gSCRMAX2;
  }
}



void G4GoudsmitSaundersonTable::LoadMSCData() {
  char* path = getenv("G4LEDATA");
  if (!path) {
    G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
		FatalException,
		"Environment variable G4LEDATA not defined");
    return;
  }
  // 
  fGSMSCAngularDistributions1.resize(gLAMBNUM*gQNUM1);
  for (G4int il=0; il<gLAMBNUM; ++il) {
    char fname[512];
    sprintf(fname,"%s/msc_GS/GSGrid_1/gsDistr_%d",path,il);
    std::ifstream infile(fname,std::ios::in);
    if (!infile.is_open()) {
      char msgc[512];
      sprintf(msgc,"Cannot open file: %s .",fname);
      G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
	 	  FatalException, msgc);
      return;
    }
    for (G4int iq=0; iq<gQNUM1; ++iq) {
      GSMSCAngularDtr *gsd = new GSMSCAngularDtr();
      infile >> gsd->fNumData;
      gsd->fQScale  = 1.0;
      gsd->fUValues = new G4double[gsd->fNumData]();
      gsd->fParamA  = new G4double[gsd->fNumData]();
      gsd->fParamB  = new G4double[gsd->fNumData]();
      G4double ddummy;
      infile >> ddummy; infile >> ddummy;
      for (G4int i=0; i<gsd->fNumData; ++i) {
        infile >> gsd->fUValues[i];
        infile >> gsd->fParamA[i];
        infile >> gsd->fParamB[i];
      }
      fGSMSCAngularDistributions1[il*gQNUM1+iq] = gsd;
    }
    infile.close();
  }
  //
  // second grid
  fGSMSCAngularDistributions2.resize(gLAMBNUM*gQNUM2);
  for (G4int il=0; il<gLAMBNUM; ++il) {
    G4double lambda = G4Exp(fLogLambda0+fLogDeltaLambda*il);
    char fname[512];
    sprintf(fname,"%s/msc_GS/GSGrid_2/gsDistr_%d",path,il);
    std::ifstream infile(fname,std::ios::in);
    if (!infile.is_open()) {
      char msgc[512];
      sprintf(msgc,"Cannot open file: %s .",fname);
      G4Exception("G4GoudsmitSaundersonTable::LoadMSCData()","em0006",
	 	  FatalException, msgc);
      return;
    }
    for (G4int iq=0; iq<gQNUM2; ++iq) {
      G4double qDelta = fDeltaQ2; // indicates that = fDeltaQ2
      G4int numData;
      infile >> numData;
      if (numData>0) {
        G4double qval   = gQMIN2+iq*fDeltaQ2;
        if (iq<gQNUM2-1) {
          G4double qvalPlus = qval+fDeltaQ2;
          G4double g1   = qvalPlus/lambda; // it's G1 that should be <1.0
          if (g1>1.0) {
            qvalPlus = 1.0*lambda;
            qDelta   = qvalPlus-qval;
          }
        }
        GSMSCAngularDtr *gsd = new GSMSCAngularDtr();
        gsd->fNumData = numData;
        gsd->fQScale  = fDeltaQ2/qDelta;  // fDeltaQ2/qDelta : remain* this => probOfUpper bin
        gsd->fUValues = new G4double[gsd->fNumData]();
        gsd->fParamA  = new G4double[gsd->fNumData]();
        gsd->fParamB  = new G4double[gsd->fNumData]();
        double ddummy;
        infile >> ddummy; infile >> ddummy;
        for (G4int i=0; i<gsd->fNumData; ++i) {
          infile >> gsd->fUValues[i];
          infile >> gsd->fParamA[i];
          infile >> gsd->fParamB[i];
        }
        fGSMSCAngularDistributions2[il*gQNUM2+iq] = gsd;
      } else {
        fGSMSCAngularDistributions2[il*gQNUM2+iq] = nullptr;
      }
    }
    infile.close();
  }
}


////////////////////////////////////////////////////////////////////////////////
// compute material dependent Moliere MSC parameters at initialisation
void G4GoudsmitSaundersonTable::InitMoliereMSCParams() {
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

