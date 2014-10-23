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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLNaturalIsotopicDistributions.cc
 * \brief Classes that stores isotopic abundances
 *
 * \date 21st October 2012
 * \author Davide Mancusi
 */

#include "G4INCLNaturalIsotopicDistributions.hh"
#include "G4INCLRandom.hh"
#include "G4INCLLogger.hh"
// #include <cassert>
#include <utility>
#include <iostream>

namespace G4INCL {

  Isotope::Isotope(const G4int A, const G4double abundance) :
    theA(A),
    theAbundance(abundance)
  {}

  IsotopicDistribution::IsotopicDistribution(IsotopeVector const &aVector) :
    theIsotopes(aVector)
  {
    G4double previousAbundance = 0.;
    // Cumulate the abundances
    for(IsotopeIter i=theIsotopes.begin(), e=theIsotopes.end(); i!=e; ++i) {
      i->theAbundance += previousAbundance;
      previousAbundance = i->theAbundance;
    }
    // Normalize the abundances to 1
    const G4double normalisation = 1./theIsotopes.back().theAbundance;
    for(IsotopeIter i=theIsotopes.begin(), e=theIsotopes.end(); i!=e; ++i)
      i->theAbundance *= normalisation;
  }

  G4int IsotopicDistribution::drawRandomIsotope() const {
    const G4double r = Random::shoot();
    for(unsigned int i=0; i<theIsotopes.size()-1; ++i) {
      if(r<=theIsotopes.at(i).theAbundance)
        return theIsotopes.at(i).theA;
    }
    return theIsotopes.back().theA;
  }

  IsotopeVector const &IsotopicDistribution::getIsotopes() const {
    return theIsotopes;
  }

  IsotopicDistribution const &NaturalIsotopicDistributions::getIsotopicDistribution(G4int Z) const {
    std::map<G4int, IsotopicDistribution>::const_iterator i = theDistributions.find(Z);
    if(i!=theDistributions.end())
      return i->second;
    else {
      INCL_FATAL("Requested natural isotopic distribution for synthetic element Z = " << Z << '\n');
      return theDistributions.begin()->second;
    }
  }

  G4int NaturalIsotopicDistributions::drawRandomIsotope(G4int Z) const {
    return getIsotopicDistribution(Z).drawRandomIsotope();
  }

  namespace {
    std::pair<G4int, Isotope> theRawDistributions[] = {
      std::pair<G4int, Isotope>(1, Isotope(1, 99.985)),
      std::pair<G4int, Isotope>(1, Isotope(2, 0.015)),
      std::pair<G4int, Isotope>(2, Isotope(3, 0.000137)),
      std::pair<G4int, Isotope>(2, Isotope(4, 99.999863)),
      std::pair<G4int, Isotope>(3, Isotope(6, 7.5)),
      std::pair<G4int, Isotope>(3, Isotope(7, 92.5)),
      std::pair<G4int, Isotope>(4, Isotope(9, 100.0)),
      std::pair<G4int, Isotope>(5, Isotope(10, 19.9)),
      std::pair<G4int, Isotope>(5, Isotope(11, 80.1)),
      std::pair<G4int, Isotope>(6, Isotope(12, 98.90)),
      std::pair<G4int, Isotope>(6, Isotope(13, 1.10)),
      std::pair<G4int, Isotope>(7, Isotope(14, 99.634)),
      std::pair<G4int, Isotope>(7, Isotope(15, 0.366)),
      std::pair<G4int, Isotope>(8, Isotope(16, 99.762)),
      std::pair<G4int, Isotope>(8, Isotope(17, 0.038)),
      std::pair<G4int, Isotope>(8, Isotope(18, 0.200)),
      std::pair<G4int, Isotope>(9, Isotope(19, 100.0)),
      std::pair<G4int, Isotope>(10, Isotope(20, 90.48)),
      std::pair<G4int, Isotope>(10, Isotope(21, 0.27)),
      std::pair<G4int, Isotope>(10, Isotope(22, 9.25)),
      std::pair<G4int, Isotope>(11, Isotope(23, 100.0)),
      std::pair<G4int, Isotope>(12, Isotope(24, 78.99)),
      std::pair<G4int, Isotope>(12, Isotope(25, 10.00)),
      std::pair<G4int, Isotope>(12, Isotope(26, 11.01)),
      std::pair<G4int, Isotope>(13, Isotope(27, 100.0)),
      std::pair<G4int, Isotope>(14, Isotope(28, 92.23)),
      std::pair<G4int, Isotope>(14, Isotope(29, 4.67)),
      std::pair<G4int, Isotope>(14, Isotope(30, 3.10)),
      std::pair<G4int, Isotope>(15, Isotope(31, 100.0)),
      std::pair<G4int, Isotope>(16, Isotope(32, 95.02)),
      std::pair<G4int, Isotope>(16, Isotope(33, 0.75)),
      std::pair<G4int, Isotope>(16, Isotope(34, 4.21)),
      std::pair<G4int, Isotope>(16, Isotope(36, 0.02)),
      std::pair<G4int, Isotope>(17, Isotope(35, 75.77)),
      std::pair<G4int, Isotope>(17, Isotope(37, 24.23)),
      std::pair<G4int, Isotope>(18, Isotope(36, 0.337)),
      std::pair<G4int, Isotope>(18, Isotope(38, 0.063)),
      std::pair<G4int, Isotope>(18, Isotope(40, 99.600)),
      std::pair<G4int, Isotope>(19, Isotope(39, 93.2581)),
      std::pair<G4int, Isotope>(19, Isotope(40, 0.0117)),
      std::pair<G4int, Isotope>(19, Isotope(41, 6.7302)),
      std::pair<G4int, Isotope>(20, Isotope(40, 96.941)),
      std::pair<G4int, Isotope>(20, Isotope(42, 0.647)),
      std::pair<G4int, Isotope>(20, Isotope(43, 0.135)),
      std::pair<G4int, Isotope>(20, Isotope(44, 2.086)),
      std::pair<G4int, Isotope>(20, Isotope(46, 0.004)),
      std::pair<G4int, Isotope>(20, Isotope(48, 0.187)),
      std::pair<G4int, Isotope>(21, Isotope(45, 100.0)),
      std::pair<G4int, Isotope>(22, Isotope(46, 8.0)),
      std::pair<G4int, Isotope>(22, Isotope(47, 7.3)),
      std::pair<G4int, Isotope>(22, Isotope(48, 73.8)),
      std::pair<G4int, Isotope>(22, Isotope(49, 5.5)),
      std::pair<G4int, Isotope>(22, Isotope(50, 5.4)),
      std::pair<G4int, Isotope>(23, Isotope(50, 0.250)),
      std::pair<G4int, Isotope>(23, Isotope(51, 99.750)),
      std::pair<G4int, Isotope>(24, Isotope(50, 4.345)),
      std::pair<G4int, Isotope>(24, Isotope(52, 83.789)),
      std::pair<G4int, Isotope>(24, Isotope(53, 9.501)),
      std::pair<G4int, Isotope>(24, Isotope(54, 2.365)),
      std::pair<G4int, Isotope>(25, Isotope(55, 100.0)),
      std::pair<G4int, Isotope>(26, Isotope(54, 5.8)),
      std::pair<G4int, Isotope>(26, Isotope(56, 91.72)),
      std::pair<G4int, Isotope>(26, Isotope(57, 2.2)),
      std::pair<G4int, Isotope>(26, Isotope(58, 0.28)),
      std::pair<G4int, Isotope>(27, Isotope(59, 100.0)),
      std::pair<G4int, Isotope>(28, Isotope(58, 68.077)),
      std::pair<G4int, Isotope>(28, Isotope(60, 26.223)),
      std::pair<G4int, Isotope>(28, Isotope(61, 1.140)),
      std::pair<G4int, Isotope>(28, Isotope(62, 3.634)),
      std::pair<G4int, Isotope>(28, Isotope(64, 0.926)),
      std::pair<G4int, Isotope>(29, Isotope(63, 69.17)),
      std::pair<G4int, Isotope>(29, Isotope(65, 30.83)),
      std::pair<G4int, Isotope>(30, Isotope(64, 48.6)),
      std::pair<G4int, Isotope>(30, Isotope(66, 27.9)),
      std::pair<G4int, Isotope>(30, Isotope(67, 4.1)),
      std::pair<G4int, Isotope>(30, Isotope(68, 18.8)),
      std::pair<G4int, Isotope>(30, Isotope(70, 0.6)),
      std::pair<G4int, Isotope>(31, Isotope(69, 60.108)),
      std::pair<G4int, Isotope>(31, Isotope(71, 39.892)),
      std::pair<G4int, Isotope>(32, Isotope(70, 21.23)),
      std::pair<G4int, Isotope>(32, Isotope(72, 27.66)),
      std::pair<G4int, Isotope>(32, Isotope(73, 7.73)),
      std::pair<G4int, Isotope>(32, Isotope(74, 35.94)),
      std::pair<G4int, Isotope>(32, Isotope(76, 7.44)),
      std::pair<G4int, Isotope>(33, Isotope(75, 100.0)),
      std::pair<G4int, Isotope>(34, Isotope(74, 0.89)),
      std::pair<G4int, Isotope>(34, Isotope(76, 9.36)),
      std::pair<G4int, Isotope>(34, Isotope(77, 7.63)),
      std::pair<G4int, Isotope>(34, Isotope(78, 23.78)),
      std::pair<G4int, Isotope>(34, Isotope(80, 49.61)),
      std::pair<G4int, Isotope>(34, Isotope(82, 8.73)),
      std::pair<G4int, Isotope>(35, Isotope(79, 50.69)),
      std::pair<G4int, Isotope>(35, Isotope(81, 49.31)),
      std::pair<G4int, Isotope>(36, Isotope(78, 0.35)),
      std::pair<G4int, Isotope>(36, Isotope(80, 2.25)),
      std::pair<G4int, Isotope>(36, Isotope(82, 11.6)),
      std::pair<G4int, Isotope>(36, Isotope(83, 11.5)),
      std::pair<G4int, Isotope>(36, Isotope(84, 57.0)),
      std::pair<G4int, Isotope>(36, Isotope(86, 17.3)),
      std::pair<G4int, Isotope>(37, Isotope(85, 72.165)),
      std::pair<G4int, Isotope>(37, Isotope(87, 27.835)),
      std::pair<G4int, Isotope>(38, Isotope(84, 0.56)),
      std::pair<G4int, Isotope>(38, Isotope(86, 9.86)),
      std::pair<G4int, Isotope>(38, Isotope(87, 7.00)),
      std::pair<G4int, Isotope>(38, Isotope(88, 82.58)),
      std::pair<G4int, Isotope>(39, Isotope(89, 100.0)),
      std::pair<G4int, Isotope>(40, Isotope(90, 51.45)),
      std::pair<G4int, Isotope>(40, Isotope(91, 11.22)),
      std::pair<G4int, Isotope>(40, Isotope(92, 17.15)),
      std::pair<G4int, Isotope>(40, Isotope(94, 17.38)),
      std::pair<G4int, Isotope>(40, Isotope(96, 2.80)),
      std::pair<G4int, Isotope>(41, Isotope(93, 100.0)),
      std::pair<G4int, Isotope>(42, Isotope(92, 14.84)),
      std::pair<G4int, Isotope>(42, Isotope(94, 9.25)),
      std::pair<G4int, Isotope>(42, Isotope(95, 15.92)),
      std::pair<G4int, Isotope>(42, Isotope(96, 16.68)),
      std::pair<G4int, Isotope>(42, Isotope(97, 9.55)),
      std::pair<G4int, Isotope>(42, Isotope(98, 24.13)),
      std::pair<G4int, Isotope>(42, Isotope(100, 9.63)),
      std::pair<G4int, Isotope>(44, Isotope(96, 5.52)),
      std::pair<G4int, Isotope>(44, Isotope(98, 1.88)),
      std::pair<G4int, Isotope>(44, Isotope(99, 12.7)),
      std::pair<G4int, Isotope>(44, Isotope(100, 12.6)),
      std::pair<G4int, Isotope>(44, Isotope(101, 17.0)),
      std::pair<G4int, Isotope>(44, Isotope(102, 31.6)),
      std::pair<G4int, Isotope>(44, Isotope(104, 18.7)),
      std::pair<G4int, Isotope>(45, Isotope(103, 100.0)),
      std::pair<G4int, Isotope>(46, Isotope(102, 1.02)),
      std::pair<G4int, Isotope>(46, Isotope(104, 11.14)),
      std::pair<G4int, Isotope>(46, Isotope(105, 22.33)),
      std::pair<G4int, Isotope>(46, Isotope(106, 27.33)),
      std::pair<G4int, Isotope>(46, Isotope(108, 26.46)),
      std::pair<G4int, Isotope>(46, Isotope(110, 11.72)),
      std::pair<G4int, Isotope>(47, Isotope(107, 51.839)),
      std::pair<G4int, Isotope>(47, Isotope(109, 48.161)),
      std::pair<G4int, Isotope>(48, Isotope(106, 1.25)),
      std::pair<G4int, Isotope>(48, Isotope(108, 0.89)),
      std::pair<G4int, Isotope>(48, Isotope(110, 12.49)),
      std::pair<G4int, Isotope>(48, Isotope(111, 12.80)),
      std::pair<G4int, Isotope>(48, Isotope(112, 24.13)),
      std::pair<G4int, Isotope>(48, Isotope(113, 12.22)),
      std::pair<G4int, Isotope>(48, Isotope(114, 28.73)),
      std::pair<G4int, Isotope>(48, Isotope(116, 7.49)),
      std::pair<G4int, Isotope>(49, Isotope(113, 4.3)),
      std::pair<G4int, Isotope>(49, Isotope(115, 95.7)),
      std::pair<G4int, Isotope>(50, Isotope(112, 0.97)),
      std::pair<G4int, Isotope>(50, Isotope(114, 0.65)),
      std::pair<G4int, Isotope>(50, Isotope(115, 0.34)),
      std::pair<G4int, Isotope>(50, Isotope(116, 14.53)),
      std::pair<G4int, Isotope>(50, Isotope(117, 7.68)),
      std::pair<G4int, Isotope>(50, Isotope(118, 24.23)),
      std::pair<G4int, Isotope>(50, Isotope(119, 8.59)),
      std::pair<G4int, Isotope>(50, Isotope(120, 32.59)),
      std::pair<G4int, Isotope>(50, Isotope(122, 4.63)),
      std::pair<G4int, Isotope>(50, Isotope(124, 5.79)),
      std::pair<G4int, Isotope>(51, Isotope(121, 57.36)),
      std::pair<G4int, Isotope>(51, Isotope(123, 42.64)),
      std::pair<G4int, Isotope>(52, Isotope(120, 0.096)),
      std::pair<G4int, Isotope>(52, Isotope(122, 2.603)),
      std::pair<G4int, Isotope>(52, Isotope(123, 0.908)),
      std::pair<G4int, Isotope>(52, Isotope(124, 4.816)),
      std::pair<G4int, Isotope>(52, Isotope(125, 7.139)),
      std::pair<G4int, Isotope>(52, Isotope(126, 18.95)),
      std::pair<G4int, Isotope>(52, Isotope(128, 31.69)),
      std::pair<G4int, Isotope>(52, Isotope(130, 33.80)),
      std::pair<G4int, Isotope>(53, Isotope(127, 100.0)),
      std::pair<G4int, Isotope>(54, Isotope(124, 0.10)),
      std::pair<G4int, Isotope>(54, Isotope(126, 0.09)),
      std::pair<G4int, Isotope>(54, Isotope(128, 1.91)),
      std::pair<G4int, Isotope>(54, Isotope(129, 26.4)),
      std::pair<G4int, Isotope>(54, Isotope(130, 4.1)),
      std::pair<G4int, Isotope>(54, Isotope(131, 21.2)),
      std::pair<G4int, Isotope>(54, Isotope(132, 26.9)),
      std::pair<G4int, Isotope>(54, Isotope(134, 10.4)),
      std::pair<G4int, Isotope>(54, Isotope(136, 8.9)),
      std::pair<G4int, Isotope>(55, Isotope(133, 100.0)),
      std::pair<G4int, Isotope>(56, Isotope(130, 0.106)),
      std::pair<G4int, Isotope>(56, Isotope(132, 0.101)),
      std::pair<G4int, Isotope>(56, Isotope(134, 2.417)),
      std::pair<G4int, Isotope>(56, Isotope(135, 6.592)),
      std::pair<G4int, Isotope>(56, Isotope(136, 7.854)),
      std::pair<G4int, Isotope>(56, Isotope(137, 11.23)),
      std::pair<G4int, Isotope>(56, Isotope(138, 71.70)),
      std::pair<G4int, Isotope>(57, Isotope(138, 0.0902)),
      std::pair<G4int, Isotope>(57, Isotope(139, 99.9098)),
      std::pair<G4int, Isotope>(58, Isotope(136, 0.19)),
      std::pair<G4int, Isotope>(58, Isotope(138, 0.25)),
      std::pair<G4int, Isotope>(58, Isotope(140, 88.48)),
      std::pair<G4int, Isotope>(58, Isotope(142, 11.08)),
      std::pair<G4int, Isotope>(59, Isotope(141, 100.0)),
      std::pair<G4int, Isotope>(60, Isotope(142, 27.13)),
      std::pair<G4int, Isotope>(60, Isotope(143, 12.18)),
      std::pair<G4int, Isotope>(60, Isotope(144, 23.80)),
      std::pair<G4int, Isotope>(60, Isotope(145, 8.30)),
      std::pair<G4int, Isotope>(60, Isotope(146, 17.19)),
      std::pair<G4int, Isotope>(60, Isotope(148, 5.76)),
      std::pair<G4int, Isotope>(60, Isotope(150, 5.64)),
      std::pair<G4int, Isotope>(62, Isotope(144, 3.1)),
      std::pair<G4int, Isotope>(62, Isotope(147, 15.0)),
      std::pair<G4int, Isotope>(62, Isotope(148, 11.3)),
      std::pair<G4int, Isotope>(62, Isotope(149, 13.8)),
      std::pair<G4int, Isotope>(62, Isotope(150, 7.4)),
      std::pair<G4int, Isotope>(62, Isotope(152, 26.7)),
      std::pair<G4int, Isotope>(62, Isotope(154, 22.7)),
      std::pair<G4int, Isotope>(63, Isotope(151, 47.8)),
      std::pair<G4int, Isotope>(63, Isotope(153, 52.2)),
      std::pair<G4int, Isotope>(64, Isotope(152, 0.20)),
      std::pair<G4int, Isotope>(64, Isotope(154, 2.18)),
      std::pair<G4int, Isotope>(64, Isotope(155, 14.80)),
      std::pair<G4int, Isotope>(64, Isotope(156, 20.47)),
      std::pair<G4int, Isotope>(64, Isotope(157, 15.65)),
      std::pair<G4int, Isotope>(64, Isotope(158, 24.84)),
      std::pair<G4int, Isotope>(64, Isotope(160, 21.86)),
      std::pair<G4int, Isotope>(65, Isotope(159, 100.0)),
      std::pair<G4int, Isotope>(66, Isotope(156, 0.06)),
      std::pair<G4int, Isotope>(66, Isotope(158, 0.10)),
      std::pair<G4int, Isotope>(66, Isotope(160, 2.34)),
      std::pair<G4int, Isotope>(66, Isotope(161, 18.9)),
      std::pair<G4int, Isotope>(66, Isotope(162, 25.5)),
      std::pair<G4int, Isotope>(66, Isotope(163, 24.9)),
      std::pair<G4int, Isotope>(66, Isotope(164, 28.2)),
      std::pair<G4int, Isotope>(67, Isotope(165, 100.0)),
      std::pair<G4int, Isotope>(68, Isotope(162, 0.14)),
      std::pair<G4int, Isotope>(68, Isotope(164, 1.61)),
      std::pair<G4int, Isotope>(68, Isotope(166, 33.6)),
      std::pair<G4int, Isotope>(68, Isotope(167, 22.95)),
      std::pair<G4int, Isotope>(68, Isotope(168, 26.8)),
      std::pair<G4int, Isotope>(68, Isotope(170, 14.9)),
      std::pair<G4int, Isotope>(69, Isotope(169, 100.0)),
      std::pair<G4int, Isotope>(70, Isotope(168, 0.13)),
      std::pair<G4int, Isotope>(70, Isotope(170, 3.05)),
      std::pair<G4int, Isotope>(70, Isotope(171, 14.3)),
      std::pair<G4int, Isotope>(70, Isotope(172, 21.9)),
      std::pair<G4int, Isotope>(70, Isotope(173, 16.12)),
      std::pair<G4int, Isotope>(70, Isotope(174, 31.8)),
      std::pair<G4int, Isotope>(70, Isotope(176, 12.7)),
      std::pair<G4int, Isotope>(71, Isotope(175, 97.41)),
      std::pair<G4int, Isotope>(71, Isotope(176, 2.59)),
      std::pair<G4int, Isotope>(72, Isotope(174, 0.162)),
      std::pair<G4int, Isotope>(72, Isotope(176, 5.206)),
      std::pair<G4int, Isotope>(72, Isotope(177, 18.606)),
      std::pair<G4int, Isotope>(72, Isotope(178, 27.297)),
      std::pair<G4int, Isotope>(72, Isotope(179, 13.629)),
      std::pair<G4int, Isotope>(72, Isotope(180, 35.100)),
      std::pair<G4int, Isotope>(73, Isotope(180, 0.012)),
      std::pair<G4int, Isotope>(73, Isotope(181, 99.988)),
      std::pair<G4int, Isotope>(74, Isotope(180, 0.13)),
      std::pair<G4int, Isotope>(74, Isotope(182, 26.3)),
      std::pair<G4int, Isotope>(74, Isotope(183, 14.3)),
      std::pair<G4int, Isotope>(74, Isotope(184, 30.67)),
      std::pair<G4int, Isotope>(74, Isotope(186, 28.6)),
      std::pair<G4int, Isotope>(75, Isotope(185, 37.40)),
      std::pair<G4int, Isotope>(75, Isotope(187, 62.60)),
      std::pair<G4int, Isotope>(76, Isotope(184, 0.02)),
      std::pair<G4int, Isotope>(76, Isotope(186, 1.58)),
      std::pair<G4int, Isotope>(76, Isotope(187, 1.6)),
      std::pair<G4int, Isotope>(76, Isotope(188, 13.3)),
      std::pair<G4int, Isotope>(76, Isotope(189, 16.1)),
      std::pair<G4int, Isotope>(76, Isotope(190, 26.4)),
      std::pair<G4int, Isotope>(76, Isotope(192, 41.0)),
      std::pair<G4int, Isotope>(77, Isotope(191, 37.3)),
      std::pair<G4int, Isotope>(77, Isotope(193, 62.7)),
      std::pair<G4int, Isotope>(78, Isotope(190, 0.01)),
      std::pair<G4int, Isotope>(78, Isotope(192, 0.79)),
      std::pair<G4int, Isotope>(78, Isotope(194, 32.9)),
      std::pair<G4int, Isotope>(78, Isotope(195, 33.8)),
      std::pair<G4int, Isotope>(78, Isotope(196, 25.3)),
      std::pair<G4int, Isotope>(78, Isotope(198, 7.2)),
      std::pair<G4int, Isotope>(79, Isotope(197, 100.0)),
      std::pair<G4int, Isotope>(80, Isotope(196, 0.15)),
      std::pair<G4int, Isotope>(80, Isotope(198, 9.97)),
      std::pair<G4int, Isotope>(80, Isotope(199, 16.87)),
      std::pair<G4int, Isotope>(80, Isotope(200, 23.10)),
      std::pair<G4int, Isotope>(80, Isotope(201, 13.18)),
      std::pair<G4int, Isotope>(80, Isotope(202, 29.86)),
      std::pair<G4int, Isotope>(80, Isotope(204, 6.87)),
      std::pair<G4int, Isotope>(81, Isotope(203, 29.524)),
      std::pair<G4int, Isotope>(81, Isotope(205, 70.476)),
      std::pair<G4int, Isotope>(82, Isotope(204, 1.4)),
      std::pair<G4int, Isotope>(82, Isotope(206, 24.1)),
      std::pair<G4int, Isotope>(82, Isotope(207, 22.1)),
      std::pair<G4int, Isotope>(82, Isotope(208, 52.4)),
      std::pair<G4int, Isotope>(83, Isotope(209, 100.0)),
      std::pair<G4int, Isotope>(90, Isotope(232, 100.0)),
      std::pair<G4int, Isotope>(92, Isotope(234, 0.0055)),
      std::pair<G4int, Isotope>(92, Isotope(235, 0.7200)),
      std::pair<G4int, Isotope>(92, Isotope(238, 99.2745))
    };

    // Cool hack to get the size of an array in C++
    template<typename T, ::std::size_t N> ::std::size_t sizeOfArray(const T(&)[ N ] ) {
      return N;
    }
  }

  NaturalIsotopicDistributions::NaturalIsotopicDistributions() {
    G4int oldZ = -1;
    IsotopeVector aVector;
    for(unsigned int i=0; i<sizeOfArray(theRawDistributions); ++i) {
      std::pair<G4int, Isotope> const &aPair = theRawDistributions[i];
      if(aPair.first == oldZ) {
        aVector.push_back(aPair.second);
      } else {
        if(oldZ!=-1)
          theDistributions.insert(std::pair<G4int, IsotopicDistribution>(oldZ, IsotopicDistribution(aVector)));
        oldZ = aPair.first;
        aVector.clear();
        aVector.push_back(aPair.second);
      }
    }
    // last element
    theDistributions.insert(std::pair<G4int, IsotopicDistribution>(oldZ, IsotopicDistribution(aVector)));
  }

}

