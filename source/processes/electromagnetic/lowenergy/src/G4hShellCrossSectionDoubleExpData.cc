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
// Author: Simona Saliceti (simona.saliceti@ge.infn.it)
//
// History:
// -----------
// 22 Apr 2004  First committed to cvs
//
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionDoubleExpData.cc,v 1.5 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4hShellCrossSectionDoubleExpData.hh"

G4hShellCrossSectionDoubleExpData::G4hShellCrossSectionDoubleExpData()
{

  FillVectorValuesEnergy();
  FillVectorValuesPar1();
  FillVectorValuesPar2(); 
  FillParameterMapEnergy();
  FillParameterMapPar1(); 
  FillParameterMapPar2();

}

G4hShellCrossSectionDoubleExpData::~G4hShellCrossSectionDoubleExpData()
{ 

  std::map< G4int,std::vector<G4double>*,std::less<G4int> >::iterator pos1;
  for (pos1 = parameterMapEnergy.begin(); pos1 != parameterMapEnergy.end(); ++pos1)
    {
      std::vector<G4double>* dataSet1 = (*pos1).second;
      delete dataSet1;
    }

  std::map< G4int,std::vector<G4double>*,std::less<G4int> >::iterator pos2;
  for (pos2 = parameterMapPar1.begin(); pos2 != parameterMapPar1.end(); ++pos2)
    {
      std::vector<G4double>* dataSet2 = (*pos2).second;
      delete dataSet2;
    }

  std::map< G4int,std::vector<G4double>*,std::less<G4int> >::iterator pos3;
  for (pos3 = parameterMapPar2.begin(); pos3 != parameterMapPar2.end(); ++pos3)
    {
      std::vector<G4double>* dataSet3 = (*pos3).second;
      delete dataSet3;
    }
}

inline void G4hShellCrossSectionDoubleExpData::InitializeVectorEnergy(std::vector<G4double> &vectEnergy, G4double value) const
{
  vectEnergy.push_back(value);
}

void G4hShellCrossSectionDoubleExpData::FillVectorValuesEnergy()
{
  InitializeVectorEnergy(energy6C,  0.35);
  InitializeVectorEnergy(energy7N,  0.6);
  InitializeVectorEnergy(energy8O,  0.7);
  InitializeVectorEnergy(energy9F,  1.0);
  InitializeVectorEnergy(energy10Ne,1.0);
  InitializeVectorEnergy(energy11Na,1.0);
  InitializeVectorEnergy(energy12Mg,2.0);
  InitializeVectorEnergy(energy13Al,2.0);
  InitializeVectorEnergy(energy14Si,3.0);
  InitializeVectorEnergy(energy15P, 3.0);
  InitializeVectorEnergy(energy16S, 5.0);
  InitializeVectorEnergy(energy17Cl,5.0);
  InitializeVectorEnergy(energy18Ar,5.0);
  InitializeVectorEnergy(energy19K, 5.0);
  InitializeVectorEnergy(energy20Ca,5.0);
  InitializeVectorEnergy(energy21Sc,5.0);
  InitializeVectorEnergy(energy22Ti,6.0);
  InitializeVectorEnergy(energy23V, 6.0);
  InitializeVectorEnergy(energy24Cr,6.0);
  InitializeVectorEnergy(energy25Mn,6.0);
  InitializeVectorEnergy(energy26Fe,3.0);
  InitializeVectorEnergy(energy27Co,3.5);
  InitializeVectorEnergy(energy28Ni,3.5);
  InitializeVectorEnergy(energy29Cu,4.0);
  InitializeVectorEnergy(energy30Zn,4.0);
  InitializeVectorEnergy(energy31Ga,4.5);
  InitializeVectorEnergy(energy32Ge,5.0);
  InitializeVectorEnergy(energy33As,5.0);
  InitializeVectorEnergy(energy34Se,5.5);
  InitializeVectorEnergy(energy35Br,6.0);
  InitializeVectorEnergy(energy36Kr,6.0);
  InitializeVectorEnergy(energy37Rb,6.0);
  InitializeVectorEnergy(energy38Sr,7.0);
  InitializeVectorEnergy(energy39Y, 7.0);
  InitializeVectorEnergy(energy40Zr,8.0);
  InitializeVectorEnergy(energy41Nb,8.0);
  InitializeVectorEnergy(energy42Mo,8.0);
  InitializeVectorEnergy(energy43Tc,9.0);
  InitializeVectorEnergy(energy44Ru,9.0);
  InitializeVectorEnergy(energy45Rh,9.0);
  InitializeVectorEnergy(energy46Pd,9.0);
  InitializeVectorEnergy(energy47Ag,10.);
  InitializeVectorEnergy(energy48Cd,10.);
  InitializeVectorEnergy(energy49In,11.);
  InitializeVectorEnergy(energy50Sn,12.);
  InitializeVectorEnergy(energy51Sb,12.);
  InitializeVectorEnergy(energy52Te,12.);
  InitializeVectorEnergy(energy53I, 14.);
  InitializeVectorEnergy(energy54Xe,14.);
  InitializeVectorEnergy(energy55Cs,16.);
  InitializeVectorEnergy(energy56Ba,16.);
  InitializeVectorEnergy(energy57La,16.);
  InitializeVectorEnergy(energy58Ce,18.);
  InitializeVectorEnergy(energy59Pr,18.);
  InitializeVectorEnergy(energy60Nd,20.);
  InitializeVectorEnergy(energy61Pm,20.);
  InitializeVectorEnergy(energy62Sm,20.);
  InitializeVectorEnergy(energy63Eu,23.);
  InitializeVectorEnergy(energy64Gd,23.);
  InitializeVectorEnergy(energy65Tb,23.);
  InitializeVectorEnergy(energy66Dy,23.);
  InitializeVectorEnergy(energy67Ho,23.);
  InitializeVectorEnergy(energy68Er,23.);
  InitializeVectorEnergy(energy69Tm,26.);
  InitializeVectorEnergy(energy70Yb,30.);
  InitializeVectorEnergy(energy71Lu,30.);
  InitializeVectorEnergy(energy72Hf,30.);
  InitializeVectorEnergy(energy73Ta,26.);
  InitializeVectorEnergy(energy74W, 30.);
  InitializeVectorEnergy(energy75Re,30.);
  InitializeVectorEnergy(energy76Os,30.);
  InitializeVectorEnergy(energy77Ir,30.);
  InitializeVectorEnergy(energy78Pt,35.);
  InitializeVectorEnergy(energy79Au,35.);
  InitializeVectorEnergy(energy80Hg,35.);
  InitializeVectorEnergy(energy81Tl,35.);
  InitializeVectorEnergy(energy82Pb,40.);
  InitializeVectorEnergy(energy83Bi,45.);
  InitializeVectorEnergy(energy84Po,45.);
  InitializeVectorEnergy(energy85At,45.);
  InitializeVectorEnergy(energy86Rn,45.);
  InitializeVectorEnergy(energy87Fr,45.);
  InitializeVectorEnergy(energy88Ra,45.);
  InitializeVectorEnergy(energy89Ac,50.);
  InitializeVectorEnergy(energy90Th,50.);
  InitializeVectorEnergy(energy91Pa,50.);
  InitializeVectorEnergy(energy92U, 55.);
}                                                                              

inline void G4hShellCrossSectionDoubleExpData::InitializeVectorPar1(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3) const
{
  vect.push_back(value1);
  vect.push_back(value2);
  vect.push_back(value3);
}

// First parameters

void G4hShellCrossSectionDoubleExpData::FillVectorValuesPar1()
{ 
  InitializeVectorPar1(parlow6C,  6.028754       , 22.76324      , 8.217967     );
  InitializeVectorPar1(parlow7N,  5.477011       , 18.81141      , 4.858727     );
  InitializeVectorPar1(parlow8O,  5.425050       , 16.80088      , 3.910342     );
  InitializeVectorPar1(parlow9F,  5.017089       , 14.56733      , 2.689939     );
  InitializeVectorPar1(parlow10Ne,4.668802       , 12.75720      , 1.897390     );
  InitializeVectorPar1(parlow11Na,4.705349       , 11.54904      , 1.700995     );
  InitializeVectorPar1(parlow12Mg,4.263983       , 10.00729      , 1.117699     );
  InitializeVectorPar1(parlow13Al,4.324701       , 9.009748      , 1.013421     );
  InitializeVectorPar1(parlow14Si,4.090029       , 7.939616      , 0.7610449    );
  InitializeVectorPar1(parlow15P, 4.098163       , 7.072802      , 0.6808610    );
  InitializeVectorPar1(parlow16S, 3.632073       , 6.268209      , 0.4380426    );
  InitializeVectorPar1(parlow17Cl,3.769954       , 5.545212      , 0.4345047    );
  InitializeVectorPar1(parlow18Ar,3.803460       , 4.852357      , 0.4018530    );
  InitializeVectorPar1(parlow19K, 3.818340       , 4.215003      , 0.3707851    );
  InitializeVectorPar1(parlow20Ca,3.767600       , 3.621560      , 0.3257078    );
  InitializeVectorPar1(parlow21Sc,3.851859       , 3.047405      , 0.3256020    );
  InitializeVectorPar1(parlow22Ti,3.739032       , 2.540971      , 0.2726244    );
  InitializeVectorPar1(parlow23V, 3.762692       , 2.028288      , 0.2595049    );
  InitializeVectorPar1(parlow24Cr,3.830856       , 1.517501      , 0.2585256    );
  InitializeVectorPar1(parlow25Mn,3.889344       , 1.032337      , 0.2568166    );
  InitializeVectorPar1(parlow26Fe,4.610710686    , 0.275389472   , 2.764554447  );
  InitializeVectorPar1(parlow27Co,3.55117935     , 0.228178229   , 2.702047259  );
  InitializeVectorPar1(parlow28Ni,2.334937582    , 0.232156421   , 2.738149196  );
  InitializeVectorPar1(parlow29Cu,1.915920665    , 0.182589984   , 2.664081065  );
  InitializeVectorPar1(parlow30Zn,1.258082776    , 0.190338228   , 2.710680405  );
  InitializeVectorPar1(parlow31Ga,1.090692834    , 0.140413681   , 2.624294005  );
  InitializeVectorPar1(parlow32Ge,0.868211657    , 0.115779844   , 2.587823968  );
  InitializeVectorPar1(parlow33As,0.588895163    , 0.121427557   , 2.630013394  );
  InitializeVectorPar1(parlow34Se,0.464947416    , 0.103777355   , 2.60655728   );
  InitializeVectorPar1(parlow35Br,0.372685678    , 0.087649176   , 2.580702817  );
  InitializeVectorPar1(parlow36Kr,0.261618542    , 0.071657004   , 2.628062571  );
  InitializeVectorPar1(parlow37Rb,0.203966062    , 0.063471976   , 2.621118231  );
  InitializeVectorPar1(parlow38Sr,0.151309868    , 0.062650983   , 2.635670955  );
  InitializeVectorPar1(parlow39Y, 0.130831069    , 0.046546373   , 2.596438847  );
  InitializeVectorPar1(parlow40Zr,0.113465174    , 0.034961204   , 2.557081599  );
  InitializeVectorPar1(parlow41Nb,0.089305327    , 0.031631918   , 2.560378036  );
  InitializeVectorPar1(parlow42Mo,0.066940214    , 0.032387014   , 2.580578224  );
  InitializeVectorPar1(parlow43Tc,0.054005087    , 0.028625811   , 2.578181335  );
  InitializeVectorPar1(parlow44Ru,0.047435043    , 0.021858688   , 2.544955569  );
  InitializeVectorPar1(parlow45Rh,0.037036289    , 0.021226071   , 2.557147432  );
  InitializeVectorPar1(parlow46Pd,0.031108855    , 0.017082016   , 2.560638611  );
  InitializeVectorPar1(parlow47Ag,0.030039299    , 0.010687233   , 2.501552766  );
  InitializeVectorPar1(parlow48Cd,0.025896549    , 0.008577291   , 2.485085985  );
  InitializeVectorPar1(parlow49In,0.021622343    , 0.007503402   , 2.479363645  );
  InitializeVectorPar1(parlow50Sn,0.019658368    , 0.005372974   , 2.448647025  );
  InitializeVectorPar1(parlow51Sb,0.016726542    , 0.004540107   , 2.440841502  );
  InitializeVectorPar1(parlow52Te,0.01210377     , 0.005793226   , 2.487112403  );
  InitializeVectorPar1(parlow53I, 0.012068341    , 0.003340949   , 2.430469903  );
  InitializeVectorPar1(parlow54Xe,0.011309929    , 0.002247705   , 2.39653105   );
  InitializeVectorPar1(parlow55Cs,0.010322464    , 0.001584264   , 2.372751551  );
  InitializeVectorPar1(parlow56Ba,0.00891645     , 0.001663984   , 2.365104652  );
  InitializeVectorPar1(parlow57La,0.008600326    , 0.000973706   , 2.330302922  );
  InitializeVectorPar1(parlow58Ce,0.00852455     , 0.000554424   , 2.286062759  );
  InitializeVectorPar1(parlow59Pr,0.006464573    , 0.000751099   , 2.320217715  );
  InitializeVectorPar1(parlow60Nd,0.005772215    , 0.000604636   , 2.308998547  );
  InitializeVectorPar1(parlow61Pm,0.005439364    , 0.000395253   , 2.284991617  );
  InitializeVectorPar1(parlow62Sm,0.004140155    , 0.000572248   , 2.319931149  );
  InitializeVectorPar1(parlow63Eu,0.004391028    , 0.00026105    , 2.263027478  );
  InitializeVectorPar1(parlow64Gd,0.004108655    , 0.00018564    , 2.242594759  );
  InitializeVectorPar1(parlow65Tb,0.003335626    , 0.000214331   , 2.262351386  );
  InitializeVectorPar1(parlow66Dy,0.003772419844 , 1.491354113E-4, 2.182513765  );
  InitializeVectorPar1(parlow67Ho,0.003171261268 , 1.610171391E-4, 2.19445375   );
  InitializeVectorPar1(parlow68Er,0.003511231198 , 6.431951513E-5, 2.131309482  );
  InitializeVectorPar1(parlow69Tm,0.002951671258 , 7.185528546E-5, 2.14407165   );
  InitializeVectorPar1(parlow70Yb,0.003090003974 , 3.232352598E-5, 2.10040042   );
  InitializeVectorPar1(parlow71Lu,0.002504560425 , 4.55286561E-5 , 2.122593025  );
  InitializeVectorPar1(parlow72Hf,0.002634273859 , 2.025018154E-5, 2.078762159  );
  InitializeVectorPar1(parlow73Ta,0.002246920523 , 2.530957413E-5, 2.087036821  );
  InitializeVectorPar1(parlow74W, 0.002138140746 , 1.685058881E-5, 2.071640226  );
  InitializeVectorPar1(parlow75Re,0.002053720518 , 1.1699924E-5  , 2.052411395  );
  InitializeVectorPar1(parlow76Os,0.001609160281 , 1.98304501E-5 , 2.07932865   );
  InitializeVectorPar1(parlow77Ir,0.001451780286 , 1.89647459E-5 , 2.077028977  );
  InitializeVectorPar1(parlow78Pt,0.1484928E-2   , 0.1015788E-4  , 2.043982     );
  InitializeVectorPar1(parlow79Au,0.1401423E-2   , 0.7618721E-5  , 2.031677     );
  InitializeVectorPar1(parlow80Hg,0.1279221E-2   , 0.7684558E-5  , 2.026922     );
  InitializeVectorPar1(parlow81Tl,0.001116529691 , 8.665290032E-6, 2.035250841  );
  InitializeVectorPar1(parlow82Pb,0.00123546161  , 3.264194012E-6, 1.984513455  );
  InitializeVectorPar1(parlow83Bi,0.1168456E-2   , 0.2565728E-5  , 1.974040     );
  InitializeVectorPar1(parlow84Po,0.001161495281 , 1.517549601E-6, 1.951280331  );
  InitializeVectorPar1(parlow85At,0.9574954E-3   , 0.2619334E-5  , 1.974164     );
  InitializeVectorPar1(parlow86Rn,0.9172934E-3   , 0.9124897E-6  , 1.966786     );
  InitializeVectorPar1(parlow87Fr,0.7900983E-3   , 0.1254855E-5  , 1.978866     );
  InitializeVectorPar1(parlow88Ra,0.7093558E-3   , 0.1436426E-5  , 1.981631     );
  InitializeVectorPar1(parlow89Ac,0.0007906194724, 4.455647522E-7, 1.934485393  );
  InitializeVectorPar1(parlow90Th,0.0007280666063, 4.668165706E-7, 1.931837161  );
  InitializeVectorPar1(parlow91Pa,0.0006494787919, 5.572898409E-7, 1.937018718  );
  InitializeVectorPar1(parlow92U, 0.0006712361299, 2.734932326E-7, 1.908976604  );
}

inline void G4hShellCrossSectionDoubleExpData::InitializeVectorPar2(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3, G4double value4, G4double value5) const
{
  vect.push_back(value1);
  vect.push_back(value2);
  vect.push_back(value3);
  vect.push_back(value4);
  vect.push_back(value5);
}

// Second parameters

void G4hShellCrossSectionDoubleExpData::FillVectorValuesPar2()
{
  InitializeVectorPar2(parhigh6C,  0.2614505E+08,  0.1960089    ,  -0.9781201,    0.0,  	0.0);
  InitializeVectorPar2(parhigh7N,  0.1469968E+08,  0.9198590E-01,  -0.9035943,    0.0,  	0.0);
  InitializeVectorPar2(parhigh8O,  9015526.0    ,  0.4561001E-01,  -0.8435532,    0.0,  	0.0);
  InitializeVectorPar2(parhigh9F,  5878802.0    ,  0.2042609E-01,  -0.8003597,    0.0,  	0.0);
  InitializeVectorPar2(parhigh10Ne,3951196.0    ,  0.8858837E-02,  -0.7637176,    0.0,  	0.0);
  InitializeVectorPar2(parhigh11Na,2838478.0    ,  0.3098775E-02,  -0.7451933,    0.0,  	0.0);
  InitializeVectorPar2(parhigh12Mg,2322929.0    ,  0.6889491E-03,  -0.7595924,    0.0,  	0.0);
  InitializeVectorPar2(parhigh13Al,1739866.0    ,  0.1930171E-03,  -0.7445110,    0.0,  	0.0);
  InitializeVectorPar2(parhigh14Si,1433611.0    ,  0.3515504E-04,  -0.7487031,    0.0,  	0.0);
  InitializeVectorPar2(parhigh15P, 1116443.0    ,  0.7842201E-05,  -0.7376342,    0.0,  	0.0);
  InitializeVectorPar2(parhigh16S, 1032940.0    ,  0.7188930E-06,  -0.7601197,    0.0,  	0.0);
  InitializeVectorPar2(parhigh17Cl,804101.7     ,  0.1521286E-06,  -0.7416714,    0.0,  	0.0);
  InitializeVectorPar2(parhigh18Ar,637702.8     ,  0.2779553E-07,  -0.7268317,    0.0,  	0.0);
  InitializeVectorPar2(parhigh19K, 519937.5     ,  0.4463855E-08,  -0.7158139,    0.0,  	0.0);
  InitializeVectorPar2(parhigh20Ca,436219.8     ,  0.5878926E-09,  -0.7102222,    0.0,  	0.0);
  InitializeVectorPar2(parhigh21Sc,368703.6     ,  0.7282178E-10,  -0.7046522,    0.0,  	0.0);
  InitializeVectorPar2(parhigh22Ti,397138.5     ,  0.1330592E-11,  -0.7196041,    0.0,  	0.0);
  InitializeVectorPar2(parhigh23V, 274427.9     ,  0.6204048E-12,  -0.6963667,    0.0,  	0.0);
  InitializeVectorPar2(parhigh24Cr,228919.3     ,  0.7995359E-13,  -0.6851181,    0.0,  	0.0);
  InitializeVectorPar2(parhigh25Mn,199919.3     ,  0.4059536E-14,  -0.6801181,    0.0,  	0.0);
  InitializeVectorPar2(parhigh26Fe,5888.956103  ,  -9427.905394 ,  4519.130051 ,  -537.9587092, 0.0);
  InitializeVectorPar2(parhigh27Co,5320.748659  ,  -8223.484184 ,  3805.325189 ,  -441.4870379, 0.0);
  InitializeVectorPar2(parhigh28Ni,4781.319646  ,  -7174.783041 ,  3219.076628 ,  -364.8043497, 0.0);
  InitializeVectorPar2(parhigh29Cu,5036.391497  ,  -7072.396111 ,  2999.984178 ,  -330.4071107, 0.0);
  InitializeVectorPar2(parhigh30Zn,4379.092796  ,  -6046.707901 ,  2513.270363 ,  -271.5536118, 0.0);
  InitializeVectorPar2(parhigh31Ga,4517.118137  ,  -5908.263246 ,  2344.263972 ,  -247.0690309, 0.0);
  InitializeVectorPar2(parhigh32Ge,3842.848695  ,  -4981.810605 ,  1950.19123  ,  -202.2213944, 0.0);
  InitializeVectorPar2(parhigh33As,3936.923647  ,  -4860.36463  ,  1825.209324 ,  -184.9650549, 0.0);
  InitializeVectorPar2(parhigh34Se,3405.652981  ,  -4159.784027 ,  1540.917232 ,  -153.778945 , 0.0);
  InitializeVectorPar2(parhigh35Br,3482.563108  ,  -4065.531877 ,  1450.303751 ,  -141.6935145, 0.0);
  InitializeVectorPar2(parhigh36Kr,27792.4379   ,  2.49277E-26  ,  -0.500883086,  0.0,          0.0);
  InitializeVectorPar2(parhigh37Rb,15982.48534  ,  3.73258E-25  ,  -0.416741096,  0.0,          0.0);
  InitializeVectorPar2(parhigh38Sr,17889.43232  ,  1.13451E-27  ,  -0.454847146,  0.0,          0.0);
  InitializeVectorPar2(parhigh39Y, 12451.29445  ,  1.53532E-27  ,  -0.404893344,  0.0,          0.0);
  InitializeVectorPar2(parhigh40Zr,10206.19067  ,  1.77653E-28  ,  -0.383987982,  0.0,          0.0);
  InitializeVectorPar2(parhigh41Nb,11017.03479  ,  8.86594E-31  ,  -0.414463124,  0.0,          0.0);
  InitializeVectorPar2(parhigh42Mo,8535.867912  ,  3.16226E-31  ,  -0.383744046,  0.0,          0.0);
  InitializeVectorPar2(parhigh43Tc,9820.134131  ,  3.36129E-34  ,  -0.423630175,  0.0,          0.0);
  InitializeVectorPar2(parhigh44Ru,7669.947859  ,  1.20762E-34  ,  -0.39423283 ,  0.0,          0.0);
  InitializeVectorPar2(parhigh45Rh,8187.067923  ,  3.36422E-37  ,  -0.420048481,  0.0,          0.0);
  InitializeVectorPar2(parhigh46Pd,6449.690254  ,  1.07461E-37  ,  -0.388615201,  0.0,          0.0);
  InitializeVectorPar2(parhigh47Ag,7131.825038  ,  1.13958E-40  ,  -0.419291952,  0.0,          0.0);
  InitializeVectorPar2(parhigh48Cd,5622.898379  ,  3.91721E-41  ,  -0.39138804 ,  0.0,          0.0);
  InitializeVectorPar2(parhigh49In,6154.018663  ,  3.57242E-44  ,  -0.419285141,  0.0,          0.0);
  InitializeVectorPar2(parhigh50Sn,5308.904577  ,  1.83756E-45  ,  -0.405840662,  0.0,          0.0);
  InitializeVectorPar2(parhigh51Sb,4205.562408  ,  7.25156E-46  ,  -0.378710249,  0.0,          0.0);
  InitializeVectorPar2(parhigh52Te,4590.865265  ,  4.9712E-49   ,  -0.405349183,  0.0,          0.0);
  InitializeVectorPar2(parhigh53I, 3930.781204  ,  4.05142E-50  ,  -0.390685876,  0.0,          0.0);
  InitializeVectorPar2(parhigh54Xe,3130.401154  ,  1.74282E-50  ,  -0.364347309,  0.0,          0.0);
  InitializeVectorPar2(parhigh55Cs,3744.705387  ,  5.90267E-55  ,  -0.404336908,  0.0,          0.0);
  InitializeVectorPar2(parhigh56Ba,677.6831684  ,  1.8374E-42   ,  -0.129113534,  0.0,          0.0);
  InitializeVectorPar2(parhigh57La,717.0531832  ,  9.51565E-45  ,  -0.15227809 ,  0.0,          0.0);
  InitializeVectorPar2(parhigh58Ce,641.4062303  ,  3.69483E-46  ,  -0.144632417,  0.0,          0.0);
  InitializeVectorPar2(parhigh59Pr,507.1392883  ,  3.50854E-46  ,  -0.117584038,  0.0,          0.0);
  InitializeVectorPar2(parhigh60Nd,590.0643943  ,  1.09726E-49  ,  -0.15428397 ,  0.0,          0.0);
  InitializeVectorPar2(parhigh61Pm,487.3334007  ,  2.75498E-50  ,  -0.133698368,  0.0,          0.0);
  InitializeVectorPar2(parhigh62Sm,358.3518729  ,  2.70531E-49  ,  -0.095301068,  0.0,          0.0);
  InitializeVectorPar2(parhigh63Eu,329.3455874  ,  5.55138E-51  ,  -0.091729198,  0.0,          0.0);
  InitializeVectorPar2(parhigh64Gd,413.9439156  ,  1.23792E-55  ,  -0.139710531,  0.0,          0.0);
  InitializeVectorPar2(parhigh65Tb,336.3177626  ,  1.01038E-55  ,  -0.116827533,  0.0,          0.0);
  InitializeVectorPar2(parhigh66Dy,-1424.412332 ,  1495.300993  ,  -579.4404515,  97.27901953,  -5.836831374 );
  InitializeVectorPar2(parhigh67Ho,-1517.735661 ,  1565.770524  ,  -596.9215188,  98.77020094,  -5.861740123 );
  InitializeVectorPar2(parhigh68Er,-1474.139215 ,  1509.677298  ,  -571.276182 ,  93.84192307,  -5.533405798 );
  InitializeVectorPar2(parhigh69Tm,-1468.335    ,  1486.750     ,  -556.4769   ,  90.48160   ,  -5.288727    );
  InitializeVectorPar2(parhigh70Yb,-1224.051265 ,  1244.798125  ,  -467.4078327,  76.09587249,  -4.438193762 );
  InitializeVectorPar2(parhigh71Lu,-1297.779    ,  1295.939     ,  -478.6122   ,  76.80235   ,  -4.430208    );
  InitializeVectorPar2(parhigh72Hf,-1343.224    ,  1320.527     ,  -480.7212   ,  76.16011   ,  -4.348449    );
  InitializeVectorPar2(parhigh73Ta,-1121.596    ,  1113.986     ,  -408.6746   ,  65.06258   ,  -3.719345    );
  InitializeVectorPar2(parhigh74W, -1132.138    ,  1109.449     ,  -402.0285   ,  63.30364   ,  -3.586027    );
  InitializeVectorPar2(parhigh75Re,-1099.159467 ,  1069.010051  ,  -384.5257059,  60.12188988,  -3.384193421 );
  InitializeVectorPar2(parhigh76Os,-1190.387975 ,  1140.34171   ,  -404.1551887,  62.32091335,  -3.468141645 );
  InitializeVectorPar2(parhigh77Ir,-1094.70162  ,  1044.989734  ,  -368.8959923,  56.62618718,  -3.134201442 );
  InitializeVectorPar2(parhigh78Pt,-1163.104391 ,  1092.439757  ,  -379.9958732,  57.58306885,  -3.155663589 );
  InitializeVectorPar2(parhigh79Au,-1144.500    ,  1068.235     ,  -369.2752   ,  55.62967   ,  -3.033460    );
  InitializeVectorPar2(parhigh80Hg,-1003.592    ,  937.0178     ,  -323.7748   ,  48.69880   ,  -2.646195    );
  InitializeVectorPar2(parhigh81Tl,-907.5804    ,  850.2096     ,  -294.3571   ,  44.29993   ,  -2.405155    );
  InitializeVectorPar2(parhigh82Pb,-1015.548    ,  929.8707     ,  -315.4740   ,  46.66766   ,  -2.501258    );
  InitializeVectorPar2(parhigh83Bi,-1076.950    ,  971.7385     ,  -325.2844   ,  47.55627   ,  -2.525684    );
  InitializeVectorPar2(parhigh84Po,-1021.285    ,  917.5267     ,  -305.7521   ,  44.49191   ,  -2.351658    );
  InitializeVectorPar2(parhigh85At,-967.7470    ,  866.0588     ,  -287.4283   ,  41.64792   ,  -2.191694    );
  InitializeVectorPar2(parhigh86Rn,-338.6518    ,  311.5210     ,  -105.8054   ,  15.52788   ,  -0.8066356   );
  InitializeVectorPar2(parhigh87Fr,-295.5938    ,  271.6075     ,  -92.01641   ,  13.43905   ,  -0.6914321   );
  InitializeVectorPar2(parhigh88Ra,-348.9885768 ,  317.6478799  ,  -106.8672925,  15.57672306,  -0.8090625926);
  InitializeVectorPar2(parhigh89Ac,-141.772905  ,  134.4449818  ,  -46.38678324,  6.754305544,  -0.3309070686);
  InitializeVectorPar2(parhigh90Th,-265.8311235 ,  240.8739384  ,  -80.49877194,  11.60745927,  -0.5910321851);
  InitializeVectorPar2(parhigh91Pa,-170.928272  ,  157.4297802  ,  -53.12415098,  7.648148663,  -0.379425533 );
  InitializeVectorPar2(parhigh92U, -344.2548753 ,  304.5108607  ,  -99.78215019,  14.21641181,  -0.7269520491);
}

void G4hShellCrossSectionDoubleExpData::FillParameterMapEnergy()
{
  parameterMapEnergy [6] =  &energy6C ;
  parameterMapEnergy [7] =  &energy7N ;
  parameterMapEnergy [8] =  &energy8O ;
  parameterMapEnergy [9] =  &energy9F ;
  parameterMapEnergy [10] = &energy10Ne;
  parameterMapEnergy [11] = &energy11Na;
  parameterMapEnergy [12] = &energy12Mg;
  parameterMapEnergy [13] = &energy13Al;
  parameterMapEnergy [14] = &energy14Si;
  parameterMapEnergy [15] = &energy15P ;
  parameterMapEnergy [16] = &energy16S ;
  parameterMapEnergy [17] = &energy17Cl;
  parameterMapEnergy [18] = &energy18Ar;
  parameterMapEnergy [19] = &energy19K ;
  parameterMapEnergy [20] = &energy20Ca;
  parameterMapEnergy [21] = &energy21Sc;
  parameterMapEnergy [22] = &energy22Ti;
  parameterMapEnergy [23] = &energy23V ;
  parameterMapEnergy [24] = &energy24Cr;
  parameterMapEnergy [25] = &energy25Mn;
  parameterMapEnergy [26] = &energy26Fe;
  parameterMapEnergy [27] = &energy27Co;
  parameterMapEnergy [28] = &energy28Ni;
  parameterMapEnergy [29] = &energy29Cu;
  parameterMapEnergy [30] = &energy30Zn;
  parameterMapEnergy [31] = &energy31Ga;
  parameterMapEnergy [32] = &energy32Ge;
  parameterMapEnergy [33] = &energy33As;
  parameterMapEnergy [34] = &energy34Se;
  parameterMapEnergy [35] = &energy35Br;
  parameterMapEnergy [36] = &energy36Kr;
  parameterMapEnergy [37] = &energy37Rb;
  parameterMapEnergy [38] = &energy38Sr;
  parameterMapEnergy [39] = &energy39Y ;
  parameterMapEnergy [40] = &energy40Zr;
  parameterMapEnergy [41] = &energy41Nb;
  parameterMapEnergy [42] = &energy42Mo;
  parameterMapEnergy [43] = &energy43Tc;
  parameterMapEnergy [44] = &energy44Ru;
  parameterMapEnergy [45] = &energy45Rh;
  parameterMapEnergy [46] = &energy46Pd;
  parameterMapEnergy [47] = &energy47Ag;
  parameterMapEnergy [48] = &energy48Cd;
  parameterMapEnergy [49] = &energy49In;
  parameterMapEnergy [50] = &energy50Sn;
  parameterMapEnergy [51] = &energy51Sb;
  parameterMapEnergy [52] = &energy52Te;
  parameterMapEnergy [53] = &energy53I ;
  parameterMapEnergy [54] = &energy54Xe;
  parameterMapEnergy [55] = &energy55Cs; 
  parameterMapEnergy [56] = &energy56Ba;
  parameterMapEnergy [57] = &energy57La;
  parameterMapEnergy [58] = &energy58Ce;
  parameterMapEnergy [59] = &energy59Pr;
  parameterMapEnergy [60] = &energy60Nd;
  parameterMapEnergy [61] = &energy61Pm;
  parameterMapEnergy [62] = &energy62Sm;
  parameterMapEnergy [63] = &energy63Eu;
  parameterMapEnergy [64] = &energy64Gd;
  parameterMapEnergy [65] = &energy65Tb;
  parameterMapEnergy [66] = &energy66Dy;
  parameterMapEnergy [67] = &energy67Ho;
  parameterMapEnergy [68] = &energy68Er;
  parameterMapEnergy [69] = &energy69Tm;
  parameterMapEnergy [70] = &energy70Yb;
  parameterMapEnergy [71] = &energy71Lu;
  parameterMapEnergy [72] = &energy72Hf;
  parameterMapEnergy [73] = &energy73Ta;
  parameterMapEnergy [74] = &energy74W ;
  parameterMapEnergy [75] = &energy75Re;
  parameterMapEnergy [76] = &energy76Os;
  parameterMapEnergy [77] = &energy77Ir;
  parameterMapEnergy [78] = &energy78Pt;
  parameterMapEnergy [79] = &energy79Au;
  parameterMapEnergy [80] = &energy80Hg;
  parameterMapEnergy [81] = &energy81Tl;
  parameterMapEnergy [82] = &energy82Pb;
  parameterMapEnergy [83] = &energy83Bi;
  parameterMapEnergy [84] = &energy84Po;
  parameterMapEnergy [85] = &energy85At;
  parameterMapEnergy [86] = &energy86Rn;
  parameterMapEnergy [87] = &energy87Fr;
  parameterMapEnergy [88] = &energy88Ra;
  parameterMapEnergy [89] = &energy89Ac;
  parameterMapEnergy [90] = &energy90Th;
  parameterMapEnergy [91] = &energy91Pa;
  parameterMapEnergy [92] = &energy92U ;
}

void G4hShellCrossSectionDoubleExpData::FillParameterMapPar1()
{
  parameterMapPar1 [6] =  &parlow6C ;
  parameterMapPar1 [7] =  &parlow7N ;
  parameterMapPar1 [8] =  &parlow8O ;
  parameterMapPar1 [9] =  &parlow9F ;
  parameterMapPar1 [10] = &parlow10Ne;
  parameterMapPar1 [11] = &parlow11Na;
  parameterMapPar1 [12] = &parlow12Mg;
  parameterMapPar1 [13] = &parlow13Al;
  parameterMapPar1 [14] = &parlow14Si;
  parameterMapPar1 [15] = &parlow15P ;
  parameterMapPar1 [16] = &parlow16S ;
  parameterMapPar1 [17] = &parlow17Cl;
  parameterMapPar1 [18] = &parlow18Ar;
  parameterMapPar1 [19] = &parlow19K ;
  parameterMapPar1 [20] = &parlow20Ca;
  parameterMapPar1 [21] = &parlow21Sc;
  parameterMapPar1 [22] = &parlow22Ti;
  parameterMapPar1 [23] = &parlow23V ;
  parameterMapPar1 [24] = &parlow24Cr;
  parameterMapPar1 [25] = &parlow25Mn;
  parameterMapPar1 [26] = &parlow26Fe;
  parameterMapPar1 [27] = &parlow27Co;
  parameterMapPar1 [28] = &parlow28Ni;
  parameterMapPar1 [29] = &parlow29Cu;
  parameterMapPar1 [30] = &parlow30Zn;
  parameterMapPar1 [31] = &parlow31Ga;
  parameterMapPar1 [32] = &parlow32Ge;
  parameterMapPar1 [33] = &parlow33As;
  parameterMapPar1 [34] = &parlow34Se;
  parameterMapPar1 [35] = &parlow35Br;
  parameterMapPar1 [36] = &parlow36Kr;
  parameterMapPar1 [37] = &parlow37Rb;
  parameterMapPar1 [38] = &parlow38Sr;
  parameterMapPar1 [39] = &parlow39Y ;
  parameterMapPar1 [40] = &parlow40Zr;
  parameterMapPar1 [41] = &parlow41Nb;
  parameterMapPar1 [42] = &parlow42Mo;
  parameterMapPar1 [43] = &parlow43Tc;
  parameterMapPar1 [44] = &parlow44Ru;
  parameterMapPar1 [45] = &parlow45Rh;
  parameterMapPar1 [46] = &parlow46Pd;
  parameterMapPar1 [47] = &parlow47Ag;
  parameterMapPar1 [48] = &parlow48Cd;
  parameterMapPar1 [49] = &parlow49In;
  parameterMapPar1 [50] = &parlow50Sn;
  parameterMapPar1 [51] = &parlow51Sb;
  parameterMapPar1 [52] = &parlow52Te;
  parameterMapPar1 [53] = &parlow53I ;
  parameterMapPar1 [54] = &parlow54Xe;
  parameterMapPar1 [55] = &parlow55Cs; 
  parameterMapPar1 [56] = &parlow56Ba;
  parameterMapPar1 [57] = &parlow57La;
  parameterMapPar1 [58] = &parlow58Ce;
  parameterMapPar1 [59] = &parlow59Pr;
  parameterMapPar1 [60] = &parlow60Nd;
  parameterMapPar1 [61] = &parlow61Pm;
  parameterMapPar1 [62] = &parlow62Sm;
  parameterMapPar1 [63] = &parlow63Eu;
  parameterMapPar1 [64] = &parlow64Gd;
  parameterMapPar1 [65] = &parlow65Tb;
  parameterMapPar1 [66] = &parlow66Dy;
  parameterMapPar1 [67] = &parlow67Ho;
  parameterMapPar1 [68] = &parlow68Er;
  parameterMapPar1 [69] = &parlow69Tm;
  parameterMapPar1 [70] = &parlow70Yb;
  parameterMapPar1 [71] = &parlow71Lu;
  parameterMapPar1 [72] = &parlow72Hf;
  parameterMapPar1 [73] = &parlow73Ta;
  parameterMapPar1 [74] = &parlow74W ;
  parameterMapPar1 [75] = &parlow75Re;
  parameterMapPar1 [76] = &parlow76Os;
  parameterMapPar1 [77] = &parlow77Ir;
  parameterMapPar1 [78] = &parlow78Pt;
  parameterMapPar1 [79] = &parlow79Au;
  parameterMapPar1 [80] = &parlow80Hg;
  parameterMapPar1 [81] = &parlow81Tl;
  parameterMapPar1 [82] = &parlow82Pb;
  parameterMapPar1 [83] = &parlow83Bi;
  parameterMapPar1 [84] = &parlow84Po;
  parameterMapPar1 [85] = &parlow85At;
  parameterMapPar1 [86] = &parlow86Rn;
  parameterMapPar1 [87] = &parlow87Fr;
  parameterMapPar1 [88] = &parlow88Ra;
  parameterMapPar1 [89] = &parlow89Ac;
  parameterMapPar1 [90] = &parlow90Th;
  parameterMapPar1 [91] = &parlow91Pa;
  parameterMapPar1 [92] = &parlow92U ;
}

void G4hShellCrossSectionDoubleExpData::FillParameterMapPar2()
{
  parameterMapPar2 [6] =  &parhigh6C ;
  parameterMapPar2 [7] =  &parhigh7N ;
  parameterMapPar2 [8] =  &parhigh8O ;
  parameterMapPar2 [9] =  &parhigh9F ;
  parameterMapPar2 [10] = &parhigh10Ne;
  parameterMapPar2 [11] = &parhigh11Na;
  parameterMapPar2 [12] = &parhigh12Mg;
  parameterMapPar2 [13] = &parhigh13Al;
  parameterMapPar2 [14] = &parhigh14Si;
  parameterMapPar2 [15] = &parhigh15P ;
  parameterMapPar2 [16] = &parhigh16S ;
  parameterMapPar2 [17] = &parhigh17Cl;
  parameterMapPar2 [18] = &parhigh18Ar;
  parameterMapPar2 [19] = &parhigh19K ;
  parameterMapPar2 [20] = &parhigh20Ca;
  parameterMapPar2 [21] = &parhigh21Sc;
  parameterMapPar2 [22] = &parhigh22Ti;
  parameterMapPar2 [23] = &parhigh23V ;
  parameterMapPar2 [24] = &parhigh24Cr;
  parameterMapPar2 [25] = &parhigh25Mn;
  parameterMapPar2 [26] = &parhigh26Fe;
  parameterMapPar2 [27] = &parhigh27Co;
  parameterMapPar2 [28] = &parhigh28Ni;
  parameterMapPar2 [29] = &parhigh29Cu;
  parameterMapPar2 [30] = &parhigh30Zn;
  parameterMapPar2 [31] = &parhigh31Ga;
  parameterMapPar2 [32] = &parhigh32Ge;
  parameterMapPar2 [33] = &parhigh33As;
  parameterMapPar2 [34] = &parhigh34Se;
  parameterMapPar2 [35] = &parhigh35Br;
  parameterMapPar2 [36] = &parhigh36Kr;
  parameterMapPar2 [37] = &parhigh37Rb;
  parameterMapPar2 [38] = &parhigh38Sr;
  parameterMapPar2 [39] = &parhigh39Y ;
  parameterMapPar2 [40] = &parhigh40Zr;
  parameterMapPar2 [41] = &parhigh41Nb;
  parameterMapPar2 [42] = &parhigh42Mo;
  parameterMapPar2 [43] = &parhigh43Tc;
  parameterMapPar2 [44] = &parhigh44Ru;
  parameterMapPar2 [45] = &parhigh45Rh;
  parameterMapPar2 [46] = &parhigh46Pd;
  parameterMapPar2 [47] = &parhigh47Ag;
  parameterMapPar2 [48] = &parhigh48Cd;
  parameterMapPar2 [49] = &parhigh49In;
  parameterMapPar2 [50] = &parhigh50Sn;
  parameterMapPar2 [51] = &parhigh51Sb;
  parameterMapPar2 [52] = &parhigh52Te;
  parameterMapPar2 [53] = &parhigh53I ;
  parameterMapPar2 [54] = &parhigh54Xe;
  parameterMapPar2 [55] = &parhigh55Cs; 
  parameterMapPar2 [56] = &parhigh56Ba;
  parameterMapPar2 [57] = &parhigh57La;
  parameterMapPar2 [58] = &parhigh58Ce;
  parameterMapPar2 [59] = &parhigh59Pr;
  parameterMapPar2 [60] = &parhigh60Nd;
  parameterMapPar2 [61] = &parhigh61Pm;
  parameterMapPar2 [62] = &parhigh62Sm;
  parameterMapPar2 [63] = &parhigh63Eu;
  parameterMapPar2 [64] = &parhigh64Gd;
  parameterMapPar2 [65] = &parhigh65Tb;
  parameterMapPar2 [66] = &parhigh66Dy;
  parameterMapPar2 [67] = &parhigh67Ho;
  parameterMapPar2 [68] = &parhigh68Er;
  parameterMapPar2 [69] = &parhigh69Tm;
  parameterMapPar2 [70] = &parhigh70Yb;
  parameterMapPar2 [71] = &parhigh71Lu;
  parameterMapPar2 [72] = &parhigh72Hf;
  parameterMapPar2 [73] = &parhigh73Ta;
  parameterMapPar2 [74] = &parhigh74W ;
  parameterMapPar2 [75] = &parhigh75Re;
  parameterMapPar2 [76] = &parhigh76Os;
  parameterMapPar2 [77] = &parhigh77Ir;
  parameterMapPar2 [78] = &parhigh78Pt;
  parameterMapPar2 [79] = &parhigh79Au;
  parameterMapPar2 [80] = &parhigh80Hg;
  parameterMapPar2 [81] = &parhigh81Tl;
  parameterMapPar2 [82] = &parhigh82Pb;
  parameterMapPar2 [83] = &parhigh83Bi;
  parameterMapPar2 [84] = &parhigh84Po;
  parameterMapPar2 [85] = &parhigh85At;
  parameterMapPar2 [86] = &parhigh86Rn;
  parameterMapPar2 [87] = &parhigh87Fr;
  parameterMapPar2 [88] = &parhigh88Ra;
  parameterMapPar2 [89] = &parhigh89Ac;
  parameterMapPar2 [90] = &parhigh90Th;
  parameterMapPar2 [91] = &parhigh91Pa;
  parameterMapPar2 [92] = &parhigh92U ;
}

std::vector<std::vector<G4double>*> G4hShellCrossSectionDoubleExpData::GetParam(G4int Z)
{   

  std::vector<std::vector<G4double>*> energyPar1Par2;

  //   if (energyPar1Par2.size() != 0) {
  //     for (size_t i=1; i==energyPar1Par2.size(); i++)
  //       {
  // 	energyPar1Par2.erase(energyPar1Par2.begin());
  //       }
  //   }

  energyPar1Par2.push_back(parameterMapEnergy[Z]);
  energyPar1Par2.push_back(parameterMapPar1[Z]);
  energyPar1Par2.push_back(parameterMapPar2[Z]);
  return energyPar1Par2;
}
