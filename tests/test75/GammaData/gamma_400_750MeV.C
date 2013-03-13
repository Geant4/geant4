////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Cu -> pi- + X (theta_piminus = 28.4 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Cu_PiMin_28() {

  Float_t mom_668_pim_28[] = { 429., 443., 457., 472., 485.,
                               477., 492., 508., 525., 538.,
		               549., 567., 585., 605., 620.,
		               619., 639., 659. };
 
  Float_t sig_668_pim_28[] = { 0.84, 0.17, 0.80, 0.86, 1.26,
                               1.09, 1.00, 1.54, 2.05, 1.42,
		               1.65, 2.44, 1.98, 1.37, 0.76,
		               0.82, 0.36, 0.24 };
 
  Float_t emom_668_pim_28[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0 };

  Float_t esig_668_pim_28[] = { 0.63, 0.59, 0.53, 0.49, 0.63,
                                0.52, 0.49, 0.46, 0.41, 0.53,
		                0.36, 0.31, 0.25, 0.18, 0.19,
		                0.18, 0.10, 0.09 };


  TGraphErrors*  gGam668Cu_PiMin_28 = new TGraphErrors( sizeof(mom_668_pim_28)/sizeof(Float_t),
					                mom_668_pim_28,sig_668_pim_28,
							emom_668_pim_28,esig_668_pim_28 );
  gGam668Cu_PiMin_28->SetMarkerStyle(8);
  gGam668Cu_PiMin_28->SetMarkerColor(kBlue);
  return gGam668Cu_PiMin_28;

}

////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Cu -> pi+ + X (theta_piplus = 28.4 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Cu_PiPls_28() {

   Float_t mom_668_pip_28[] = { 429., 443., 457., 473., 485.,
                                476., 491., 507., 524., 538.,
				549., 566., 584., 604., 620., 
				618., 638., 658. };

   Float_t sig_668_pip_28[] = { 1.16, 1.05, 0.78, 0.97, 0.95,
                                1.40, 1.37, 1.01, 1.98, 2.02,
				1.61, 1.52, 1.47, 1.31, 0.75,
				1.32, 0.25, 0.14 };

   Float_t emom_668_pip_28[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0 };
				
   Float_t esig_668_pip_28[] = { 0.66, 0.62, 0.57, 0.52, 0.72,
                                 0.58, 0.54, 0.49, 0.45, 0.61,
				 0.31, 0.24, 0.19, 0.14, 0.13,
				 0.17, 0.07, 0.05 };

   TGraphErrors*  gGam668Cu_PiPls_28 = new TGraphErrors( sizeof(mom_668_pip_28)/sizeof(Float_t),
					                 mom_668_pip_28,sig_668_pip_28,
							 emom_668_pip_28,esig_668_pip_28 );
   gGam668Cu_PiPls_28->SetMarkerStyle(8);
   gGam668Cu_PiPls_28->SetMarkerColor(kBlue);
   return gGam668Cu_PiPls_28;

}

////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Cu -> pi- + X (theta_piminus = 44.2 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Cu_PiMin_44() {

   Float_t mom_668_pim_44[] = { 430., 443., 458., 473., 485.,
                                476., 491., 507., 524., 538.,
				549., 566., 584., 604., 619. };

   Float_t sig_668_pim_44[] = { 0.51, 1.19, 1.15, 0.79, 1.03,
                                0.81, 1.61, 0.93, 1.07, 1.03,
				0.96, 0.69, 0.48, 0.31, 0.13 };

   Float_t emom_668_pim_44[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   
   Float_t esig_668_pim_44[] = { 0.29, 0.27, 0.25, 0.22, 0.29,
                                 0.24, 0.22, 0.19, 0.15, 0.19,
				 0.15, 0.12, 0.08, 0.05, 0.06 };

   TGraphErrors*  gGam668Cu_PiMin_44 = new TGraphErrors( sizeof(mom_668_pim_44)/sizeof(Float_t),
					                 mom_668_pim_44,sig_668_pim_44,
							 emom_668_pim_44,esig_668_pim_44 );
   gGam668Cu_PiMin_44->SetMarkerStyle(8);
   gGam668Cu_PiMin_44->SetMarkerColor(kBlue);
   return gGam668Cu_PiMin_44;

}

////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Cu -> pi+ + X (theta_piplus = 44.2 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Cu_PiPls_44() {

   Float_t mom_668_pip_44[] = { 430., 444., 458., 473., 486., 
                                477., 492., 508., 525., 539.,
				549., 566., 584., 604., 620. };

   Float_t sig_668_pip_44[] = { 0.92, 1.00, 0.87, 1.13, 1.26,
                                1.15, 1.11, 0.92, 1.33, 1.19,
				1.09, 1.08, 0.79, 0.42, 0.17 };

   Float_t emom_668_pip_44[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

   Float_t esig_668_pip_44[] = { 0.30, 0.27, 0.25, 0.23, 0.30,
                                 0.24, 0.22, 0.18, 0.16, 0.20,
				 0.19, 0.15, 0.11, 0.07, 0.06 };

   TGraphErrors*  gGam668Cu_PiPls_44 = new TGraphErrors( sizeof(mom_668_pip_44)/sizeof(Float_t),
					                 mom_668_pip_44,sig_668_pip_44,
							 emom_668_pip_44,esig_668_pip_44 );
   gGam668Cu_PiPls_44->SetMarkerStyle(8);
   gGam668Cu_PiPls_44->SetMarkerColor(kBlue);
   return gGam668Cu_PiPls_44;

}

////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Pb -> pi- + X (theta_piminus = 44.2 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Pb_PiMin_44() {

   Float_t mom_668_Pb_pim_44[] = { 430., 444., 458., 473., 485.,
                                   477., 492., 508., 525., 538.,
				   550., 568., 586., 606., 621. };

   Float_t sig_668_Pb_pim_44[] = { 2.27, 1.51, 2.40, 2.38, 2.02,
                                   2.70, 2.43, 3.28, 3.07, 2.64,
				   2.60, 1.58, 1.43, 0.72, 0.20 };

   Float_t emom_668_Pb_pim_44[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

   Float_t esig_668_Pb_pim_44[] = { 0.87, 0.83, 0.77, 0.67, 0.87,
                                    0.55, 0.48, 0.44, 0.36, 0.43,
				    0.41, 0.32, 0.23, 0.14, 0.18 };

   TGraphErrors*  gGam668Pb_PiMin_44 = new TGraphErrors( sizeof(mom_668_Pb_pim_44)/sizeof(Float_t),
					                 mom_668_Pb_pim_44,sig_668_Pb_pim_44,
							 emom_668_Pb_pim_44,esig_668_Pb_pim_44 );
   gGam668Pb_PiMin_44->SetMarkerStyle(8);
   gGam668Pb_PiMin_44->SetMarkerColor(kBlue);
   return gGam668Pb_PiMin_44;

}

////////////////////////////////////////////////////////////////////////
//  668 MeV gamma + Pb -> pi+ + X (theta_piplus = 44.2 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam668Pb_PiPls_44() {

   Float_t mom_668_Pb_pip_44[] = { 429., 443., 457., 472., 485., 
                                   477., 492., 508., 525., 539.,
				   550., 567., 586., 605., 621. };

   Float_t sig_668_Pb_pip_44[] = { 2.32, 2.27, 1.55, 2.41, 1.63,
                                   1.39, 3.05, 2.82, 2.11, 2.27,
				   2.01, 1.52, 1.61, 1.26, 0.35 };

   Float_t emom_668_Pb_pip_44[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

   Float_t esig_668_Pb_pip_44[] = { 0.65, 0.59, 0.57, 0.51, 0.65, 
                                    0.55, 0.51, 0.44, 0.37, 0.48,
				    0.41, 0.36, 0.25, 0.16, 0.17 };

   TGraphErrors*  gGam668Pb_PiPls_44 = new TGraphErrors( sizeof(mom_668_Pb_pip_44)/sizeof(Float_t),
					                 mom_668_Pb_pip_44,sig_668_Pb_pip_44,
							 emom_668_Pb_pip_44,esig_668_Pb_pip_44 );
   gGam668Pb_PiPls_44->SetMarkerStyle(8);
   gGam668Pb_PiPls_44->SetMarkerColor(kBlue);
   return gGam668Pb_PiPls_44;

}


/* Not as much of interest as 668MeV data */
/*
////////////////////////////////////////////////////////////////////////
//  418 MeV gamma + Cu -> pi- + X (theta_piminus = 28.4 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam418Cu_PiMin_28() {

  Float_t mom1[] = { 250., 257., 266., 274., 281., 
                     287., 296., 306., 316., 324.,
		     336., 346., 357., 369., 379., 
		     382., 394., 407., 420., 431. };
 
  Float_t sig1[] = {-0.36, 0.89, 0.87, 0.97, 2.02,
                     0.84, 1.73, 1.50, 1.91, 1.37,
		     1.71, 1.83, 1.41, 0.77, 0.67,
		     0.69, 0.26, 0.33, 0.04, 0.07 };
 
  Float_t emom1[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  Float_t esig1[] = { 0.60, 0.59, 0.57, 0.51, 0.69, 
                      0.43, 0.43, 0.39, 0.34, 0.38, 
		      0.35, 0.28, 0.23, 0.15, 0.17,
		      0.19, 0.11, 0.11, 0.04, 0.07 };


  TGraphErrors*  gGam418Cu_PiMin_28 = new TGraphErrors( sizeof(mom1)/sizeof(Float_t),
					                mom1,sig1,emom1,esig1 );
  gGam418Cu_PiMin_28->SetMarkerStyle(8);
  gGam418Cu_PiMin_28->SetMarkerColor(kBlue);
  return gGam418Cu_PiMin_28;
}

////////////////////////////////////////////////////////////////////////
//  468 MeV gamma + Cu -> pi- + X (theta_piminus = 28.4 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam468Cu_PiMin_28() {

  Float_t mom2[] = { 288., 297., 306., 317., 325.,
                     336., 346., 357., 369., 379.,
		     382., 394., 407., 420., 429.,
		     431., 443., 457., 472., 485.};
 
  Float_t sig2[] = { 2.27, 0.72, 0.82, 1.07, 1.58,
                     1.26, 1.61, 1.37, 2.03, 2.09,
		     1.86, 1.84, 1.33, 0.89, 0.79,
		     0.47, 0.39, 0.12, 0.05, 0.00 };
 
  Float_t emom2[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  Float_t esig2[] = { 0.55, 0.52, 0.49, 0.46, 0.55,
                      0.41, 0.38, 0.32, 0.25, 0.34,
		      0.31, 0.23, 0.20, 0.13, 0.15,
		      0.13, 0.90, 0.04, 0.02, 0.02 };


  TGraphErrors*  gGam468Cu_PiMin_28 = new TGraphErrors( sizeof(mom2)/sizeof(Float_t),
					                mom2,sig2,emom2,esig2 );
  gGam468Cu_PiMin_28->SetMarkerStyle(8);
  gGam468Cu_PiMin_28->SetMarkerColor(kBlue);
  return gGam468Cu_PiMin_28;

}

////////////////////////////////////////////////////////////////////////
//  518 MeV gamma + Cu -> pi- + X (theta_piminus = 28.4 degrees)
//
//  Data from: K. Baba et al., Nucl. Phys. A322, 349 (1979) 
////////////////////////////////////////////////////////////////////////

TGraphErrors* gGam518Cu_PiMin_28() {

  Float_t mom3[] = { 336., 346., 358., 369., 379.,
                     382., 394., 407., 420., 431.,
		     429., 442., 456., 472., 484.,
		     477., 492., 508., 525., 539. };
 
  Float_t sig3[] = { 0.59, 0.70, 1.31, 1.24, 1.42, 
                     1.57, 1.56, 1.47, 2.09, 1.89,
		     1.96, 1.62, 1.26, 0.88, 0.36,
		     0.74, 0.43, 0.15, 0.00, 0.03 };
 
  Float_t emom3[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  Float_t esig3[] = { 0.38, 0.37, 0.33, 0.30, 0.42,
                      0.34, 0.30, 0.26, 0.22, 0.26,
		      0.22, 0.17, 0.11, 0.09, 0.07, 
		      0.13, 0.09, 0.04, 0.02, 0.03 };


  TGraphErrors*  gGam518Cu_PiMin_28 = new TGraphErrors( sizeof(mom3)/sizeof(Float_t),
					                mom3,sig3,emom3,esig3 );
  gGam518Cu_PiMin_28->SetMarkerStyle(8);
  gGam518Cu_PiMin_28->SetMarkerColor(kBlue);
  return gGam518Cu_PiMin_28;

}
*/
