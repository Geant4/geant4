{
  cout << endl << "::: Welcome to ROOT" << " v " << gROOT->GetVersionInt() << endl;
  cout << "::: * Aliases defined in rootalias.C" << endl; 
  cout << ":::    - List them with q<TAB> " << endl; 
  cout << "::: * Visualisation style defined in rootlogon.C" << endl;
  cout << ":::    - Add line gROOT->SetStyle(\"clearRetro\"); to your script to use it" << endl;


  TStyle *hipStyle= new TStyle("clearRetro","HIP plots style for publications");

  // use plain black on white colors

  hipStyle->SetFrameBorderMode(0);
  hipStyle->SetCanvasBorderMode(0);
  hipStyle->SetPadBorderMode(0);
  hipStyle->SetPadBorderSize(0);
  hipStyle->SetPadColor(0);
  hipStyle->SetCanvasColor(0);
  hipStyle->SetTitleColor(0);
  hipStyle->SetStatColor(0);
  hipStyle->SetFillColor(0);

  hipStyle->SetTextSize(0.05);
  // use bold lines 
  hipStyle->SetHistLineWidth(2);
  hipStyle->SetLineWidth(2);

hipStyle->SetPadBottomMargin(0.10);
  hipStyle->SetTitleBorderSize(0);
  hipStyle->SetPadLeftMargin(0.20);


// hipStyle->SetTitleAlign(01); //fix :::
/*
  // set the paper & margin sizes
// hipStyle->SetPaperSize(20,26);
  hipStyle->SetPadTopMargin(-0.30);
//  hipStyle->SetPadRightMargin(0.10);
//  hipStyle->SetPadBottomMargin(0.10);
//  hipStyle->SetPadLeftMargin(0.10);

  // use large Times-Roman fonts
  //hipStyle->SetTextFont(132);

  //*-*  =============================================================
    //*-*   Font ID       X11                       Win32 TTF       lfItalic  lfWeight  x 10
	//*-*        1 : times-medium-i-normal      "Times New Roman"      1           4
	//*-*        2 : times-bold-r-normal        "Times New Roman"      0           7
	//*-*        3 : times-bold-i-normal        "Times New Roman"      1           7
	//*-*        4 : helvetica-medium-r-normal  "Arial"                0           4
	//*-*        5 : helvetica-medium-o-normal  "Arial"                1           4
	//*-*        6 : helvetica-bold-r-normal    "Arial"                0           7
	//*-*        7 : helvetica-bold-o-normal    "Arial"                1           7
	//*-*        8 : courier-medium-r-normal    "Courier New"          0           4
	//*-*        9 : courier-medium-o-normal    "Courier New"          1           4
	//*-*       10 : courier-bold-r-normal      "Courier New"          0           7
	//*-*       11 : courier-bold-o-normal      "Courier New"          1           7
	//*-*       12 : symbol-medium-r-normal     "Symbol"               0           6
	//*-*       13 : times-medium-r-normal      "Times New Roman"      0           4
	//*-*       14 :                            "Wingdings"            0           4
	//*-*
	const int kHipFont=102;  // Preferred fonts in orde 102 (courier typewriter),42, 52, 62, 82
  //hipStyle->SetTextFont(52);   //42 OK
//  hipStyle->SetTextFont(kHipFont);   //42 OK

//  hipStyle->SetLabelFont(kHipFont,"x");
//  hipStyle->SetLabelFont(kHipFont,"y");
 // hipStyle->SetLabelFont(kHipFont,"z");

//  hipStyle->SetTitleFont(kHipFont);
//  hipStyle->SetTitleFont(kHipFont,"X");
//  hipStyle->SetTitleFont(kHipFont,"Y");
//  hipStyle->SetTitleFont(kHipFont,"Z");

  hipStyle->SetTitleBorderSize(0);
  hipStyle->SetTitleBorderSize(0);
  hipStyle->SetTitleFillColor(0);
  //hipStyle->SetTitleAlign(2); //fix :::

  hipStyle->SetLabelSize(0.03,"x");
  hipStyle->SetTitleSize(0.03,"x");
  hipStyle->SetLabelSize(0.03,"y");
  hipStyle->SetTitleSize(0.03,"y");
  hipStyle->SetLabelSize(0.03,"z");
  hipStyle->SetTitleSize(0.03,"z");

  // markers

  //hipStyle->SetMarkerStyle(7);  // strong dot
  hipStyle->SetMarkerStyle(2);  // cross +

  //hipStyle->SetMarkerStyle(26); // open triangle up
  //hipStyle->SetMarkerStyle(22); // triangle up 

  //hipStyle->SetMarkerStyle(25); // open square
  //hipStyle->SetMarkerStyle(21); // full square 

  //hipStyle->SetMarkerStyle(24); // open circle
  //hipStyle->SetMarkerStyle(20); // full circle

  hipStyle->SetMarkerStyle(25);



  hipStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  hipStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  hipStyle->SetOptTitle(0); // if you want title add hipStyle->SetOptTitle(1); to your script
  hipStyle->SetOptStat(0);
  hipStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  hipStyle->SetPadTickX(1);
  hipStyle->SetPadTickY(1);

  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111111);
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  */
gROOT->SetStyle("clearRetro");
}

