// is in use in high energy cases, new Plot.C files: 800, 1200, 1500, 1600, 3000 MeV, for all elements

{
Double_t x1[3000] = {.0};
Double_t y1[3000] = {.0};
Double_t errx1[3000] = {.0};
Double_t erry1[3000] = {.0};

string file45 = "test45.txt";                
ifstream  data_file(file45.c_str(),  ios::in );
if(data_file.fail()) cout << "Error opening data file " << file45 << endl;
else                 cout << "data file " << file45 << " is opened" << endl;

TString linetype;
Double_t angle;

for (i = 0; i < npl; i++) {
  data_file >> linetype >> bin >> angle;
  if(linetype == "END") break;
  cout << bin << " " << angle << endl;
  for (j=0; j<bin; j++) {
    data_file >> x1[j] >> errx1[j] >> y1[j] >> erry1[j];
    cout  << x1[j] << " " << y1[j] << endl;
  }
  cout << "***" << endl;
  cout << i << " " << bin << endl;
  gr[i] = new TGraphErrors(bin, x1, y1, errx1, erry1);
}
}
