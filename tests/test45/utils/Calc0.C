{
Double_t yy[300];
Double_t x11[300];
Double_t y11[300];
Double_t errx11[300];
Double_t erry11[300];
Double_t erry[3000] = {.0};
Double_t x22[300];
Double_t y22[300];
Double_t errx22[300];
Double_t erry22[300];
Double_t x33[300];
Double_t y33[300];
Double_t errx33[300];
Double_t erry33[300];
Double_t x44[300];
Double_t y44[300];
Double_t errx44[300];
Double_t erry44[300];
Double_t x55[300];
Double_t y55[300];
Double_t errx55[300];
Double_t erry55[300];
Double_t x66[300];
Double_t y66[300];
Double_t errx66[300];
Double_t erry66[300];
string linetype;
Double_t angle;
Int_t l;

TH1F *h;
Int_t j1 = 11; 
Int_t j2 = 41; 
Double_t z1, z2, z3, z4;

string file45 = "test45.txt";                
ifstream  data_file(file45.c_str(),  ios::in );
if (data_file.fail()) cout << "Error opening data file " << file45 << endl;
else {cout << "data file " << file45 << " is opened" << endl;}

for (i = 0; i < npl; i++) {
  data_file >> linetype >> bin >> angle;
  if(linetype == "END") {cout << "Buy " << endl; break;}
 // cout << bin << " " << angle << endl;

  if (i==0){
  for (j=0; j<bin; j++) {
    data_file >> x11[j] >> errx11[j] >> y11[j] >> erry11[j];
   // cout  << x11[j] << " " << y11[j] << endl;
  }
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }

  if (i==1){
  for (j=0; j<bin; j++) {
    data_file >> x22[j] >> errx22[j] >> y22[j] >> erry22[j];
   // cout  << x22[j] << " " << y22[j] << endl;
 }
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }

  if (i==2){
  for (j=0; j<bin; j++) {
    data_file >> x33[j] >> errx33[j] >> y33[j] >> erry33[j];
   // cout  << x33[j] << " " << y33[j] << endl;
}
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }

  if (i==3){
  for (j=0; j<bin; j++) {
    data_file >> x44[j] >> errx44[j] >> y44[j] >> erry44[j];
   // cout  << x44[j] << " " << y44[j] << endl;
 }
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }

  if (i==4){
  for (j=0; j<bin; j++) {
    data_file >> x55[j] >> errx55[j] >> y55[j] >> erry55[j];
   // cout  << x55[j] << " " << y55[j] << endl;
 }
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }

  if (i==5){
  for (j=0; j<bin; j++) {
    data_file >> x66[j] >> errx66[j] >> y66[j] >> erry66[j];
   // cout  << x66[j] << " " << y66[j] << endl;
 }
 // cout << "***" << endl;
 // cout << i << " " << bin << endl;
  }
 // cout << "Buy... " << endl;
}

}

