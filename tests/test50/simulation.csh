$G4WORKDIR/bin/Linux-g++/test50 < gamma_Pb.in> piombo_low.out
mv Test50_output.txt piombo_low.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_st_Al.in > alluminio_st_SP.out
mv Test50_output.txt SP_st_Al.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_st_Au.in > oro_st_SP.out
mv Test50_output.txt SP_st_Au.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_st_Pb.in > piombo_st_SP.out
mv Test50_output.txt SP_st_Pb.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_low_Al.in > alluminio_low_SP.out
mv Test50_output.txt SP_low_Al.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_low_Au.in > oro_low_SP.out
mv Test50_output.txt SP_low_Au.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_low_Pb.in > piombo_low_SP.out
mv Test50_output.txt SP_low_Pb.txt
$G4WORKDIR/bin/Linux-g++/test50 < back_Al_low.in > Al_low_back_electron.out
mv Test50_output.txt Al_back_electron_low.txt
$G4WORKDIR/bin/Linux-g++/test50 < back_Al_st.in > Al_st_back_electron.out
mv Test50_output.txt Al_back_electron_st.txt
$G4WORKDIR/bin/Linux-g++/test50 < trans_Al_st_positron.in > Al_st_trans_e+.out
mv Test50_output.txt Al_trans_positron_st.txt
$G4WORKDIR/bin/Linux-g++/test50 < trans_Al_st.in > Al_st_trans_e-.out
mv Test50_output.txt Al_trans_electron_st.txt
$G4WORKDIR/bin/Linux-g++/test50 < SP_st_Al_p.in > alluminio_st_SP_p.out
mv Test50_output.txt SP_st_Al_proton.txt




