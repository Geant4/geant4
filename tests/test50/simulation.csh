$G4WORKDIR/bin/Linux-g++/test50 < gamma_Pb.in> piombo_low.out
mv Test50_output.txt piombo_low.txt
$G4WORKDIR/bin/Linux-g++/test50 < gamma_Al.in >alluminio_low.out
mv Test50_output.txt allumionio_low.txt
$G4WORKDIR/bin/Linux-g++/test50 < gamma_Au.in > oro_low.out
mv Test50_output.txt oro_low.txt
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


