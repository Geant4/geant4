#example: 
#setenv INPUT ~/g4/geant4alpha/prototype/particle+matter/processes/hadronic/test/AIX-AFS/KaonMinusAbsorptionAtRest.out
#awk -fG4KaonMinusA_pionmomentum.awk $INPUT > piminus.vec
#
begin{i=0; pim-0; pip=0; piz=0; lam=0; sigm=0; sigp=0; sigz=0}
{
  if($2="Name") i++
  if($1!="!") {
   if($1=="pi-")  pim++
   if($1=="pi+")  pip++
   if($1=="pi0")  piz++
   if($1=="lambda") lam++
   if($1=="sigma-") sigm++
   if($1=="sigma+") sigp++
   if($1=="sigma0") sigz++
   
   px=$2
   py=$3
   pz=$4
   e=$5
   if($1=="pi-")  printf("%6.1f\n", sqrt(px*px+py*py+pz*pz))
  } 
}

END {
#printf("  pi-  pi+  pi0  lam sigm sigp sigz \n")
#printf("----------------------------------\n")
#printf("%5d%5d%5d%5d%5d%5d%5d\n",pim,pip,piz,lam,sigm,sigp,sigz)
#printf("%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f%5.2f\n",pim/i,pip/i,piz/i,lam/i,sigm/i,sigp/i,sigz/i)
#for(x=1;x<=j;x++)
#printf("%6.1f    ",ks3[x])
#printf("\n")
#for(x=1;x<=j;x++)
#printf("%6.1f    ",ks4[x])
#printf("\n")
}

