#Lizard 3.6.6


hm.selectStore("brachytherapy.hbk")
hm.ls()
#spectrum of initial particles      
hEn=hm.load1D("20") 
hplot(hEn)
#if you want to print the results:pl.psPrint("Brachytherapy.ps") 
h2=hm.load2D("10")
hplot(h2)
pl.reset()
