
# tell HistoManger where to find the histograms:
hm.selectStore("fluotesthisto.hbook")

# ... and load them into memory:
histoGamDet     = hm.load1D(10)
histoGamDetPre  = hm.load1D(20)
histoGamLeavSam = hm.load1D(30)
histoEleLeavSam = hm.load1D(40)
histoGamLST     = hm.load1D(50)
histoGamLSP     = hm.load1D(60)
histoGamBornSam = hm.load1D(70) 
histoEleBornSam = hm.load1D(80)

# set plotter to 4*2 zones
pl.zone(4,2)
# ... and plot the histograms
hplot(histoGamDet)
hplot(histoGamDetPre)
hplot(histoGamLeavSam)
hplot(histoEleLeavSam)
hplot(histoGamLST)
hplot(histoGamLSP)
hplot(histoGamBornSam)
hplot(histoEleBornSam)

pl.psPrint("histos.ps")

# prompt user for <return>
wait()

# set plotter to 1*1 zones
pl.zone(1,1)

#titoli:
pl.xAxisOption ("title","energy of the electrons generated in the sample")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "red")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","EleBorn")


hplot(histoEleBornSam)

pl.psPrint("hEleBornSam.ps")

wait()
#titoli:
pl.xAxisOption ("title","energy of the gammas generated in the sample")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "blue")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","GamBorn")


hplot(histoGamBornSam)

pl.psPrint("hGamBornSam.ps")

#titoli:
pl.xAxisOption ("title","energy deposit in the detector")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "black")


#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","EnergyDep")


hplot(histoGamDet)

pl.psPrint("histoGamDet.ps")









