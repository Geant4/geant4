
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
histoSpectrum   = hm.load1D(90)

# set plotter to 4*2 zones
pl.zone(5,2)
# ... and plot the histograms
hplot(histoGamDet)
hplot(histoGamDetPre)
hplot(histoGamLeavSam)
hplot(histoEleLeavSam)
hplot(histoGamLST)
hplot(histoGamLSP)
hplot(histoGamBornSam)
hplot(histoEleBornSam)
hplot(histoSpectrum)

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

wait()
pl.psPrint("hEleBornSam.ps")



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
wait()
pl.psPrint("hGamBornSam.ps")

#titoli:
pl.xAxisOption ("title","energy deposit in the detector")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "green")


#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","EnergyDep")


hplot(histoGamDet)
wait()
pl.psPrint("histoGamDet.ps")


pl.xAxisOption ("title","spectrum of the incident photons")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "yellow")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","M-Flare")


hplot(histoSpectrum)
wait()
pl.psPrint("hSpectrum.ps")


pl.xAxisOption ("title","spectrum of the photons at the detector")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "blue")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","photonsDet")


hplot(histoGamDetPre)
wait()
pl.psPrint("hGamDetPre.ps")

pl.xAxisOption ("title","spectrum of the electrons leaving the sample")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "red")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","electronsLeav")



hplot(histoEleLeavSam)
wait()
pl.psPrint("hEleLeav.ps")

pl.xAxisOption ("title","spectrum of the photons leaving the sample")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.textStyle("fontsize","10.")
#colore dell'istogramma
pl.dataStyle("linecolor", "red")

#disegna il reticolo anche sopra e a destra
#pl.zoneOption("mirrorAxis","yes") 

#mette una legenda per collegare il colore all'istogramma
pl.dataOption("legend","photonsLeav")

hplot(histoGamLeavSam)
wait()
pl.psPrint("hGamLeav.ps")



pl.xAxisOption ("title","history of photons")
pl.xAxisOption ("label","E (keV)")
pl.yAxisOption("title","counts")

pl.dataStyle("linecolor", "yellow")

pl.dataOption("legend","photonsBorn")

hplot(histoGamBornSam)

pl.dataStyle("linecolor", "red")

pl.dataOption("legend","photonsLeav")

v1=vm.from1D(histoGamLeavSam)

pl.overlay(v1)

pl.dataStyle("linecolor", "blue")

pl.dataOption("legend","photonsDet")

v2=vm.from1D(histoGamDetPre)

pl.overlay(v2)

wait()
pl.psPrint("history.ps")