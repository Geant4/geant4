hm.selectStore("xrayfluo.his")
histoGamDet     = hm.load1D(1)

# prepare to fit the histo h1
# creating a vector from a histogram
v1=vm.from1D(histoGamDet)

vm.list()

# create a new fitter:
fit=Fitter()
# perform fit using the histogram		
fit.setData(v1)	
v1.toAscii("EnergyDep.dat")  	
# set the model to "gaussian"   
fit.setModel("G")	
# define (and init) parameters for fit function		
fit.setParameter("amp"  ,653.)
amp = fit.fitParameter ("amp")
amp.setStart(653.)
amp.setStep(1.)
amp.setBounds(500.,700.)
# ... for these the order is relevant
fit.setParameter("mean" ,1.72)
mean = fit.fitParameter ("mean")
mean.setStart(1.72)
mean.setStep(0.01)
mean.setBounds(1.62,1.82)

# ... as there is no "intelligent parsing" possible :-(
fit.setParameter("sigma",.1)
sigma = fit.fitParameter ("sigma")
sigma.setStart(0.1)
sigma.setStep(0.005)
sigma.setBounds(0.01,0.3)
# perform fit using the histogram
fit.chiSquareFit()	
fit.printResult()		
# get vector wiht fitted function for overlay

vfit=fit.fittedVector(vm)			

# comment this if you want to have a look at the resource file for NAG
shell("rm -f e04ucc.r")

#pl.zoneOption("mirrorAxis","yes")
pl.xAxisOption("title","energy release in the detector")
pl.yAxisOption("label","counts")
# set error bars (full version)
#pl.dataOption("Representation","Error") 
# plot data with overlayed fit-result
pl.textStyle("fontsize","10.")
pl.dataStyle("linecolor", "green")
pl.dataOption("legend","fit")
#pl.dataStyle("lineshape","dashdot")
pl.plot(vfit)

wait()
pl.psPrint("fitEnegyDep.ps")	
pl.reset()

pl.textStyle("fontsize","10.")
pl.dataStyle("linecolor", "red")
pl.dataOption("legend","energyDeposit")
pl.plot(v1)
wait()

pl.dataStyle("linecolor", "green")
pl.dataOption("legend","fit")
pl.dataStyle("lineshape","dashdot")
pl.overlay(vfit)

pl.psPrint("fitEnegyOverl.ps")		
pl.reset()
del fit


