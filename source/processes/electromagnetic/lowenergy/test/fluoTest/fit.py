
# prepare to fit the histo h1
# creating a vector from a histogram
v1=vm.from1D(histoGamDet)

vm.list()

# create a new fitter:
fit=Fitter()
# perform fit using the histogram		
fit.setData(v1)	  	
# set the model to "gaussian"   
fit.setModel("G")	
# define (and init) parameters for fit function		
fit.setParameter("amp"  ,767)	
# ... for these the order is relevant
fit.setParameter("mean" ,8.04)

# ... as there is no "intelligent parsing" possible :-(
fit.setParameter("sigma",.1)	
# perform fit using the histogram
fit.chiSquareFit()	
fit.printResult()		
# get vector wiht fitted function for overlay

vfit=fit.fittedVector()			

# comment this if you want to have a look at the resource file for NAG
shell("rm -f e04ucc.r")

pl.zoneOption("mirrorAxis","yes")
#pl.xAxisOption("XAxisTitle","x-title")
#pl.yAxisOption("YAxisLabel","y-label")
#pl.dataOption("Representation","Error") # set error bars (full version)
pl.plot(v1,vfit)			# plot data with overlayed fit-result
pl.reset()
del fit


