# Examples
Here we provide examples for running CaTS (without Opticks). The examples consist of a gdml file to define the detector configuration. 

## homogeneous Crystal Calorimeter
Here we just read out the energy deposited in a huge Crystal.

    time ./CaTS -g homogeneous_pbwo_crystalcal.gdml  -pl 'FTFP_BERT'  -m pip_IO_Calo.mac
    time ./readCalorimeterHits CaloHits_Run0.root CaloHistos.root CalorimeterVolume

The created histogram can then be looked at using ROOT

![alt text](https://github.com/hanswenzel/CaTS/blob/master/images/CaloHistos.png)

## Dual Read out Crystal Calorimeter

Here in addtion to the deposited energy in the calorimeter cell we readout the number of Cerenkov photons that have been produced.

    ./CaTS -g crystalcal_pbwo.gdml -pl 'FTFP_BERT+OPTICAL'  -m pip_IO_DR.mac 
    ./readDRCalorimeterHits  DRCaloHits_Run0.root DRCaloHhstos.root CalorimeterVolume

## multiple scattering at a thin Diamond layer

In this example we register the momentum of a particle as it undergaos multiple scatering in a thin Diamond layer. 

    ./CaTS -g DiamondTarget.gdml -pl 'FTFP_BERT' -m msc.mac
    ./readMscHits  Msc_Run0.root ms_histos.root 'volTarget'


The created histogram can then be looked at using ROOT

![alt text](https://github.com/hanswenzel/CaTS/blob/master/images/mscHistos.png)
