# DRS4Analysis (see also https://github.com/StrahinjaLukic
Reading the raw file of DRS4 and analysis of waveforms of muons/pions 

It reads the waveforms from the raw file (WaveForm::ProcessFile()) and put their parameters into the ROOT TTree (like arrival time, FWHM, base line, amplitude...)

Experiment setup consists of two trigger scintillator counters (S3 and S4) in front, the anticoincidence counter (not in the acquisition) then the 30cm iron cube absorber, PM1 and PM2 connected on the same scintillator and the 30cm iron cube backing absorber

      |     pion beam 12 - 60 GeV
      |
      V   
     ___   S3
     ___   S4
   ___ ___ anticoincidence veto (large scintillator with a hole to let beam pass, if it detects something - probably shower 
______________                    from upstream)
|             |
|    iron     |
|  absorber   |
|_____________|
  PM1 __ PM2
______________
|             |
|    iron     |
|  absorber   |
|_____________|

trigger for DSR4 is S3 and S4 and not veto


