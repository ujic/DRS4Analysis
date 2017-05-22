// DA.cpp
/// *******************
/// Data Analisys of Fast Hadronic Calorimeter - FHCal
/// Created 27/4/2017
/// by Predrag Ujic
/// *****

/// 



// 24.4.2017. USHORT changed to SSHORT for amplitudes, except trigger_cell  and date structure!!!!!!

// TO DO:

/** eventID doesn't start allways from zero, software divided files starts from int(10000) !!!!!!!! **/

/// to introduce spike rejection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// to introduce deconvolution all that in the 27.4.2017. folder where I'll form a github conection
			// - first find unit signals
			// then using their positions and estimated amplitudes make a fit
/// use baseLine RMS as a criterium for a rejection of the baseLine
// to put pdf as one canvas and all hist on that Canvas (with hist titles as EventID)
// do the fit part fit 

// maybe partial FFT of the part of histogram after the peak value, where the rippling probably starts
// change 


/// Ch4 (S4) is a refferent one, as it has the best time resolution

/// TriggerHeight kako definisati ???
/// delay to be well known

/// PM1 and PM2 starts ~40 ns, S3 and S4 starts ~30 ns, for the baseLine it calculates until 35. ns and 25. ns, respectively

///      ChPM1 and ChPM2 MUST BE DEFINED AS CHANNEL WHERE PM1 AND PM2 ARE CONNECTED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/DA.h"

using namespace DRS4_data;

int main (int argc, const char * argv[]) {
	
TH1::AddDirectory(kFALSE);
	
	float RawTimeArr[1024];
	float RawVoltArr[1024];
	DRS4_data::Observables* WFonline;
	
   char filename[256];
   char FWHMhistname[256];
   char* chTemp;
   char chTemp2[20];
   date dateStamp;
   char cchh;
   float trgIN;
   float delayIN;
     
   	
   	
	if (argc > 3) {
		strcpy(filename, argv[1]);
		trgIN = strtof(argv[2], &chTemp);
		delayIN = strtof(argv[3], &chTemp);
	}
    else {
      printf("Usage: ./DA <filename> <TriggerHeight(in V)> <delay(in ns)> \n");
      return 0;
    }
    
    //cout<<"delayIN:"<<delayIN<<endl;
 
// reading the binary file....
	
	DataHandler->setTriggerHeight(trgIN);
	DataHandler->setDelay(delayIN);
	
	//DataHandler->AmplitudeTimeCorrection(filename); // this should make amplitude time correction of two trigger scintillator
	
	DataHandler->ProcessFile(filename);
	
	cout<<"DA.cxx"<<flush<<endl;
	
	//dateStamp = DataHandler->getDateStamp();
	
// ANALYSIS ....




exit(EXIT_SUCCESS);
                
}
    

	     
     
