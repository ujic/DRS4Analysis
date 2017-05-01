// WaveProcessor.cxx
/// *******************
/// Data Analisys of Fast Hadronic Calorimeter - FastHCal
/// Created 27/4/2017
/// by Predrag Ujic
/// *****


#include "../include/WaveProcessor.h"

using namespace DRS4_data;

WaveProcessor::WaveProcessor() :
    triggerHeight(0), delay(0), eventID(0), baseLineAVG(BASELINE), baseLineCNT(0),
    aligned(false), No_of_Ch(0),
    TempShapeCh1(NULL),
    TempShapeCh2(NULL),
    TempShapeCh3(NULL),
    TempShapeCh4(NULL),
    RawTempShape(NULL)
{}

WaveProcessor::~WaveProcessor(){
/*	
	free_matrix(TimeBinWidth,1,No_of_Ch,0,1023);
	free_matrix(BinVoltage,1,No_of_Ch,0,1023);
	free_matrix(TimeBinVoltage,1,No_of_Ch,0,1023);
*/
/*
	delete TimeShapeCh1;
	if (No_of_Ch>1) delete TimeShapeCh2; 
	if (No_of_Ch>2) delete TimeShapeCh3;
	if (TimeShapeCh4) delete TimeShapeCh4; // should work if TimeShapeCh4 is defined i.e. is not a null pointer any more
*/
}

void WaveProcessor::InitializeAnalysisTree(){
	
	
	paramTree = new TTree("paramTree","Tree with parameters of the waveforms");
	//PM1	
	paramTree -> Branch("arrivalTimePM1",&WFParamPM1.arrivalTime,"arrivalTimePM1/F");
	paramTree -> Branch("EtotPM1",&WFParamPM1.Etot,"EtotPM1/F");
	paramTree -> Branch("T90PM1",&WFParamPM1.T90,"T90PM1/F");
	paramTree -> Branch("T70PM1",&WFParamPM1.T70,"T70PM1/F");
	paramTree -> Branch("T50PM1",&WFParamPM1.T50,"T50PM1/F");
	paramTree -> Branch("maxValPM1",&WFParamPM1.maxVal,"maxValPM1/F");
	paramTree -> Branch("baseLinePM1",&WFParamPM1.baseLine,"baseLinePM1/F");
	paramTree -> Branch("baseLineRMSPM1",&WFParamPM1.baseLineRMS,"baseLineRMSPM1/F");	
	paramTree -> Branch("Eof10nsPM1",&WFParamPM1.Eof10ns,"Eof10nsPM1/F");
	paramTree -> Branch("FWHMPM1",&WFParamPM1.FWHM,"FWHMPM1/F");
	paramTree -> Branch("FW10pcntMPM1",&WFParamPM1.FW10pcntM,"FWH10pcntPM1/F");	
	//PM2
	paramTree -> Branch("arrivalTimePM2",&WFParamPM2.arrivalTime,"arrivalTimePM2/F");
	paramTree -> Branch("EtotPM2",&WFParamPM2.Etot,"EtotPM2/F");
	paramTree -> Branch("T90PM2",&WFParamPM2.T90,"T90PM2/F");
	paramTree -> Branch("T70PM2",&WFParamPM2.T70,"T70PM2/F");
	paramTree -> Branch("T50PM2",&WFParamPM2.T50,"T50PM2/F");
	paramTree -> Branch("maxValPM2",&WFParamPM2.maxVal,"maxValPM2/F");
	paramTree -> Branch("baseLinePM2",&WFParamPM2.baseLine,"baseLinePM2/F");
	paramTree -> Branch("baseLineRMSPM2",&WFParamPM2.baseLineRMS,"baseLineRMSPM2/F");
	paramTree -> Branch("Eof10nsPM2",&WFParamPM2.Eof10ns,"Eof10nsPM2/F");
	paramTree -> Branch("FWHMPM2",&WFParamPM2.FWHM,"FWHMPM2/F");
	paramTree -> Branch("FW10pcntMPM2",&WFParamPM2.FW10pcntM,"FWH10pcntPM2/F");	
	//PM the rest - later
	
	
}



////////////////////////// set and get of arrays of DRS4 /////////////////////////////////////////////////////////////////////

void WaveProcessor::set_time_calibration(int ch, int bin, float value){
	aligned = false; // any change of time induce that the alignment is no more valid
    TimeBinWidth[ch][bin]=value;
}
float WaveProcessor::get_time_calibration(int ch, int bin ) const {return TimeBinWidth[ch][bin];}

void WaveProcessor::set_bin_time_n_voltage(int ch, int bin, SSHORT voltage, SSHORT range, USHORT trigger_Cell){
	float aux_f=0; 
	aligned=false; // any change of time induce that the alignment is no more valid
	
	//BinVoltage[ch][bin]=(-((float)voltage)/65536. + ((float)range)/1000.); // this way the amplitude is inverted (positive in our case)
	BinVoltage[ch][bin]=(-((float)voltage)/10. + ((float)range)/1000.); // because the values writen in the DAQ file are 0.1 mV, not just USHORT to be divided by 65536
	for (int k=0; k<bin ; k++)	aux_f += TimeBinWidth[ch][((k+trigger_Cell) % 1024)];
	TimeBinVoltage[ch][bin]=aux_f;
}



/////////////////////////////////// allignCells0 //////////////////////////////////////////////////////////////////////////

void WaveProcessor::allignCells0(USHORT trigger_cell){ // align cell #0 of all channels
	int trigCell, t1, t2, dt, chn, i;
	
	trigCell = (int) trigger_cell;
	if(aligned) { 					
		cout<<"Already aligned !!! \n Exiting..."<<endl;
		exit(EXIT_FAILURE);
		}
	
	aligned=true;
	
	if(No_of_Ch==4) t1 = TimeBinVoltage[1][(1024-trigCell) % 1024]; // ch1 is a referent chanel
	else {
		cout<<"Alignement on given channel not possible \n See WaveProcessor::allignCells0 function   \n Exiting..."<<endl;
		exit(EXIT_FAILURE);
	}

	for (chn=1 ; chn<=3 ; chn++) {
		t2 = TimeBinVoltage[chn][(1024-trigCell) % 1024];
		dt = t1 - t2;
		for (i=0 ; i<1024 ; i++) TimeBinVoltage[chn][i] += dt;
	}
}

void WaveProcessor::AmplitudeTimeCorrection(char* filename){ // this function read the file as the FileProcess function does
	
	//TGraphErrors* greS3 = new TGraphErrors(); 
	//TGraphErrors* greS4 = new TGraphErrors();
	//TGraphErrors* greS3S4 = new TGraphErrors();
	int ChS3(CHS3), ChS4(CHS4);
	Float_t S3Time[20000], S3Amplitude[20000], S3AmpErr[20000], S3TimeErr[20000], S4Time[20000], S4Amplitude[20000], S4AmpErr[20000], S4TimeErr[20000];
	int eventNo, eventS3(0), eventS4(0);// not all event will be writen in TGraph
	
	ifstream DAfile (filename, ios::in|ios::binary);
	char word[4], small_word[2];
	
	UnitPosAmpl triggScint, triggScint2;

	int ChCounter(1);
	int i, j;

	SSHORT aux_SSHORT;
	USHORT trigCell;
	string string1;
      /*
    WaveformParam WFParamCh1; // for one counter of the scintillator counter
	WaveformParam WFParamCh2; // for the other one
	WaveformParam WFParamCh3; 
	WaveformParam WFParamCh4; 
	
	TH1F* FWHMCh1 = new TH1F("FWHMCh1", "FWHM of S1", 100, 4, 15);
	TH1F* FWHMCh2 = new TH1F("FWHMCh2", "FWHM of S2", 100, 4, 15);
	TH1F* FWHMCh3 = new TH1F("FWHMCh3", "FWHM of S3", 100, 4, 15);
	TH1F* FWHMCh4 = new TH1F("FWHMCh4", "FWHM of S4", 100, 4, 15);
	*/
	
	string1=filename;
	
if(DAfile.is_open())
	{    
		for (i=0; i<4; i++) DAfile.read(word, 4); // first four words are allways file header, time header, board serial No and Channel1 header
		string1="";
		for (i=0; i<4; i++) string1+=word[i];
		if (DEBUG) cout<<"Ch. header: "<<string1<<endl<<flush;
		string1="";
		while (string1!="EHDR"){ // prvi put kad naidje na EHDR, znaci da vise nema kanala za citanje kalibracije
			for (i=0; i<1024; i++){ //after that allways followed by the time calibration of the CH1
				DAfile.read(word, 4);
							
				if(DEBUG) cout<<"Calib: ChCounter="<< ChCounter<<", bin="<<i<<", value="<<convertChtoF(word)<<endl<<flush;
				
				 set_time_calibration(ChCounter, i, convertChtoF(word));
			}
			
			DAfile.read(word, 4);
			string1="";
			for (i=0; i<4; i++) string1+=word[i];
			if (string1=="C002") ChCounter=2;
			if (string1=="C003") ChCounter=3;
			if (string1=="C004") ChCounter=4;
		if (DEBUG) cout<<"Next ch. header: "<<string1<<endl<<flush;
			//cin>>cchh;
		}
		
		SetNoOfChannels(ChCounter);
		
				string1="";
		// ocitano je da "EHDR" i ide se na ocitavanje dogadjaja
		while (!DAfile.eof()) {//this loop until there is no more of any event i.e end of the file
			DAfile.read(word, 4);
			eventNo=convertChtoUint(word); // first word after "EHDR"
			//cout<<"Event: "<<eventID<<endl;
			if ((float)((int)((float)eventNo/500.))==(float)eventNo/500.) cout << "Event: "<<eventNo<<endl;

			DAfile.read(small_word,2);  
			dateStmp.year=convertChtoUSHORT(small_word);  
			DAfile.read(small_word,2);
			dateStmp.month=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.day=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.hour=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.minute=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.second=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.milisecond=convertChtoUSHORT(small_word);
			
			DAfile.read(small_word,2);
			range=convertChtoSSHORT(small_word); 
			if (DEBUG) cout<<" range = "<<((float)range)/1000.<<endl;

			DAfile.read(word, 4); // this reads the board serial number, not needed for us
			DAfile.read(small_word,2); // this is always "T#" marking that the next small_word is Number of first readout cell

			DAfile.read(small_word,2);
			trigCell=convertChtoUSHORT(small_word);
			
			for (j=1; j<=ChCounter; j++) {
				DAfile.read(word, 4);
				string1="";
				for (i=0; i<3; i++) string1+=word[i];
				if (string1!="C00") {
					cout<<"string1:"<<string1<<endl;
					cout<<"Error, C00x expected. Exiting..."<<endl;
					exit(EXIT_FAILURE);
				}
				DAfile.read(word, 4);
				scaller = convertChtoF(word);

				
				for (i=0; i<1024; i++){ 
					DAfile.read(small_word,2);
					//aux_USHORT=convertChtoUSHORT(small_word);
					aux_SSHORT = convertChtoSSHORT(small_word);

					 set_bin_time_n_voltage(j, i, aux_SSHORT, range, trigCell);
				}//end for i
			} //end for j
			
			if (DEBUG) cout<<"Aligning channels ..."<<eventID<<endl;
			
			 allignCells0(trigCell); // align cell #0 of all channels
			
			if (DEBUG) cout<<"Create and fill histograms of an event: "<<eventID<<endl;
						 CreateTempHistograms(); // <slower procedure>
			
			triggScint = give_time_amplitude(ChS3);
//			if (DEBUG2) cout<<"amplitude="<<triggScint.amplitude<<", startPos="<<triggScint.startPosition<<endl;
			if (triggScint.amplitude>=0) { 	// amplitude is writen as -1 if its smaller then 4. mV above baseLine, 
											// or if the fit is bad
				eventS3++;
				S3Amplitude[eventS3] 	= triggScint.amplitude; 
				S3Time[eventS3]			= triggScint.startPosition;
				S3AmpErr[eventS3]		= triggScint.amplitude_err;
				S3TimeErr[eventS3]		= 0.1; 
			}
			


			//write to graph
			triggScint2 = give_time_amplitude(ChS4);	
			if (triggScint2.startPosition > 40.) cout<<triggScint2.startPosition<<", event:"<<eventNo<<endl;
			
			if (triggScint2.amplitude>=0) { 	// amplitude is writen as -1 if its smaller then 4. mV above baseLine, 
											// or if the fit is bad
				eventS4++;
				S4Amplitude[eventS4] 	= triggScint2.amplitude; 
				S4Time[eventS4]			= triggScint2.startPosition;
				S4AmpErr[eventS4]		= triggScint2.amplitude_err; 
				S4TimeErr[eventS4]		= 0.1;
			}
			
		
			
			

			
			
			
			
						DAfile.read(word, 4); // if not eof, then it must be a new event
			string1="";
			for (i=0; i<4; i++) string1+=word[i];
			if ((string1!="EHDR")&&(!DAfile.eof())){
				
				cout<<"Error, EHDR expected instead of"<< string1 <<endl;
				exit(EXIT_FAILURE);
			}
			
			 DeleteTempHistograms();	
		
		}// while !eof, go to read new event
				
} // end of DAfile.is_open

TGraphErrors* greS3 = new TGraphErrors(eventS3, S3Amplitude, S3Time, S3AmpErr, S3TimeErr);
TGraphErrors* greS4 = new TGraphErrors(eventS4, S4Amplitude, S4Time, S4AmpErr, S4TimeErr);
			TF1 *fitfun = new TF1("fitfun", "[0]+[1]/x", 0., 2000.);
			fitfun->SetParameter(0, 35.);
			fitfun->SetParameter(1, 100.);
			//rcfun->SetRange(3., 125.);
			//TFitResultPtr pfit = 
			greS3->Fit(fitfun, "S EX0");
			double a0S3 = fitfun->GetParameter(0);
			double a1S3 = fitfun->GetParameter(1);
			greS4->Fit(fitfun, "S EX0");
			double a0S4 = fitfun->GetParameter(0);
			double a1S4 = fitfun->GetParameter(1);
			
			//cout<<a01<<" "<<a11<<endl;

	TCanvas *canv = new TCanvas("canv", "greS#",800,400);
	canv->Divide(2,1,0.05,0.05);
	canv->cd(1);
	greS3->Draw("ap");
	canv->cd(2);
	greS4->Draw("*ap");
	canv->SaveAs("greS.pdf");
	
						
DAfile.close();						

}

void WaveProcessor::ProcessFile(char* filename){
	
	int Ch_PM1(CHPM1), Ch_PM2(CHPM2);
	ifstream DAfile (filename, ios::in|ios::binary);
	char word[4], small_word[2];
	char filenameROOT[256];
	int ChCounter(1);
	int i, j;
	date dateStmp;
	SSHORT aux_SSHORT;
	USHORT trigCell;
	string string1;
      
    WaveformParam WFParamCh1; // for one counter of the scintillator counter
	WaveformParam WFParamCh2; // for the other one
	WaveformParam WFParamCh3; 
	WaveformParam WFParamCh4; 
	
	TH1F* FWHMCh1 = new TH1F("FWHMCh1", "FWHM of S1", 100, 4, 15);
	TH1F* FWHMCh2 = new TH1F("FWHMCh2", "FWHM of S2", 100, 4, 15);
	TH1F* FWHMCh3 = new TH1F("FWHMCh3", "FWHM of S3", 100, 4, 15);
	TH1F* FWHMCh4 = new TH1F("FWHMCh4", "FWHM of S4", 100, 4, 15);
	
	string1="root/";
	string1.append(filename);
	string1.erase(string1.end()-4, string1.end());
	string1.append("_trig=");
	string1.append(to_string(triggerHeight));
	string1.append("_delay=");
	string1.append(to_string(delay));
	string1.append("_analys.root");
	
	std::size_t lengthS = string1.copy(filenameROOT,string1.length(),0);
	filenameROOT[lengthS]='\0';
	
	cout<<"Open " << filenameROOT <<endl;
	f = new TFile(filenameROOT, "recreate");
	
	InitializeAnalysisTree();
	
if(DAfile.is_open())
	{    
		for (i=0; i<4; i++) DAfile.read(word, 4); // first four words are allways file header, time header, board serial No and Channel1 header
		string1="";
		for (i=0; i<4; i++) string1+=word[i];
		if (DEBUG) cout<<"Ch. header: "<<string1<<endl<<flush;
		string1="";
		while (string1!="EHDR"){ // prvi put kad naidje na EHDR, znaci da vise nema kanala za citanje kalibracije
			for (i=0; i<1024; i++){ //after that allways followed by the time calibration of the CH1
				DAfile.read(word, 4);
							
				if(DEBUG) cout<<"Calib: ChCounter="<< ChCounter<<", bin="<<i<<", value="<<convertChtoF(word)<<endl<<flush;
				
				 set_time_calibration(ChCounter, i, convertChtoF(word));
			}
			
			DAfile.read(word, 4);
			string1="";
			for (i=0; i<4; i++) string1+=word[i];
			if (string1=="C002") ChCounter=2;
			if (string1=="C003") ChCounter=3;
			if (string1=="C004") ChCounter=4;
		if (DEBUG) cout<<"Next ch. header: "<<string1<<endl<<flush;
			//cin>>cchh;
		}
				 SetNoOfChannels(ChCounter);
		
				string1="";
		// ocitano je da "EHDR" i ide se na ocitavanje dogadjaja
		while (!DAfile.eof()) {//this loop until there is no more of any event i.e end of the file
			DAfile.read(word, 4);
			eventID=convertChtoUint(word); // first word after "EHDR"
			//cout<<"Event: "<<eventID<<endl;
			
			 setEventID(eventID);
			 
			 if ((float)((int)((float)eventID/500.))==(float)eventID/500.) cout << "Event: "<<eventID<<endl;
			
			DAfile.read(small_word,2);  
			dateStmp.year=convertChtoUSHORT(small_word);  
			DAfile.read(small_word,2);
			dateStmp.month=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.day=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.hour=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.minute=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.second=convertChtoUSHORT(small_word);
			DAfile.read(small_word,2);
			dateStmp.milisecond=convertChtoUSHORT(small_word);
			
			DAfile.read(small_word,2);
			range=convertChtoSSHORT(small_word); 
			if (DEBUG) cout<<" range = "<<((float)range)/1000.<<endl;
			
			setDateStamp(dateStmp);

			DAfile.read(word, 4); // this reads the board serial number, not needed for us
			DAfile.read(small_word,2); // this is always "T#" marking that the next small_word is Number of first readout cell

			DAfile.read(small_word,2);
			trigCell=convertChtoUSHORT(small_word);
			
			for (j=1; j<=ChCounter; j++) {
				DAfile.read(word, 4);
				string1="";
				for (i=0; i<3; i++) string1+=word[i];
				if (string1!="C00") {
					cout<<"string1:"<<string1<<endl;
					cout<<"Error, C00x expected. Exiting..."<<endl;
					exit(EXIT_FAILURE);
				}
				DAfile.read(word, 4);
				scaller = convertChtoF(word);

				
				for (i=0; i<1024; i++){ 
					DAfile.read(small_word,2);
					//aux_USHORT=convertChtoUSHORT(small_word);
					aux_SSHORT = convertChtoSSHORT(small_word);

					 set_bin_time_n_voltage(j, i, aux_SSHORT, range, trigCell);
				}//end for i
			} //end for j
			
			if (DEBUG) cout<<"Aligning channels ..."<<eventID<<endl;
			
			 allignCells0(trigCell); // align cell #0 of all channels
			
			if (DEBUG) cout<<"Create and fill histograms of an event: "<<eventID<<endl;
						 CreateTempHistograms(); // <slower procedure>
						
						
			if(eventID==1)  FilterFFTofCurrentHist(1);
			
			//cin>>cchh;

			if (DEBUG) cout<<"give_waveform of event: "<<eventID<<endl;
			
			

			if (ChCounter>=Ch_PM1) WFParamPM1 =  give_waveform_parameters(Ch_PM1); 
			if (ChCounter>=Ch_PM2) WFParamPM2 =  give_waveform_parameters(Ch_PM2);
			if (DEBUG) cout<<"Filling the tree..."<<flush<<endl;
			paramTree -> Fill();
			if (DEBUG) cout<<"The tree filled."<<flush<<endl;			

			if (eventID<20)
			{ 
				PrintCurrentHist(3);
				PrintCurrentHist(4);
			}
			
			
			DAfile.read(word, 4); // if not eof, then it must be a new event
			string1="";
			for (i=0; i<4; i++) string1+=word[i];
			if ((string1!="EHDR")&&(!DAfile.eof())){
				
				cout<<"Error, EHDR expected instead of"<< string1 <<endl;
				exit(EXIT_FAILURE);
			}
			
			 DeleteTempHistograms();

	}// while !eof, go to read new event
				
} // end of DAfile.is_open

paramTree -> Write();
f->Close();
DAfile.close();	

	
}

//////////////////////// CreateHistograms() ///////////////////////////////////////////////////////////////////////////////////////

void WaveProcessor::CreateTempHistograms(){
	//Float_t time_aux[1024];
	Int_t pomi;
//	if (DEBUG2) cout<<"CreateHistograms"<<endl;
	if((!aligned)&&(No_of_Ch>1)) { 					
		cout<<" Chanels not aligned !!! \n Exiting..."<<endl;
		exit(EXIT_FAILURE);
		}
	for (int j=1; j<=No_of_Ch; j++){
		for (int i=0; i<1024; i++) {
			//time_aux[i] = TimeBinVoltage [j][i];
			if (DEBUG3) cout<<"Allocation time_aux["<<i<<"]="<<TimeBinVoltage [j][i]<<endl;
		}
		
		// temporary histograms for time shape, must be deleted after each event because they are rebinned each time
		switch (j) {
			case 1: TempShapeCh1 = new TH1F("TempShapeCh1", "Time Shape of channel 1", 1023, TimeBinVoltage [j]); break;// 1024 bins (1023 in definition), width in ns
			case 2: TempShapeCh2 = new TH1F("TempShapeCh2", "Time Shape of channel 2", 1023, TimeBinVoltage [j]); break;
			case 3: TempShapeCh3 = new TH1F("TempShapeCh3", "Time Shape of channel 3", 1023, TimeBinVoltage [j]); break;
			case 4: TempShapeCh4 = new TH1F("TempShapeCh4", "Time Shape of channel 4", 1023, TimeBinVoltage [j]); break;
		}
		
		if (DEBUG3) cout<<"AllHist, bin 512 at: "<<	TempShapeCh1->GetXaxis()->GetBinCenter(512)<<endl;
		
	}
	
	if (DEBUG3) cout<<"FillHistograms, No of Channels: "<<No_of_Ch<<endl;
	
	for (int j=1; j<=No_of_Ch; j++){
	//if(DEBUG) cout <<"j in FillHistograms: "<<j<<endl;
	for (int i=0; i<1024; i++){
		//if(DEBUG) cout <<"i in FillHistograms: "<<i<<endl;
		if (DEBUG3) cout<<"BinVoltage["<<j<<"]["<<i<<"]="<<BinVoltage[j][i]<<endl;
		if (DEBUG3) cout<<"TimeBinVoltage["<<j<<"]["<<i<<"]="<<TimeBinVoltage[j][i]<<endl;
		
		switch (j) {
			case 1: TempShapeCh1->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
			case 2: TempShapeCh2->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
			case 3: TempShapeCh3->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
			case 4: TempShapeCh4->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
		}
		
		if (DEBUG3) pomi = TempShapeCh1->FindBin(TimeBinVoltage[j][i]);
		if (DEBUG3) cout<<"Filling the bin No: "<< pomi <<endl;
		if (DEBUG3) cout<<"Left bin edge:"<<TempShapeCh1->GetXaxis()->GetBinLowEdge(pomi)<<", Right Edge: "<<(TempShapeCh1->GetXaxis()->GetBinLowEdge(pomi)+TempShapeCh1->GetXaxis()->GetBinWidth(pomi))<<endl;
		if (DEBUG3) cout<<"TempShapeCh1(TimeBinVoltage)="<<TempShapeCh1->GetBinContent(TempShapeCh1->FindBin(TimeBinVoltage[j][i]))<<flush<<endl<<endl;
		
		}// end loop i	
	}// end loop j
	
	
	
}

///////////////////////////////////// DeleteHistograms /////////////////////////////////////////////////////////////////

void WaveProcessor::DeleteTempHistograms(){
	TempShapeCh1->Delete();
	if (TempShapeCh2) TempShapeCh2->Delete(); 
	if (No_of_Ch>2) delete TempShapeCh3;
	if (No_of_Ch>3) delete TempShapeCh4; // should work if TimeShapeCh4 is defined i.e. is not a null pointer any more

}
	

/////////////////////////////////////  GetHistogram ///////////////////////////////////////////////////////////////////////////////////

TH1F* WaveProcessor::GetTempHist(int Ch) const {
	if ((Ch=1)&&(No_of_Ch>0)) return TempShapeCh1;
	if ((Ch=2)&&(No_of_Ch>1)) return TempShapeCh2;
	if ((Ch=3)&&(No_of_Ch>2)) return TempShapeCh3;
	if ((Ch=4)&&(No_of_Ch>3)) return TempShapeCh4;
	
	cout<<"Channel "<<Ch<<" is not defined GetTempHist. Exiting..."<<endl;
	exit(EXIT_FAILURE);
}

void WaveProcessor::PrintCurrentHist(int ch) const {
	TH1F* tmp;
	char histfilename[60];
	switch (ch) {
		case 1: tmp = TempShapeCh1; break;
		case 2: tmp = TempShapeCh2; break;
		case 3: tmp = TempShapeCh3; break;
		case 4: tmp = TempShapeCh4; break;
	}
	if (DEBUG) cout<<" PRINTCURRENTHIST"<<endl;
	

	if (DEBUG) cout<<"tmp(10ns)inPrint="<<tmp->GetBinContent(tmp->FindBin(10.))<<endl;
	
	sprintf(histfilename, "TempShapeCh%i_event%i_d%im%iy%i_%i:%i:%i::%i.pdf", ch,  eventID, dateStmp.day, dateStmp.month, dateStmp.year, dateStmp.hour, dateStmp.minute, dateStmp.second, dateStmp.milisecond);
	TCanvas *canvTempShape = new TCanvas("TempShape", histfilename,1);
	if (DEBUG) cout<<"Canvas defined"<<endl<<flush;
	

	tmp->Draw("hist");
	

	if (DEBUG) cout<<"Histogram Draw"<<endl<<flush;
	canvTempShape->SaveAs(histfilename);

}

TH1* WaveProcessor::FilterFFTofCurrentHist(int ch){ // to finish it if needed
	TH1F* tmp;
	char histfilename[60];
	int n(1024), i ;
	switch (ch) {
		case 1: tmp = TempShapeCh1; break;
		case 2: tmp = TempShapeCh2; break;
		case 3: tmp = TempShapeCh3; break;
		case 4: tmp = TempShapeCh4; break;
	}
	sprintf(histfilename, "FFTCh%i_event%i_d%im%iy%i_%i:%i:%i::%i.root", ch,  eventID, dateStmp.day, dateStmp.month, dateStmp.year, dateStmp.hour, dateStmp.minute, dateStmp.second, dateStmp.milisecond);
	TCanvas *canvFFT = new TCanvas("FFT", histfilename,1);
	
	TH1 *histMagnitude =0;
	TVirtualFFT::SetTransform(0);
	histMagnitude = tmp->FFT(histMagnitude, "MAG");

   histMagnitude->SetTitle("Magnitude of the 1st transform");
   histMagnitude->Draw();
	canvFFT->SaveAs(histfilename);
	
	TH1 *histPhase = 0;
   histPhase = tmp->FFT(histPhase, "PH");
	
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	   
	Double_t *re_full = new Double_t[n];
   Double_t *im_full = new Double_t[n];
   fft->GetPointsComplex(re_full,im_full);
   
   // Frequency filter (the noise should be checked on the histograms and then apply the filter)
   
   for (i=90;i<130;i++){ // there was a bump at this position, apparently coming from the noise of 10 ns period seen on the waveform
	   re_full[i]=25;
	   im_full[i]=25;
	
   }
	

 //inverse transform:
   TVirtualFFT *fft_inverse = TVirtualFFT::FFT(1, &n, "C2R M K");
   fft_inverse->SetPointsComplex(re_full,im_full);
   fft_inverse->Transform();
   TH1 *hFiltered = 0;

   hFiltered = TH1::TransformHisto(fft_inverse,hFiltered,"Re");
//   hFiltered->SetTitle("inverse transform filtered");
//   hFiltered->Draw();
   
//	canvFFT->SaveAs("inverse.root");

return hFiltered;
	
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// ANALYSIS //////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Float_t WaveProcessor::GetFWHM(int Ch, Float_t BL) {return WaveProcessor::GetFWHM(Ch, BL, 50.);} // FWHM is at 50 % 

Float_t WaveProcessor::GetFWHM(int Ch, Float_t BL, Float_t height) { // height [%] is the height at which the width is measured FWHM, height = 50 

	TH1F* AnalysisHist;
	int leftFWHM(0), rightFWHM(0), maxBin(0); 
	int i(0);
	Float_t FWHM;
	
	switch (Ch) {
		case 1: AnalysisHist = TempShapeCh1; break;
		case 2: AnalysisHist = TempShapeCh2; break;
		case 3: AnalysisHist = TempShapeCh3; break;
		case 4: AnalysisHist = TempShapeCh4; break;
	}

	
	maxBin = AnalysisHist->GetMaximumBin();
	
	//cout<<"maxVal="<<AnalysisHist->GetBinContent(maxBin)<<", BL="<<BL<<endl;
	
	height = ((AnalysisHist->GetBinContent(maxBin))-BL)/(100./height); // height in [%]
	

	
	while ((leftFWHM==0)||(rightFWHM==0)){
		i++;
		if ((maxBin-i<0)||(maxBin+i>1024)) break; // couldn't find FWHM, will return 0 
		if (leftFWHM==0) leftFWHM = ((AnalysisHist->GetBinContent(maxBin-i))<height)?(maxBin-i):leftFWHM ;
		if (rightFWHM==0) rightFWHM= ((AnalysisHist->GetBinContent(maxBin+i))<height)?(maxBin+i):rightFWHM ;
	}	
	FWHM = AnalysisHist->GetBinCenter(rightFWHM) - AnalysisHist->GetBinCenter(leftFWHM);
	
	return FWHM;

}

UnitPosAmpl WaveProcessor::give_time_amplitude(int Ch) {
	TH1F* AnalysisHist;
	UnitPosAmpl output;
	Float_t baseLine;
	Float_t time_of_maximum;
		
		switch (Ch) {
		case 1: AnalysisHist = TempShapeCh1; break;
		case 2: AnalysisHist = TempShapeCh2; break;
		case 3: AnalysisHist = TempShapeCh3; break;
		case 4: AnalysisHist = TempShapeCh4; break;
	}
	//if (DEBUG2)	cout<<"Ch:"<<Ch<<endl;
/// check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
/// 35. may be wrong for S3 and S4	-  IT should be 25. ns
	baseLine = AnalysisHist->Integral(0, AnalysisHist->FindBin(25.-delay), "width") / (25.-delay) ;
	//if (DEBUG2)	cout<<"baseLine:"<<baseLine<<endl;	
	
	if ((AnalysisHist->GetMaximum()-baseLine)<4.) {
		output.startPosition = -1.;
		output.amplitude = -1; 
		return output; // ignore event - could be false event - too risky
	}
	
	time_of_maximum = AnalysisHist->GetXaxis()->GetBinCenter(AnalysisHist->GetMaximumBin());
	
	output.startPosition = AnalysisHist->GetXaxis()->GetBinCenter(AnalysisHist->FindFirstBinAbove(triggerHeight+baseLine));
	
	//if (DEBUG2)	cout<<"starPosition:"<<output.startPosition<<", time_of_maximum="<<time_of_maximum<<endl;
	
	TF1* gauss = new TF1("gauss","gaus",(time_of_maximum - 4.), (time_of_maximum + 4.));
	
	AnalysisHist->Fit(gauss,"NQ","", (time_of_maximum - 4.), (time_of_maximum + 4.)) ; 
	
	//if (DEBUG2) cout<<"Fit done.."<<flush<<endl;
	
	/*
	TCanvas *canv = new TCanvas("canv", "greS#",1);

	AnalysisHist->Draw("");
	canv->SaveAs("AnalysisHist.pdf");
	
	exit(EXIT_SUCCESS);
	*/
	
	

	Double_t Chi2 = gauss->GetChisquare();
	if (Chi2<1.){
		output.amplitude = gauss->GetParameter(0);
		output.amplitude_err = gauss->GetParError(0);
	}
	else output.amplitude = -1.; // this will indicate that the fit was bad (here Chi2 is <1 if the fit is OK
	
	if ((output.startPosition<20.)||(output.startPosition>40.)) output.amplitude=-1.; // this also indicates that something is wrong
	
	if (DEBUG2) {
				Double_t Chi2 = gauss->GetChisquare();
		if ((output.amplitude_err>1))  cout<<": output.amplitude="<<output.amplitude<<", output.amplitudeErr="<<output.amplitude_err<<", Chi2="<<Chi2<<flush<<endl;
	}
	//Double_t chi2=func1->GetChisquare(); 
		

		
return output;

}

///////////// give_waveform_parameters /////////////////////

WaveformParam WaveProcessor::give_waveform_parameters(int Ch) {

WaveformParam output; //how many parameters to be returned

	TH1F* AnalysisHist;
	//char histfilename[60];
	int i;
	
	switch (Ch) {
		case 1: AnalysisHist = TempShapeCh1; break;
		case 2: AnalysisHist = TempShapeCh2; break;
		case 3: AnalysisHist = TempShapeCh3; break;
		case 4: AnalysisHist = TempShapeCh4; break;
	}
	
	if (DEBUG) cout<<"TempShapeCh1:"<<TempShapeCh1<<", AnalysisHist:"<<AnalysisHist<<endl;
	
	// DRS4 writes everything what happened 40 ns before the trigger, when delay=0, otherwise the delay is substracted from 40 ns
	// I took 35 ns, just to have a safe distance.
	output.baseLine = AnalysisHist->Integral(0, AnalysisHist->FindBin(35.-delay), "width") / (35.-delay) ; 
	output.baseLineRMS = CalcHistRMS(AnalysisHist, 1, AnalysisHist->FindBin(35.-delay));

	int ArrivalTimeBin = AnalysisHist->FindFirstBinAbove(triggerHeight+output.baseLine); // bin should be transformed to ns according to axis
	
	if (DEBUG) cout<<"ArrivalTimeBin="<<ArrivalTimeBin<<endl;
	
	if (ArrivalTimeBin==-1) { 
		//cout<<"In give_waveform_parameters, couldn't find the signal. Empty histogram or the threshold is too high! returning 0..."<<flush<<endl; 
			cout<<"Null Event ID: "<<eventID<<endl;
			output.arrivalTime = 0.; // indicates that something is wrong
		//return output; 
		} //exit(1);}
	output.arrivalTime = AnalysisHist->GetXaxis()->GetBinCenter(ArrivalTimeBin);
	output.FWHM = GetFWHM(Ch, output.baseLine);
	output.FW10pcntM = GetFWHM(Ch, output.baseLine, 10);
	
	// Get Amplitude from the fit


	output.Etot = AnalysisHist->Integral(ArrivalTimeBin, AnalysisHist->GetXaxis()->GetNbins(), "width"); 
	
	if (DEBUG) cout<<"Etot="<<output.Etot<<endl;
	

	/*
	 * RMS of the baseLine will be used as a rejection condition
	output->hist->GetXaxis()->SetRange(0, output->hist->FindBin(35.-delay));
	output->Value(baseLineRMS) = output->hist->GetRMS();
	output->hist->GetXaxis()->UnZoom();
	*/
	
/*	
	if (output.baseLine != -100) {
		if (baseLineCNT==0) baseLineAVG=0; // baseLineAVG is set to BASELINE (4 mV), just in case if the first event is baseLine = -100
											// if the first event is OK, then no need to use this arbitrary value
		baseLineAVG = (baseLineAVG*((float)baseLineCNT) + output.baseLine)/((float)(baseLineCNT+1));
		baseLineCNT++;
	}
	else output.baseLine=baseLineAVG;
	*/
	
	output.Etot -= output.baseLine*(AnalysisHist->GetXaxis()->GetBinUpEdge(1023) - output.arrivalTime);//(AnalysisHist->GetBinCenter(1023)-(35.-delay));
	
	output.maxVal  = AnalysisHist->GetMaximum() - output.baseLine;	

	TH1F *tmpHist = (TH1F*) AnalysisHist->Clone(); // baseLine should be subtracted in order to get time of 90% energy deposition of real signal (baseLine excluded)
	for (i=0; i<1024; i++){
		tmpHist->SetBinContent(i, (tmpHist->GetBinContent(i)-output.baseLine)); 
	}
	TH1* hcumul = tmpHist->GetCumulative();
	hcumul->Scale(1/tmpHist->Integral()); // normalize cumulative histogram to 1
	output.T90 = AnalysisHist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.9))-output.arrivalTime; 
	output.T70 = AnalysisHist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.7))-output.arrivalTime; 
	output.T50 = AnalysisHist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.5))-output.arrivalTime;	
	
	output.Eof10ns = AnalysisHist->Integral(ArrivalTimeBin, AnalysisHist->FindBin(output.arrivalTime + 10.), "width") 
	- output.baseLine*10. ; // Integral of first 10 ns of signal
	
		if (DEBUG) cout<<"tmpHist->Delete()"<<flush<<endl;

	tmpHist->Delete();

	return output;
}

Observables *WaveProcessor::ProcessOnline( Float_t* RawTimeArr,
    Float_t* RawVoltArr,
    Int_t RawArrLength)
{
  return ProcessOnline( RawTimeArr, RawVoltArr, RawArrLength, triggerHeight, delay);
}

Observables* WaveProcessor::ProcessOnline( Float_t* RawTimeArr,
                                                  Float_t* RawVoltArr,
                                                  Int_t RawArrLength,
                                                  float threshold,
                                                  float trigDelay)
{
  Observables *output = new Observables;
	int i;
	static float historyLen = 185;

	output->hist = new TH1F("RawTempShape", "RawTempShape", RawArrLength - 1, RawTimeArr); // nonequidistand histogram
	
	//for(i=0; i<RawArrLength; i++) RawTempShape -> SetBinContent((i+RawTrigCell)%RawArrLength, RawVoltArr[i]);
	for (i=0; i<RawArrLength; i++) output->hist -> Fill(RawTimeArr[i], -RawVoltArr[i]);
	
  int ArrivalTimeBin = output->hist->FindFirstBinAbove(threshold); // bin should be transformed to ns according to axis
	if (ArrivalTimeBin==-1) return output; // I guess will be 0 and the event ignored
	
	// baseLine must be calculated first.
	output->Value(baseLine) = output->hist->Integral(0, output->hist->FindBin(historyLen-trigDelay), "width") / (historyLen-trigDelay) ;

	output->Value(baseLineRMS) = CalcHistRMS(output->hist, 1, output->hist->FindBin(historyLen-trigDelay));

	// baseline subtraction
	TH1F *tmpHist = (TH1F*) output->hist->Clone(); // baseLine should be subtracted in order to get time of 90% energy deposition of real signal (baseLine excluded)
	for (i=0; i<RawArrLength; i++) tmpHist->SetBinContent(i, (tmpHist->GetBinContent(i)-output->Value(baseLine)));

	TH1* hcumul = tmpHist->GetCumulative();
	hcumul->Scale(1/tmpHist->Integral()); // normalize cumulative histogram to 1
	
	output->Value(arrivalTime) = output->hist->GetXaxis()->GetBinCenter(ArrivalTimeBin);
	output->Value(eTot) = output->hist->Integral(ArrivalTimeBin, output->hist->GetXaxis()->GetNbins(), "width")
						- output->Value(baseLine)*(output->hist->GetXaxis()->GetBinUpEdge(1023)- output->Value(arrivalTime));
	output->Value(maxVal) = output->hist->GetMaximum() - output->Value(baseLine);
	output->Value(dt90) = output->hist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.9))-output->Value(arrivalTime);
	output->Value(dt70) = output->hist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.7))-output->Value(arrivalTime);
	output->Value(dt50) = output->hist->GetXaxis()->GetBinCenter(hcumul->FindFirstBinAbove(0.5))-output->Value(arrivalTime);
	
	output->Value(ePrompt) = output->hist->Integral(ArrivalTimeBin, output->hist->FindBin(output->Value(arrivalTime) + 10.), "width")
	- output->Value(baseLine)*10. ; // Integral of first 10 ns of signal

	delete hcumul;
	delete tmpHist;

	return output;
}

float WaveProcessor::CalcHistRMS(const TH1F *hist, int first, int last){
	TH1F *RMShist;
	Float_t width;
	
	Double_t max = hist->GetMaximum();
	Double_t min = hist->GetMaximum();
	
	RMShist = new TH1F("RMShist", "RMShist", 100, min, max);
	
	for(int i=first; i<=last; i++){
		width = hist->GetXaxis()->GetBinWidth(i);
		//widthSum+=width;
		RMShist->Fill(hist->GetBinContent(i), 1/width); // longer integral measurement carries less of information about the time variation
	}
	//RMShist->Scale(1/widthSum); // not necessary
	
	float rms = RMShist->GetRMS();
	delete RMShist;
	return rms;
	
}




