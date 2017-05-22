#include "../include/WaveProcessor.h"

void WaveProcessor::ProcessFile(char* filename){
	
	int Ch_PM1(CHPM1), Ch_PM2(CHPM2), Ch_S3(CHS3), Ch_S4(CHS4);
	ifstream DAfile (filename, ios::in|ios::binary);
	char word[4], small_word[2];
	char filenameROOT[256];
	int ChCounter(1);
	int i, j, cnttmp(0);
	date dateStmp;
	SSHORT aux_SSHORT;
	USHORT trigCell;
	string string1;
	float t3,t4;
	Units_in_PeakPM1 = new UnitPosAmpl[41]; // unit responses in signal with their amplitudes and time positions 
	Units_in_PeakPM2 = new UnitPosAmpl[41];	
	//int first(FIRST), last(LAST);
	int chi2PM1cnt(0);
	Float_t NUnitchi2PM1sum(0), Chi2reject(CHI2REJECT);
	Chi2reject/=10.;
	
 
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

	//string1="root/";
	//string1.append(filename);
	string1="";
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
	InitializeUnitParameters(); // form a UnitShape histogram containint the shape of the Unit response 
								// GetUnitHistogramRaw subtracted by GetBaseLineHistogramRaw
	
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
	
			DAfile.read(word, 4);
			eventStart=eventID=convertChtoUint(word); // first word after "EHDR"
			//cout<<"Event: "<<eventID<<endl;
			
			 setEventID(eventID);			
				
		// ocitano je da "EHDR" i ide se na ocitavanje dogadjaja
		while (!DAfile.eof()) {//this loop until there is no more of any event i.e end of the file

			 
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
				set_bin_time_n_voltage(j, 1024, 0, 0, 0); // this is only to set upper limit of the last bin
			} //end for j
			
			if (DEBUG) cout<<"Aligning channels ..."<<eventID<<endl;
			
			 allignCells0(trigCell); // align cell #0 of all channels

			//RemoveSpikes(1.5, 2); // (threshold, width in bins)
			
			if (DEBUG) cout<<"Create and fill histograms of an event: "<<eventID<<endl;

			// CreateTempHistograms() will set bins 0,1, 1023,1022,1021 and 1020 same as bin 3 ... , because these bin have spikes-gliches
			 CreateTempHistograms(); // <slower procedure>
			 
			 
						
						
			if(eventID==1)  FilterFFTofCurrentHist(1);
			
			//cin>>cchh;

			if (DEBUG) cout<<"give_waveform of event: "<<eventID<<endl;
			
			

			if (ChCounter>=Ch_PM1) WFParamPM1 =  give_waveform_parameters(Ch_PM1); 
			if (ChCounter>=Ch_PM2) WFParamPM2 =  give_waveform_parameters(Ch_PM2);
			if (ChCounter>=CHS3) WFParamS3 =  give_waveform_parameters(CHS3); 
			if (ChCounter>=CHS4) WFParamS4 =  give_waveform_parameters(CHS4);
			
			t3 = ArrivalTime (GetTempHist(CHS3), getTriggerHeight(), WFParamS3.baseLine, 3., 0.4); // not necessary to have a same fraction (FRACTION) as PM1 and PM2
			t4 = ArrivalTime (GetTempHist(CHS4), getTriggerHeight(), WFParamS4.baseLine, 3., 0.4);
//cout<< "WFParamPM1.FWHM="<< WFParamPM1.FWHM <<endl;
		if (eventID<20) PrintCurrentHist(1);
			
			TimeRef=(t4+t3)/2;
			
			// if FWHM < . it doesn't qualify to be a response, also saturated peaks don't qualify
			if ((WFParamPM1.FWHM>6.)&&(WFParamPM1.maxVal<490.)) {
				UnitResponseFinder(GetTempHist(1), TimeRef, &NUnitsPM1, Units_in_PeakPM1);// Units_in_Peak holds pas and ampl, NUnits - how many are find
				if (DEBUG5) cout<<"Calling: Units_in_PeakPM1[0].amplitude="<<Units_in_PeakPM1[0].amplitude<<endl;
				if (DEBUG5) cout<<"Calling: Units_in_PeakPM1[0].position="<<Units_in_PeakPM1[0].position<<endl;
				NUnitchi2PM1=FitUnitResponses(TempShapeCh1, UnitShape, NUnitsPM1, Units_in_PeakPM1, 
				triggerHeight, eventID-eventStart, filename, TimeRef, UnitParameters->maxVal); // there is ~5 bins per ns
				NUnitchi2PM1sum+=NUnitchi2PM1;
								//cout<<"chi2(inPF)="<<NUnitchi2PM1<<endl;
				chi2PM1cnt++; // counts how many chi2 is calculated (for average chi2)
										

			}
			
			

			/// see how to compile Units_in_Peak
			
						//NUnitchi2PM2=FitUnitResponses(TempShapeCh2, NUnitsPM2, Units_in_PeakPM2);
			// the original values of Unit_in_Peak are rewriten by the fitted values
			if (DEBUG5) cout<<"NUnitsPM1="<<NUnitsPM1<<endl;
			for (i=0; i<=NUnitsPM1; i++){
				if (DEBUG)cout<<i<<": position="<<Units_in_PeakPM1[i].position<<", amplitude="<<Units_in_PeakPM1[i].amplitude<<flush<<endl;
				if (DEBUG) cout<<UnRespDistrPM1<<endl;
				if (NUnitchi2PM1<Chi2reject) UnRespDistrPM1->Fill((Double_t)(Units_in_PeakPM1[i].position), (Double_t)(Units_in_PeakPM1[i].amplitude));
			//UnRespDistrPM1->SetBinContent(1, 0.2);
			}
			if (DEBUG) cout<<"here"<<endl;
			WFParamPM1.arrivalTimeCorrected = WFParamPM1.arrivalTime2 - TimeRef + 30.; 
			if (WFParamPM1.arrivalTime2==-1.) WFParamPM1.arrivalTimeCorrected=-1.;
			WFParamPM2.arrivalTimeCorrected = WFParamPM2.arrivalTime2 - TimeRef + 30.;
			if (WFParamPM2.arrivalTime2==-1.) WFParamPM2.arrivalTimeCorrected=-1.; 			
			// 30 ns added just to put peaks on the Total histograms, where they are usually
			
			
			if (DEBUG) cout<<"Filling the tree..."<<flush<<endl;
			paramTree -> Fill();
			if (DEBUG) cout<<"The tree filled."<<endl;
			//cout<<"ok"<<endl;
			FillTotHistograms(WFParamPM1.arrivalTimeCorrected, WFParamPM2.arrivalTimeCorrected, WFParamS3.arrivalTime2, WFParamS4.arrivalTime2);
			if (DEBUG) cout<<"Filling FillTotHistogramsNonAlligned..."<<flush<<endl;
			FillTotHistogramsNonAlligned(WFParamPM1.arrivalTime2, WFParamPM2.arrivalTime2, WFParamS3.arrivalTime2, WFParamS4.arrivalTime2);
			if (DEBUG) cout<<"Total histograms filled."<<flush<<endl;
			


			if ((WFParamPM1.FW10pcntM>30.)&&(WFParamPM1.maxVal>15.)&&(cnttmp<5))
			{ 
				cnttmp++;
				PrintCurrentHist(1);
				PrintCurrentHist(4);
				cout<<"baseLine="<<WFParamPM1.baseLine<<", baseLineRMS="<<WFParamPM1.baseLineRMS<<", maxVal="<<WFParamPM1.maxVal
				<<", FW10pcntM="<<WFParamPM1.FW10pcntM<<endl;
			}
			if (DEBUG) cout<<"check if EHDR"<<endl;
			
			DAfile.read(word, 4); // if not eof, then it must be a new event
			string1="";
			for (i=0; i<4; i++) string1+=word[i];
			if ((string1!="EHDR")&&(!DAfile.eof())){
				
				cout<<"Error, EHDR expected instead of"<< string1 <<endl;
				exit(EXIT_FAILURE);
			}
			if (DEBUG) cout<<"Deleting Temp histograms"<<endl;
			 DeleteTempHistograms();
			 
			//if(eventID>eventStart+3000) break;
			if(!DAfile.eof()){ // don't read if end of file, it will erase eventID
				DAfile.read(word, 4);
				eventID=convertChtoUint(word); // first word after "EHDR"
				setEventID(eventID);
			}

	}// while !eof, go to read new event
	

				
} // end of DAfile.is_open

if (DEBUG) cout<<"Scaling TotShapeHistograms"<<endl;

	TotShapeNA[Ch_PM1]->Scale(1/(float)(GoodArrivalCNT[Ch_PM1]));
	TotShapeNA[Ch_PM2]->Scale(1/(float)(GoodArrivalCNT[Ch_PM2]));
	TotShapeNA[Ch_S3]->Scale(1/(float)(eventID-eventStart));
	TotShapeNA[Ch_S4]->Scale(1/(float)(eventID-eventStart));	
	

for (i=1;i<=4;i++){
	cout<<"eventID-eventStart="<<eventID-eventStart<<", NullEventCNT["<<i<<"]="<<NullEventCNT[i]<<", GoodArrival["<<i<<"]="<<GoodArrivalCNT[i]<<endl;
	if (NullEventCNT[i]!=0) NullEventShape[i]->Scale(1/(float)NullEventCNT[i]);
	if (NullEventCNT[i]!=0) NullEventShapeNA[i]->Scale(1/(float)NullEventCNT[i]);
}

for (i=1;i<=4;i++){
	cout<<"Write into histograms"<<endl;
	TotShape[i]->Write();
	TotShapeNA[i]->Write();
	NullEventShape[i]->Write();
	NullEventShapeNA[i]->Write();
}

UnRespDistrPM1->Write();

cout<<"Null Event % = "<<((float)NullEventCNT[1]/2+NullEventCNT[2]/2)/(float)(eventID-eventStart)*100<<endl; // Null event for PM1 and PM2
cout<<"Chi2="<<NUnitchi2PM1sum/chi2PM1cnt<<endl;
paramTree -> Write();

cout<<"Tree writen..., closing root file..."<<flush<<endl;

f->Close();

cout<<"root file closed, closing dat file..."<<flush<<endl;
DAfile.close();	
cout<<"closed dat file..."<<flush<<endl;

delete Units_in_PeakPM1;
delete Units_in_PeakPM2;
	
}
