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
{
	TotShape[0] = new TH1F("Stack", "Stack", NBINS, 0, 200);
	TotShape[1] = new TH1F("TotShapeCh1", "Time Shape of channel 1", NBINS, 0, 200); 
	TotShape[2] = new TH1F("TotShapeCh2", "Time Shape of channel 2", NBINS, 0, 200); 
	TotShape[3] = new TH1F("TotShapeCh3", "Time Shape of channel 3", NBINS, 0, 200); 
	TotShape[4] = new TH1F("TotShapeCh4", "Time Shape of channel 4", NBINS, 0, 200); 
	}

WaveProcessor::~WaveProcessor(){
// 
	delete TotShape[1]; // should work if TimeShapeCh4 is defined i.e. is not a null pointer any more
	delete TotShape[2];
	delete TotShape[3];
	delete TotShape[4];
}

void WaveProcessor::InitializeAnalysisTree(){
	
	
	paramTree = new TTree("paramTree","Tree with parameters of the waveforms");
	//PM1	
	paramTree -> Branch("arrivalTimePM1",&WFParamPM1.arrivalTime,"arrivalTimePM1/F");
	paramTree -> Branch("arrivalTimeRawPM1",&WFParamPM1.arrivalTimeRaw,"arrivalTimeRawPM1/F");
	paramTree -> Branch("arrivalTime2PM1",&WFParamPM1.arrivalTime2,"arrivalTime2PM1/F");
	paramTree -> Branch("arrivalTimeCorrectedPM1",&WFParamPM1.arrivalTimeCorrected,"arrivalTimeCorrectedPM1/F");	
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
	paramTree -> Branch("arrivalTimeRawPM2",&WFParamPM2.arrivalTimeRaw,"arrivalTimeRawPM2/F");
	paramTree -> Branch("arrivalTime2PM2",&WFParamPM2.arrivalTime2,"arrivalTime2PM2/F");
	paramTree -> Branch("arrivalTimeCorrectedPM2",&WFParamPM2.arrivalTimeCorrected,"arrivalTimeCorrectedPM2/F");
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
	//tref
	paramTree -> Branch("TimeRef", &TimeRef, "TimeRef/F");
	
	
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
	
	if(No_of_Ch==4) t1 = TimeBinVoltage[4][(1024-trigCell) % 1024]; // ch1 is a referent chanel
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



void WaveProcessor::ProcessFile(char* filename){
	
	int Ch_PM1(CHPM1), Ch_PM2(CHPM2);
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
			
			t3 = ArrivalTime (GetTempHist(CHS3), getTriggerHeght(), WFParamS3.baseLine, 3., 0.4);
			t4 = ArrivalTime (GetTempHist(CHS4), getTriggerHeght(), WFParamS4.baseLine, 3., 0.4);
			
			TimeRef=(t4+t3)/2;
			
			WFParamPM1.arrivalTimeCorrected = WFParamPM1.arrivalTime - TimeRef + 30.; 
			WFParamPM2.arrivalTimeCorrected = WFParamPM2.arrivalTime - TimeRef + 30.; 			
			// 30 ns added just to put peaks on the Total histograms, where they are usually
			
			if (DEBUG) cout<<"Filling the tree..."<<flush<<endl;
			paramTree -> Fill();
			if (DEBUG) cout<<"The tree filled."<<flush<<endl;
			FillTotHistograms(WFParamPM1.arrivalTime, WFParamPM2.arrivalTime, WFParamS3.arrivalTime, WFParamS4.arrivalTime);
	
			
		

			if ((WFParamPM1.FW10pcntM>30.)&&(WFParamPM1.maxVal>15.)&&(cnttmp<20))
			{ 
				cnttmp++;
				PrintCurrentHist(1);
				PrintCurrentHist(4);
				cout<<"baseLine="<<WFParamPM1.baseLine<<", baseLineRMS="<<WFParamPM1.baseLineRMS<<", maxVal="<<WFParamPM1.maxVal
				<<", FW10pcntM="<<WFParamPM1.FW10pcntM<<endl;
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

TotShape[1]->Write();
TotShape[2]->Write();
TotShape[3]->Write();
TotShape[4]->Write();
paramTree -> Write();

cout<<"Tree writen..., closing root file..."<<flush<<endl;

f->Close();

cout<<"root file closed, closing dat file..."<<flush<<endl;
DAfile.close();	
cout<<"closed dat file..."<<flush<<endl;
	
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
		// it has to be done since the bins' widths are different, after the "rotation" every bin will not come in the coresponding
		// bin, if their widths wouldn't be rotated as well
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
		
		// the bins 0,1, 1023, 1022, 1021, 1020 have very high spikes too often 
		BinVoltage[j][0] = BinVoltage[j][3];
		BinVoltage[j][1] = BinVoltage[j][4];
		for (int k=0; k<4; k++){
			BinVoltage[j][1023-k] = BinVoltage[j][1019-k];
		}
			
		//if(DEBUG) cout <<"j in FillHistograms: "<<j<<endl;
		for (int i=0; i<1024; i++){
			//if(DEBUG) cout <<"i in FillHistograms: "<<i<<endl;
			if (DEBUG3) cout<<"BinVoltage["<<j<<"]["<<i<<"]="<<BinVoltage[j][i]<<endl;
			if (DEBUG3) cout<<"TimeBinVoltage["<<j<<"]["<<i<<"]="<<TimeBinVoltage[j][i]<<endl;

	/*		
			switch (j) {
				case 1: TempShapeCh1->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
				case 2: TempShapeCh2->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
				case 3: TempShapeCh3->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
				case 4: TempShapeCh4->Fill(TimeBinVoltage[j][i], BinVoltage[j][i]); break;
			}
	*/
	// Filling is done by the SetbinContent, because it's faster and now it's OK since the widths are also rotated...
	// May, 3 2017
			switch (j) {
				case 1: TempShapeCh1->SetBinContent(i, BinVoltage[j][i]); break;
				case 2: TempShapeCh2->SetBinContent(i, BinVoltage[j][i]); break;
				case 3: TempShapeCh3->SetBinContent(i, BinVoltage[j][i]); break;
				case 4: TempShapeCh4->SetBinContent(i, BinVoltage[j][i]); break;
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
	// canvFFT->SaveAs(histfilename); // to check the graph
	
	TH1 *histPhase = 0;
   histPhase = tmp->FFT(histPhase, "PH");
	
	TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
	   
	Double_t *re_full = new Double_t[n];
   Double_t *im_full = new Double_t[n];
   fft->GetPointsComplex(re_full,im_full);
   

   // Frequency filter (the noise should be checked on the histograms and then apply the filter)

      // this should be re-made, but if found to be necessary   
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
	
	//if ((AnalysisHist->GetBinCenter(maxBin)>70)||(AnalysisHist->GetBinCenter(maxBin)<20)) return 200. ; // this means that something is wrong - probably a spike
	
	while ((leftFWHM==0)||(rightFWHM==0)){
		i++;
		if ((maxBin-i<0)||(maxBin+i>1024)) break; // couldn't find FWHM, will return 0 
		if (leftFWHM==0) leftFWHM = ((AnalysisHist->GetBinContent(maxBin-i))<(height+BL))?(maxBin-i):leftFWHM ;
		if (rightFWHM==0) rightFWHM= ((AnalysisHist->GetBinContent(maxBin+i))<(height+BL))?(maxBin+i):rightFWHM ;
	}	
	FWHM = AnalysisHist->GetBinCenter(rightFWHM) - AnalysisHist->GetBinCenter(leftFWHM);
	
	return FWHM;

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
	output.arrivalTimeRaw = AnalysisHist->GetXaxis()->GetBinCenter(ArrivalTimeBin);
	
	// this line is for constant fraction discrimination :
	output.arrivalTime = ArrivalTime(AnalysisHist, triggerHeight, output.baseLine, 3., 0.4); // risetime 3. ns, by eye, fraction 0.2

	output.arrivalTime2 = ArrivalTime2(AnalysisHist, output.baseLine, 0.3); // last nuber is a fraction at 

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

float WaveProcessor::MeanAndRMS(const TH1F *hist, int first, int last, float &mean, float &rms){

  if ( !hist ) return -1.;
  if ( first < 1 || first > hist->GetNbinsX() || last < 1 || last > hist->GetNbinsX() ) {
    return -1;
  }
  if ( last < first ) {
    int temp = first;
    first = last;
    last = temp;
  }

  float sum = 0, sumSq = 0;

  for(int i=first; i<=last; i++){
    float c = hist->GetBinContent(i);
    sum += c;
    sumSq += c*c;
  }

  int nBins = (last - first + 1);
  mean = sum / nBins;
  float meanSq = sumSq / nBins;
  rms = sqrt( meanSq - mean*mean );
  return mean;

}


float WaveProcessor::ArrivalTime(TH1F* hist, float threshold, float baseline,
                                  float risetime, float fraction) // risetime=3ns, saw by eye as interrupt of arrival time of low amplitude signals
{

  if ( !hist ) return -1.;

  int maxBin = hist->GetMaximumBin();
  float maxVal = hist->GetBinContent(maxBin) - baseline;

  if (maxVal < threshold) return -1.;

  float tt1 = 0;
  int nBint1 = 0;

  for (int ibin=3; ibin<=maxBin; ibin++) {

    if (hist->GetBinContent(ibin) - baseline > threshold) {
      tt1 = hist->GetBinLowEdge(ibin);
      nBint1 = ibin;
      break;
    }
  }

  // Return simple threshold crossing point (uncomment for debugging)
 // return tt1;

  if (0.9*maxVal < threshold) {
    return tt1 - risetime;
  }

  if (maxVal*fraction < threshold) {

    for (int ibin=nBint1; ibin<=maxBin; ibin++) {

      if (hist->GetBinContent(ibin) - baseline > 0.9*maxVal) {
        return hist->GetBinLowEdge(ibin) - risetime;
      }
    }

  } // maxVal < 2*threshold

  /*** Constant fraction for tall signals ***/

  float endfit = 0;
  int endFitBin = 0;
  float startfit = 0;
  int startFitBin = 0;

  // Fraction how far to integrate in both directions from the crossing point
  float fspan = 0.5;

  if (maxVal*fraction*(1.-fspan) < 2*threshold) {
    startfit = tt1;
    startFitBin = nBint1;
    for (int ibin=nBint1; ibin<=maxBin; ibin++) {
      if (hist->GetBinContent(ibin) - baseline > 2*fraction*maxVal-threshold) {
        endfit = hist->GetBinLowEdge(ibin+1);
        endFitBin = ibin;
        break;
      }
    }
  }
  else { // Threshold is below (1.-fspan) of const fraction
    for (int ibin=nBint1; ibin<=maxBin; ibin++) {
      if (hist->GetBinContent(ibin) - baseline > maxVal*fraction*(1.-fspan)) {
        startfit = hist->GetBinLowEdge(ibin);
        startFitBin = ibin;
        break;
      }
    }
    for (int ibin=nBint1; ibin<=maxBin; ibin++) {
      if (hist->GetBinContent(ibin) - baseline > maxVal*fraction*(1.+fspan)) {
        endfit = hist->GetBinLowEdge(ibin+1);
        endFitBin = ibin;
        break;
      }
    }
  } // Selection startfit and endfit

  assert (startFitBin <= endFitBin);
  if (startFitBin == endFitBin) {
    startFitBin--;
    endFitBin++;
  }

  if (startFitBin+1 == endFitBin) {
    endFitBin++;
  }

  float sumt=0, sumw=0;

  for (int ibin=startFitBin; ibin<endFitBin; ibin++) {
//    float w = exp( -fabs(hist->GetBinContent(ibin)/maxVal-fraction)/fspan );
    // cusp weighting
    float w = 1. - sqrt(fabs(hist->GetBinContent(ibin)/maxVal-fraction)/fspan);
    sumt += hist->GetBinLowEdge(ibin) * w;
    sumw += w;
  }

 //  t0 = v0/maxVal; // Replacing time by actual fraction (for debugging)
  float t0 = sumt / sumw;

  return t0 ;
}


void WaveProcessor::RemoveSpikes(float threshold, short spikeWidth)
  {
	  int kNumberOfBins(NBINS);
     int spikePos[kNumberOfBins];
     memset(spikePos, 0, sizeof(spikePos));

     const unsigned nChan = 4;
     int sp[nChan][10];
     int rsp[10];
     int n_sp[nChan], n_rsp;
     int  nNeighbor, nSymmetric;


     memset(sp, 0, sizeof(sp));
     memset(n_sp, 0, sizeof(n_sp));
     memset(rsp, 0, sizeof(rsp));
     n_rsp = 0;


     /* find spikes with a high-pass filter */
     for (unsigned iChan=0 ; iChan<nChan ; iChan++) {

        for (unsigned ibin=0 ; ibin<kNumberOfBins-1 ; ibin++) {

           float diff = - ( BinVoltage[iChan][ibin] ) / 2;
           for (unsigned spikeBin=1; spikeBin<=spikeWidth; spikeBin++) {
             diff += (BinVoltage[iChan][(ibin+spikeBin) % kNumberOfBins]) / spikeWidth;
           }
           diff -= ( BinVoltage[iChan][(ibin+spikeWidth+1) % kNumberOfBins] ) / 2;

           float slope = ( BinVoltage[iChan][(ibin+spikeWidth+1) % kNumberOfBins] )
                       - ( BinVoltage[iChan][ibin] ) ;

           if (diff > threshold && diff > slope) {
             n_sp[iChan]++;
             sp[iChan][n_sp[iChan]] = ibin;
             spikePos[ibin]++;
           }
        } // Loop over bins
     } // Loop over chans

     /* go through all spikes and look for neighbors */
     for (unsigned iChan=0 ; iChan<nChan ; iChan++) {
        for (unsigned ispike =0 ; ispike<n_sp[iChan] ; ispike++) {

           /* check if there is a spike at the same position in other channels */
           nNeighbor=0;
           for (unsigned jChan=0 ; jChan<nChan ; jChan++) {
              if (iChan != jChan) {
                 for (unsigned lspike=0 ; lspike<n_sp[jChan] ; lspike++)
                    if ( sp[iChan][ispike] == sp[jChan][lspike] )
                    {
                       nNeighbor++;
                       break;
                    }
              }
           }


           /* if at least two matching spikes, treat this as a real spike */
           if (nNeighbor >= 2) {
              // Check if this spike is already registered as real
              unsigned jspike;
              for (jspike=0 ; jspike<n_rsp ; jspike++)
                 if (rsp[jspike] == sp[iChan][ispike])
                    break;
              // If not registered, register
              if (n_rsp < 100 && jspike == n_rsp)
                 rsp[n_rsp++] = sp[iChan][ispike];
           }
        }
     } // End search for neighbors

     if (n_rsp > 10) {
       std::cout << "WARNING: More than 10 spikes in event!\n";
     }


     /* Correct spikes */
     for (unsigned ispike=0 ; ispike<n_rsp ; ispike++) {
        for (unsigned iChan=0 ; iChan<nChan ; iChan++) {
          /* remove single spike */
          float x = BinVoltage[iChan][rsp[ispike]];
          float y = BinVoltage[iChan][(rsp[ispike]+spikeWidth+1) % kNumberOfBins];

          double slope = static_cast<double>(y-x)/(spikeWidth+1);
          for (unsigned spikeBin=1; spikeBin<=spikeWidth; spikeBin++) {
            BinVoltage[iChan][(rsp[ispike]+spikeBin) % kNumberOfBins] = (x + spikeBin*slope);
          }
        } // Loop over iChan
     } // Loop over ispike


     /* find spikes at cell #0 and #1023*/
     for (unsigned iChan=0 ; iChan<nChan ; iChan++) {

        float diff = 0;
        for (unsigned spikeBin=0; spikeBin<spikeWidth; spikeBin++) {
          diff += (BinVoltage[iChan][spikeBin]) / spikeWidth;
        }
        diff -= (BinVoltage[iChan][spikeWidth]);

        // Correct immediately. False spikes have low impact here.
        if ( fabs(diff) > threshold) {
          for (unsigned spikeBin=0; spikeBin<spikeWidth; spikeBin++) {
             BinVoltage[iChan][spikeBin] = BinVoltage[iChan][spikeWidth];
          }
        }

        diff = 0;
        for (unsigned spikeBin=kNumberOfBins-spikeWidth; spikeBin<kNumberOfBins; spikeBin++) {
          diff += (BinVoltage[iChan][spikeBin]) / spikeWidth;
        }
        diff -= (BinVoltage[iChan][kNumberOfBins-spikeWidth-1]);

        // Correct immediately. False spikes have low impact here.
        if (fabs(diff) > threshold) {
          for (unsigned spikeBin=kNumberOfBins-spikeWidth; spikeBin<kNumberOfBins; spikeBin++) {
             BinVoltage[iChan][spikeBin] = BinVoltage[iChan][kNumberOfBins-spikeWidth-1];
          }
        }
     }

} // RemoveSpikes()

void WaveProcessor::FillTotHistograms(Float_t ArTm1, Float_t ArTm2, Float_t ArTm3, Float_t ArTm4){
Float_t BinContent;
int i, j, ArBin, Nbins(NBINS);

	for (j=1; j<=No_of_Ch; j++){
		for (i=0; i<Nbins; i++){
			switch (j) {
				case 1: ArBin = TempShapeCh1->FindBin(ArTm1); BinContent=TempShapeCh1->GetBinContent((ArBin+i)%(Nbins)); break;
				case 2: ArBin = TempShapeCh2->FindBin(ArTm2); BinContent=TempShapeCh2->GetBinContent((ArBin+i)%(Nbins)); break;
				case 3: ArBin = TempShapeCh3->FindBin(ArTm3); BinContent=TempShapeCh3->GetBinContent((ArBin+i)%(Nbins)); break;
				case 4: ArBin = TempShapeCh4->FindBin(ArTm4); BinContent=TempShapeCh4->GetBinContent((ArBin+i)%(Nbins)); break;
			}
			
			if (DEBUG) cout<<"j:"<<j<<", i:"<<i<<", ArBin="<<ArBin<<", BinContent="<<BinContent<<endl;
		TotShape[j]->SetBinContent((70+i)%Nbins, TotShape[j]->GetBinContent((70+i)%Nbins)+BinContent); 
		}
	}

}

float WaveProcessor::ArrivalTime2(TH1F* hist, float baseLine, float fraction){
	float tArrival;
	int i;
	Int_t maxBin = hist->GetMaximumBin();
	Float_t maxVal = hist->GetMaximum()-baseLine;

	Float_t thrsBin = hist->FindFirstBinAbove(triggerHeight+baseLine);
	if ((maxVal*fraction)>triggerHeight) tArrival = hist->GetXaxis()->GetBinLowEdge(hist->FindFirstBinAbove(maxVal*fraction+baseLine));
	else {
		for (i=0; i<20; i++){ // 10 bins-> 2 ns, the 10% of amplitude cannot be more than 2 ns before thrs
			if (hist->GetBinContent(thrsBin-i)<(maxVal*fraction+baseLine)) 
				{tArrival = hist->GetXaxis()->GetBinUpEdge(thrsBin-i); i=100;}
			}
		if (i!=100) { 
			// couldn't find the fraction of maxVal (probably because the base line is higher just before the peak)
			// thus extrapolation to fraction*maxVal ..
			tArrival = hist->GetXaxis()->GetBinLowEdge((int)((maxBin-thrsBin)/(maxVal-triggerHeight)*(0.1*maxVal)+0.5)); 
		}
			
	}
return tArrival;
}


