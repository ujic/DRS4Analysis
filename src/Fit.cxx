#include "../include/WaveProcessor.h"
#include "../include/Fit.h"
/*
Double_t fitValEst(Double_t *x, Double_t *par){
	return WaveProcessor::Bridge(x, par);
}
	
Double_t WaveProcessor::Bridge(Double_t *x, Double_t *par) {
	return ftotal(x, par);
}	
*/


Float_t FitUnitResponses(TH1F *TempShape, TH1F* UnitShape, int NUnits, UnitPosAmpl* Unit_in_Peak, const float triggerHeight, int eventNo, char* filename, Float_t TRef, float UnitAvEn){
	/*** if the signal amplitude is less than E(mip) it has difficulties to fit it ***/
	/*** thus, I force to fit to only one UnitResponse, beacause it is better to have it's precise timing ***/
	/*** than to have some questionable form-structure ***/
	
	Float_t chi2;
	int i;
	bool drawFlag(false), minimalAmplFlag(false);
	char histfilename[70];
	float delay(DELAY);
	float fract(FRACTION);
	int first(FIRST), last(LAST);
	assert(NUnits<=40);
	string stringPriv;
	//freeUR = FREE_FIT_UR;
/*
	extern int NUnitsBridge;
extern TH1F *TempShapeBridge;
 extern TH1F *UnitShapeBridge;
 extern Float_t refArrTimeShape;
*/
	//PeakDerivatives *PeakParamTmpSh;
	fract/=100.; // because FRACTION is in %
	Float_t baseLineEnd(BASELINEEND);
	m_gStyle = new TStyle();
	//baseLineTS = TempShape->Integral(0, TempShape->FindBin(baseLineEnd-delay), "width") / (baseLineEnd-delay);
	baseLineUS = UnitShape->Integral(0, UnitShape->FindBin(baseLineEnd-delay), "width") / (baseLineEnd-delay);	 
	NUnitsBridge=NUnits;
	TempShapeBridge = (TH1F*)TempShape->Clone();
	UnitShapeBridge = (TH1F*)UnitShape->Clone();
	BaseLineShape = WaveProcessor::GetBaseLineHistogramRaw();
	baseLineBL = BaseLineShape->Integral(0, BaseLineShape->FindBin(baseLineEnd-delay), "width") / (baseLineEnd-delay);
	if (DEBUG5) cout<<"baseLineBL="<<baseLineBL<<", baseLineTS="<<baseLineTS<<endl;
	//TempShapeBridge->Add(BaseLineShape, -1); <- making problems on many basis
	
	for (i=1; i<=1024; i++) TempShapeBridge->SetBinContent(i, TempShapeBridge->GetBinContent(i)-BaseLineShape->GetBinContent(i));
	
	baseLineTS=TempShapeBridge->Integral(0, TempShapeBridge->FindBin(baseLineEnd-delay), "width") / (baseLineEnd-delay);
	//TF1 *background = new TF1("background", "-1*baseLineTS", 0, 200);
	//TempShapeBridge->Add(background, 1); 
	
	
	if (DEBUG5) cout<<"tirgHeight="<<triggerHeight<<", fract="<<fract<<endl;
	refArrTimeShape = WaveProcessor::ArrivalTime2(TempShape, baseLineTS, fract, triggerHeight);
	refArrUnitShape = WaveProcessor::ArrivalTime2(UnitShape, baseLineUS, fract, triggerHeight);
	
	NUnitsBridge=NUnits;
	
	stringPriv="";
	stringPriv.append(filename);
	stringPriv.erase(stringPriv.end()-4, stringPriv.end());
	stringPriv.append("_evt");	
	stringPriv.append(std::to_string(eventNo));
	stringPriv.append(".pdf");
		
	std::size_t lengthS = stringPriv.copy(histfilename, stringPriv.length(),0);
	histfilename[lengthS]='\0';
	TCanvas *c = new TCanvas("fit", histfilename, 1);
	
	minimalAmplFlag=(TempShapeBridge->GetMaximum()/UnitShapeBridge->GetMaximum())<1.;
	if (minimalAmplFlag) {
		NUnitsBridge=0;
		freeUR=1;
		last=Unit_in_Peak[0].position+TRef+10.;
		//cout<<"last="<<last<<endl;
	}
	
	TF1 *func = new TF1("fit", ftotal, first, last, (NUnitsBridge+freeUR)*2+1); 
	//(each found peak is two param (position, amplitude) + additional peak (2 param) + background (1 param); 
	
	
	if (DEBUG5) cout<<"NUnits="<<NUnits<<endl;
	if (!minimalAmplFlag){
		for (i=0;i<=(NUnitsBridge-1);i++){
			if (DEBUG5) cout<<"Initial Unit_in_Peak["<<i<<"].amplitude="<<Unit_in_Peak[i].amplitude*UnitAvEn<<endl;
			if (DEBUG5) cout<<"Initial Unit_in_Peak["<<i<<"].position="<<Unit_in_Peak[i].position+TRef<<endl;
			if (Unit_in_Peak[i].position>5.) drawFlag = true;
			func->SetParameter(i*2, Unit_in_Peak[i].amplitude);
			func->SetParLimits(i*2, Unit_in_Peak[i].amplitude*0.8,Unit_in_Peak[i].amplitude*1.2);
			func->SetParameter(i*2+1, Unit_in_Peak[i].position);//+TempShapeBridge->GetBinCenter(refArrTimeShape));
			func->SetParLimits(i*2+1, (Unit_in_Peak[i].position-4.),(Unit_in_Peak[i].position+4.));
		}
		func->SetParameter(NUnitsBridge*2, 0.3*TempShape->GetMaximum()); //small amplitude for late component
		func->SetParameter(NUnitsBridge*2+1, Unit_in_Peak[i].position+5.); // search around 15 ns after first arrival
		//func->SetParameter(NUnits*2+2, 0.3); //small amplitude for late component
		//func->SetParameter(NUnits*2+3, Unit_in_Peak[i].position); // search around 15 ns after first arrival	

	}
	else {
		func->SetParameter(0, Unit_in_Peak[0].amplitude); //use just initial founding, but NO LIMITS
		func->SetParameter(1, Unit_in_Peak[1].position-2.); // 
	}

	func->SetParameter((NUnitsBridge+freeUR)*2, 1.); // baseLine ~1 mV
	
	if (DEBUG5) cout<<TempShapeBridge->GetBinCenter(refArrTimeShape)<<endl;
	if (DEBUG5) cout<<UnitShapeBridge->GetBinCenter(refArrUnitShape)<<endl;
	if (DEBUG5) cout<<TRef<<endl;
	
	m_gStyle->SetOptStat(1111);
	m_gStyle->SetOptFit(1111);
	//cout<<"EventNo="<<eventNo<<endl;
	//TempShapeBridge->Draw("");
    TempShapeBridge->Fit("fit","RQ");

    //TempShapeBridge->Fit("fit","RMEQ");     // repeat to improve fit
    

   
    //UnitShapeBridge->Draw("same");
    //BaseLineShape->Draw("same");
	for (i=0;i<=(NUnitsBridge+freeUR-1);i++){// FreeUR
		Unit_in_Peak[i].amplitude=func->GetParameter(i*2);
		Unit_in_Peak[i].position= func->GetParameter(i*2+1);

		if (DEBUG5) cout<<"Fitted Unit_in_Peak["<<i<<"].amplitude="<< Unit_in_Peak[i].amplitude<<endl;//*UnitAvEn<<endl;
		if (DEBUG5) cout<<"Fitted Unit_in_Peak["<<i<<"].position ="<<  Unit_in_Peak[i].position<<endl;//+TRef <<endl;
	}

	chi2=func -> GetChisquare();
	//cout<<"chi2/ndf="<<chi2/(float)((last-first)*5-NUnitsBridge*2+4+1)<<endl;
	chi2=chi2/(float)((last-first)*5-(NUnitsBridge+freeUR)*2-1);


   TPaveText *pave = new TPaveText(.5, .74, .9, .88, "NDC");
   pave->SetTextFont(42);
   pave->SetTextSize(.05);
   pave->SetBorderSize(0);
   pave->SetFillStyle(4001);
   pave->SetFillColor(kWhite);
   pave->AddText(Form("#chi^{2}/N_{df} = %4.02f", chi2));
   pave->Draw();
	
	if (DEBUG) cout<<"Histogram Draw"<<endl<<flush;
	if (drawFlag||(eventNo<50)) c->SaveAs(histfilename);
	
    delete c;
	delete func;
	delete pave;
	TempShapeBridge->Delete();
	UnitShapeBridge->Delete();
	//cout<<"chi2 after="<<chi2<<endl;	
	return chi2;
   
}
Double_t ftotal(Double_t *x, Double_t *par) {
// remember the true position depends on the fraction
/// the refferent point must be the crossBin of the Temp Shape
/// also take care when a UnitResponses is coming after the Tempshape you have to add some values to at the begining
/// find which bin corresponds to crossBin, that one should be the refferent
// you need to match the positions of the histogram and the position of the UnitShape 

	Double_t value(0);
	Int_t bin(0);
	Double_t xx = x[0];
	/*
	extern int NUnitsBridge;
	extern TH1F *TempShapeBridge;
 extern TH1F *UnitShapeBridge;
extern Float_t refArrTimeShape;
*/
	//if (DEBUG5) cout<<"NUnits in ftotal = "<<NUnitsBridge;
	for (int i=0; i<(NUnitsBridge+freeUR); i++){
		bin = UnitShapeBridge->GetXaxis()->FindBin(xx-par[1+i*2]); //+TempShapeBridge->GetBinCenter(refArrTimeShape));
						//-UnitShapeBridge->GetBinCenter(refArrUnitShape));
		//if (DEBUG5) cout<<"bin="<<bin<<endl;
		value += par[i*2]*UnitShapeBridge->GetBinContent(bin);
	}
	value += par[(NUnitsBridge+freeUR)*2];
	//if (DEBUG5) cout<<" value="<<value;
	return value;
}


