// WaveProcessor.h
/// *******************
/// Data Analisys of Fast Hadronic Calorimeter - FastHCal
/// Created 27/4/2017
/// by Predrag Ujic
/// *****


#ifndef WAVEPROCESSOR_H
#define WAVEPROCESSOR_H

#define DEBUG 0
#define DEBUG2 0 // 
#define DEBUG3 0 // CreateHistograms debug (
#define CHPM1 1 // first scintillator
#define CHPM2 2 // second
#define CHS3 3 // ch of first trigger
#define CHS4 4 // ch of second trigger
#define BASELINE 4 // ad-hoc baseline in mV
#define NBINS 1024 // usually for DRS4 it's 1024 ...



#include <string>
#include <fstream>
#include <assert.h>


// root

#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TVirtualFFT.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResultPtr.h"


#include "observables.h"




using namespace std;

typedef unsigned short USHORT;
typedef signed short SSHORT; 





struct date {
    USHORT year;
    USHORT month;
    USHORT day;
    USHORT hour;
    USHORT minute;
    USHORT second;
    USHORT milisecond;
    };

struct WaveformParam {
    Float_t arrivalTime;
    Float_t Etot;
    Float_t T90;
    Float_t T70;
    Float_t T50;
    Float_t maxVal;
    Float_t baseLine;
    Float_t baseLineRMS;
    Float_t Eof10ns;
    Float_t FWHM;
    Float_t FW10pcntM; // full width at 10 % of Maximum
};

float convertChtoF (char*) ; // converts ch aray of 4 to the float
USHORT convertChtoUSHORT (char*) ; // converts ch aray of 2 to the unsigned short
unsigned int convertChtoUint (char*); // converts ch aray of 4 to the unsigned int
SSHORT convertChtoSSHORT (char*) ; // converts ch aray of 2 to the unsigned short

// The Unit respose is a response to high energy muon, i.e it's a kind of Green function of the system
// since the amplitudes of the response to muons are not allways the same
// The hadronic shower is constituted of many superimposed Green functions and our goal is to find their positions   
struct UnitPosAmpl{
	Float_t startPosition; // start position of the Unit response 
	Float_t amplitude;// amplitude
	Float_t amplitude_err;
};

class WaveProcessor {

// TimeBinWidth, BinVoltage are showing to which channel they belong to 
// BinVoltage[1][321]CH1, bin 321 CH STARTS WITH 1 !!!!! to corespond to real channel, mark 0 is ignored
// bin starts as usual - form 0

    public:

    WaveProcessor();   //constructor
    ~WaveProcessor();  //destructor

    // setget DRS4 waveform analysis parameters
    void setTriggerHeight(float value) {triggerHeight=value;}
    float getTriggerHeght() const {return triggerHeight;}
    void setDelay(float value) {delay = value;}
    float getDelay() const {return delay;}
    void setEventID(unsigned int value) {eventID = value;}
    int getEventID() const {return eventID;}
    void setDateStamp(date value) {dateStmp=value;}
    date getDateStamp() const {return dateStmp;}
    void SetNoOfChannels(int Ch) {No_of_Ch = Ch;}
    int GetNoOfChannels() const {return No_of_Ch;}
    
    
    void InitializeAnalysisTree();
    void set_time_calibration (int,int,float); // CH, bin, value
    float get_time_calibration (int,int) const;
    void set_bin_time_n_voltage (int, int, SSHORT, SSHORT, USHORT);


    void allignCells0(unsigned short); // align cell #0 of all channels
    void CreateTempHistograms(); // fill the histograms with the colected waveforms of one event, one histogram per channel
    void DeleteTempHistograms(); // they must be deleted at the end of each event
    TH1F* GetTotHist(int Ch) const {return TotShape[Ch];}; // return the histogram of given chanel
    TH1F* GetTempHist(int) const;
    void PrintCurrentHist(int) const; // print pdf of temporary histogram (TempShape) of given channel
    void FillTotHistograms(Float_t ArTm1, Float_t ArTm2, Float_t ArTm3, Float_t ArTm4);

// analysis methods

    WaveformParam give_waveform_parameters(int);  // gives time and amplitude of a given channel
    DRS4_data::Observables* ProcessOnline(Float_t* , Float_t* , Int_t);
    static DRS4_data::Observables* ProcessOnline(Float_t* time, Float_t* amplitude, Int_t length, float threshold, float trigDelay);

    Float_t GetFWHM(int, Float_t);
    Float_t GetFWHM(int, Float_t, Float_t); // third float is in percent, 50% for FWHM, 10 % - width at 10
    
    TH1* FilterFFTofCurrentHist(int);
    void ProcessFile(char*);
    void AmplitudeTimeCorrection(char*);
    UnitPosAmpl give_time_amplitude(int);
    
    private:										// PRIVATE:
       
    static float CalcHistRMS(const TH1F*, int, int );
    static float MeanAndRMS(const TH1F*, int first, int last, float &mean, float &rms);
    static float ArrivalTime(TH1F*, float threshold, float baseline,
                                    float risetime, float fraction = 0.2);
    void RemoveSpikes(float, short); // removes spikes of given threshold (float) in mV and given width (in bins) - usually 2
                                    
    
    TFile* f;
	TTree* paramTree;
    
    float triggerHeight;
    float delay;
    float range;
    float scaller;
    unsigned int eventID;
    date dateStmp;
    float baseLineAVG;
    int baseLineCNT;
    bool aligned; // flag that the 0 cells of chanels are aligned
    int No_of_Ch;
    
    WaveformParam WFParamPM1, WFParamPM2, WFParamS3, WFParamS4;
        
    float TimeBinWidth[5][1024]; // this is the time width of given bin according to the calibration
    float BinVoltage[5][1024]; // this is the voltage of the given bin, for us (-0.5 V, +0.5 V)
    float TimeBinVoltage[5][1024]; // this is the time in ns between the trigger and the given bin

    // should contain the shape of only one event
    TH1F* TempShapeCh1;
    TH1F* TempShapeCh2;
    TH1F* TempShapeCh3;     
    TH1F* TempShapeCh4;
    
    TH1F* TotShape[5];
    
    TH1F* RawTempShape;
  
    
    
};
#endif //WAVEPROCESSOR_H
