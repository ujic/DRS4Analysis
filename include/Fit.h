#ifndef FIT_H
#define FIT_H

#include <TStyle.h>
#include <TPaveText.h>

int NUnitsBridge;
int freeUR(FREE_FIT_UR);
Float_t refArrTimeShape; // this is to tell the fit function where is the refference bin
Float_t refArrUnitShape;
TH1F *UnitShapeBridge;
TH1F *TempShapeBridge;
Float_t baseLineTS;
Float_t baseLineUS;
TH1F *BaseLineShape = WaveProcessor::GetBaseLineHistogramRaw();
Float_t	baseLineBL;

TStyle* m_gStyle;

#endif
