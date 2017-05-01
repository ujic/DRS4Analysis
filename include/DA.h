// DA.h
/// *******************
/// Data Analisys of Fast Hadronic Calorimeter - FastHCal
/// Created 27/4/2017
/// by Predrag Ujic
/// *****


#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>

//root

#include <TSystem.h>

//DA
#include "WaveProcessor.h"

WaveProcessor* DataHandler = new WaveProcessor();

float aux_f;
unsigned int eventSerNo;
//USHORT year, month, day, hour, minute, second, milisecond, aux_USHORT;
//SSHORT range, aux_SSHORT;
float scaller;

int i, j, k;













