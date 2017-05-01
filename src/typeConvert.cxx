#include "../include/WaveProcessor.h"

float convertChtoF (char* char_array) {
	float the_float;
	memcpy(&the_float, char_array, sizeof the_float);
	return the_float;
}
		
USHORT convertChtoUSHORT (char* char_array) {
	USHORT shorty;
	memcpy(&shorty, char_array, sizeof shorty);
	return shorty;
}

unsigned int convertChtoUint (char* char_array) {
	unsigned int Uint;
	memcpy(&Uint, char_array, sizeof Uint);
	return Uint;
}

SSHORT convertChtoSSHORT (char* char_array) {
	SSHORT shorty;
	memcpy(&shorty, char_array, sizeof shorty);
	return shorty;
}
