
#ifndef ImageProcessingMsg_H
#define ImageProcessingMsg_H

typedef struct {
    double beaconID[1000]		//!< current reference catalog ID of visible beacons in image
    double beaconPL[1000][2]		//!< current pixel/line coordinate of visible beacons 
    double sigma_BN[3];               	//!< current MRPs (inertial)
    uint64_t timeImage			//!< time of image [ns]  				
}ImageProcessingMsg;
#endif

