/*
 * lightSim.h
 *
 *  Created on: Apr 25, 2017
 *      Author: jp
 */

#ifndef LIGHTSIM_H_
#define LIGHTSIM_H_

#include <cstring>
#include <iostream>
#include <vector>

#include "simMessages/lightSimMsg.h"
#include "_GeneralModuleFiles/sys_model.h"
#include "simMessages/spicePlanetStateSimMsg.h"
#include "simMessages/scPlusStatesSimMsg.h"
#include "simMessages/atmoPropsSimMsg.h"

/*! \addtogroup SimModelGroup
 * @{
 */

//! @brief Container for the properties of a simple illuminated celestial body. */
typedef struct {
    double visMagnitude;            //!< Visual Magnitude
    double RA;                    	//!< deg   Right Ascension in HAEJ2000
    double Dec;                		//!< deg Declination in HAEJ2000
}lightSimProperties;


class LightingSimulation: public SysModel {
public:

	LightingSimulation(std::string in1);
	~LightingSimulation();

	// BSK format
    void SelfInit();
    void CrossInit();
    void UpdateState(uint64_t CurrentSimNanos);
    void WriteOutputMessages(uint64_t CurrentClock);
    bool ReadInputs();

	void setCBodyParam();
	void setCamParam();
    void updateLighting(uint64_t CurrentClock);
	void addSpacecraftToModel(std::string);
	void updateRelativePos(SpicePlanetStateSimMsg& cbodyState, SCPlusStatesSimMsg& scState);

    std::vector<std::string> scStateInMsgNames;		//!< Vector of the spacecraft position/velocity message names
    std::vector<uint64_t> scStateInMsgIds;
    std::vector<SCPlusStatesSimMsg> scStates;

    std::vector<std::string> LightSimOutMsgNames; 	//!< Vector of strings containing light sim output message names
    std::vector<uint64_t> LightSimMsgIds;

    std::string cbodyPosInMsgName;					//!< Message name for the planet's SPICE position message
    int64_t cbodyPosInMsgId;
    SpicePlanetStateSimMsg bodyState;

    std::string cb_name;

	/* needs to be pulled from outside module, units in km */
	double cb_albedo;
	double cb_radius;
	double cb_position[3];
	double obs_position[3];
	double dcm_BN_obs[3][3];

	/* locally generated or passed variables */
	double cb_distance;							//HAEJ2000 km
	double obs2cb_position[3];					//HAEJ2000 km
	double obs2cb_position_body[3];				//Camera Body km
	double cb2obs_distance;						//km

	/* variables used by other modules */
	double vis_magnitude;
	double RaDec[2];							//HAEJ200 RA and Dec deg

private:

	void mapSphere();
	void lumos();
	void project2CamView();
	void getVMagnitude();
	void getRaDec();
	void printMessage();
	bool checkFoV(double pos[3]);

	std::vector<lightSimMsg> lSimOutBuffer; //!< -- Message buffer for density messages
	uint64_t OutputBufferCount;	//!< number of output buffers for messaging system

	double FLUX_ref;
	double FLUX_ref_distance;
	int npts_spheremap;
	int npts_visible;
	const int max_spheremap = 1000;
	double spheremap_latlong[1000][3];
	double spheremap_ijk[1000][3];
	double spheremap_ijk_visible[1000][3];
	double spheremap_albedo[1000];
	double spheremap_albedo_visible[1000];
	double facet_area;
	int res_lat;
	int res_long;

	// variables pulled from other messages
	double vert_fov;
	double horiz_fov;

};

#endif /* LIGHTSIM_H_ */
