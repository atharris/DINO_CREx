/*
 * lightSim.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: jp
 */

#include "lightSim.h"
#include "utilities/linearAlgebra.h"
#include "utilities/astroConstants.h"
#include "utilities/rigidBodyKinematics.h"
#include "architecture/messaging/system_messaging.h"
#include "simMessages/lightSimMsg.h"
#include "architecture/messaging/blank_storage.h"
#include "architecture/system_model/sys_model_task.h"
#include "dynamics/_GeneralModuleFiles/stateData.h"
#include "_GeneralModuleFiles/sys_model.h"
#include "environment/spice/spice_interface.h"
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
using namespace std;

LightingSimulation::LightingSimulation(std::string cbname) {

	this->cb_name = cbname;
	setCBodyParam();
	setCamParam();
}

LightingSimulation::~LightingSimulation() {
	return;
}

void LightingSimulation::addSpacecraftToModel(std::string tmpScMsgName){
  std::string tmpLSimMsgName;
  this->scStateInMsgNames.push_back(tmpLSimMsgName);
  tmpLSimMsgName = "light_sim"+ std::to_string(this->scStateInMsgNames.size()-1)+"_data";
  this->LightSimOutMsgNames.push_back(tmpLSimMsgName);
  return;
}

/*! SelfInit for this method creates a seperate light sim message for each of the spacecraft
that were added using AddSpacecraftToModel.
*/
void LightingSimulation::SelfInit()
{
    std::string tmpLightSimName;
    uint64_t tmpLightSimId;
    //! Begin method steps
    std::vector<std::string>::iterator it;
    std::vector<std::string>::iterator nameIt;

    for(it = this->LightSimOutMsgNames.begin(); it!=this->LightSimOutMsgNames.end(); it++){
    	tmpLightSimId = SystemMessaging::GetInstance()->CreateNewMessage(*it, sizeof(lightSimMsg), this->OutputBufferCount, "lightSimMsg", moduleID);

      this->LightSimMsgIds.push_back(tmpLightSimId);
    }

    return;
}

void LightingSimulation::setCBodyParam() {

	/* placeholder --- update with parameters pulled from BSK or custom database file*/
	if (this->cb_name == "Mercury"){
		this->cb_radius = 2439700;		// km
		this->cb_albedo = 0.142;
	}
}

void LightingSimulation::setCamParam() {

	// PLACEHOLDER - INSERT FLUX CALCULATION FROM SORCE DATA BASED on Camera wavelength limits
	FLUX_ref = SOLAR_FLUX_EARTH;
	FLUX_ref_distance = AU;

	FLUX_ref = 5.63898324842E20;  // W/km^2
	FLUX_ref_distance =  695700;  // km 

	// PLACEHOLDER - Pull camera attitude from Cam_model message (or store relative cam-sc attitude and pull sc attitude)
	//SCPlusStatesSimMsg

	// PLACEHOLDER - Pull camera FoV limits from Cam_model message
    horiz_fov = 90;
    vert_fov = 90;
}

void LightingSimulation::CrossInit()
{
  // assume SpicePlanetStateSimMsg has celestial bodies
  this->cbodyPosInMsgId = SystemMessaging::GetInstance()->subscribeToMessage(
		  this->cbodyPosInMsgName, sizeof(SpicePlanetStateSimMsg), moduleID);

  std::vector<std::string>::iterator it;
  for(it = scStateInMsgNames.begin(); it!=scStateInMsgNames.end(); it++){
    this->scStateInMsgIds.push_back(SystemMessaging::GetInstance()->subscribeToMessage(*it, sizeof(SCPlusStatesSimMsg), moduleID));
  }
  return;
}

bool LightingSimulation::ReadInputs()
{
  SCPlusStatesSimMsg tmpState;
  //! Begin method steps
  SingleMessageHeader localHeader;
  memset(&this->bodyState, 0x0, sizeof(SpicePlanetStateSimMsg));
  memset(&tmpState, 0x0, sizeof(SCPlusStatesSimMsg));
  scStates.clear();
  if(scStateInMsgIds[0] >= 0)
  {
    //! Iterate over spacecraft message ids
    std::vector<uint64_t>::iterator it;
    for(it = scStateInMsgIds.begin(); it!= scStateInMsgIds.end(); it++){
      SystemMessaging::GetInstance()->ReadMessage(*it, &localHeader,
                                                  sizeof(SCPlusStatesSimMsg),
                                                  reinterpret_cast<uint8_t*>(&tmpState),
                                                  moduleID);
    this->scStates.push_back(tmpState);

    //assuming only one s/c at the moment
    v3Copy(tmpState.r_BN_N, this->obs_position);
    v3Scale(1/1000, this->obs_position, this->obs_position);

    //convert mrp to DCM and find net camera dcm_BN
    double dcm_BN[3][3];
    double dcm_SB[3][3];
    MRP2C(tmpState.sigma_BN, dcm_BN);
    m33Transpose(tmpState.dcm_BS, dcm_SB);
    m33MultM33(dcm_SB, dcm_BN, this->dcm_BN_obs);

      }
  }

  if(cbodyPosInMsgId >= 0)
  {
      SystemMessaging::GetInstance()->ReadMessage(this->cbodyPosInMsgId , &localHeader,
                                                  sizeof(SpicePlanetStateSimMsg), reinterpret_cast<uint8_t*>(&this->bodyState), moduleID);
  }
  return(true);
}

void LightingSimulation::updateRelativePos(SpicePlanetStateSimMsg& cbodyState, SCPlusStatesSimMsg& scState)
{
    if(cbodyPosInMsgId >= 0)
    {
      v3Subtract(cbodyState.PositionVector, scState.r_BN_N, this->obs2cb_position);
      m33MultV3(this->dcm_BN_obs, this->obs2cb_position, this->obs2cb_position_body);
    }
    //std::cout<<"Relative Pos: "<<this->relativePos <<std::endl;
    return;
}


bool LightingSimulation::checkFoV(double pos[3]){

	// pos needs to be in observer body frame
	// Need to update to new camera body frame (currently +i to right, +j forward, +k updward)!	

	double pos_norm;
	double az;
	double el;
	bool chk;

	pos_norm = v3Norm(pos);

	az = atan2(pos[0], pos[1])*R2D;
	el = asin(pos[2]/pos_norm)*R2D;

	cout << az << "  ,  " << el << endl;

	if (abs(az) < horiz_fov && abs(el) < vert_fov){
		chk = true;
	}
	else
		chk = false;

	return chk;

}

void LightingSimulation::mapSphere() {

	
	/* placeholder -- update with pulling from database file or variable resolution maps */
	double delta_lat;
	double delta_long_equator;
	double delta_long_distance;

	// placeholder -- update with variable resolution based on distance to obs-CB distance
	this-> res_lat = 30;
	this-> res_long = 30;

	delta_lat = 180/res_lat;
	delta_long_equator = 180/res_long;

	// Find arclength between longitude points - map longitude with arclength constant at all latitudes
	delta_long_distance = cb_radius * delta_long_equator*D2R;

	int npts = 0;
	int ind_c = 0;
	for (int kk = 1 ; kk <= res_lat; kk++){

		// compute latitude and longitude for facet vertices
		double current_lat = -90.0 + (kk-.5)*delta_lat;

		if  (current_lat <= 90) {

			// horizontal semi-circle radius at current latitude
			double current_radius  = abs(cos(current_lat*D2R)*cb_radius);
			// number of latitude points required at current latitude to keep horizontal spacing constant
			double current_nlong = floor((M_PI*current_radius)/delta_long_distance);
			// longitude spacing at current latitutde (deg)
			double current_delta_long = (delta_long_distance / current_radius)*R2D;
			// calculate initial longitude value (make surface mapping points centered horizontally)
			double long_i = (180 - (current_nlong*current_delta_long))/2;

			// calculate center facet locations
			double current_radius_c  = abs(cos(current_lat*D2R)*cb_radius);
			double current_nlong_c = floor((M_PI*current_radius_c)/delta_long_distance);
			double current_delta_long_c = (delta_long_distance / current_radius_c)*R2D;

			if (current_lat < 90) {

				//cout << "Delta Longitude  " <<  current_delta_long_c << "  Latitude " << current_lat << endl;

				for (int kk = 0; kk < current_nlong_c; kk++){
					spheremap_latlong[ind_c][0] = current_lat;

					if ( kk == 0 ){
						if (current_nlong_c == 1){
							spheremap_latlong[ind_c][1] = 90;
						}
						else {
							spheremap_latlong[ind_c][1] = long_i + .5*current_delta_long_c;
						}
					}
					else {
						spheremap_latlong[ind_c][1] = long_i + (kk+.5)*current_delta_long_c;
					}
					ind_c++;
					npts++;
				}
			}
		}
		else{
			for (int nn = npts; nn < max_spheremap ; nn++){
				spheremap_latlong[ind_c][0] = nan("");
				spheremap_latlong[ind_c][1] = nan("");
				ind_c++;
			}
		}
	}

	// populate empty values with NaN
	for (int nn = npts; nn < max_spheremap ; nn++){
		spheremap_latlong[ind_c][0] = nan("");
		spheremap_latlong[ind_c][1] = nan("");
		ind_c++;
	}

	// check values
	bool chk = true;
	if (chk == true) {
		for (int mm = 0; mm < npts+5; mm++){
			cout << "Facet Center Lat Long, " << spheremap_latlong[mm][0] << ",  " << spheremap_latlong[mm][1] << endl;
		}
	}

	// Calculate facet area - all facets will be equal totalling semi-spherical surface area
	this->facet_area = (4*M_PI*pow(cb_radius,2)) / npts;
	
	// cout << "Facet Area, "<< this->facet_area << " # Pts, " << npts << endl;
}

void LightingSimulation::lumos() {

	double ibody_helio[3];
	double ijnormbody_helio[3];
	double kbody_helio[3];
	double jbody_helio[3];
	double dcm_NB[3][3];
	double e_sun_body[3] = {1,0,0};

	// Transform CB position from helio-inertial to helio-body frame {+i from CB to sun, +j normal i in helio-inertial IJ plane, +k normal to IJ plane}

    // unit vector to sun (+i)
    v3Normalize((this->cb_position), ibody_helio);
    v3Scale(-1, ibody_helio, ibody_helio);

    // calculate other unit vector directions (+j, +k)
    v3Set(-ibody_helio[1], ibody_helio[0], 0, ijnormbody_helio);
    v3Normalize(ijnormbody_helio, ijnormbody_helio);
    v3Cross(ibody_helio, ijnormbody_helio, kbody_helio);
    v3Normalize(kbody_helio, kbody_helio);
    v3Cross(kbody_helio, ibody_helio, jbody_helio);

    // DCM from helio-body to an inertia frame
    m33Set(ibody_helio[0], jbody_helio[0], kbody_helio[0],
    		ibody_helio[1], jbody_helio[1], kbody_helio[1],
			ibody_helio[2], jbody_helio[2], kbody_helio[2], dcm_NB);
    m33PrintScreen("Helio NB matrix: ", dcm_NB);

    // generate normal unit vector of surface normal for each surface points in body frame then convert to helio frame
    for (int aa = 0; aa < max_spheremap; aa++){
    	if (!isnan(spheremap_latlong[aa][0])){
			double current_lat = spheremap_latlong[aa][0];
			double current_long = spheremap_latlong[aa][1];

			double zbody = sin(current_lat*D2R);
			double ybody = cos(current_lat*D2R)*cos(current_long*D2R);
			double xbody = cos(current_lat*D2R)*sin(current_long*D2R);

			double enorm_body[3];
			v3Set(xbody, ybody, zbody, enorm_body);

			//cout << "Body Vector: "<< "\t" << xbody << "\t" << ybody << "\t" << zbody<< endl;

			// transform into helio-inertial coordinates
			double enorm_helio[3];
			m33MultV3(dcm_NB, enorm_body, enorm_helio);
			this->spheremap_ijk[aa][0] = enorm_helio[0];
			this->spheremap_ijk[aa][1] = enorm_helio[1];
			this->spheremap_ijk[aa][2] = enorm_helio[2];

			// cout << this->spheremap_ijk[aa][0] << ",  " << this->spheremap_ijk[aa][1]<< ",  " << this->spheremap_ijk[aa][2] << endl;

			// albedo due to solar phase angle and due to flux decay
			double phase_solar = acos(v3Dot(enorm_body, e_sun_body))*R2D;
			double flux_decay_sun2cb = pow((FLUX_ref_distance/cb_distance),2);

			// Note: will have to move this if multiple observers are to be considered
			double flux_decay_cb2obs = pow((cb_radius/cb2obs_distance),2);

			spheremap_albedo[aa] = cb_albedo * cos(phase_solar*D2R) * flux_decay_sun2cb * FLUX_ref;
			// cout << spheremap_albedo[aa] << endl;
    	}
    	else {
    		this->spheremap_ijk[aa][0] = nan("");
    		this->spheremap_ijk[aa][1] = nan("");
    		this->spheremap_ijk[aa][2] = nan("");
    		this->spheremap_albedo[aa] = nan("");
    	}
    }
}

void LightingSimulation::project2CamView() {

	// compute position of celestial body relative to observer in observer body coordinates
	double e_obs2cb_helio[3];
	double e_obs2cb_body[3];
	double pt_body[3];

	v3Normalize(obs2cb_position, e_obs2cb_helio);
	v3Normalize(obs2cb_position_body, e_obs2cb_body);

	// enter visible surface points into spheremap_ijk_visible and spheremap_albedo_visible
	int ind_vis = 0;
	for (int cc = 0 ; cc < max_spheremap; cc++){

		if ( !isnan(spheremap_ijk[cc][0])){
			double current_ijk[3] = {spheremap_ijk[cc][0], spheremap_ijk[cc][1], spheremap_ijk[cc][2]};

			// if phase angle is > 90 degree for obs2cb vector and surface normal, then surface point is visible, store current values in visible arrays
			if (abs(acos(v3Dot(e_obs2cb_helio, current_ijk))*R2D) > 90){

				double pt_position_body[3];
				double current_ijk_scaled[3];
				double current_ijk_scaled_body[3];

				v3Scale(cb_radius, current_ijk, current_ijk_scaled);
				m33MultV3(dcm_BN_obs, current_ijk_scaled, current_ijk_scaled_body);
				v3Add(obs2cb_position_body, current_ijk_scaled_body, pt_position_body);
				v3PrintScreen("Pt position body" , pt_position_body);

				if (checkFoV(pt_position_body) == true) {

					v3Set(spheremap_ijk[cc][0], spheremap_ijk[cc][1], spheremap_ijk[cc][2], spheremap_ijk_visible[ind_vis]);
					spheremap_albedo_visible[ind_vis] = spheremap_albedo[cc];
					ind_vis++;
				}
			}
		}
	}
	this->npts_visible = ind_vis;
}

void LightingSimulation::getVMagnitude() {
	/* calculate the visual magnitude of a body using the sun as a reference */

	double cb_lum = 0;
	double e_obs2cb_body[3];
	double obs2sun_distance;
	double flux_sun_obs;
	double flux_cb_obs;
	double m_sun = -26.73;

	// assume LOS vector sufficient for now
	v3Normalize(obs2cb_position_body, e_obs2cb_body);

	// calculate distance to sun for magnitude comparison
	obs2sun_distance = v3Norm(obs_position);

	for( int aa = 0; aa < npts_visible; aa++) {

		//double facet_area_effective;
		double facet_lum;

		// m33MultV3(dcm_BN_obs, current_ijk, current_ijk_body);
		// facet_area_effective = facet_area*abs(v3Dot(e_obs2cb_body, current_ijk_body));

		// compute luminance from facet and add to sum
		facet_lum = facet_area * spheremap_albedo_visible[aa];
		// cout << facet_area_effective << "  ,  "  << spheremap_albedo_visible[aa] << endl;
		cb_lum += facet_lum;
		cout << cb_lum << "  , " << facet_area << "  ,  " << spheremap_albedo_visible[aa] << endl;
	}

	// calculate sun's irradiance at observer
	flux_sun_obs = FLUX_ref* pow((FLUX_ref_distance/obs2sun_distance),2);
	flux_cb_obs = cb_lum/ (2*M_PI* pow(cb2obs_distance,2));

	// calculate visual magnitude using sun's value at 1AU as a reference
	this->vis_magnitude = m_sun - 2.5*log10(flux_cb_obs / flux_sun_obs);
}

void LightingSimulation::getRaDec(){

	// RA: from (+I) in helio IJ plane [0 360] (assumed to be 1st point of Aries in ecliptic plane)
	// Dec: angle upwards from helio IJ plane [-90 90]

	if (this->cb_position[0] == 0){
		this->RaDec[0] = 0;
	}
	else if (cb_position[0]<0){
		this->RaDec[0] = 360- R2D*acos(cb_position[1]/cb_position[0]);
	}
	else
		this->RaDec[0] = R2D*acos(cb_position[1]/cb_position[0]);

	this->RaDec[1] = R2D*asin(cb_position[2]/cb_distance);
}

void LightingSimulation::printMessage() {

	cout << "CB Distance: " << this->cb_distance << endl;
	cout << this->cb_name << endl;
	cout << "Visual Magnitude: " << this->vis_magnitude << endl;
	cout << "Ra Dec in Inertial Frame: " << this->RaDec[0] << ", " << this->RaDec[1] << endl;
}

void LightingSimulation::WriteOutputMessages(uint64_t CurrentClock)
{
    lightSimMsg tmpLSim;
    std::vector<uint64_t>::iterator it;
    std::vector<lightSimMsg>::iterator lSimIt;

    //How Does this Work????
    lSimIt = lSimOutBuffer.begin();
    for(it = LightSimMsgIds.begin(); it!= LightSimMsgIds.end(); it++, lSimIt++){
      tmpLSim = *lSimIt;
      //std::cout<<"WriteMsg: "<<tmpAtmo.neutralDensity<<std::endl;
      SystemMessaging::GetInstance()->WriteMessage(*it,
                                                  CurrentClock,
                                                  sizeof(lightSimMsg),
                                                  reinterpret_cast<uint8_t*>(&tmpLSim),
                                                  moduleID);
    }
}

void LightingSimulation::updateLighting(uint64_t CurrentSimNanos) {

	std::vector<SCPlusStatesSimMsg>::iterator it;
	lightSimMsg tmpData;
	uint64_t lSimInd = 0;
	this->lSimOutBuffer.clear();

	for(it = scStates.begin(); it != scStates.end(); it++, lSimInd++){

		this->updateRelativePos(this->bodyState, *it);

		if (this->checkFoV(this->obs2cb_position_body)){

			this->lumos();
			this->project2CamView();
			this->getVMagnitude();
			this->getRaDec();

			tmpData.visMagnitude = this->vis_magnitude;
			tmpData.RA = this->RaDec[0];
			tmpData.Dec = this->RaDec[1];

			this->lSimOutBuffer.push_back(tmpData);
		}
	}
}

void LightingSimulation::UpdateState(uint64_t CurrentSimNanos)
{
    if(this->ReadInputs())
    {
        updateLighting(CurrentSimNanos*1.0E-9);
    }
    WriteOutputMessages(CurrentSimNanos);

}

