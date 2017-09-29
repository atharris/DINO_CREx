/*
 ISC License

 Copyright (c) 2016-2017, Autonomous Vehicle Systems Lab, University of Colorado at Boulder

 Permission to use, copy, modify, and/or distribute this software for any
 purpose with or without fee is hereby granted, provided that the above
 copyright notice and this permission notice appear in all copies.

 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

 */

#ifndef SpiceInterface_H
#define SpiceInterface_H

#include "dynamicEffector.h"
#include "_GeneralModuleFiles/sys_model.h"
#include "utilities/linearAlgebra.h"
#include "simMessages/spiceTimeSimMsg.h"	
#include "architecture/messaging/system_messaging.h"
#include <vector>
#include <Eigen/Dense>
#include "simMessages/spicePlanetStateSimMsg.h"


//!@brief Container for spice body data
/*! This class is designed to hold all of the information for a Spice celestial
 body.  The nominal use-case has it initialized at the python level and
 attached to dynamics using the AddGravityBody method.
 */
class SpiceBodyData
{
public:
    // Default constructor
    SpiceBodyData();
    ~SpiceBodyData();
    
    void initBody(uint64_t moduleID); 		//!<        Method to initialize the gravity body
    void loadEphemeris(uint64_t moduleID); 	//!< Command to load the ephemeris data
    
public:
    double ephemTime;               //!< [s]      Ephemeris time for the body in question
    double ephIntTime;              //!< [s]      Integration time associated with the ephem data
    double radEquator;              //!< [m]      Equatorial radius for the body
    double geo_albedo; 		    //!< [-]      Geometric albedo for body
    SpicePlanetStateSimMsg localPlanet;	//!< [-]      Class storage of ephemeris info from scheduled portion
    SingleMessageHeader localHeader;	//!  [-]      Header information for ephemeris storage
    std::string bodyInMsgName;      //!<          Gravitational body name
    std::string outputMsgName;      //!<          Ephemeris information relative to display frame
    std::string planetEphemName;    //!<          Ephemeris name for the planet
    int64_t outputMsgID;            //!<          ID for output message data
    int64_t bodyMsgID;              //!<          ID for ephemeris data message
};


/*! \addtogroup SimModelGroup
 *  This group is used to model parts of the vehicle and the surrounding environment
 *  in the simulation system.  All components/dynamics/environment models are a
 *  part of this group.
 * @{
 */
//! The SPICE interface class gets time and planetary body information from the JPL ephemeris library. Also, allows for message forming and writing
class SpiceBeacon: public SysModel {
public:

    SpiceBeacon();
    ~SpiceBeacon();
    void SelfInit();
    void CrossInit();
    void UpdateState(uint64_t CurrentSimNanos);
    //void linkInStates(DynParamManager& statesIn);  	 // may be unecessary
    //void registerProperties(DynParamManager& statesIn);  // may be unecessary 
    void updateInertialPosAndVel();


private:
    void writeOutputMessages(uint64_t currentSimNanos);
    
public:
    // Originally Gravity Effector
    std::string vehiclePositionStateName;          //! [-] Name of the vehicle position state
    std::string vehicleVelocityStateName;          //! [-] Name of the vehicle position state
    std::string systemTimeCorrPropName;            //! [-] Name of the correlation between times
    std::vector<SpiceBodyData*> spiceBodies;       //! [-] Vector of Spice bodies
    std::string inertialPositionPropName;          //! [-] Name of the inertial position property
    std::string inertialVelocityPropName;          //! [-] Name of the inertial velocity property

private:
    // originally GravityEffector
    StateData *posState;                            //! [-] Position state of the vehicle
    StateData *velState;                            //! [-] Position state of the vehicle
    StateData *hubSigma;                            //! [-] sigmaBN for the hub
    //Eigen::MatrixXd *gravProperty;                  //! [-] g_N property for output
    Eigen::MatrixXd *timeCorr;                      //! [-] Time correlation property
    int64_t centralBodyOutMsgId;                //! [-] Id for the central body spice data output message
    std::string centralBodyOutMsgName;              //! [-] Unique name for the central body spice data output message
    Eigen::MatrixXd *inertialPositionProperty;             //! [m] r_N inertial position relative to system spice zeroBase/refBase coordinate frame, property for output.
    Eigen::MatrixXd *inertialVelocityProperty;             //! [m/s] v_N inertial velocity relative to system spice zeroBase/refBase coordinate frame, property for output.
};





/*! @} */

#endif
