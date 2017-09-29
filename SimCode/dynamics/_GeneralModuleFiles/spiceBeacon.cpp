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


#include "spiceBeacon.h"
#include "simFswInterfaceMessages/macroDefinitions.h"
#include "utilities/avsEigenMRP.h"


/*--------------------------------------------------------------------------------------------------*/
// SpiceBodyData implementation

/*!
 @brief Use this constructor to use the class as the old structure. Should be deprecated soon.
 */
SpiceBodyData::SpiceBodyData()
{
    this->ephemTime = 0;               //!< [s]      Ephemeris time for the body in question
    this->ephIntTime = 0;              //!< [s]      Integration time associated with the ephem data
    this->radEquator = 0;              //!< [m]      Equatorial radius for the body
    this->geo_albedo = 0;	       //!< [-]      Geometric albedo
    
    // Default these values to zero just in case they don't get populated
    this->localPlanet.J2000Current = 0.0;
    v3SetZero(this->localPlanet.PositionVector);
    v3SetZero(this->localPlanet.VelocityVector);
    m33SetIdentity(this->localPlanet.J20002Pfix);
    m33SetZero(this->localPlanet.J20002Pfix_dot);
}

/*!
 @brief Destructor.
 */
SpiceBodyData::~SpiceBodyData()
{
    return;
}

void SpiceBodyData::initBody(uint64_t moduleID)
{
    this->bodyMsgID = SystemMessaging::GetInstance()->subscribeToMessage(
                    this->bodyInMsgName, sizeof(SpicePlanetStateSimMsg), moduleID);

    // unclear if this is correct below
    //this->radEquator = spherFound ? this->spherHarm.radEquator : this->radEquator;
    //this->geo_albedo = spherFound ? this->spherHarm.radEquator : this->radEquator;
}


void SpiceBodyData::loadEphemeris(uint64_t moduleID)
{
    SystemMessaging::GetInstance()->ReadMessage(this->bodyMsgID, &this->localHeader,
        sizeof(SpicePlanetStateSimMsg), reinterpret_cast<uint8_t *>(&this->localPlanet));
}


SpiceBeacon::SpiceBeacon()
{
    this->vehiclePositionStateName = "hubPosition";
    this->vehicleVelocityStateName = "hubVelocity";
    this->systemTimeCorrPropName = "systemTime";
    this->centralBodyOutMsgName = "central_body_spice";
    this->inertialPositionPropName = "r_BN_N";
    this->inertialVelocityPropName = "v_BN_N";
    return;
}

SpiceBeacon::~SpiceBeacon()
{
    return;
}

void SpiceBeacon::SelfInit()
{
    return;
}

void SpiceBeacon::CrossInit()
{
    //! Begin method steps
    //! - For each spice body in the data vector, find message ID
    //! - If message ID is not found, alert the user and disable message
    std::vector<SpiceBodyData *>::iterator it;
    for(it = this->spiceBodies.begin(); it != this->spiceBodies.end(); it++)
    {
        (*it)->initBody(this->moduleID);
    }
}

void SpiceBeacon::UpdateState(uint64_t CurrentSimNanos)
{
    //! Begin method steps
    //! - For each gravity body in the data vector, find message ID
    //! - If message ID is not found, alert the user and disable message
    std::vector<SpiceBodyData *>::iterator it;
    //for(it = this->spiceBodies.begin(); it != this->spiceBodies.end(); it++)
    //{
    //    (*it)->loadEphemeris(this->moduleID);
    //    if((*it)->isCentralBody)
    //    {
    //        this->centralBody = (*it);
    //   }
    //}
    this->writeOutputMessages(CurrentSimNanos);
}

void SpiceBeacon::writeOutputMessages(uint64_t currentSimNanos)
{
    //if (this->centralBodyOutMsgId > 0) {
    //    SystemMessaging::GetInstance()->WriteMessage(this->centralBodyOutMsgId, currentSimNanos, sizeof(SpicePlanetStateSimMsg), reinterpret_cast<uint8_t*> (&this->centralBody->localPlanet), this->moduleID);
    //}
}

//void SpiceBeacon::registerProperties(DynParamManager& statesIn)
//{
    //Eigen::Vector3d gravInit;
    //gravInit.fill(0.0);
    //this->gravProperty = statesIn.createProperty(this->vehicleGravityPropName, gravInit);
//    this->inertialPositionProperty = statesIn.createProperty(this->inertialPositionPropName, gravInit);
//    this->inertialVelocityProperty = statesIn.createProperty(this->inertialVelocityPropName, gravInit);
//}

void SpiceBeacon::updateInertialPosAndVel()
{
    // Here we explicitly update the system inertial spacecraft position
    // in the spice reference frame if we are computing dynamics
    // relative to a central body
    *this->inertialPositionProperty = this->posState->getState();
    *this->inertialVelocityProperty = this->velState->getState();
   
}

