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

#ifndef OPNAV_CAMERA_H
#define OPNAV_CAMERA_H

#include <vector>
//#include "_GeneralModuleFiles/sys_model.h"
//#include "dynamics/spacecraftPlus/spacecraftPlus.h"
//#include "environment/spice/spice_interface.h"
#include <math.h>
#include "_GeneralModuleFiles/sys_model.h"
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h> 
#include <cmath>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>


using namespace std;
using Record = std::vector<std::string>;
using Records = std::vector<Record>;

struct OpnavMessageStruct {
    double pixel_out [ 2000 ];
    double line_out [ 2000 ];
    double mag_out [ 2000 ];       
};


class OpnavCamera {
public:
    OpnavCamera();
    ~OpnavCamera();
    double deg2rad(double degrees);                       //!< [-] Right Ascension of Stars taken from database
    int updateState();                                  //!< [-] Right Ascension of Stars taken from database
    double camera_frame2focal_plane(double c2,double c3);

    vector<double>       RA_stars_only;                             //!< [-] Right Ascension of Stars taken from database
    vector<double>       DE_stars_only;                             //!< [-] Declination of Stars taken from database
    vector<double>       n1_stars_only;                             //!< [-] Declination of Stars taken from database
    vector<double>       n2_stars_only;                             //!< [-] Declination of Stars taken from database
    vector<double>       n3_stars_only;                             //!< [-] Declination of Stars taken from database
    vector<double>       vismag_stars_only;
    vector<double>       RA;                             //!< [-] Right Ascension of Stars taken from database
    vector<double>       DE;                             //!< [-] Declination of Stars taken from database
    vector<double>       n1;                             //!< [-] Inertial coordinate 1 for star unit vectors
    vector<double>       n2;                             //!< [-] Inertial coordinate 2 for star unit vectors
    vector<double>       n3;                             //!< [-] Inertial coordinate 3 for star unit vectors
    vector<double>       c1;                             //!< [-] Camera Frame coordinate 1 for star unit vectors
    vector<double>       c2;                             //!< [-] Camera Frame coordinate 2 for star unit vectors
    vector<double>       c3;                             //!< [-] Camera Frame coordinate 3 for star unit vectors
    vector<double>       all_mags;
    double               RA_float;                       //!< [-] Right Ascension of Stars taken from database
    double               DE_float;                       //!< [-] Right Ascension of Stars taken from database
    double               vismag_float;                    //!< [-] Right Ascension of Stars taken from database
    double               alpha;                          //!< [-] Right Ascension of Stars taken from database
    double               beta;                           //!< [-] Right Ascension of Stars taken from database
    double               gamma;                          //!< [-] Right Ascension of Stars taken from database
    double               alpha_detector;
    double               beta_detector;
    double               a;
    double               b;
    double               c;
    double               f;
    double               gamma_detector;
    double               c2_min;   //!< [-] Right Ascension of Stars taken from database
    double               c2_max;   //!< [-] Right Ascension of Stars taken from database  
    double               c3_min;    //!< [-] Right Ascension of Stars taken from database 
    double               c3_max;    //!< [-] Right Ascension of Stars taken from database 

    std::map<char,double>  FOV_dimensions;
    std::map<char,double>  calculate_FOV(double f, double a, double b);
    double alpha_resolution;
    double beta_resolution;
    double pixel_out [ 2000 ];
    double line_out [ 2000 ];
    double mag_out [ 2000 ];  
    struct OpnavMessageStruct OpnavMessage;

private:
    double               phi;                            //!< [-] Right Ascension of Stars taken from database
    double               theta;                          //!< [-] Right Ascension of Stars taken from database 
    double               vismag_tmp;
    double               line_temp;
    double               pixel_temp;
    double               n1_temp;                        //!< [-] Right Ascension of Stars taken from database  
    double               n2_temp;                        //!< [-] Right Ascension of Stars taken from database
    double               n3_temp;                        //!< [-] Right Ascension of Stars taken from database  
    double               c1_temp;                        //!< [-] Right Ascension of Stars taken from database
    double               c2_temp;                        //!< [-] Right Ascension of Stars taken from database
    double               c3_temp;                        //!< [-] Right Ascension of Stars taken from database
    vector<double>       vismag;                          //!< [-] Right Ascension of Stars taken from database
    vector<double>       line;
    vector<double>       pixel;

};




#endif
