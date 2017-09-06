
//
//  main.cpp
//  matt_work
//
//  Created by Matt Muszynski on 4/13/17.
//  Copyright © 2017 Matt Muszynski. All rights reserved.
//

#include "opnav_camera.h"

//time_t start = time(0);

OpnavCamera::OpnavCamera () {
    this->alpha = deg2rad(85); //from BSK
    this->beta = deg2rad(0); //from BSK
    this->gamma = deg2rad(72); //from BSK
    this->f = 3; //user input
    this->a = 1.5; //user input
    this->b = 2.5; //user input
    this->alpha_resolution = 512; //user input
    this->beta_resolution = 1024; //user input
}

OpnavCamera::~OpnavCamera () {
    return;
}

std::map<char,double> OpnavCamera::calculate_FOV(double f, double a, double b) {

    //This script is documented in the DINO C-REx CoDR Presentation.
    c = sqrt(a*a + b*b);

    //angular distance of diagonal of FOV
    gamma_detector = atan2(c/2,f);
    
    alpha_detector = sqrt(gamma_detector*a*a/c/c);
    beta_detector = sqrt(gamma_detector*b*b/c/c);
    FOV_dimensions['g'] = gamma_detector;
    FOV_dimensions['a'] = alpha_detector;
    FOV_dimensions['b'] = beta_detector;

    return FOV_dimensions;
}
double OpnavCamera::deg2rad(double degrees) {
    double radians = degrees*2.*M_PI/360.;
    return radians;
}

double OpnavCamera::camera_frame2focal_plane(double c2, double c3) {
    return 1;
}


int OpnavCamera::updateState() {
    FOV_dimensions = OpnavCamera::calculate_FOV(this->f,this->a,this->b);
    alpha_detector = FOV_dimensions['a'];
    beta_detector = FOV_dimensions['b'];
    c2_min = sin(-beta_detector/2);
    c2_max = sin(beta_detector/2);
    c3_min = sin(-alpha_detector/2);
    c3_max = sin(alpha_detector/2);
    
    
    //Add Stars to RA and DE arrays

    if (RA_stars_only.size() == 0) { 
    //if you haven't already loaded the star catalog, read tycho_small.csv and create
    //RA_stars_only and DE_stars_only
    std::string line;
    std::ifstream myfile ("/Users/Matt/Basilisk/SimCode/sensors/opnav_camera/tycho_small.csv");
    
    
        if (myfile.is_open())  {
            while ( getline (myfile,line) ) {
                std::string buf; // Have a buffer string
                std::stringstream ss(line); // Insert the string into a stream
                vector<string> tokens; // Create vector to hold our words
                while (ss >> buf)
                    tokens.push_back(buf);
                RA_float = ::atof(tokens[0].c_str());
                DE_float = ::atof(tokens[1].c_str());
                vismag_float = ::atof(tokens[2].c_str());
                RA_stars_only.insert(RA_stars_only.end(), 1, RA_float);
                DE_stars_only.insert(DE_stars_only.end(), 1, DE_float);
                vismag_stars_only.insert(vismag_stars_only.end(), 1, vismag_float);

                phi = deg2rad(90.-DE_float);
                theta = deg2rad(RA_float);
                n1_temp = sin(phi)*cos(theta);
                n2_temp = sin(phi)*sin(theta);
                n3_temp = cos(phi);

                n1_stars_only.insert(n1_stars_only.end(), 1, n1_temp);
                n2_stars_only.insert(n2_stars_only.end(), 1, n2_temp);
                n3_stars_only.insert(n3_stars_only.end(), 1, n3_temp);
                }
            myfile.close();
        } else {
            cout << "Unable to open tycho_small.csv\n";
        }
    } 

    RA = RA_stars_only;
    DE = DE_stars_only;
    n1 = n1_stars_only;
    n2 = n2_stars_only;
    n3 = n3_stars_only;
    vismag = vismag_stars_only;

    //add beacons
    std::vector<double>  Beacon_RA = {85, 232.5, 187 };
    std::vector<double>  Beacon_DE = {0, 80, -59 };
    std::vector<double>  Beacon_Mag = {11, 12, 13 };

    for (int i=0; i<Beacon_RA.size(); ++i) {
            RA_float = Beacon_RA[i];
            DE_float = Beacon_DE[i];
            vismag_float = Beacon_Mag[i];

            RA.insert(RA.begin(), 1, RA_float);
            DE.insert(DE.begin(), 1, DE_float);
            vismag.insert(vismag.begin(), 1, vismag_float);

            phi = deg2rad(90.-DE_float);
            theta = deg2rad(RA_float);
            n1_temp = sin(phi)*cos(theta);
            n2_temp = sin(phi)*sin(theta);
            n3_temp = cos(phi);
            n1.insert(n1.begin(), 1, n1_temp);
            n2.insert(n2.begin(), 1, n2_temp);
            n3.insert(n3.begin(), 1, n3_temp);
    }

    

    //calculate which stars are in FOV and compute their camera frame coordinates
    for (int i=0; i<RA.size(); ++i) {
            n1_temp = n1[i];
            n2_temp = n2[i];
            n3_temp = n3[i];
            vismag_tmp = vismag[i];
            //Compute c3_temp
            c3_temp =
                n1_temp*(
                         cos(gamma)*sin(beta)*cos(alpha) +
                         sin(gamma)*sin(alpha)
                         ) +
                n2_temp*(
                         cos(gamma)*sin(beta)*sin(alpha) -
                         sin(gamma)*cos(alpha)
                         ) +
                n3_temp*(
                         cos(gamma)*cos(beta)
                         );
            
            //if c3_temp is in range, compute c2_temp and check if it is in range
            if ( c3_temp > c3_min && c3_temp < c3_max ) {
                //Compute c2_temp
                c2_temp =
                    n1_temp*(
                         sin(gamma)*sin(beta)*cos(alpha) -
                         cos(gamma)*sin(alpha)
                         ) +
                    n2_temp*(
                         sin(gamma)*sin(beta)*sin(alpha) +
                         cos(gamma)*cos(alpha)
                         ) +
                    n3_temp*(
                         sin(gamma)*cos(beta)
                         );

            
                //if c2_temp is in range, compute c1_temp and add all to arrays
                if ( c2_temp > c2_min && c2_temp < c2_max ) {
                    c1_temp =
                        n1_temp*cos(beta)*cos(alpha) +
                        n2_temp*cos(beta)*sin(alpha) -
                        n3_temp*sin(beta);
                    
                    if (c1_temp > 0) {
                        line_temp = alpha_resolution*(c3_temp +
                                    sin(alpha_detector/2))/
                                    (2*sin(alpha_detector/2));
                        pixel_temp = -beta_resolution*(c2_temp -
                                    sin(beta_detector/2))/
                                    (2*sin(beta_detector/2));
                        c1.insert(c1.end(),1,c1_temp);
                        c2.insert(c2.end(),1,c2_temp);
                        c3.insert(c3.end(),1,c3_temp);
                        vismag.insert(vismag.end(), 1, vismag_float);
                        line.insert(line.end(),1, line_temp);
                        pixel.insert(pixel.end(),1, pixel_temp);
                        all_mags.insert(all_mags.end(),1,vismag_tmp);
                    }
                }
            }
        }

    for (int i=0; i<2000; ++i) {
        if (i < pixel.size()) {
            OpnavMessage.pixel_out[i] = pixel[i];
            OpnavMessage.line_out[i] = line[i];
            OpnavMessage.mag_out[i] = all_mags[i];
        } else {
            OpnavMessage.pixel_out[i] = NAN;
            OpnavMessage.line_out[i] = NAN;
            OpnavMessage.mag_out[i] = NAN;
        }
        
        
    }
    return 0;
}



int main()
    {
        OpnavCamera this_is_a_camera;
        //updateState();
        
        this_is_a_camera.updateState();


    }


    /*typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::milliseconds milliseconds;
    Clock::time_point t0 = Clock::now();
    time_t end2 = time(0);
    double elapsed2 = difftime(end2, start) * 1000.0;
    std::cout << elapsed2 << std::endl;
    t0 = Clock::now();*/
