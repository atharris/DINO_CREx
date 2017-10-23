
% inertial position vectors of spacecraft and object
state_spacecraft = [1000,0,0]';
state_object     = [1200,1000,400]';

%% camera parameters

% resolution
resolution = [2048,512];

% focal length constant (millimeters)
FoL = 100;

% pixel and line center off set
p0 = resolution(1) / 2;
l0 = resolution(2) / 2;

% pixel width
pixel_width = 2.5;

% pixel height
pixel_height = 10;

% x to pixel constant
Kx = 1/pixel_width;

% y to line constant
Ky = 1/pixel_height;

% x to pixel direction
Dx = 1;

% y to line direction
Dy = 1;

%% camera angles

% define right ascension (RA)
RA = 0;

% define declination (dec)
dec = 0;

% define camera twist (twist)
twist = pi/2;

%%

% equation 5.9
% get the apparent direction vector
spacecraft_to_object = state_object - state_spacecraft;

% equation 4.1
% create rotation matrix for inertial to camera
T_TV_I = rot3( RA ) * rot2( dec ) * rot1( twist );

% equation 5.13
% create inertial unit vector and rotate to camera
A_hat_I  = spacecraft_to_object / norm( spacecraft_to_object );
A_hat_TV = T_TV_I * A_hat_I;

% equation 5.14
% convert from camera unit vector to field of view millimeters
x = FoL / A_hat_TV( 3 ) * A_hat_TV(1);
y = FoL / A_hat_TV( 3 ) * A_hat_TV(2);
% x = FoL  * A_hat_TV(1);
% y = FoL  * A_hat_TV(2);

% equation 5.17
% convert from millimeter camera view to pixel and line
p = Kx * Dx * x + p0;% + pixel_width/2;
l = Ky * Dy * y + l0;% + pixel_height/2;

disp( p )
disp( l )

