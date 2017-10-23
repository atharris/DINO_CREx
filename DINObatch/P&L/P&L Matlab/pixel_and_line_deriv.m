syms x_sc y_sc z_sc x_ob y_ob z_ob FoL RA dec twist Kx Ky Dx Dy p0 l0
% inertial position vectors of spacecraft and object
state_spacecraft = [x_sc, y_sc, z_sc]';
state_object     = [x_ob, y_ob, z_ob]';

% equation 5.9
% get the apparent direction vector
spacecraft_to_object = state_object - state_spacecraft;

% equation 4.1
% create rotation matrix for inertial to camera
T_TV_I = rot3( RA ) * rot2( dec ) * rot1( twist );

% equation 5.13
% create inertial unit vector and rotate to camera
A_hat_I  = spacecraft_to_object / norm( spacecraft_to_object );

A_hat_I_deriv_input = diff(A_hat_I, x_sc);

A_hat_TV = T_TV_I * A_hat_I;

% equation 5.14
% convert from camera unit vector to field of view millimeters
x = Kx * Dx * FoL / A_hat_TV( 3 ) * A_hat_TV(1) + p0;
y = Ky * Dy * FoL / A_hat_TV( 3 ) * A_hat_TV(2) + l0;

x_y = [x;y];

x_y_deriv_input = diff( x_y, x_sc );

A_hat_I_deriv = matlabFunction( A_hat_I_deriv_input,...
                               'file','A_hat_I_deriv');

x_y_deriv     = matlabFunction( x_y_deriv_input,...
                               'file','x_y_deriv');

%% numerical stuff

state_spacecraft = [1000,0,0]';
state_object     = [1200,1000,400]';

x_sc = state_spacecraft(1);
y_sc = state_spacecraft(2);
z_sc = state_spacecraft(3);

x_ob = state_object(1);
y_ob = state_object(2);
z_ob = state_object(3);

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
dec = pi/4;

% define camera twist (twist)
twist = pi/2;

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

%% derivatives

dAdX_sym   = A_hat_I_deriv(x_ob,x_sc,y_ob,y_sc,z_ob,z_sc);

dx_ydX_sym = x_y_deriv(Dx,Dy,FoL,Kx,Ky,RA,dec,twist,x_ob,x_sc,y_ob,y_sc,z_ob,z_sc);

A_hat_norm = sqrt( (x_ob-x_sc)^2 + (y_ob-y_sc)^2 + (z_ob-z_sc)^2 );

dAdX_alg   = [ (x_ob - x_sc)^2 / A_hat_norm^3 - 1 / A_hat_norm; ...
               (x_ob - x_sc) * (y_ob - y_sc) / A_hat_norm^3;...
               (x_ob - x_sc) * (z_ob - z_sc) / A_hat_norm^3 ];

dA_TVdX_alg =  [ T_TV_I(1,1) * dAdX_alg(1) + T_TV_I(1,2) * dAdX_alg(2) + T_TV_I(1,3) * dAdX_alg(3);...
                 T_TV_I(2,1) * dAdX_alg(1) + T_TV_I(2,2) * dAdX_alg(2) + T_TV_I(2,3) * dAdX_alg(3);...
                 T_TV_I(3,1) * dAdX_alg(1) + T_TV_I(3,2) * dAdX_alg(2) + T_TV_I(3,3) * dAdX_alg(3) ];

dx_ydX_alg  =   FoL * [Kx * Dx; Ky*Dy].*...
             [ -1 / A_hat_TV(3)^2 * dA_TVdX_alg(3) * A_hat_TV(1) + ...
                                dA_TVdX_alg(1) / A_hat_TV(3);...
               -1 / A_hat_TV(3)^2 * dA_TVdX_alg(3) * A_hat_TV(2) + ...
                                dA_TVdX_alg(2) / A_hat_TV(3);
                       ];

dx_ydX_sym - dx_ydX_alg 
dx_ydX_alg
