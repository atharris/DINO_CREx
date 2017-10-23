function rot1_mat = rot1( angle )

rot1_mat = [1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)];

end