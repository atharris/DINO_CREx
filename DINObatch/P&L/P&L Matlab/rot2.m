function rot2_mat = rot2( angle )

rot2_mat = [cos(angle),0,-sin(angle);0,1,0;sin(angle),0,cos(angle)];

end