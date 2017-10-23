function MM_deriv_z_input = MM_deriv_z(Dx,Dy,FoL,Kx,Ky,RA,dec,twist,x_ob,x_sc,y_ob,y_sc,z_ob,z_sc)
%MM_DERIV_Z
%    MM_DERIV_Z_INPUT = MM_DERIV_Z(DX,DY,FOL,KX,KY,RA,DEC,TWIST,X_OB,X_SC,Y_OB,Y_SC,Z_OB,Z_SC)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    23-Oct-2017 15:52:32

t2 = conj(x_ob);
t3 = conj(x_sc);
t4 = t2-t3;
t5 = abs(t4);
t12 = conj(y_ob);
t13 = conj(y_sc);
t14 = t12-t13;
t6 = abs(t14);
t8 = conj(z_ob);
t9 = conj(z_sc);
t10 = t8-t9;
t7 = abs(t10);
t11 = t5.^2;
t15 = t6.^2;
t16 = t7.^2;
t17 = t11+t15+t16;
t18 = 1.0./sqrt(t17);
t19 = cos(dec);
t20 = sin(twist);
t21 = sin(dec);
t22 = cos(twist);
t23 = sin(RA);
t24 = cos(RA);
t25 = sign(t10);
t26 = t20.*t23;
t27 = t21.*t22.*t24;
t28 = t26+t27;
t29 = 1.0./t17.^(3.0./2.0);
t30 = t22.*t23;
t31 = t30-t20.*t21.*t24;
t32 = t4.*t18.*t21;
t33 = t10.*t18.*t19.*t22;
t34 = t14.*t18.*t19.*t20;
t35 = t32+t33+t34;
t36 = 1.0./t35;
t37 = t20.*t24;
t42 = t21.*t22.*t23;
t38 = t37-t42;
t39 = t22.*t24;
t40 = t20.*t21.*t23;
t41 = t39+t40;
t43 = 1.0./t35.^2;
t44 = t4.*t7.*t21.*t25.*t29;
t45 = t7.*t10.*t19.*t22.*t25.*t29;
t46 = t7.*t14.*t19.*t20.*t25.*t29;
t47 = t44+t45+t46-t18.*t19.*t22;
MM_deriv_z_input = [Dx.*FoL.*Kx.*t36.*(t18.*t28-t7.*t10.*t25.*t28.*t29+t7.*t14.*t25.*t29.*t31+t4.*t7.*t19.*t24.*t25.*t29)-Dx.*FoL.*Kx.*t43.*t47.*(-t10.*t18.*t28+t14.*t18.*t31+t4.*t18.*t19.*t24);Dy.*FoL.*Ky.*t36.*(t18.*t38-t7.*t10.*t25.*t29.*t38+t7.*t14.*t25.*t29.*t41-t4.*t7.*t19.*t23.*t25.*t29)+Dy.*FoL.*Ky.*t43.*t47.*(t10.*t18.*t38-t14.*t18.*t41+t4.*t18.*t19.*t23)];
