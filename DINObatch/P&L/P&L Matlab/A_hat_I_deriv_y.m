function A_hat_I_deriv_y_input = A_hat_I_deriv_y(x_ob,x_sc,y_ob,y_sc,z_ob,z_sc)
%A_HAT_I_DERIV_Y
%    A_HAT_I_DERIV_Y_INPUT = A_HAT_I_DERIV_Y(X_OB,X_SC,Y_OB,Y_SC,Z_OB,Z_SC)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    23-Oct-2017 15:52:28

t2 = conj(y_ob);
t3 = conj(y_sc);
t4 = t2-t3;
t5 = conj(x_ob);
t6 = conj(x_sc);
t7 = t5-t6;
t8 = abs(t7);
t9 = abs(t4);
t13 = conj(z_ob);
t14 = conj(z_sc);
t15 = t13-t14;
t10 = abs(t15);
t11 = t8.^2;
t12 = t9.^2;
t16 = t10.^2;
t17 = t11+t12+t16;
t18 = sign(t4);
t19 = 1.0./t17.^(3.0./2.0);
A_hat_I_deriv_y_input = [t7.*t9.*t18.*t19;-1.0./sqrt(t17)+t4.*t9.*t18.*t19;t9.*t15.*t18.*t19];
