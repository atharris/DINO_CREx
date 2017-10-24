function A_hat_I_deriv_x_input = A_hat_I_deriv_x(x_ob,x_sc,y_ob,y_sc,z_ob,z_sc)
%A_HAT_I_DERIV_X
%    A_HAT_I_DERIV_X_INPUT = A_HAT_I_DERIV_X(X_OB,X_SC,Y_OB,Y_SC,Z_OB,Z_SC)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    23-Oct-2017 16:57:06

t5 = conj(x_ob);
t6 = conj(x_sc);
t7 = t5-t6;
t2 = abs(t7);
t9 = conj(y_ob);
t10 = conj(y_sc);
t11 = t9-t10;
t3 = abs(t11);
t13 = conj(z_ob);
t14 = conj(z_sc);
t15 = t13-t14;
t4 = abs(t15);
t8 = t2.^2;
t12 = t3.^2;
t16 = t4.^2;
t17 = t8+t12+t16;
t18 = sign(t7);
t19 = 1.0./t17.^(3.0./2.0);
A_hat_I_deriv_x_input = [-1.0./sqrt(t17)+t2.*t7.*t18.*t19;t2.*t11.*t18.*t19;t2.*t15.*t18.*t19];
