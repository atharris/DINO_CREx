#! /usr/bin/env python3
import pdb
from numpy import array, arange, pi, average, log10, exp, sqrt
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'dependencies')
from scipy.optimize import curve_fit
from em import planck
from util import interpolate_lambda_dependent
import camera
from constants import r_sun, ly
from lmfit import minimize, Minimizer, Parameters, Parameter,report_fit

a = '3511.4 0.000 4494.7 3532.8 0.007 4525.9 3554.2 0.015 4557.2 3575.5 0.029 4588.5 3596.9 0.050 4619.7 3618.3 0.074 4651.0 3639.7 0.102 4682.3 3661.1 0.134 4713.6 3682.5 0.169 4744.8 3703.9 0.208 4776.1 3725.3 0.249 4807.4 3746.7 0.292 4838.6 3768.1 0.338 4869.9 3789.4 0.387 4901.2 3810.8 0.436 4932.4 3832.2 0.486 4963.7 3853.6 0.536 4995.0 3875.0 0.584 5026.3 3896.4 0.630 5057.5 3917.8 0.673 5088.8 3939.2 0.713 5120.1 3960.6 0.748 5151.3 3982.0 0.779 5182.6 4003.3 0.806 5213.9 4024.7 0.828 5245.2 4046.1 0.846 5276.4 4067.5 0.862 5307.7 4088.9 0.874 5339.0 4110.3 0.885 5370.2 4131.7 0.896 5401.5 4153.1 0.905 5432.8 4174.5 0.915 5464.0 4195.9 0.926 5495.3 4217.2 0.937 5526.6 4238.6 0.948 5557.9 4260.0 0.959 5589.1 4281.4 0.970 5620.4 4302.8 0.980 5651.7 4324.2 0.989 5682.9 4345.6 0.996 5714.2 4367.0 1.000 5745.5 4388.4 0.996 5776.8 4409.8 0.983 5808.0 4431.2 0.957 5839.3 4452.5 0.916 5870.6 4473.9 0.862 5901.8 4495.3 0.797 5933.1 4516.7 0.724 5964.4 4538.1 0.645 5995.7 4559.5 0.565 6026.9 4580.9 0.485 6058.2 4602.3 0.411 6089.5 4623.7 0.343 6120.7 4645.1 0.283 6152.0 4666.4 0.236 6183.3 4687.8 0.195 6214.5 4709.2 0.164 6245.8 4730.6 0.139 6277.1 4752.0 0.120 6308.4 4773.4 0.105 6339.6 4794.8 0.094 6370.9 4816.2 0.085 6402.2 4837.6 0.076 6433.4 4859.0 0.068 6464.7 4880.3 0.059 6496.0 4901.7 0.050 6527.3 4923.1 0.042 6558.5 4944.5 0.033 6589.8 4965.9 0.025 6621.1 4987.3 0.017 6652.3 5008.7 0.011 6683.6 5030.1 0.006 6714.9 5051.5 0.003 6746.2 5072.9 0.001 6777.4'
a = array([float(i) for i in a.split()])

b = '0.000 0.014 0.042 0.094 0.175 0.280 0.399 0.523 0.638 0.734 0.807 0.863 0.901 0.930 0.953 0.970 0.983 0.993 0.998 1.000 0.998 0.992 0.982 0.970 0.955 0.938 0.920 0.901 0.881 0.859 0.837 0.813 0.789 0.764 0.739 0.713 0.687 0.659 0.631 0.602 0.573 0.544 0.515 0.486 0.457 0.427 0.398 0.370 0.343 0.317 0.291 0.267 0.245 0.225 0.207 0.189 0.172 0.157 0.143 0.131 0.120 0.110 0.100 0.091 0.082 0.073 0.065 0.058 0.050 0.042 0.033 0.024 0.016 0.008'
b = array([float(i) for i in b.split()])

c = '2942.0 0.000 3566.0 0.000 4718.2 0.000 5945.3 2959.1 0.012 3591.6 0.005 4758.6 0.017 5990.7 2976.2 0.024 3617.3 0.010 4799.1 0.041 6036.0 2993.3 0.037 3642.9 0.016 4839.5 0.078 6081.3 3010.4 0.051 3668.6 0.027 4880.0 0.134 6126.7 3027.5 0.066 3694.2 0.044 4920.4 0.211 6172.0 3044.6 0.082 3719.9 0.064 4960.8 0.311 6217.3 3061.7 0.098 3745.5 0.092 5001.3 0.429 6262.6 3078.8 0.115 3771.2 0.148 5041.7 0.553 6308.0 3095.9 0.133 3796.8 0.232 5082.2 0.675 6353.3 3113.0 0.152 3822.5 0.325 5122.6 0.788 6398.6 3130.1 0.172 3848.1 0.419 5163.0 0.879 6444.0 3147.2 0.192 3873.8 0.510 5203.5 0.943 6489.3 3164.3 0.212 3899.4 0.599 5243.9 0.983 6534.6 3181.5 0.234 3925.0 0.688 5284.4 1.000 6580.0 3198.6 0.256 3950.7 0.763 5324.8 0.993 6625.3 3215.7 0.279 3976.3 0.810 5365.3 0.975 6670.6 3232.8 0.302 4002.0 0.839 5405.7 0.949 6716.0 3249.9 0.325 4027.6 0.869 5446.1 0.915 6761.3 3267.0 0.349 4053.3 0.895 5486.6 0.877 6806.6 3284.1 0.373 4078.9 0.916 5527.0 0.835 6851.9 3301.2 0.398 4104.6 0.936 5567.5 0.790 6897.3 3318.3 0.423 4130.2 0.954 5607.9 0.742 6942.6 3335.4 0.449 4155.9 0.967 5648.4 0.694 6987.9 3352.5 0.475 4181.5 0.977 5688.8 0.648 7033.3 3369.6 0.501 4207.2 0.985 5729.2 0.605 7078.6 3386.7 0.528 4232.8 0.992 5769.7 0.565 7123.9 3403.8 0.556 4258.5 0.996 5810.1 0.525 7169.3 3420.9 0.583 4284.1 0.999 5850.6 0.486 7214.6 3438.0 0.611 4309.8 1.000 5891.0 0.446 7259.9 3455.1 0.639 4335.4 0.998 5931.4 0.407 7305.3 3472.2 0.667 4361.1 0.991 5971.9 0.370 7350.6 3489.3 0.695 4386.7 0.981 6012.3 0.334 7395.9 3506.4 0.722 4412.3 0.970 6052.8 0.301 7441.2 3523.5 0.749 4438.0 0.954 6093.2 0.270 7486.6 3540.6 0.775 4463.6 0.933 6133.7 0.240 7531.9 3557.8 0.800 4489.3 0.910 6174.1 0.214 7577.2 3574.9 0.824 4514.9 0.887 6214.5 0.189 7622.6 3592.0 0.847 4540.6 0.865 6255.0 0.167 7667.9 3609.1 0.870 4566.2 0.843 6295.4 0.146 7713.2 3626.2 0.892 4591.9 0.821 6335.9 0.128 7758.6 3643.3 0.913 4617.5 0.797 6376.3 0.111 7803.9 3660.4 0.933 4643.2 0.772 6416.7 0.096 7849.2 3677.5 0.952 4668.8 0.746 6457.2 0.082 7894.5 3694.6 0.970 4694.5 0.719 6497.6 0.069 7939.9 3711.7 0.985 4720.1 0.691 6538.1 0.058 7985.2 3728.8 0.996 4745.8 0.660 6578.5 0.047 8030.5 3745.9 1.000 4771.4 0.628 6619.0 0.038 8075.9 3763.0 0.997 4797.1 0.595 6659.4 0.030 8121.2 3780.1 0.985 4822.7 0.564 6699.8 0.024 8166.5 3797.2 0.962 4848.3 0.533 6740.3 0.019 8211.9 3814.3 0.929 4874.0 0.502 6780.7 0.015 8257.2 3831.4 0.885 4899.6 0.471 6821.2 0.012 8302.5 3848.5 0.830 4925.3 0.440 6861.6 0.010 8347.9 3865.6 0.767 4950.9 0.410 6902.1 0.008 8393.2 3882.7 0.697 4976.6 0.378 6942.5 0.007 8438.5 3899.8 0.625 5002.2 0.347 6982.9 0.006 8483.8 3917.0 0.552 5027.9 0.316 7023.4 0.006 8529.2 3934.1 0.481 5053.5 0.284 7063.8 0.005 8574.5 3951.2 0.413 5079.2 0.253 7104.3 0.004 8619.8 3968.3 0.349 5104.8 0.226 7144.7 0.004 8665.2 3985.4 0.290 5130.5 0.204 7185.1 0.003 8710.5 4002.5 0.238 5156.1 0.182 7225.6 0.003 8755.8 4019.6 0.191 5181.8 0.160 7266.0 0.003 8801.2 4036.7 0.151 5207.4 0.138 7306.5 0.002 8846.5 4053.8 0.117 5233.1 0.116 7346.9 0.002 8891.8 4070.9 0.089 5258.7 0.094 7387.4 0.002 8937.2 4088.0 0.066 5284.4 0.074 7427.8 0.001 8982.5 4105.1 0.047 5310.0 0.057 7468.2 0.001 9027.8 4122.2 0.033 5335.6 0.041 7508.7 0.001 9073.1 4139.3 0.022 5361.3 0.026 7549.1 0.001 9118.5 4156.4 0.013 5386.9 0.016 7589.6 0.000 9163.8 4173.5 0.008 5412.6 0.010 7630.0 0.000 9209.1 4190.6 0.003 5438.2 0.005 7670.4 0.000 9254.5'
c = array([float(i) for i in c.split()])

d = '5507.0 0.000 7028.9 5553.0 0.106 7058.8 5598.9 0.241 7088.7 5644.8 0.431 7118.6 5690.7 0.631 7148.5 5736.6 0.769 7178.4 5782.5 0.853 7208.3 5828.4 0.910 7238.2 5874.4 0.948 7268.1 5920.3 0.975 7297.9 5966.2 0.989 7327.8 6012.1 0.997 7357.7 6058.0 1.000 7387.6 6103.9 0.999 7417.5 6149.9 0.997 7447.4 6195.8 0.992 7477.3 6241.7 0.986 7507.2 6287.6 0.979 7537.1 6333.5 0.970 7567.0 6379.4 0.960 7596.9 6425.3 0.948 7626.7 6471.3 0.934 7656.6 6517.2 0.918 7686.5 6563.1 0.900 7716.4 6609.0 0.884 7746.3 6654.9 0.867 7776.2 6700.8 0.849 7806.1 6746.7 0.828 7836.0 6792.7 0.805 7865.9 6838.6 0.780 7895.8 6884.5 0.755 7925.7 6930.4 0.730 7955.6 6976.3 0.704 7985.4 7022.2 0.678 8015.3 7068.2 0.652 8045.2 7114.1 0.626 8075.1 7160.0 0.601 8105.0 7205.9 0.576 8134.9 7251.8 0.549 8164.8 7297.7 0.521 8194.7 7343.6 0.492 8224.6 7389.6 0.463 8254.5 7435.5 0.434 8284.4 7481.4 0.407 8314.3 7527.3 0.379 8344.1 7573.2 0.352 8374.0 7619.1 0.325 8403.9 7665.0 0.299 8433.8 7711.0 0.274 8463.7 7756.9 0.251 8493.6 7802.8 0.229 8523.5 7848.7 0.208 8553.4 7894.6 0.188 8583.3 7940.5 0.167 8613.2 7986.4 0.148 8643.1 8032.4 0.130 8673.0 8078.3 0.114 8702.8 8124.2 0.101 8732.7 8170.1 0.090 8762.6 8216.0 0.080 8792.5 8261.9 0.071 8822.4 8307.9 0.063 8852.3 8353.8 0.055 8882.2 8399.7 0.047 8912.1 8445.6 0.040 8942.0 8491.5 0.033 8971.9 8537.4 0.027 9001.8 8583.3 0.021 9031.7 8629.3 0.017 9061.5 8675.2 0.014 9091.4 8721.1 0.012 9121.3 8767.0 0.009 9151.2 8812.9 0.006 9181.1 8858.8 0.003 9211.0'
d = array([float(i) for i in d.split()])

e = '0.000 0.025 0.051 0.083 0.138 0.210 0.285 0.366 0.452 0.540 0.621 0.691 0.754 0.812 0.858 0.890 0.919 0.944 0.962 0.976 0.989 0.995 0.998 0.999 1.000 0.999 0.998 0.997 0.996 0.995 0.995 0.995 0.997 0.998 0.999 0.998 0.997 0.995 0.990 0.985 0.979 0.974 0.967 0.961 0.952 0.942 0.930 0.918 0.905 0.892 0.877 0.859 0.837 0.815 0.784 0.743 0.698 0.649 0.593 0.530 0.467 0.404 0.343 0.282 0.226 0.178 0.134 0.094 0.064 0.043 0.024 0.010 0.005 0.003'
e = array([float(i) for i in e.split()])
B_T_lambda = a[0::3]/10
S_B_T = a[1::3]
V_T_lambda = a[2::3]/10
S_V_T = b

V_J_lambda = c[4::7]/10
S_V_J = c[5::7]
B_J_lambda = c[2::7]/10
S_B_J = c[3::7]

T = 9602

lambda_set_nm = arange(100,3010,10)*1e-9
bb = planck(T,lambda_set_nm)

bb_S_V_T = planck(T,V_T_lambda*1e-9)*S_V_T
dlambda_S_V_T = average(V_T_lambda[1:]-V_T_lambda[0:-1])

bb_S_B_T = planck(T,B_T_lambda*1e-9)*S_B_T
dlambda_S_B_T = average(B_T_lambda[1:]-B_T_lambda[0:-1])

bb_S_V_J = planck(T,V_J_lambda*1e-9)*S_V_J
dlambda_S_V_J = average(V_J_lambda[1:]-V_J_lambda[0:-1])

bb_S_B_J = planck(T,B_J_lambda*1e-9)*S_B_J
dlambda_S_B_J = average(B_J_lambda[1:]-B_J_lambda[0:-1])

F_V_T_vega = sum(bb_S_V_T)*dlambda_S_V_T
F_B_T_vega = sum(bb_S_B_T)*dlambda_S_B_T
F_V_J_vega = sum(bb_S_V_J)*dlambda_S_V_J
F_B_J_vega = sum(bb_S_B_J)*dlambda_S_B_J

if 1:
	plt.plot(V_T_lambda,planck(T,V_T_lambda*1e-9))
	plt.plot(B_T_lambda,planck(T,B_T_lambda*1e-9))
	plt.plot(B_T_lambda,bb_S_B_T)
	plt.plot(V_T_lambda,bb_S_V_T)
	plt.show()

pdb.set_trace()
BVJ = []
BVT = []
T_arr = []
from numpy import logspace 
for T in logspace(3,4):
	bb_S_V_T = planck(T,V_T_lambda*1e-9)*S_V_T
	bb_S_B_T = planck(T,B_T_lambda*1e-9)*S_B_T
	bb_S_V_J = planck(T,V_J_lambda*1e-9)*S_V_J
	bb_S_B_J = planck(T,B_J_lambda*1e-9)*S_B_J

	F_V_T = sum(bb_S_V_T)*dlambda_S_V_T
	F_B_T = sum(bb_S_B_T)*dlambda_S_B_T
	F_V_J = sum(bb_S_V_J)*dlambda_S_V_J
	F_B_J = sum(bb_S_B_J)*dlambda_S_B_J

	BVJ.append(-2.5*log10((F_B_J*F_V_J_vega)/(F_V_J*F_B_J_vega)))
	BVT.append(-2.5*log10((F_B_T*F_V_T_vega)/(F_V_T*F_B_T_vega)))
	T_arr.append(float(T))
	if 0:
		plt.plot(V_T_lambda,planck(T,V_T_lambda*1e-9),color='black',label='Blackbody Curve')
		plt.plot(B_T_lambda,planck(T,B_T_lambda*1e-9),color='black')
		plt.plot(B_T_lambda,bb_S_B_T,color='blue',label='Tcyho Blue')
		plt.plot(V_T_lambda,bb_S_V_T,color='green',label='Tcyho Visible')
		plt.legend
		plt.show()

BVJ = array(BVJ)
BVT = array(BVT)
T_arr = array(T_arr)

def T_model(BV,a,b,c,f):
	return a*(f + c*BV)**-b

def BV_model(a,b,c,T):
	return a/(T + b) + c

# plt.plot(BVJ,T_arr, label='Johnson')
plt.plot(BVT,T_arr)
plt.plot(BVT,T_arr)
plt.xlabel('Color Index')
plt.ylabel('Temperature (K)')
plt.title('Color Index and Temperature Relation')
plt.legend()


# define objective function: returns the array to be minimized
def fcn2min(params, BV, T):
    a = params['a']
    b = params['b']
    c = params['c']
    f = params['f']
    model = T_model(BV,a,b,c,f)
    return model - T

# create a set of Parameters
params = Parameters()
params.add('a', value= 1,min=0)
params.add('b', value= 1, min=0)
params.add('c', value= 1, max=1)
params.add('f', value= 1, min=-min(BVJ))

minner = Minimizer(fcn2min, params, fcn_args=(BVJ, T_arr))
result = minner.minimize()
pdb.set_trace()
final = T_arr + result.residual
report_fit(result)

a_J = result.params['a'].value
b_J = result.params['b'].value
c_J = result.params['c'].value
f_J = result.params['f'].value

# try to plot results
plt.figure()
plt.plot(BVJ, T_arr, 'k.',label='Derived Blackbody Relation')
# plt.plot(BVJ, final, 'r')
plt.plot(BVJ, T_model(BVJ,a_J,b_J,c_J,f_J),label='Johnson LM Fit')
plt.title('Emprical Fit for (B-V)_ J and T')
plt.xlabel('Johnson Color Index')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()

params = Parameters()
params.add('a', value= 1,min=0)
params.add('b', value= 1, min=0)
params.add('c', value= 1, max=1)
params.add('f', value= 1, min=-min(BVT))

minner = Minimizer(fcn2min, params, fcn_args=(BVT, T_arr))
result = minner.minimize()
pdb.set_trace()
final = T_arr + result.residual
report_fit(result)

a_T = result.params['a'].value
b_T = result.params['b'].value
c_T = result.params['c'].value
f_T = result.params['f'].value


# try to plot results
plt.figure()
plt.plot(BVT, T_arr, 'k.',label='Derived Blackbody Relation')
plt.plot(BVT, final, 'r')
plt.plot(BVT, T_model(BVT,a_T,b_T,c_T,f_T))
plt.plot(BVT, T_model(BVT,a_T,b_T,c_T,f_T),label='Tycho LM Fit')
plt.title('Emprical Fit for (B-V)_ T and T')
plt.xlabel('Tycho Color Index')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()


pdb.set_trace()

import sqlite3
db = 'db/tycho.db'
conn = sqlite3.connect(db)
c = conn.cursor()
d = conn.cursor()

#COMMENT OUT THE UPDATE SINCE IT'S BEEN DONE ALREADY
if 1:
	select_string = "SELECT id, BTmag, VTmag, published_temperature from tycho_data"
	c.execute(select_string)
	for row in c:
		if row[1] != None and row[2] != None:
			BVT = row[1] - row[2]
			T = T_model(BVT,a_T,b_T,c_T,f_T)
			update_string = "UPDATE tycho_data SET computed_temperature=\'" + \
			str(T) + "\' where id=\'" + str(row[0]) + "'"
			d.execute(update_string)
			print(update_string)
		elif row[3] != None:
			update_string = "UPDATE tycho_data SET computed_temperature=\'" + \
			str(row[3]) + "\' where id=\'" + str(row[0]) + "'"
			d.execute(update_string)
			print(update_string)
		else:
			BVT = 0
			T = T_model(BVT,a_T,b_T,c_T,f_T)
			update_string = "UPDATE tycho_data SET computed_temperature=\'" + \
			str(T) + "\' where id=\'" + str(row[0]) + "'"
			d.execute(update_string)
			print(update_string)
	conn.commit()


select_string = "SELECT name, computed_temperature, published_temperature from tycho_data"
c.execute(select_string)

plt.figure()
plt.plot([3000,13000],[3000,13000],'--',color='black',label="Perfect Correlation")
computed_temps = []
published_temps = []
for row in c:
	if row[2] != None and row[1] != '':
		print(row)
		computed_temperature = row[1]
		published_temperature = row[2]
		computed_temps.append(computed_temperature)
		published_temps.append(published_temperature)
		err = abs(published_temperature-computed_temperature)/published_temperature
		if err > 0.2: 
			color = 'red'
		elif err > 0.1:
			color = 'orange'
		else:
			color = 'green'
		plt.plot(computed_temperature,published_temperature,'.',color=color)
plt.xlim([2500,14000])
plt.ylim([2500,58000])

plt.xlabel('Computed Temperature')
plt.ylabel('Published Temperature')
plt.title('All Stars Studied')
plt.plot([-100,-100],'.g',label='Error < 10%')
plt.plot([-100,-100],'.',color='orange',label='10% < Error < 20%')
plt.plot([-100,-100],'.r',label='Error > 20%')
plt.legend()


select_string = "SELECT name, computed_temperature, published_temperature from tycho_data"
c.execute(select_string)

plt.figure()
plt.plot([3000,13000],[3000,13000],'--',color='black',label="Perfect Correlation")

for row in c:
	if row[2] != None and row[1] != '':
		print(row)
		computed_temperature = row[1]
		published_temperature = row[2]
		err = abs(published_temperature-computed_temperature)/published_temperature
		if err > 0.2: 
			color = 'red'
		elif err > 0.1:
			color = 'orange'
		else:
			color = 'green'
		plt.plot(computed_temperature,published_temperature,'.',color=color)

plt.xlabel('Computed Temperature')
plt.ylabel('Published Temperature')
plt.xlim([2500,12500])
plt.ylim([2500,12500])
plt.title('Linear Regime: 2500K to 12500K')
plt.plot([-100,-100],'.g',label='Error < 10%')
plt.plot([-100,-100],'.',color='orange',label='10% < Error < 20%')
plt.plot([-100,-100],'.r',label='Error > 20%')
plt.legend()
plt.show()
