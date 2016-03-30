function utility(x)

    lambda = zeros(x[4:6]);
    lambda[1] = x[4];
    lambda[2] = x[5];
    lambda[3] = x[6];

    # Begin parts produced by Maple
    t1 = lambda[2];
  	t2 = c[3];
  	t4 = lambda[1];
  	t6 = x[1];
  	t7 = 0.1e1 / t6;
  	t9 = c[2];
  	u = atan(t1 * t2 / t4 * t7 / t9);
  	t13 = t6 ^ 2;
  	t16 = x[3];
  	t18 = exp(-beta * R * t16);
  	t20 = sqrt(rho0 * t18);
  	t22 = 0.10e2 * t13 * t6 * t20;
  	t23 = Sm * rho0;
  	t26 = cos(u);
  	t32 = x[2];
  	t33 = sin(t32);
  	t35 = 0.1e1 + t16;
  	t36 = t35 ^ 2;
  	t37 = 0.1e1 / t36;
  	t43 = sin(u);
  	t47 = cos(t32);
  	t49 = 0.1e1 / R;
  	h = t22 + t4 * (-t23 * t18 * t13 * (c[1] - t9 * t26) / 0.2e1 - g * t33 * t37) + t1 * (t23 * t18 * t6 * t2 * t43 / 0.2e1 + t6 * t47 * t49 / t35 - g * t47 * t7 * t37) + lambda[3] * t6 * t33 * t49;
  	phi = t22;
    return phi,h,u
end
