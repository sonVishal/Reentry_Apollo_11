# Calculation of starting values for the Lagrange multipliers and
# the control during reentry
# by Stoer / Bulirsch 1973

# Parts of the function created on Maple

function multiplier_start(t,x)
    p = x[4:5];
    lambda = zeros(x[4:6]);

    # Begin parts produced by Maple
    t5 = atan((-1000 * p[2] + 1000 * t));
  	u = -p[1] * t5;
  	lambda[1] = -1;
  	t7 = tan(u);
  	t8 = x[1];
  	t10 = c[2];
  	t11 = c[3];
  	t12 = 0.1e1 / t11;
  	lambda[2] = -(t7 * t8 * t10 * t12);
  	t15 = t8 ^ 2;
  	t18 = x[3];
  	t20 = exp(-beta * R * t18);
  	t21 = rho0 * t20;
  	t22 = sqrt(t21);
  	t23 = t15 * t8 * t22;
  	t24 = t11 * R;
  	t25 = cos(u);
  	t26 = t24 * t25;
  	t29 = t23 * t11;
  	t30 = R * t18;
  	t34 = t18 ^ 2;
  	t35 = R * t34;
  	t40 = t20 * t15;
  	t41 = Sm * rho0 * t40;
  	t42 = c[1];
  	t55 = x[2];
  	t56 = sin(t55);
  	t60 = t10 * Sm;
  	t64 = t60 * t21;
  	t65 = t15 * t11;
  	t71 = sin(u);
  	t72 = t71 * t10;
  	t73 = cos(t55);
  	t74 = t15 * t73;
  	t84 = 0.20e2 * t23 * t26 + 0.40e2 * t29 * t30 * t25 + 0.20e2 * t29 * t35 * t25 + t41 * t24 * t42 * t25 + 0.2e1 * t41 * t24 * t42 * t18 * t25 + t41 * t24 * t42 * t34 * t25 + 0.2e1 * g * t56 * t26 - t60 * rho0 * t40 * t24 - 0.2e1 * t64 * t65 * t30 - t64 * t65 * t35 - 0.2e1 * t72 * t74 - 0.2e1 * t72 * t74 * t18 + 0.2e1 * t72 * g * t73 * R;
  	lambda[3] = -(t84 / t25 / t8 / t56 * t12 / (0.1e1 + 0.2e1 * t18 + t34) / 0.2e1);
    return lambda, u
end
