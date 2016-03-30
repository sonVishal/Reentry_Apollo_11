# Right side of the ODE system for the auxiliary boundary value problem
# with the given control 'u' during reentry for the problem by Stoer / Bulirsch 1973

# Parts of the function created on Maple

function aux_f(t,x,dx)
    if dir == "forward"
        aux_f_(t,x,dx);
    elseif dir == "backward"
    	tmp = zeros(dx);
        aux_f_(1-t,x,tmp);
        dx[:] = -tmp;
    end

    return nothing
end

function aux_f_(t,x,dx)

	p = x[4:5];
	T = x[6];

    # Begin parts produced by Maple

	t5 = atan((-1000 * p[2] + 1000 * t));
	u = -p[1] * t5;
	t7 = Sm * rho0;
	t9 = x[3];
	t11 = exp(-beta * R * t9);
	t12 = x[1];
	t13 = t12 ^ 2;
	t17 = cos(u);
	t23 = x[2];
	t24 = sin(t23);
	t26 = 0.1e1 + t9;
	t27 = t26 ^ 2;
	t28 = 0.1e1 / t27;
	dx[1] = T * (-t7 * t11 * t13 * (c[1] - c[2] * t17) / 0.2e1 - g * t24 * t28);
	t34 = sin(u);
	t38 = cos(t23);
	t40 = 0.1e1 / R;
	dx[2] = T * (t7 * t11 * t12 * c[3] * t34 / 0.2e1 + t12 * t38 * t40 / t26 - g * t38 / t12 * t28);
	dx[3] = T * t12 * t24 * t40;
	dx[4] = 0.0e0;
	dx[5] = 0.0e0;
	dx[6] = 0.0e0;

    return nothing
end
