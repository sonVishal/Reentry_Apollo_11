# Solution of the reentry problem for the Apollo capsule
# by Stoer / Bulirsch 1973

# Import required modules
using ODEInterface
@ODEInterface.import_huge

include("reentry_f.jl");
include("reentry_bc.jl");
include("aux_f.jl");
include("aux_bc.jl");
include("showSolution.jl");
include("multiplier_start.jl");

# Parameters
R    = 209.0;
beta = 4.26;
rho0 = 2.704e-3;
g    = 3.2172e-4;
Sm   = 53200.0;
c    = [1.174, 0.9, 0.6];

# Boundary Conditions
v0     = 0.36;
v1     = 0.27;
gamma0 = -8.1*pi/180;
gamma1 = 0.0;
h0     = 4/209;
h1     = 2.5/209;

# bvpsol IVP solver
odesolver = odex;

# Tolerance
tol = 1e-8;

function reentry(T0)

	# Plotting parameters
	global plotVar = Array{Gadfly.Plot}(8);

	# Labels and titles for plotting
	xAxis = "rel. maneuver time";
	yAxis = ["[10<sup>5</sup> ft/sec]","[Degrees]","[10<sup>5</sup> ft]",
		"[Degrees]","","","",""];
	titles = ["Velocity (ν)","Flight-path angle (γ)",
		"Height (ξ)","Control parameter (u)","adj. Var. λ<sub>ν</sub>",
		"adj. Var. λ<sub>γ</sub>","adj. Var. λ<sub>ξ</sub>",
		"Hamilton-func"];

	# Set up empty plots
	for i = 1:8
		plotVar[i] = plot(Guide.xlabel(xAxis,orientation=:horizontal),Guide.ylabel(yAxis[i]),
		Guide.title(titles[i]),Guide.xticks(ticks=[0:0.2:1;]));
	end

	println("§§ Solution with T0 = $T0 §§");

  J_opt = -Inf;
  flag = 1;

  p = [1.0, 0.5];

  # Forward Shooting
  println("** Solution of the auxiliary problem, 1st attempt: Forward Shooting...");
  global dir = "forward";
  th  = [0.0, 1.0];

  # Initial values of states = boundary values
  x = Matrix(6,2)

  x[1:3,1] = [v0; gamma0; h0];
  x[1:3,2] = [v1; gamma1; h1];

  # Starting values of adjoint variables for the control during the maneuver
  x[4:6,1] = [p[1]; p[2]; T0];
  x[4:6,2] = x[4:6,1];

  # Solver options
	# IVP Solver options
	ivpopt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
	OPT_RTOL => tol, OPT_ATOL => tol);

	# BVP solver options
	opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
	OPT_MAXSTEPS => 100, OPT_RTOL => tol, OPT_BVPCLASS => 3,
	OPT_SOLMETHOD => 0, OPT_IVPOPT => ivpopt);

  tic();
  (_,yh,retcode,_) = bvpsol(aux_f,aux_bc,th,x,odesolver,opt);
  elapsedTime = toq();

  println(@sprintf "    CPU-time Boundary Value Solver = %3.1f sec" elapsedTime);

  if retcode > 0
      println("    # Newton-iterations = $retcode")
  end

  p[1] = yh[4,1];   # Amplitude parameter for control u
  p[2] = yh[5,1];   # Switching parameter for control u
  T1   = yh[6,1];   # Maneuver time

  if retcode < 0
      println("-- Termination of the BVP solver with error code $retcode")
      flag = 0;
  elseif p[1] <= 0 || p[1] >= 1 || p[2] <= 0 || p[2] >= 1
      println("-- Solution is unusable: Incorrect control profile")
      flag = 0;
  end

  if flag == 0
      # Backward Shooting
      println("** Solution of the auxiliary problem, 2nd attempt: Backward Shooting...");
      dir = "backward";

      x = x[:,end:-1:1];

      # Solver options
			# IVP Solver options
			ivpopt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
			OPT_RTOL => tol, OPT_ATOL => tol);

			# BVP solver options
			opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
			OPT_MAXSTEPS => 100, OPT_RTOL => tol, OPT_BVPCLASS => 3,
			OPT_SOLMETHOD => 0, OPT_IVPOPT => ivpopt);

      tic();
      (_,yh,retcode,_) = bvpsol(aux_f,aux_bc,th,x,odesolver,opt);
      elapsedTime = toq();

      println(@sprintf "    CPU-time Boundary Value Solver = %3.1f sec" elapsedTime);

      if retcode > 0
          println("    # Newton-iterations = $retcode")
      end

      p[1] = yh[4,1];   # Amplitude parameter for the control
      p[2] = yh[5,1];   # Switching point for the control
      T1   = yh[6,1];   # Maneuver time

      if retcode < 0
          println("-- Termination of the BVP solver with error code $retcode")
          return nothing
      elseif (p[1] <= 0 || p[1] >= pi/2 || p[2] <= 0 || p[2] >= 1)
          println("-- Solution is unusable: Incorrect control profile")
          return nothing
      end
  end

  # Generating the initial values for the optimal control problem
  println("** Calculation of starting trajectories for the optimal control problem...")

  K = 4;
  t_msm = zeros(K);
  t_msm[:] = linspace(0,1,K);
  t = zeros(K);
  t[:] = t_msm;

  if dir == "backward"
      t[:] = 1-t_msm[end:-1:1];
  end

  # Solver options
	# IVP Solver options
	ivpopt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
	OPT_RTOL => tol, OPT_ATOL => tol);

  tic();
  (t_0,x,retcode,_) = odecall(odesolver,aux_f,t,yh[:,1],ivpopt);
  elapsedTime = toq();

  println(@sprintf "    CPU-time Boundary Value Solver = %3.1f sec" elapsedTime);

  x = x';
  if dir == "backward"
      t = 1-t[end:-1:1];
      x = x[:,end:-1:1];
  end

  # Restore separated linear boundary conditions
  if norm(x[1:3,1]-[v0; gamma0; h0])/norm([v0; gamma0; h0]) > 0.1 ||
  	norm(x[1:3,end]-[v1; gamma1; h1])/norm([v0; gamma0; h0]) > 0.1
      println("-- Solution of the auxiliary problem unusable: more than 10% deviation in the boundary values")
      return nothing
  else
      println("++ Solution of the auxiliary problem will be used as the starting trajectory")
      println(@sprintf "** Time of the auxiliary maneuver T = %5.2f sec" T1)
  end

  x[1:3,1] = [v0; gamma0; h0];
  x[1:3,end] = [v1; gamma1; h1];

  lambda = zeros(3,K);
  u = zeros(K);

  for i=1:K
      (lambda[:,i], u[i]) = multiplier_start(t[i],x[:,i]);
  end

  x = [x[1:3,:]; lambda; T1*ones(K)'];
  x[isinf(x)] = 0.0;

  # Solution of the optimal control problem with the multiple shooting method
  println("** Solution of the optimal control problem...")

  # Calling bvpsol
  # Solver options
	# IVP Solver options
	ivpopt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
	OPT_RTOL => tol, OPT_ATOL => tol);

	# BVP solver options
	opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
	OPT_MAXSTEPS => 100, OPT_RTOL => tol, OPT_BVPCLASS => 3,
	OPT_SOLMETHOD => 0, OPT_IVPOPT => ivpopt);

  tic();
  (_,y,retcode,_) = bvpsol(reentry_f,reentry_bc,t_msm,x,odesolver,opt);
  elapsedTime = toq();

  println(@sprintf "    CPU-time Boundary Value Solver = %3.1f sec" elapsedTime);

  if retcode > 0
      println("    # Newton-Iterations = $retcode");
      T = y[7,1];
      println(@sprintf "** Manuever time T = %5.2f" T);

      # Plot the optimal solution
      J_opt = showSolution(t_msm,y,T,false,false,colorant"#00FF00","");

      # Plot 0th iteration of the Multiple shooting method
      showSolution(t_0,x,T0,false,false,colorant"#FF0000","");

      # Plot solution of auxiliary problem
      J_aux = showSolution(th,yh,T1,true,true,colorant"#0000FF",string("T0_",string(convert(Int64,T0)),".pdf"));

      saving = (J_aux-J_opt)/J_aux*100;

      println(@sprintf "** Improved the cost functional by %3.1f%% compared to the auxiliary problem!" saving);
  else
      println("-- Program terminated with error code $retcode");
  end

  println("The optimal heating J = $J_opt");

  return nothing
end
