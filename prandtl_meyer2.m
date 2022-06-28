%% PRANDTL-MEYER SOLVER
%{
Written by: Liam Foley
Contact:
        s1605969@ed.ac.uk
        liam.foley@salachy.co.uk
Paper: Design and Simulation of Passive Control Surfaces on a Supersonic Sounding Rocket

Purpose: A solver to determine the value of the flow properties following
an expansion fan, using the Prandtl-Meyer relations as defined in Anderson
(1991), Modern Compressible Flow with Historical Perspective
Inputs: 
    thetdad = The angle which the fan is turning through (Degrees)
    Ma1 = The pre-shock Mach Number (Mach)
    p1 = The pre-shock pressure (Pa)
Outputs: 
    Ma2 = The post-shock Mach Number (Mach)
    p2 = The post-shock pressure (Pa)
%}

%% MAIN CODE
function [Ma2, p2] = prandtl_meyer2(thetad, Ma1, p1)

theta = thetad*pi/180; %Convert degrees to radians
gamma=1.4; 
nmax=100; %Define maximum allowed iterations of any loop
err=0.0001; %Define desired error
%% Bisection Method
t_bi=tic; %Start timer for bisection method
xlow=1.1; %Define upper and lower limits
xhigh=3;
Kx=(Ma1^2)-1; %Define M^2-1 for Prandtl-Meyer relation
v1=sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*Kx))-atan(sqrt(Kx));
v2 = v1 + theta; %Calculate Prandtl-Meyer relationship for initial value
funclow = v2 - testfan(xlow); %Calculate initial function values at limits
funchigh = v2 - testfan(xhigh);
if funclow*funchigh>0 
    error('xlow and xhigh are on the same side of a root');
end
bi_it=1;%Start iteration counter
bi_app_err=1;%Set initial value of error
while bi_app_err>=err && bi_it<nmax%While loop continues running until error is within acceptable limits
    bi_root_val = (xlow + xhigh)/2; %Calculate new midpoint
    midval = double(v2 - testfan(bi_root_val)); %Calculate value of function at midpoint
    if midval == 0.0 %If the value is 0, solution is solved
        e = 0;
        break 
    end
    if funclow*midval < 0 %If two points multiplied together are on opposite sides of the root
        xhigh=bi_root_val; %then when they're multiplied together, the value will be negative, and as such
    else                  %the higher value is moved lower.
        xlow=bi_root_val; %However, if this isn't true, then the low point and the midpoint are on the same side,
    end                   %and then the low value is moved to the midpoint to be closer to the root.
    bi_it=bi_it+1; %Iteration counter is stepped
    bi_app_err = (xhigh-xlow)/2; %New error calculated as half the difference between the low and high points.   
end

%% Result 
Ma2 = bi_root_val; %Calculate post-shock values of Mach and Pressure
c1=1+((gamma-1)/2)*Ma1^2;
c2=1+((gamma-1)/2)*Ma2^2;
p2=p1/((c2/c1)^(gamma/(gamma-1)));


    function vtest = testfan(MachTest) %Function to test the value determined by the bisection solver
        K = MachTest^2 -1;
        vtest = sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*K))-atan(sqrt(K));
    end
end