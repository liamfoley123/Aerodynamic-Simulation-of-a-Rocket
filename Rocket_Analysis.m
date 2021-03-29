%% ROCKET ANALYSIS
%{
Written by: Liam Foley
Contact:
        s1605969@ed.ac.uk
        liam.foley@salachy.co.uk
Paper: Design and Simulation of Passive Control Surfaces on a Supersonic Sounding Rocket

Purpose: An overall function which brings together the other elements to
determine a final value of the properties of interest for the rocket,
notably a value for the location of the Center of Pressure, along with a
value for the total axial force on the rocket. 
Inputs: None.
Outputs: 
    Node = A struct which contains the values at every location along the nose
    and body
    X_Tot = Analytically calculated distance of the Center of Pressure from the
    nose tip (m)
    TotalRocketDrag = Total Axial drag on the rocket (N)
    CNA_Tot = Total value of the Neutral Coefficient Derivative

J.Prac = Practical Methods for Calculating
J.Prog = Program for Calculating
%}

%% MAIN CODE
clearvars -except RocketBodyPressure

%% Fitting Functions for conical shock data
CNaCFit=CNaC; %Syntax: CNaFit(Angle, Mach Number)
MaCFit = MaC;
CpCFit = CpC;

%% Definition of Rocket
Rocket = [3];

%% Call setup function to define variables of various parts of the rocket
[Flowstate, BodDims, FinDims] = RocketSetup(Rocket);


%% Call solver functions
% Outputs = [Drag Coefficient, Normal Coefficient Derivative, Center of
% Pressure distance from front]

[FinCD, FinCNA,FinCPX,FinTotDrag] = SolveFin(FinDims, Flowstate);
[BodCD, BodCNA, BodCPX, Node] = SolveNose(CNaCFit, MaCFit, CpCFit,BodDims, Flowstate);

% Define Beta
Beta = sqrt(Flowstate.M^2 - 1);

%Define Tau for interference drag coefficient calculations
tau = (FinDims.SPAN + BodDims.D/2)/(BodDims.D/2);
invtau = 1/tau;

% Interference factor: Tail in presence of body (J.Prac 3-96)
Ktb = (2/(pi*(1-invtau)^2)) * ( (1 + invtau^4) * (0.5 * atan(0.5 * (tau-invtau)))...
    - invtau^2 * ( (tau-invtau) + 2*atan(invtau) ));

% Center of pressure independent of interference flow, find Center of
% pressure from nose tip
Xtb = BodDims.L_T + abs(FinCPX);

% Interference factor: body in presence of tail (J.Prac 3-98)
Kbt = ( (1 - invtau^2)^2 / (1-invtau)^2  ) - Ktb;

% Body center of pressure based on approx fit of J.Prac Figure 3-11
BDbyCR = (Beta*BodDims.D)/FinDims.CR;
if BDbyCR >=1
    XbyCR = 0.67;
else
    XbyCR = 1.1*0.17*BDbyCR + 0.5;
end
Xbt = FinDims.CR * XbyCR + BodDims.L_T;

% Calculation of CNa values
CNAtb = (BodDims.NFin/2) * FinCNA * Ktb; 
CNAbt = BodCNA*Kbt;

% Calculation of total CNa (J.Prac 3-106)
CNA_Tot = BodCNA + CNAtb + CNAbt;

% Calculation of total Center of Pressure (J.Prac 3-107)
X_Tot = (abs(BodCPX)*BodCNA + Xtb*CNAtb + Xbt * CNAbt)/CNA_Tot;

%% Skin Friction Drag
% Find rocket Reynolds Number
Re = (Flowstate.rho*Flowstate.U*BodDims.L_Rocket)/Flowstate.mu;

% Assume flow is turbulent and use Hana's formula for skin frinction
% coefficient (J.Prac 4-4)
Cf = 1/((3.46*log(Re)-5.6)^2);

% For supersonic compressible flow, modify skin friction according to a
% curve fit of Figure 4-3 (J.Prac 4-12)
k=0.1;
Cfc = Cf/( (1+k * Flowstate.M^2)^0.58);

% Find drag coefficient of the fins (J.Prac 4-15)
CDft = 2*BodDims.NFin*Cfc * (FinDims.A_planform/FinDims.REFA);

% Find skin friction drag coefficient of body (J.Prac 4-16)
CDfb = (1+(0.5/BodDims.fB))*Cfc * (BodDims.Wet_Area/FinDims.REFA);

%% Pressure Drag
DeltaCD = 1.214 - (0.502/Flowstate.M^2)+(0.1095/Flowstate.M^4)...
    +(0.0231/Flowstate.M^6); % (J.Prac 4-19)
CDlt = 2 * BodDims.NFin * ((FinDims.SPAN*FinDims.RL)/(FinDims.REFA))...
    * cos(FinDims.ZETAL)^2 * DeltaCD; % (J.Prac 4-22)


CDBt = (BodDims.NFin * (1 - 0.52 * Flowstate.M^-1.19) * (FinDims.A_Base/FinDims.REFA))...
    /((1-0.18 * Cfc * (FinDims.Troot/FinDims.HRoot)^2)*Flowstate.M^2); % (J.Prac 4-30)

%% Overall Tail Drag
CDtt = CDft + CDfb + CDlt + CDBt; % (J.Prac 4-38)

%% Body Drag
%BodCD = BodCD/(Flowstate.q_inf*BodDims.Nose_Wet_Area); %Find CD based on actual 
% body pressure drag from nose solver

CDbStar = (BodDims.Base_Area/FinDims.REFA) * (0.185 + 1.15*(FinDims.HRoot/FinDims.CR)); %(J.Prac 4-51)

Mcr = (0.892)/(sqrt(CDbStar)); % J.Prac 4-56
if Flowstate.M<Mcr
    CDb = CDbStar * (0.88 + 0.12*exp(-3.58*(Flowstate.M-1))); %J.Prac 4-54
else
    CDb = (0.7*(BodDims.Base_Area/FinDims.REFA))/Flowstate.M^2; %J.Prac 4-55
end

CDTB = CDb + BodCD; %J.Prac 4-57

%% Total Drag
Rocket_CD = CDTB + CDtt; %J.Prac 4-58
TotalRocketDrag = BodCD + FinTotDrag * 3;
%% Plot Data
X = [0, Node.nodepos];
P = [Flowstate.p, Node.pTan];

RocketBodyPressure = readtable('Rocket Body Pressure.csv', 'VariableNamingRule', 'preserve');

RocketBodyPressure2 = readtable('Rocket Body Pressure Run2.csv', 'VariableNamingRule', 'preserve');
ExperimentalBodyPressure = readtable('Experimental Pressure Distribution.csv', 'VariableNamingRule', 'preserve');
RocketCrossSection = readtable('Rocket Cross Section.csv', 'VariableNamingRule', 'preserve');
OpenrocketCoP = 1.91;
SimmedCoP = 1.6704;
CoMStart = 1.65;
CoMBurnout = 1.13;

figure(1)
plot(X/1000, P-Flowstate.p, 'LineWidth', 3);
hold on
grid on
plot(RocketBodyPressure{6:665,1}, RocketBodyPressure{6:665,2}, 'kx','Markersize', 6)
title('Comparison of analytical and simulated pressure distribution')
xlabel('Distance from rocket tip (m)')
ylabel('Relative pressure (Pa)')
legend('Analytical Results','Simulated Data')
set(gcf,'color','w');
xlim([0 3])
hold off

figure(2)
plot(RocketBodyPressure{6:633,1}, RocketBodyPressure{6:633,2}, 'k-')
hold on
plot(RocketBodyPressure2{7:662,1}, RocketBodyPressure2{7:662,2}, 'b-')
grid on
title('Comparison of pressure distribution at varying refinement levels')
xlabel('Distance from rocket tip (m)')
ylabel('Relative pressure (Pa)')
legend('Initial Results','Further Refinement Results')
set(gcf,'color','w');
xlim([0 0.5])
hold off

figure(3)
plot(RocketBodyPressure{6:633,1}/0.1524, (RocketBodyPressure{6:633,2}+26325)/26325)
hold on
grid on
plot(ExperimentalBodyPressure{:,1}, ExperimentalBodyPressure{:,4},'kx', 'MarkerSize', 15)
title('Comparison of simulated and experimental data')
xlabel('Body Calibers from nose tip (X/D)')
ylabel('Wall Pressure/Freestream Pressure')
legend('Simulated Results (Mach 1.8, Nose Fineness Ratio 3.6)','Experimental Results (Mach 2.0, Nose Fineness Ratio 3.0')
set(gcf,'color','w');
xlim([0 6])
hold off

figure(4)
plot(RocketCrossSection{:,1}, RocketCrossSection{:,2},'r-')
hold on
plot(RocketCrossSection{:,1}, -RocketCrossSection{:,2},'r-')
grid on
plot(OpenrocketCoP, 0, 'kx', 'MarkerSize', 15)
plot(SimmedCoP, 0, 'bx', 'MarkerSize', 15)
plot(X_Tot/1000, 0, 'rx', 'MarkerSize', 15)
plot(CoMStart, 0, 'gx', 'MarkerSize', 15)
plot(CoMBurnout, 0, 'cx', 'MarkerSize', 15)
set(gcf,'color','w');
hold off