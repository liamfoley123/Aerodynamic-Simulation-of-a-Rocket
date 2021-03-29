%% ROCKET SETUP
%{
Written by: Liam Foley
Contact:
        s1605969@ed.ac.uk
        liam.foley@salachy.co.uk
Paper: Design and Simulation of Passive Control Surfaces on a Supersonic Sounding Rocket

Purpose: 
Inputs: 
    Rocket = A single variable to check which of the rocket setups to use.
Outputs: 
    Flowstate, BodDims, FinDims = Structs which contain the variables
    of the rocket geometry, along with the flow properties of the rocket.
%}

%% MAIN CODE
function [Flowstate, BodDims, FinDims] = RocketSetup(Rocket)

in2mm = 1; %0.0254;
Mach2ms = 300;

if Rocket == [1] %Aerobee 550
    Flowstate.M = 3; %freestream mach Number
    Flowstate.U = Flowstate.M * Mach2ms; %freestream velocity (m/s)
    Flowstate.p = 100000; %freestream Pressure (Pa)
    Flowstate.rho = 1.225; %freestream density (kg/m3)
    Flowstate.T = 300; %freestream temperature (K)
    Flowstate.mu = 1.81e-5; %freestream dynamic viscosity (Pa/s)
    Flowstate.q_inf = 0.5 * Flowstate.rho*Flowstate.U^2;


    %% Fin Dimensions
    FinDims.CR = 44.25 * in2mm;
    FinDims.LLR = 15 * in2mm;
    FinDims.LTR = 0;
    FinDims.SweepL = 45; %leading edge sweep angle in deg
    FinDims.SweepT = 15; %trailing edge sweep angle in deg
    FinDims.Sweep1 = 45; %Region 1 boundary sweep angle  in deg
    FinDims.Sweep2 = 15; %Region 2 boundary sweep angle  in deg
    FinDims.ZETAL = 1.5; %leading edge half angle in deg
    FinDims.ZETAT = 0; %Trailing edge half angle in deg
    FinDims.RL = 0.25; %Leading edge radius
    FinDims.SPAN = 34.2 * in2mm; %Fin span
    FinDims.REFA = 381 * in2mm^2; %Fin Reference Area
    FinDims.REFL = 22 * in2mm; %Fin Reference Length
    FinDims.CT = FinDims.CR - FinDims.SPAN * tan(FinDims.SweepL * pi/180) - ...
        FinDims.SPAN * tan(FinDims.SweepT * pi/180);%tip chord
    FinDims.A_planform = 0.5 * (FinDims.CR + FinDims.CT)*FinDims.SPAN; %planform fin area
    FinDims.Troot = tand(FinDims.ZETAL) * FinDims.LLR;
    FinDims.A_Base = FinDims.SPAN * FinDims.Troot;
    FinDims.HRoot = FinDims.Troot;

    %% Nose/Body Dimensions
    BodDims.D = 22;
    BodDims.L_Rocket = 469.62;
    BodDims.L_Nose = 103.86;
    BodDims.L_T = BodDims.L_Rocket - FinDims.CR;
    BodDims.NFin = 4;
    BodDims.Base_Area = pi*(BodDims.D/2)^2;
    R_Ogive = ((BodDims.D/2)^2 + BodDims.L_Nose^2)/(BodDims.D);
    BodDims.Nose_Wet_Area = 2*pi*R_Ogive *( ( BodDims.D/2-R_Ogive) * asin(BodDims.L_Nose/R_Ogive)+BodDims.L_Nose);
    Bod_Wet_Area = pi*(BodDims.D/2)^2 * (BodDims.L_Rocket - BodDims.L_Nose);
    BodDims.Wet_Area = BodDims.Nose_Wet_Area + Bod_Wet_Area;
    BodDims.fB = BodDims.L_Rocket/BodDims.D;
    BodDims.BoatTailAngle = 0;
    BodDims.L_Boattail = 0;
    
elseif Rocket == [3] %Darwin II
    
    Flowstate.M = 1.8; %freestream mach Number
    Flowstate.U = Flowstate.M * Mach2ms; %freestream velocity (m/s)
    Flowstate.p = 26500; %freestream Pressure (Pa)
    Flowstate.rho = 0.414; %freestream density (kg/m3)
    Flowstate.T = 223.3; %freestream temperature (K)
    Flowstate.mu = 1.81e-5; %freestream dynamic viscosity (Pa/s)
    Flowstate.q_inf = 0.5 * Flowstate.rho*Flowstate.U^2;


    %% Fin Dimensions
    FinDims.CR = 400 ;
    FinDims.LLR = 25 ;
    FinDims.LTR = 0;
    FinDims.SweepL = 30; %leading edge sweep angle in deg
    FinDims.SweepT = 0; %trailing edge sweep angle in deg
    FinDims.Sweep1 = FinDims.SweepL; %Region 1 boundary sweep angle  in deg
    FinDims.Sweep2 = FinDims.SweepT; %Region 2 boundary sweep angle  in deg
    FinDims.ZETAL = 4.3; %leading edge half angle in deg
    FinDims.ZETAT = 0; %Trailing edge half angle in deg
    FinDims.RL = 0; %Leading edge radius
    FinDims.SPAN = 199 ; %Fin span
    FinDims.REFL = 220 ; %Fin Reference Length
    FinDims.CT = FinDims.CR - FinDims.SPAN * tan(FinDims.SweepL * pi/180) - ...
    FinDims.SPAN * tan(FinDims.SweepT * pi/180);%tip chord
    FinDims.A_planform = 0.5 * (FinDims.CR + FinDims.CT)*FinDims.SPAN; %planform fin area
    FinDims.REFA = FinDims.A_planform ; %Fin Reference Area
    FinDims.Troot = tand(FinDims.ZETAL) * FinDims.LLR;
    FinDims.REFA = FinDims.A_planform; %Fin Reference Area    
    FinDims.REFL = (2/3) * (FinDims.CR + FinDims.CT - FinDims.CR*FinDims.CT)/(FinDims.CR + FinDims.CT);
    FinDims.A_Base = FinDims.SPAN * FinDims.Troot; 
    FinDims.HRoot = FinDims.Troot;

    %% Nose/Body Dimensions
    BodDims.D = 152.4; %Diameter 
    BodDims.L_Rocket = 2870; %Length (Including Nose)
    BodDims.L_Nose = 550; %Length of nose
    BodDims.L_T = BodDims.L_Rocket - FinDims.CR; %Length of body influenced by tail
    BodDims.NFin = 3; %Number of fins
    BodDims.Base_Area = pi*(BodDims.D/2)^2; %Area of base (assumes no aft cone)
    R_Ogive = ((BodDims.D/2)^2 + BodDims.L_Nose^2)/(BodDims.D); %Ogive radius of nose
    BodDims.Nose_Wet_Area = 2*pi*R_Ogive *( ( BodDims.D/2-R_Ogive) ...
            * asin(BodDims.L_Nose/R_Ogive)+BodDims.L_Nose); %Wetted area of nose
    Bod_Wet_Area = pi*(BodDims.D/2)^2 * (BodDims.L_Rocket - BodDims.L_Nose); %Wetted area of body%
    BodDims.Wet_Area = BodDims.Nose_Wet_Area + Bod_Wet_Area; %Total wetted area
    BodDims.fB = BodDims.L_Rocket/BodDims.D; %Fineness ratio of body
    BodDims.BoatTailAngle = 0; %Unused variable to define an aft cone/boattail
    BodDims.L_Boattail = 0;  %Unused variable to define an aft cone/boattail
    
end

