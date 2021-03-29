%% SOLVE FIN
%{
Written by: Liam Foley
Contact:
        s1605969@ed.ac.uk
        liam.foley@salachy.co.uk
Paper: Design and Simulation of Passive Control Surfaces on a Supersonic Sounding Rocket

Purpose: Use the Busemann Second Order Airfoil calculation method to
calculate the properties of a supersonic rocket fin. Code is based almost
verbatim on the program FIN, created in FORTRAN by Barrowman (1966) for the
same purpose, with updates to make it usable in MATLAB.
Inputs: 
    FinDims, Flowstate = Inputs from RocketSetup.m defining the shape of
    the body and the nature of the flowfield.    
A single variable to check which of the rocket setups to use.
Outputs: 
    FinCD = Value of the drag coefficient of the fin, based on the
    axial force
    BodCNA = Value of the Normal Force Derivative of the fin, for use in
    determining stability
    BoDCPX = Value of the Center of Pressure of the fin, reported as
    a number of length units behind the fin tip, current setup uses
    millimeters
    FinTotDrag = Total value of the drag force on the rocket at 0 degrees,
    used in comparing axial force between the analytical and simulated
    cases

NOTE: This code is KNOWN to be innacurate in determining the Center of
Pressure, as the distribution is significantly different from that
calculated by Barrowman in his own paper. Other properties of the fins are
calculated with good agreement, but use with caution. 

%}

%% MAIN CODE
function [FinCD, FinCNA, FinCPX,FinTotDrag] = SolveFin(FinDims, Flowstate)

%% Flow Properties (Input from setup script)
M = Flowstate.M;
T0 = Flowstate.T;
p0 = Flowstate.p;
rho0 = Flowstate.rho;

%% Fin Dimensions (Input from setup script)
CR = FinDims.CR; 
LLR = FinDims.LLR;
LTR = FinDims.LTR;
SweepL = FinDims.SweepL; %leading edge sweep angle in deg
SweepT = FinDims.SweepT; %trailing edge sweep angle in deg
Sweep1 = FinDims.Sweep1; %Region 1 boundary sweep angle  in deg
Sweep2 = FinDims.Sweep2; %Region 2 boundary sweep angle  in deg
ZETAL = FinDims.ZETAL; %leading edge half angle in deg
ZETAT = FinDims.ZETAT; %Trailing edge half angle in deg
SPAN = FinDims.SPAN;
REFA = FinDims.REFA;
REFL = FinDims.REFL;


%% Solver Setup Steps
n = 20; %number of strips along the fin
DelY = SPAN/n; %length of strip in m
Beta = sqrt(M^2 - 1);

for a = [0 1]
ALPHA = a; %angle of attack
ALPHAR = ALPHA*pi/180;

%% FIN ANGLE VALUES
LAM = [SweepL Sweep1 Sweep2 SweepT];
LAR = LAM .* pi/180;
TLAM = sin(LAR)./cos(LAR);
ZETALR = ZETAL*pi/180;
ZETATR = ZETAT*pi/180;
CZL = cos(ZETALR);
CZT = cos(ZETATR);

%% STRIP WIDTH
XNS = n;
DS = SPAN/XNS;
XNS = 2*DS/REFA;

%% MACH NUMBER
MSQ = M^2;
BETA = sqrt(MSQ-1);
DEM = -2*MSQ+4/3;
DEN = 1/BETA^7;
K1 = 2/BETA;
K2 = (1.2*MSQ^2 - 2*BETA^2)/BETA^4;
K3 = (0.4*MSQ^4 - (10.88/6)*MSQ^3 + 4*MSQ*MSQ+DEM)*DEN;
B3 = (0.04 * MSQ^4 - .32*MSQ^3 + 0.4 * MSQ*MSQ)*DEN;

%% MACH ANGLE
CMU = (TLAM(1)+BETA);
MU = SPAN * CMU;

TLIFT = 0;
TDRAG = 0;
TMMOM = 0;
TLMOM = 0;
%% CL VS ALPHA
S = DS/2;
while S < SPAN
    C = S*(TLAM(4)-TLAM(1))+CR;
    LL = S*(TLAM(2)-TLAM(1))+LLR;
    LT = S*(TLAM(4)-TLAM(3))+LTR;
    SL=LL/CZL;
    ST=LT/CZT;
    SM = C-LL-LT;
    XL = S*TLAM(1)*cos(ALPHAR);
    ETA(1) = ZETAL - ALPHA;
    ETA(2) = -ALPHA;
    ETA(3) = -ZETAT - ALPHA;
    ETA(4) = ZETAL+ALPHA;
    ETA(5) = ALPHA;
    ETA(6) = -ZETAT + ALPHA;
    ETAR = ETA.*pi/180;
    D(1) = SL* cos(ETAR(1));
    D(2) = SM* cos(ETAR(2));
    D(3) = ST* cos(ETAR(3));
    D(4) = SL* cos(ETAR(4));
    D(5) = SM* cos(ETAR(5));
    D(6) = ST* cos(ETAR(6));
    N(1) = SL* sin(ETAR(1));
    N(2) = SM* sin(ETAR(2));
    N(3) = ST* sin(ETAR(3));
    N(4) = SL* sin(ETAR(4));
    N(5) = SM* sin(ETAR(5));
    N(6) = ST* sin(ETAR(6));       
    for I=1:6
        SAVE = 1.05*(K1*ETAR(I)+K2*ETAR(I)^2+K3*ETAR(I)^3);
        if I >=4
            R(I)=SAVE - B3*ETAR(4)^3;
        elseif ETAR(1) > 0
            R(I)=SAVE- B3*ETAR(1)^3;
        else
            R(I)=SAVE;
        end
    end
%% FIN TIP MACH CONE CORRECTION
    LW=MU-S*CMU;
    SL = C-LT;
    if LW >= C
        PL= [0 0 0 0 0 0];
    elseif LW >= SL
        PL(1) = 0;
        PL(2) = 0;
        PL(3) = (C-LW)/LT;
        PL(6) = PL(3);
        PL(4) = 0;
        PL(5) = 0;
    elseif LW >= LL
        PL(1) = 0;
        PL(2) = (C-LT-LW)/SM;
        PL(3) = 1;
        PL(4) = 0;
        PL(5) = PL(2);
        PL(6) = PL(3);
    else
        PL(1) = (LL-LW)/LL;
        PL(4) = PL(1);
        PL(2) = 1;
        PL(3) = 1;
        PL(5) = 1;
        PL(6) = 1;
    end
%% STRIP COEFFICIENTS
    for I=1:6
        if any(I==1:3)
            FL(I) = -R(I)*D(I)*(1-0.5*PL(I));
        else
            FL(I) = R(I)*D(I)*(1-0.5*PL(I));
        end
        
        FD(I) = R(I)*N(I)*(1-0.5*PL(I));
    end
    XP(1) = XL;
    XP(2) = XL + D(1);
    XP(3) = XL + D(1) + D(2);
    XP(4) = XL;
    XP(5) = XL + D(4);
    XP(6) = XL + D(4) + D(5);
    
    MMOMT = 0;
    for I=1:6
        if D(I) ~= 0
            CW = 0.5*D(I)*(1-PL(I)+0.5*PL(I)^2+XP(I)*(2-PL(I))/D(I))/(1-0.5*PL(I));
        else
            CW = 1;
        end
        MMOMT = MMOMT - CW*FL(I);
    end
    LIFT = 0;
    DRAG = 0;
    
    for I=1:6
        LIFT=LIFT+FL(I);
        DRAG = DRAG+FD(I);
    end
    
    LMOM = S*LIFT;
    TLIFT=TLIFT+LIFT;
    TDRAG = TDRAG+DRAG;
    TMMOM = TMMOM+MMOMT;
    TLMOM=TLMOM+LMOM;
    S=S+DS;  
end

%% TOTAL COEFFICIENTS
count = a+1;
CL(count) = TLIFT*XNS;
CD(count) = TDRAG*XNS;
CM(count) = TMMOM*XNS/REFL;
CLM(count) = 2*TLMOM*XNS/REFL;
CPX(count) = TMMOM/(TLIFT*cos(ALPHAR));
CPS(count) = TLMOM/TLIFT;
end

FinCD = CD(1);
FinCNA = (CL(2)/(pi/180));
FinCPX = CPX(2);
FinTotDrag = CD(1)/XNS;