%% SOLVE NOSE
%{
Written by: Liam Foley
Contact:
        s1605969@ed.ac.uk
        liam.foley@salachy.co.uk
Paper: Design and Simulation of Passive Control Surfaces on a Supersonic Sounding Rocket

Purpose: Implement the equations proposed by Syvertson & Dennis (1957) to
carry out a Second Order Shock-Expansion analysis of a curvent, tangent
ogival nose-cone on a rocket. 
Inputs: 
    CNaCFit, MaCFit, CpCFit = Interpolated surfaces of experimental results
    for equivalent cones, to be used as part of the shock-expansion method.
    BodDims, Flowstate = Inputs from RocketSetup.m defining the shape of
    the body and the nature of the flowfield.
Outputs: 
    BodCD = Value of the drag coefficient of the body alone, based on the
    axial force and the forward area (NOTE: Current implementation reports
    *only* the value of the Axial Force in Newtons, the final result does
    NOT consider the frontal area, due to the needs of the report this code
    was written for)
    BodCNaB = Value of the Normal Force Derivative of the body, for use in
    determining stability
    BoDCPX = Value of the Center of Pressure of the body alone, reported as
    a number of length units behind the nose tip, current setup uses
    millimeters
    Node = A struct which contains the values at every location along the nose
    and body
%}

%% MAIN CODE
function [BodCD, BodCNaB, BodCPX, Node] = SolveNose(CNaCFit, MaCFit, CpCFit,BodDims, Flowstate)

Ma0 = Flowstate.M;
T0 = Flowstate.T;
p0 = Flowstate.p;
rho0 = Flowstate.rho;
gamma = 1.4;
q_inf = 0.5 * rho0*(Ma0*343)^2;
%% Nose Definition
D = BodDims.D;
L_Rocket = BodDims.L_Rocket;
L_Nose = BodDims.L_Nose;
L_Afterbody = L_Rocket - L_Nose;
L_Boattail = BodDims.L_Boattail;
Area = pi*(D/2)^2;

R_Ogive = ((D/2)^2 + L_Nose^2)/(D);


n = 500; %number of segments to use
d_x = L_Nose/n; %axial length of each segment (note, leads to inconsistent section length but with enough segments this is negligible)
n2 = round(L_Afterbody/d_x);

%% Main Loop
for i=1:n+n2 %For loop to define every Node
    %% Conical Section Definition
    Node(i).nodepos=i*d_x - d_x/2; %Calculates position of each node
    nodegrad = OgiveGrad(Node(i).nodepos, R_Ogive, L_Nose, D);
    Node(i).Angle = atan((d_x*nodegrad)/d_x)*(180/pi);
    Node(i).noderad = OgiveY(Node(i).nodepos, R_Ogive, L_Nose, D);
    nodestart = Node(i).nodepos - d_x/2;
    nodeend = Node(i).nodepos+d_x/2;

    if i == 1
        delta = Node(i).Angle;
    elseif i>n
        delta = Node(i-1).Angle ;
        Node(i).Angle=0;
        %Node(i).nodepos = L_Nose + L_Afternose/2
        Node(i).nodepos=i*d_x - d_x/2;       
        Node(i).noderad = D/2;
        nodestart = Node(i).nodepos-d_x/2;
        nodeend = Node(i).nodepos+d_x/2;
    else
        delta = Node(i-1).Angle - Node(i).Angle;
    end
    %% Equivalent Cone Definition (Related only to free stream Mach)
    MaC = MaCFit(Node(i).Angle, Ma0);
    CpC = CpCFit(Node(i).Angle, Ma0);
    pc = CpC * q_inf+p0;
    CNaC = CNaCFit(Node(i).Angle, Ma0);
    
    %% Main Solution Step
    if i == 1 %determine if Node is the initial surface (coneflow)
        
        Node(i).MaTan = MaC;
        Node(i).Ma3 = MaC;
        Node(i).pTan = pc;
        p1=pc;
        Node(i).p3 = pc;
        Node(i).pgrad3 = 0;
        Node(i).AlphaTan = (tand(delta)*CNaC);
        Node(i).Alpha3 = (tand(delta)*CNaC);
        Node(i).AlphaR = Node(i).AlphaTan*Node(i).noderad;
        Node(i).AlphaRX = Node(i).AlphaTan*Node(i).noderad*Node(i).nodepos;
        
    elseif i>1 && i<=n %For Nodes along the curved surface
        p1 = Node(i-1).p3;
        pgrad1 = Node(i-1).pgrad3;
        Ma1 = Node(i-1).Ma3;
        [Ma2, p2] = prandtl_meyer2(delta, Ma1, p1) ;
                
        Beta1 = (gamma*p1*Ma1^2)/(2*(Ma1^2-1));
        Beta2 = (gamma*p2*Ma2^2)/(2*(Ma2^2-1));
        
        Omega1 = (1/Ma1) * (  (1+((gamma-1)/2)*Ma1^2)  / ((gamma+1)/2))^((gamma+1)/(2*(gamma-1)));
        Omega2 = (1/Ma2) * (  (1+((gamma-1)/2)*Ma2^2)  / ((gamma+1)/2))^((gamma+1)/(2*(gamma-1)));
        
        pgrad2 = (Beta2/Node(i).noderad) * ( (Omega1/Omega2)*sind(Node(i).Angle)-sind(Node(i).Angle)) +...
            Beta2/Beta1 * Omega1/Omega2 * pgrad1;
        EtaTan = pgrad2 * (Node(i).nodepos - nodestart) / ( (pc-p2)* cosd(Node(i).Angle));
        Node(i).pTan = pc - (pc-p2)*exp(EtaTan);
        Node(i).MaTan = sqrt(((p2/Node(i).pTan)^((gamma-1)/gamma) * (1 + ((gamma-1)/2)*Ma2^2) -1)/((gamma-1)/2));
        
        Lambda1 = (2*gamma*p1) / (sin(2*asin(1/Ma1)));
        Lambda2 = (2*gamma*p2) / (sin(2*asin(1/Ma2)));
        
        Node(i).AlphaTan = (1-exp(EtaTan)) * tand(Node(i).Angle) * CNaC +...
            Lambda2/Lambda1*exp(EtaTan)*Node(i-1).Alpha3;
        
        
        Eta3 = 4*pgrad2 * (nodeend - nodestart) / ( (pc-p2)* cosd(Node(i).Angle));
        Node(i).p3 = (pc - (pc-p2)*exp(Eta3));
        Node(i).pgrad3 = ( (pc-Node(i).p3)/(pc - p2) )*pgrad2;
        Node(i).expT=exp(EtaTan);
        Node(i).Ma3 = sqrt(((p2/Node(i).p3)^((gamma-1)/gamma) * (1 + ((gamma-1)/2)*Ma2^2) -1)/((gamma-1)/2));
        Node(i).Alpha3 = (1-exp(Eta3)) * tand(Node(i).Angle) * CNaC +...
            Lambda2/Lambda1*exp(Eta3)*Node(i-1).Alpha3;
        Node(i).AlphaR = Node(i).AlphaTan*Node(i).noderad;
        Node(i).AlphaRX = Node(i).AlphaTan*Node(i).noderad*(Node(i).nodepos);
        Node(i).CPDelta =  ((Node(i).pTan)); %* sind(Node(i).Angle));
        Node(i).AxForce = Node(i).CPDelta * pi * ((Node(i).noderad/1000)^2 - (Node(i-1).noderad/1000)^2)* sind(Node(i).Angle);%((pi * Node(i).noderad^2)-(pi*Node(i-1).noderad^2));
    end
    
                
    
end

After_pc = p0;
Afterp1 = Node(n).p3;
Afterpgrad1 = Node(n).pgrad3;
AfterMa1 = Node(n).Ma3;
AfterMa2 = AfterMa1;
Afterp2 = Afterp1;

AfterBeta1 = (gamma*Afterp1*AfterMa1^2)/(2*(AfterMa1^2-1));
AfterBeta2 = (gamma*Afterp2*AfterMa2^2)/(2*(AfterMa2^2-1));

AfterOmega1 = (1/AfterMa1) * (  (1+((gamma-1)/2)*AfterMa1^2)  / ((gamma+1)/2))^((gamma+1)/(2*(gamma-1)));
AfterOmega2 = (1/AfterMa2) * (  (1+((gamma-1)/2)*AfterMa2^2)  / ((gamma+1)/2))^((gamma+1)/(2*(gamma-1)));

Afterpgrad2 = (AfterBeta2/(D/2)) * ( (AfterOmega1/AfterOmega2)*sind(Node(n).Angle)-sind(Node(n).Angle)) +...
AfterBeta2/AfterBeta1 * AfterOmega1/AfterOmega2 * Afterpgrad1;

AfterAlpha1 = Node(n).Alpha3;

for j = 1:n2
    Node(n+j).nodepos=L_Nose+j*d_x - d_x/2;
    
    
    
    EtaTan = 4.75*Afterpgrad2 * (Node(n+j).nodepos-L_Nose)/(After_pc-Afterp2)*cosd(Node(n-1).Angle);
    Node(n+j).pTan = After_pc - (After_pc - Afterp2)*exp(EtaTan);
    Node(n+j).MaTan = sqrt(((Afterp2/Node(n+j).pTan)^((gamma-1)/gamma) * (1 + ((gamma-1)/2)*AfterMa2^2) -1)/((gamma-1)/2));

    Lambda1 = (2*gamma*Afterp1) / (sin(2*asin(1/AfterMa1)));
    Lambda2 = (2*gamma*Afterp2) / (sin(2*asin(1/AfterMa2)));    
    Node(n+j).AlphaTan =  (1-exp(EtaTan)) * tand(Node(n).Angle) * CNaC + Lambda2/Lambda1*exp(EtaTan)*AfterAlpha1;
    Node(n+j).AlphaR = Node(n+j).AlphaTan*(D/2);
    Node(n+j).AlphaRX =Node(n+j).AlphaTan*(D/2)*(Node(n+j).nodepos);    

    
             Node(n+j).CPDelta =  0;
        Node(n+j).AxForce = 0;   
end

%% Integration to find CNalpha

AlphaR = [Node.AlphaR];
AlphaRX = [Node.AlphaRX];

IntegralCNA = trapz(d_x,AlphaR(1:end)) ;
IntegralCMA = trapz(d_x,AlphaRX(1:end));

BodCNaB = (2*pi/Area)*IntegralCNA;
BodCMaB = -2*pi/(Area*D)*IntegralCMA;

BodCPX = (BodCMaB/BodCNaB)*D;
SumF = [Node.AxForce];
BodCD = sum(SumF);


function y = OgiveY(x, R_Ogive, L, D) %Function to determine the radius of the node
    y = sqrt(R_Ogive^2 - (L - x)^2) + (D/2) - R_Ogive;
end
function y_diff = OgiveGrad(x, R_Ogive, L, D) %Function to determine the gradient at the node
    y_diff = (L-x) / sqrt(R_Ogive^2 - (L - x)^2);
end
end