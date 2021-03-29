function fitresult = CNaC
Angles = 0:2.5:25;
Mach_Numbers = [1.5 1.75 2 2.5 3 3.5 4 4.5 5];
Normal_Force_Derivative = [2 1.9770483 1.9336127 1.8821566 1.8280992 1.7735619 1.7187816 1.6630926 ...
    1.6057182 1.5462941 1.4852936;... Mach 1.5 Row
    2 1.9731742 1.9258447 1.8737533 1.8224578 1.7730646 1.7243971 1.6746400 ...
    1.6223425 1.5666956 1.5075064;... Mach 1.75 Row
    2 1.9690287 1.9180672 1.8661051 1.8183423 1.7742540 1.7310117 1.6859216 ...    
    1.6372641 1.5842135 1.5265736;... Mach 2 Row
    2 1.9602967 1.9036740 1.8548454 1.8161039 1.7820973 1.7472553 1.7083329 ...
    1.6638860 1.6134668 1.5571606;... Mach 2.5 Row
    2 1.9514206 1.8917127 1.8493550 1.8202950 1.7941549 1.7645092 1.7287097 ...
    1.6859461 1.6362061 1.5798551;... Mach 3 Row
    2 1.9426977 1.8826075 1.8487600 1.8282076 1.8074148 1.7806466 1.7461812 ...
    1.7037859 1.6537994 1.5968132;... %Mach 3.5 Row
    2 1.9343213 1.8763907 1.8516843 1.8377630 1.8203259 1.7949749 1.7608321 ...
    1.7181330 1.6674853 1.6096523;... %Mach 4 Row
    2 1.9264473 1.8728601 1.8568202 1.8477350 1.8322685 1.8074272 1.7730368 ...
    1.7296967 1.6782310 1.6195188;... %Mach 4.5 Row
    2 1.9191785 1.8715909 1.8632046 1.8574765 1.8430453 1.8181556 1.7832011 ...
    1.7390731 1.6867491 1.6272149]; %Mach 5 Row
%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( Angles, Mach_Numbers, Normal_Force_Derivative );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );