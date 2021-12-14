classdef CGS
% CGS constants definitions
% The constant defintitions are derived from the class SI
% Major constants are taken from the NIST physical constants database:
% http://physics.nist.gov/constants
%
% Constants are defined as class properties and are accessed CGS.name where
% name is the name of the desired constant.
%
% Fundamental Constants:
% c         Speed of light              [cm s-1]
% e         Unit charge                 [statC]    
% h         Plank Constant              [erg s] or [g cm2 s-1] (Base)
% hbar      Reduced Plank constant      [erg s] or [g cm2 s-1] (Base)
% G         Gravitiational Constant     [cm^3 g^-1 s^-2]
%
% Mechanical Constants
% g         Earth Surface gravity       [cm s-2]
%
% Atomic Constants
% u        Atomic mass unit             [g]  
% me        Electron mass               [g] 
% mp        Free proton mass            [g]
% mn        Free neutron mass           [g]    
% md        Free deuteron mass          [g]
% mt        Free triton mass            [g]
% ma        Free alpha mass             [g]
% a         Fine structure constant     [dimesionless]
% Ry        Rydberg Constant            [cm-1]
% re        Classical electron radius   [cm]
% a0        Bohr Radius                 [cm]
%
% Thermondynamic Constants
% R         Ideal Gas constant          [erg mol-1 K-1] or [g cm2 s-2 mol-1 K-1] (Base)
% Na        Avagadro's Number           [dimentionless]
% kB        Boltzmann constant for      [erg/K] or [g cm2 s-2 K-1] (Base)
%           Kelvin temperatures
% kBeV      Boltzmann constant for      [erg/eV] or [g cm2 s-2 eV-1](Base)
%           eV temperatures
% ST        Standard Temperature        [K]
% SP        Standard Pressure           [Ba] or [g cm-1 s-2] (Base)
% torr      Pressure of 1 Torr          [Ba] or [g cm-1 s-2] (Base)
% psi       Pressure of 1 psi           [Ba] or [g cm-1 s-2] (Base)
% cal       Energy in a calorie         [erg] or [g cm2 s-2] (Base)
%
% Written Swadling Jan 2017
    properties(Constant)
        c  = SI.c*1e2;              % Light Speed               [cm s-1]                                        
        e  = 0.1*CGS.c*SI.e;        % Unit charge               [statC]                     
        h  = SI.h/SI.erg;           % Plank Constant            [erg s]
        hbar = CGS.h/(2*pi);        % Reduced Plank Constatn    [erg s]
        G  = SI.G*1e3;              % Gravitiational Constant   [cm^3 g^-1 s^-2]
        % Mechanical
        g = SI.g*1e2;               % Earth Surface gravity     [cm s-2]
        % Atomic             
        u = SI.u*1e3;               % Atomic mass unit          [g]       
        me = SI.me*1e3;             % Electron mass             [g]  
        mp = SI.mp*1e3;             % Free proton mass          [g]       
        mn = SI.mn*1e3;             % Free neutron mass         [g]
        md = SI.md*1e3;             % Free deuteron mass        [g]
        mt = SI.mt*1e3;             % Free triton mass          [g]
        ma = SI.ma*1e3;             % Free alphaparticle mass   [g]
        a = SI.a;                   % Fine structure constant   [dimesionless]
        Ry = SI.Ry*1e-2;            % Rydberg Constant          [cm-1]
        a0 = SI.a0*1e2;             % Bohr radius               [cm]
        re = SI.re*1e2;             % Classical electron radius [cm]         
        % Thermodynamic                                   
        R = SI.R/SI.erg;            % Ideal Gas constant        [erg mol^-1 K-1]
        Na = SI.Na;                 % Avagadro's Number         [dimentionless]       
        kb = SI.kb/SI.erg           % Boltzman Const. (K)       [erg K^-1]
        kbeV = SI.kbeV/SI.erg;      % Boltzman Const. (eV)      [erg/eV]           
        ST =  273.15;               % Standard Temperature      [K] 
        SP = SI.SP/SI.Ba;           % Standard Pressure (Atmos) [Ba]                                  
        torr = SI.torr/SI.Ba;       % Pressure of 1 Torr        [Ba]
        psi = SI.psi/SI.Ba;         % Pressure of 1 PSI         [Ba]
        cal = SI.cal/SI.erg;        % Energy in a calorie       [erg]
        sigma = SI.sigma/...        % Stefan Boltzmann constant [erg cm-2 s-1 K-4]
                (1e4*SI.erg);
    end
end

