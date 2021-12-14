classdef SI
% SI constants definitions
% Major constants are taken from the NIST physical constants database:
% http://physics.nist.gov/constants
%
% Constants are defined as class properties and are accessed SI.name where
% name is the name of the desired constant.
%
% Fundamental Constants:
%
% c         Speed of light              [m s-1]
%           *(Defined by SI unit system)
% u0        Vacuum Permeability         [T m A-1] or [H m-1] or [V s A-1 m-1] (Base)
%           *(Defined by SI unit system)
% e0        Vacuum Permittivity         [C V-1 m-1] or [F m-1] or [A2 s4 kg-1 m-3] (Base)
%           *(Defined e0 = 1/(c^2*u0))  
% e         Unit charge                 [C] or [A s] (Base)     
% h         Plank Constant              [J s] or [m2 kg s-1] (Base)
% hbar      Reduced Plank constant      [J s] or [m2 kg s-1] (Base)
% G         Gravitiational Constant     [m^3 kg^-1 s^-2]
%
% Mechanical Constants
% g         Earth Surface gravity       [m s-2]
%
% Atomic Constants
% u        Atomic mass unit             [kg]  
% me        Electron mass               [kg] 
% mp        Free proton mass            [kg]
% mn        Free neutron mass           [kg]    
% md        Free deuteron mass          [kg]
% mt        Free triton mass            [kg]
% ma        Free alpha mass             [kg]
% a         Fine structure constant     [dimesionless]
% Ry        Rydberg Constant            [m-1]
% a0        Bohr Radus                  [m]
% re        Classical electron radius   [m]
%
% Thermondynamic Constants
% R         Ideal Gas constant          [J mol-1 K-1] or [kg m2 s-2 mol-1 K-1] (Base)
% Na        Avagadro's Number           [dimentionless]
% kB        Boltzmann constant for      [J/K] or [m^2 kg s^2 K-1] (Base)
%           Kelvin temperatures
% kbeV      Boltzmann constant for      [J/eV] or [m^2 kg s^2 eV-1] (Base)
%           eV temperatures
% ST        Standard Temperature        [K]
% SP        Standard Pressure           [Pa] or [kg m-1 s-2] (Base)
% torr      Pressure of 1 Torr          [Pa] or [kg m-1 s-2] (Base)
% psi       Pressure of 1 psi           [Pa] or [kg m-1 s-2] (Base)
% cal       Energy in a calorie         [J] or [kg m2 s-2] (Base)
% erg       Energey in an erg           [J] or [kg m2 s-2] (Base)
% 
% Radiation and field Constants
%
% sigma     Stefan Boltzmann Constant   [W m-2 K-4]
% gauss     Magnetic field equal to     [T]
%           one gauss   
%
% Written Swadling Jan 2017
properties(Constant)
        % Fundamental
        
        c  = 299792458;             % Light Speed               [m s-1]              
        u0 = 4e-7*pi;               % Vacuum Permeability       [H m^-1]
        e0 = 1/(SI.u0.*SI.c^2)      % Vaccum permittivity       [F m^-1]                                    
        e  = 1.6021766208e-19;      % Unit charge               [C]                     
        h  = 6.626070040e-34;       % Plank Constant            [J s]
        hbar = SI.h/(2*pi);         % Reduced Plank Constatn    [J s]
        G  = 6.67408e-11;           % Gravitiational Constant   [m^3 kg^-1 s^-2]
        % Mechanical
        
        g = 9.80665;                % Earth Surface gravity     [m s-2]
        % Atomic       
        u = 1.660539040e-27;        % Atomic mass unit          [kg]       
        me = 9.10938356e-31;        % Electron mass             [kg]  
        mp = 1.672621898e-27;       % Free proton mass          [kg]       
        mn = 1.674927471e-27;       % Free neutron mass         [kg]
        md = 3.343583719e-27;       % Free deuteron mass        [kg]
        mt = 5.007356665e-27;       % Free triton mass          [kg]
        ma = 6.644657230e-27;       % Free alphaparticle mass   [kg]
        a = SI.e^2/...              % Fine structure constant   [dimesionless]
            (4*pi*SI.e0*SI.hbar*SI.c);
        Ry = SI.a^2*SI.me*SI.c/...  % Rydberg Constant          [m-1]
             (2*SI.h);
        a0 = 4*pi*SI.e0* ...        % Bohr radius               [m]
             SI.hbar^2/(SI.me*SI.e^2)
        re = (4*pi*SI.e0)^-1*...    % Classical electron radius [m]
             SI.e^2/(SI.me*SI.c^2);          
        % Thermodynamic 
        
        R = 8.3144598;              % Ideal Gas constant        [J mol^-1 K-1]
        Na = 6.022140857e23;        % Avagadro's Number         [dimentionless]       
        kb = SI.R/SI.Na;            % Boltzman Const. (K)       [J K^-1]
        kbeV = SI.e;                % Boltzman Const. (eV)      [J/eV]           
        ST =  273.15;               % Standard Temperature      [K] 
        SP = 1.0133e5;              % Standard Pressure (Atmos) [Pa]                                   
        torr = 1.3332e2;            % Pressure of 1 Torr        [Pa]
        psi = 6894.75729;           % Pressure of 1 PSI         [Pa]
        Ba = 1e-1;                  % Pressure in one Barye     [Pa]
        cal = 4.186                 % Energy in a calorie       [J]
        erg = 1e-7;                 % Energy in an erg          [J]
        % Radiation and fields
        
        sigma = 2*pi^5*SI.kb^4/...  % Stefan Boltzmann Constant [W m-2 K-4]
            (15*SI.h^3*SI.c^2);
        gauss = 1e-4;               % Magnetic field equal to   [T]
                                    % 1 gauss
        Eh = 27.2.*SI.e;            % Hartree Energy energy (27.2eV)
end
    
    
end