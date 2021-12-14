function [scatteringCrossSection,formFactor,totalScatteringRatio,alphaM]  ...
    = thomsonCrossSection(measured,probe,neTotal,Te,fract,Ti,A,Z,vI,J,scatAngle,order,cgsUnits,nmLambda,superGaussian,verbose)
%% function [scatteringCrossSection,formFactor,totalScatteringRatio,alpahM]  ...
%    = thomsonCrossSection(measured,probe,neTotal,Te,fract,Ti,A,Z,vI,J,scatAngle,order,cgsUnits,nmLambda,superGaussian,verbose)
%
% This function calculates the expected Thomson scattering cross section.
% Takes account of multiple ionic species.
%
% Inputs
%
% For multi component plasmas these inputs should be supplied as arrays
% of length equal to the number of ion components, each labeled by (j).
%
%                                                        SI      CGS     nm
% measured     wavelengths to calculate                  [m]     [cm]    [nm]   or  [photon class]
% probe    probe wavelength                              [m]     [cm]    [nm]   or  [photon class]
% neTotal    Total electron density                      [m-3]   [cm-3]
% Te         electron Temperature                        [eV]
% Ti(j)      array of ion species temperatures           [eV]
% A(j)       array of ion species atomic masses          [AU]
% Z(j)       array of ion species average ionisations    []
% vI(j)      array of ion species flow velocities        [m s-1] [cm s-1]
%            (Flow velocity is defined in direction
%            of scattering k-vector)
% J          net Current density                         [A cm^-2]
% scatAngle  Scattering angle                            [deg]
% order      Defines order in beta.
%            By default 1, can be set to 2.
% cgsUnits   Sets input units to CGS                     [0 or 1]
% nmLambda   Sets wavelength units to nm                 [0 or 1]
%
%% Outputs:
%
% scatteringCrossSection (This is the cross section times the density)
%
% Spectrally resolved differential cross-section for Thomson scattering.
% d\sigma/(d\lambda d\Omega). Intergrate over laser intensity in scattering
% volume with polarisation factor, and integrate over collected solid angle
% to find measured signal.
% Units: [solid angle pathLength wavelength]
% SI:       [m^-1   m^-1 sr^-1]
% SI(nm):   [m^-1  nm^-1 sr^-1]
% CGS:      [cm^-1 cm^-1 sr^-1]
% CGS(nm):  [cm^-1 nm^-1 sr^-1]
%
%
% Form factor S(w,k) for Thomson scattering. Recast as spectral density as
% a function of wavelength rather than frequency.
%
% Units:
% SI:       [ m^-1]
% SI(nm):   [nm^-1]
% CGS:      [cm^-1]
% CGS(nm):  [nm^-1]
%
% totalScatteringRatio
%
% Spectrally integrated cross-section (fraction of light scattered)
%
% Units:
% SI:  [m^-1 sr^-1]
% CGS:  [cm^-1 sr^-1]
%
% Written by George 12/15
% Updated 2/3/15 with some minor bug fixes
% Update 2/17 With unit changes and removed unneccesary code
% Updates 7/17 with corrections to SI units of output

% Convert all inputs to doubles:

%% Default inputs:


arguments 
    % Demo mode values
    measured=300e-9:100e-12:600e-9;
    probe=532e-9;
    neTotal=1e26;
    Te=2000;
    fract=1;
    Ti=2000;
    A=9;
    Z=4;
    vI=0;
    J=0;
    scatAngle=60;
    % Default to first order
    order =1;
     % Default to SI units
    cgsUnits=0;
    % Default to normal wavlength scale
    nmLambda =0;
    superGaussian =2.0;
    verbose =0;
end
if nargin <11
    warning('demo mode')
end


neTotal =double(neTotal);
Te =double(Te);
fract =double(fract);
Ti =double(Ti);
A =double(A);
Z =double(Z);
vI =double(vI);
J =double(J);
scatAngle =double(scatAngle);


% If GCS units then convert to SI for internal calculations
if cgsUnits==1
    neTotal = neTotal.*1e6;
    vI = vI.*1e-2;
end
if ~isa(measured,'Photon')
    measured =double(measured);
    probe =double(probe);
    if nmLambda==1
        measured = Photon(measured);
        probe = Photon(probe);
    elseif cgsUnits==1
        measured = Photon(measured.*1e7);
        probe= Photon(probe.*1e7);
    else
        measured = Photon(measured.*1e9);
        probe= Photon(probe.*1e9);
    end
end

%% Clean up the user input and format data for calculation
if abs(sum(fract)-1)>1e-2
    warning('Input fractions have been normalized');
    fract = fract./sum(fract);
    
end
%% Determine populated populations, discard unpopulated inputs.
pop= fract > 0; fract=fract(pop);
Ti=Ti(pop);A=A(pop);Z=Z(pop);vI =vI(pop);

%% Calculate the "individual" population densities and associated partial electron densities
% Mean ionization state:
Zbar = sum(fract.*Z);
% Total ion density
niT = neTotal./Zbar;
% Population densities
ni =fract.*niT;
% Z weighted electron densities per population
ne = ni.*Z;

% Save Input shape
shapeLambda = size(measured.wavelength);
% Shape Input for code
measured = Photon(reshape([measured.nm],1,[]));

% Correct multi species shapes.
Ti = Ti(:)';
A = A(:)';
Z = Z(:)';
vI = vI(:)';
ne = ne(:)';


[scatteringCrossSection,formFactor,alphaM]  ...
    = coreThomsonCrossSection(measured,probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,superGaussian);
% Reshape to match input...
formFactor = reshape(formFactor,shapeLambda);
scatteringCrossSection = reshape(scatteringCrossSection,shapeLambda);

% Calculate the total Scattering Cross section:
if nargout>2 && numel(scatteringCrossSection)
    
    totalScatteringRatio = trapz(measured.nm,scatteringCrossSection);
    %[m^-1 sr-1] [per path length per solid angle]
else
    
    
    
    totalScatteringRatio = [];
end

%% Unit convertsion logic

if cgsUnits==0 && nmLambda == 0
    formFactor = formFactor.*1e9;
    %[nm^-1] -> [m-1]
    scatteringCrossSection = scatteringCrossSection.*1e9;
    %[sr^-1 m^-1 nm^-1] -> [sr-1 m-1 m-1]
end

% CGS units with nm spectral
if cgsUnits ==1 && nmLambda == 1
    %formFactor = formFactor;
    %[nm^-1] -> [nm-1]
    scatteringCrossSection = scatteringCrossSection.*1e-2;
    %[sr^-1 m^-1 nm^-1] -> %[sr^-1 cm^-1 nm^-1]
    totalScatteringRatio =totalScatteringRatio.*1e-2;
    %[sr-1 m^-1] -> [sr-1 cm^-1]
end

% SI units with nm spectral Everything is done!
if cgsUnits==1 && nmLambda==0
    formFactor = formFactor.*1e7;
    %[nm^-1] -> [cm-1]
    scatteringCrossSection = scatteringCrossSection.*1e5;
    %[sr^-1 m^-1 nm^-1] -> %[sr^-1 cm^-1 cm^-1]
    totalScatteringRatio = totalScatteringRatio.*1e-2;
    %[sr-1 m^-1] -> [sr-1 cm^-1]
end
end
