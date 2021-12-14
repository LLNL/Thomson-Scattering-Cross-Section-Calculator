function [scatteringCrossSection,formFactor,alphaM]  ...
    = coreThomsonCrossSection(measured,probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,superGaussian)
%% function [scatteringCrossSection,formFactor,totalScatteringRatio]  ...
%    = coreThomsonCrossSection(measured,probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,M)
%
% This function calculates the expected Thomson scattering cross section.
% Takes account of multiple ionic species.
%
%% Inputs
%
% For multi component plasmas these inputs should be supplied as arrays
% of length equal to the number of ion components, each labeled by (j).
%
%
% measured     wavelengths to calculate                  [nm]
% probe    probe wavelength                              [nm]
% neTotal    Total electron density                      [m-3]
% ne(j)      population electron density                 [m-3]
% Te         electron Temperature                        [eV]
% Ti(j)      array of ion species temperatures           [eV]
% A(j)       array of ion species atomic masses          [AU]
% Z(j)       array of ion species average ionisations    []
% vI(j)      array of ion species flow velocities        [m s-1]
%            (Flow velocity is defined in direction
%            of scattering k-vector)
% J          net Current density                         [A cm^-2]
% scatAngle  Scattering angle                            [deg]
% order      Defines order in beta.
%            By default 1, can be set to 2.
% superGaussian Super-gaussian distribution function (2 default)
%
%% Outputs:
%
%% scatteringCrossSection (This is the cross section times the density)
%
% Spectrally resolved differential cross-section for Thomson scattering.
% d\sigma/(d\lambda d\Omega). Intergrate over laser intensity in scattering
% volume with polarisation factor, and integrate over collected solid angle
% to find measured signal.
% Units: [solid angle pathLength wavelength]
% SI(nm):   [m^-1  nm^-1 sr^-1]
%
% a
%% Form factor S(w,k) for Thomson scattering. Recast as spectral density as
% a function of wavelength rather than frequency.
%
% Units:
% SI(nm):   [nm^-1]
%
%% totalScatteringRatio
%
% Spectrally integrated cross-section (fraction of light scattered)
%
% Units:
% SI:  [m^-1]
%
% Written by George 12/15
% Updated 2/3/15 with some minor bug fixes
% Update 2/17 With unit changes and removed unneccesary code
% Updates 7/17 with corrections to SI units of output
% Convert all inputs to doubles:

if nargin < 13
    M=2;
else
    M= superGaussian;
end

% Convert SI -> CGS for calculation
neTotal = neTotal.*1e-6;
ne = ne.*1e-6;
vI = vI.*1e2;

measured = Photon(measured);
probe = Photon(probe);


% Find the frequency shift required to observe light of frequency
% measured.omega
deltaOmega =  measured.omega - probe.omega;


%% Calculate the k-vectors before and after scattering
% Calculate electron plasma frequency at scattering volume
omegaPe = plasmaFrequency(neTotal,0,-1,1);
% Find Magnitudes of incoming and outgoing wave-vectors
kIn = sqrt(probe.omega.^2-omegaPe.^2)./CGS.c;
kOut = real(sqrt(measured.omega.^2-omegaPe.^2)./CGS.c);
kOut(measured.omega<omegaPe)=0;

% Build vectors of incoming and outgoing k-vectors, rotate to correct
% angles
kIn = repmat([0;kIn],1,length(kOut));
kOut =[zeros(1,length(kOut));kOut];
% Rotate output k vectors
kOut = [cosd(scatAngle),sind(scatAngle) ; -sind(scatAngle), cosd(scatAngle)]*kOut;
% Find specific k shift for scattering
k = kOut-kIn;
% Magnitude of k
kMagnitude = vecnorm(k);
% Unit of k
kUnit = k./repmat(kMagnitude,2,1);
kInUnit =  kIn./repmat(vecnorm(kIn),2,1);
%% Calculate alpha (electrons)
alpha = 1./(kMagnitude.*debyeLength(neTotal,Te,1));
alphaM = mean(alpha);%
%% Start with calculation of (f, chi) for electron population
% Mean electron velocity (based on weighted sum of ion velocities)
meanVi = sum(vI.*ne)./neTotal;
driftVelocity = (-J./(neTotal.*SI.e));
Ve = meanVi+driftVelocity;
Ve = kUnit.*Ve;
omegaShifted = deltaOmega-(sum(k.*Ve,1));
% Calculate the Phase velocity of scatterers
vPhase = (omegaShifted./kMagnitude); %in cm s^-1
%vPhase(vPhase==0) = 1e-2;
% Check for relativistic scattering
if any([vPhase(1) vPhase(end)] > 0.2.*CGS.c)
    if order==2 && any([vPhase(1) vPhase(end)] < 0.5.*CGS.c)
    else
        warning(['Approaching relativistic limit B ~ ' num2str(round(max([vPhase(1),vPhase(end)]./CGS.c)*1000)/1000)]);
    end
end

%% Caclulate the electron velocity distribution (fE) and electron succeptibility (chiE)
if M == 2.0
    % Calculate the EPW based on gaussian distribution function
    vTe = sqrt(2.*SI.kbeV.*Te./SI.me).*1e2;
    x = vPhase./vTe; % Ratio of phase velocity over thermal velocity
    
    % Get the dispersion function gradient - second output of this function
    [~, Zp] = plasmaDispersionFunction(x);
    % Calculate the susceptibility for the elctrons and ions
    chiE = -alpha.*alpha.*Zp./2;
    % Get distribution function:
    fE = sqrt(1/(pi*vTe^2)).*exp(-x.^2); 
    %fE =imag(-Zp)./(2.*pi.*vPhase);
else
    % Calculate the EPW based on super-gaussian distribution function
    % Calculate weighted normalized phase velocity
    g5 = gamma(5.0./M);
    g3 = gamma(3.0./M);
    g2 = gamma(2.0./M);
    x = (vPhase./thermalVelocity(Te,0,1)).*sqrt(2.*g5/(3.0.*g3));
    % Calculate the derivative plasma dispersion fuction
    [~, Zprime] = plasmaDispersionFunction(x);
    % Calulate the modified electron suceptibility for Mth order super
    % gaussian distribution function
    %{
    integrand = @(y) y./(y-x).*(exp(-abs(y).^M)-exp(-abs(x).^M+x.^2-y.^2));
    chiRe = -sqrt(pi)./2.0.*real(Zprime).*exp(x.^2-abs(x).^M) +...
        integral(integrand ,-inf, inf, 'ArrayValued', true);
    %}
    chiRe = -sqrt(pi)./2.0.*real(Zprime).*exp(x.^2-abs(x).^M) +...
        fastSuperGaussianIntegral(x,M);
    chiIm = pi.*x.*exp(-abs(x).^M);
    chiE= alpha.^2.*(M./6.0).*(g5./g3.^2).*(chiRe + 1i*chiIm);
    % Calculate the weighted electron distribution function for the Mth
    % order super gaussian distribution function
    fE =0.5.*abs(x./vPhase).*g2./g3.*gammainc(abs(x).^M,2.0./M,'upper');
end
%% Now calculate for velocity distribution and ion succeptibility
chiI = zeros(length(ne),length(chiE));
fI = chiI;
for j = 1: length(ne)
    % Calculate doppler Shift
    vIVec= kUnit.*vI(j);
    omegaShifted = deltaOmega-(sum(k.*vIVec,1));
    % Calculate the Phase velocity vector for electrons
    vPhase = (omegaShifted./kMagnitude); %in cm s^-1
    % Caluclate the ion thermal velocity and use this to calculate x, the
    % input to the PlasmaDispersionRelation
    vTi = sqrt(2.*SI.kbeV.*Ti(j)./(SI.u.*A(j))).*1e2;
    x = vPhase./vTi; % Ratio of phase velocity over thermal velocity
    % Get the dispersion function gradient - second output of this function
    [~, Zp] = plasmaDispersionFunction(x);
    % Calculate the susceptibility for the elctrons and ions
    % Includes a population weighting factor which effectively corrects
    % alpha^2 term.
    chiI(j,:) = -(ne(j)./neTotal).*(alpha.^2./2).*(Z(j)*Te/Ti(j)).*Zp;
    % Get distribution function:
    fI(j,:) = sqrt(1/(pi*vTi^2)).*exp(-x.^2); 
   % fI(j,:) = imag(-Zp)./(2.*pi.*vPhase);
end

chiITotal = sum(chiI,1);

%% Find total relative permitivity for all populations
epsilon = 1+chiE+chiITotal;

%% Calculate the form factor
%% First Order Calculation
if order ==1
    electronTerm = abs(1-(chiE./epsilon)).^2.*fE;
    ZTerm = repmat((Z.*ne./neTotal)',1,size(fI,2));
    chiEterm=repmat(abs(chiE./epsilon),size(ne,2),1);
    ionTerm = sum(ZTerm.*chiEterm.^2.*fI,1);
    % First order relativistic correction!
    relativisticCorrection  = (1+(2*(deltaOmega)/probe.omega));
    formFactor = relativisticCorrection.*(2*pi./kMagnitude).*(electronTerm + ionTerm);
    
    %% Second Order Caclulation
elseif order ==2
    % Define beta as the phase velocity of the waves over c
    beta = (deltaOmega./kMagnitude)./CGS.c;
    % Zeta is the dot product of the K and the scattering vector
    zeta =dot(kInUnit,kUnit);
    
    % Electron term
    electronTerm = fE.*(...
        ((1+deltaOmega./probe.omega).^2-beta.^2+2.*zeta.*(deltaOmega./probe.omega).*beta)...
        .*abs(1-(chiE./epsilon)).^2 ...
        +(zeta.*(deltaOmega./probe.omega).*beta)...
        .*abs(chiE./epsilon).^2 ...
        +((omegaPe./kMagnitude)./CGS.c).^2.*...
        real((epsilon-chiE)./abs(epsilon).^2) ...
        -0.5.*beta.^2);
    
    % Ion term
    ZTerm = repmat((Z.*ne./neTotal)',1,size(fI,2));
    chiTerm =repmat(...
        ((1+deltaOmega./probe.omega).^2 - beta.^2 + 2.*zeta.*(deltaOmega./probe.omega).*beta)...
        .* (abs(chiE./epsilon)).^2 ...
        -((omegaPe./kMagnitude)./CGS.c).^2 ...
        .*real(chiE./abs(epsilon).^2)...
        ,size(ne,2),1) ;
    
    ionTerm = sum(ZTerm.*chiTerm.*fI,1);
    
    % Form Factor
    formFactor =(2.*pi./(kMagnitude)).*(electronTerm + ionTerm);
end
if neTotal ==0
    formFactor = zeros(size(measured));
end

%% Correct as function of wavelength rather than angular frequency
formFactor = formFactor.*2.*pi.*CGS.c./(measured.wavelengthCM.^2);
% unit is in [cm^-1] (i.e. per wavelength in cm)

%% Calibrate absolutely
% Actual cross-section [cm^2 cm^-1 sr^-1] [area per wavelength per solid angle]
scatteringCrossSection = ((CGS.re.^2)./(2.*pi)).*formFactor;

% What I output - weighted by ne [cm^-1 cm^-1 sr^-1] [per probe path length per wavelength per solid angle]

scatteringCrossSection = neTotal.*scatteringCrossSection;

%% Calculation is perfomed in CGS.
%% These lines covert to SI and to nm spectral units

formFactor = formFactor.*1e-7;
%[cm^-1] -> [nm-1]
scatteringCrossSection = scatteringCrossSection.*1e-5;
%[sr^-1 cm^-1 cm^-1] -> %[sr^-1 m^-1 nm^-1]
  
   
end
 function y =fastSuperGaussianIntegral(x,m)
        % This is a fast interpolator for the integral in the superGaussian
        % spectrum
        
        % Uses a persistent variable to store the interpolator
        persistent interpolatorFn
        % Range for the interpolator
        mRange = [2  6];
        
        % If interpolator function has not been loaded then load / generate
        % it.
        if isempty(interpolatorFn)
            % Try to load the interpolator function
            folder = functionPath;
            file = fullfile(folder,'fastSuperGaussianIntegral.mat');
            try    
                loadedData = load(file);
                interpolatorFn = loadedData.interpolatorFn;
            catch
                % If load fails then calculate the interpolator
                % Generate a logarithmically spaced interpolation grid
                M =log10(mRange(1)):0.01:log10(mRange(2));
                M = 10.^M;
                X = [0:0.001:1 5 10 100];
                X = (10.^X)-1; 
                Y = zeros(length(X),length(M));
                for ind = 1:length(M)
                    integrand = @(y) y./(y-X).*(exp(-abs(y).^M(ind))-exp(-abs(X).^M(ind)+X.^2-y.^2));
                    Y(:,ind)= integral(integrand ,-inf, inf, 'ArrayValued', true);
                end
                % Create interpolator function
                interpolatorFn = @(m,x) interp2(M,X,Y,m,x);             
                save(file,'interpolatorFn');            
            end
        end
        
       
        if m>=mRange(1) && m<= mRange(2)
             % If in range then use interpolator
        x= abs(x);
            y= interpolatorFn(m,x);
        else
            % Else calculate answer
            integrand = @(y) y./(y-x).*(exp(-abs(y).^m)-exp(-abs(x).^m+x.^2-y.^2));
            y(:)= integral(integrand ,-inf, inf, 'ArrayValued', true);
        end
        % Shape ouput to match input
        y = reshape(y,size(x));
    end