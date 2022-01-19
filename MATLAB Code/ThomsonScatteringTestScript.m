% test script for Thomson code
% Makes a plot of a Thomson scattering spectrum

ne= 4e25 ; % electron density m^{-3}
Te=1000; % electron temperature eV
Ti=400; % ion temperature eV
Z=2; % average ionization 
A=4; % atomic mass 
vI=0; % ion velocity [m s^-1]
J=0;% plasma current [A cm^-3];
fract=1; % fractional contribution to plasma (multi spieces plasma)
lambda0 = (532); % Probe wavelength
lambdaRange = (400:0.01:700); % wavelengths to calculate the crosssection at
scatAngle = 90; % Scattering angle
order =1; % Order of calcualtion (1st or 2nd order relativistic correction)
cgsUnits = 0; % CGS unit flag (0 defaults to SI)
nmLambda =1; % nm wavelength units flag

%% Calculate the underlying spectrum
 [scatteringCrossSection,formFactor,totalScatteringRatio]  ...
    = thomsonCrossSection(lambdaRange,lambda0,ne,Te,fract,Ti,A,Z,vI,J,scatAngle,order,cgsUnits,nmLambda);

figure('Color','white');
nmLambda =1;
plot(lambdaRange,scatteringCrossSection)
xlabel('Wavelength [nm]');
ylabel('Cross Section [m^-1 nm^-1 sr^-1]')
set(gca,'YScale','linear');


%% Apply Arbitary Spectrometer Broadening
instrumentWidth = 0.5; % gaussian stanrdard deviation of instrument fn in nm.
scatteringCrossSection =gaussianBroadening(lambdaRange,scatteringCrossSection,instrumentWidth);

figure('Color','white');

plot(lambdaRange,scatteringCrossSection)
xlabel('Wavelength [nm]');
ylabel('Cross Section [m^-1 nm^-1 sr^-1]')
set(gca,'YScale','linear');
