classdef Photon
   % Class that handles optical units. Knows frequency /
   % wavelength / energy in vacuum
    properties(Access = private)
        frequencyData
        spectralWidthData
        refractiveIndexData = 1;
    end
    properties(Access = public, Dependent)
        frequency
        wavelength
        vacuumWavelength
        omega
        k
        eV
        spectralWidth
        wavelengthNM
        nm
        spectralWidthNM
        wavelengthCM   
        cm
        joule       
        kCM
        refractiveIndex
        electronDensity
    end
    
    methods
        function obj = Photon(input)
            % input is either a Photon class of wavelength in NM
            if nargin >0
                if isa(input,'Photon')
                    obj = input;
                else
                    obj.wavelengthNM = input;
                end
            end
        end
        
         function obj = set.spectralWidth(obj,spectralWidth)
            obj.spectralWidthData = spectralWidth;
        end
        function spectralWidth =  get.spectralWidth(obj)
            spectralWidth = obj.spectralWidthData;
        end
            
         function obj = set.refractiveIndex(obj,refractiveIndex)
            obj.refractiveIndexData = refractiveIndex;
        end
        function refractiveIndex =  get.refractiveIndex(obj)
            refractiveIndex = obj.refractiveIndexData;
        end
        
        function obj = set.electronDensity(obj,electronDensity)
            obj.refractiveIndexData = refractiveIndex(electronDensity,obj.vacuumWavelength);
        end
        function electronDensity =  get.electronDensity(obj)
            electronDensity =(1-obj.refractiveIndexData.^2).*obj.criticalDensity;           
        end
        function obj = set.spectralWidthNM(obj,spectralWidth)
            obj.spectralWidthData = spectralWidth.*1e-9;
        end
        function spectralWidth =  get.spectralWidthNM(obj)
            spectralWidth = obj.spectralWidth.*1e9;
        end
        
        function obj = set.frequency(obj,frequency)
            obj.frequencyData = frequency;
        end
        function frequency =  get.frequency(obj)
            frequency = obj.frequencyData;
        end
           function obj = set.omega(obj,omega)
            obj.frequencyData = omega./(2*pi);
        end
        function omega =  get.omega(obj)
            omega = obj.frequencyData.*(2*pi);
        end     
        function obj = set.vacuumWavelength(obj,vacuumWavelength)
            obj.frequencyData = SI.c./vacuumWavelength ;
        end
        function vacuumWavelength =  get.vacuumWavelength(obj)              
            vacuumWavelength = SI.c./obj.frequencyData;
        end    
        function obj = set.wavelength(obj,wavelength)
            obj.vacuumWavelength = wavelength./obj.refractiveIndexData;
        end 
        function wavelength =  get.wavelength(obj)              
            wavelength = obj.vacuumWavelength./obj.refractiveIndexData;
        end          
        function obj = set.wavelengthNM(obj,wavelengthNM)
            obj.wavelength = wavelengthNM.*1e-9;
        end
        function wavelengthNM =  get.wavelengthNM(obj)
            wavelengthNM = obj.wavelength.*1e9;
        end        
        function obj = set.wavelengthCM(obj,wavelengthCM)
            obj.wavelength = wavelengthCM.*1e-2;
        end
        function wavelengthCM =  get.wavelengthCM(obj)
            wavelengthCM = obj.wavelength.*1e2;
        end     
        function obj = set.k(obj,k)
            obj.wavelength = (2*pi)./k;
        end
        function k =  get.k(obj)
            k = (2*pi)./obj.wavelength;
        end
        function obj = set.kCM(obj,kCM)
            obj.k = kCM.*1e2;
        end
        function kCM =  get.kCM(obj)
            kCM = obj.k.*1e-2;
        end
        function obj = set.joule(obj,joule)                       
            obj.frequencyData= joule./SI.h;
        end         
        function joule =  get.joule(obj)
            joule = obj.frequencyData.*SI.h;
        end        
        function obj = set.eV(obj,eV) 
            obj.joule= eV.*SI.e;
        end         
        function eV =  get.eV(obj)
            eV = obj.joule./SI.e;
        end
        % These properties are just aliases
        function cm =  get.cm(obj)
            cm = obj.wavelengthCM;
        end     
        function obj = set.cm(obj,cm)
            obj.wavelengthCM = cm;
        end
        function nm =  get.nm(obj)
            nm = obj.wavelengthNM;
        end     
        function obj = set.nm(obj,nm)
            obj.wavelengthNM = nm;
        end
        
    end    
    methods
       
        function nc = criticalDensity(obj)
            nc= criticalDensity(obj.vacuumWavelength);           
        end
        
    end
    
end

