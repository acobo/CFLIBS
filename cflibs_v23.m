% MATLAB implementation of a calibration-free algorithm based on:
% CIUCCI, A., et al. New procedure for quantitative elemental analysis by laser-induced plasma spectroscopy. Applied spectroscopy, 1999, vol. 53, no 8, p. 960-964.
% and some other works referenced in the code
% version:
%  v23 mar2021 first version posted on github.
%
%{
MIT License
Copyright (c) 2021 Adolfo Cobo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

%% CONFIGURATION FOLDER WITH EXPERIMENTAL DATA TO LOAD

PN='C:\Users\adolf\Documents\UNICAN\Proyecto deepRAMP - LIBS\Arqueo\PatellaDepressa_Oct2012\14112019-130724_LAN.554_DualPulse_0.05';

load([ PN, '\', 'matlabData'] , '-regexp', '^(?!PN$|objWin|objDoc|objDocs|objSpe|objExp|objPul$|tmp_audio|tmp_spectra_snr|frame)\w'); % read the variables in the .MAT data file, 
spectra_is_valid(:,:,:)=1;




%% CONFIGURATION SECTION
DO_EXCEL = 0;  % export results in excel
WHICH_RATIO_TO_EXPORT = 'molar_cf_z1'; %   'lir_from_mean' = line intensity ratio from mean spectrum; 'mean_lir' = mean of ratios from individual spectra; 'molar_cf' = molar concentration from cf libs; 'molar_cf_z1' with z1 substitution ; 'weight_cf_z1' weight concentration (not molar)
DO_SSA = 0;
SSA_SUAVIZADO = 2;  % grado de suavizado SSA: 2=fuerte, 4=flojo, 1=stronger, changing dimensions M=16 gives even more filtering
QUITAR_RESINA = 0;
REARRANGE = 0; % si es mayor que cero hace un nuevo reparto de bloques de x espectros en param1. OJO!! con Ncleaningshots
DISCARD_NONVALID_SPECTRA = 1; % non valid spectra are discarded , different for CFLIBS and non CFLBS processing si vale 1, antes de procesar, se marcan como nó validos los picos demasiado flojos o saturados; =2 no se comprueba los picos de magnesio, es para muestras donde el elemento minoritario está al 0%
DISCARD_ALSO_WEAK_PEAKS = 0; % if 1, an entire spectrum is marked as not valid if ONE of the peaks in composition is below the value of 
REMOVE_BACKGROUND = 1;  % =1 SUBSTRACT the ELECTRONIC background using the first spectra with the laser not firing. software baseline correction (removing emission background) has been moved to flag REMOVE_BASELINE
REMOVE_BASELINE = 0; % =1  applied Yun2009 algorithm to remove baseline (emission background), not related to electronic background removing
BUSCAR_PIXEL = 1; % poniendo 2 no busca el pixel, con 0 si (sin warnings) y con 1 también (con warnings)
INTEGRAR = 3;    %% NOT FOR CF:  0= altura en bruto; 1= integra en bruto ; 2= altura relativa a la linea interpoladora con solo dos puntos laterales;3= altura relativa respecto a la línea de regresión con dos zonas laterales completas; 4=ajuste lorentziano respecto a linea de ajuste del background lateral; 5=ajuste lorentziano pero devuelve la altura total y no necesita el background lateral
OUTLIERS = 0.2; % porcentaje de desviación de la intensidad total para descartar un espectro, poner un valor alto para dejar todos
CORREGIR_IRRADIANCIA = 1; % 1= se dividen todos los espectros por la respuesta, igual que con CF=1
USAR_HR2000=0;  %=1 copia los espectros del hr2000 a spectra()
USAR_PIMAX3=0; %=1 se supone que hemos medido en paralelo con el pimax y está en otra variable, lo copia todo a spectra()
USAR_AVANTES=1; %=1 usamos los datos de avantes, que son directamente spectra()
REMOVE_RESONANT_LINES=0; % =1 remove resonant lines (piResonant()=1) from database (piUse=0)
REMOVE_LOWER_LEVEL_ZERO=0; % =1 remove lines with Ei (lower level) < 0.2 !!!!!! do not use until all lines have their Ei declared!!!
ONLY_CF_MEAN=1; %=1 does not process individual spectra, only mean spectra
%peakModelWidthFactor=3.0; % width of the peaks for fitting, with respect to the instrumental FWHM declared for each spectrometer
CORRECT_SELFABSORPTION = 0; % 0= no correction, 1=Praher2010, 2=Sun2009 (IRSAC)
FIRST_SHOT_FOR_AVERAGING = 0; % the first laser shot to calculate the average,useful to reach a stable plasma temperature with many shots on the same spot
LAST_SHOT_FOR_AVERAGING =9999; % the last laser shot to calculate the average,useful to select a small subset of shots with simmilare plasma parameters
AVERAGE_ONLY_FIRST=0; %if >0 then the N first valid spectra (starting at firtLaserPulse) are marked as valid, all other discarded
REMOVE_OVERLAPPED_LINES = 1; % if =1 the wavelengths of used lines (piUse=1) are sorted, and too close lines discarded
ALTERNATE_DELAY_WHICH_ONE = 0; % 0=do nothing, leave spectra() at it is; =1 difference of hotter and colder spectra; =2 leave only the HOTTER spectra; =3 leave only the COLDER spectra
HALFA_PERCENTAGE=0; % 0=only other's peaks Stark is used for Ne ; 100=only Halfa Stark is used for Ne
LINE_INTENSITY = 0; % for CF with peak fitting: 0=height over baseline, 1=area (calculated using AREA_METHOD)
AREA_METHOD = 3; % for CF, how to calculate peak area: 0=area of modeled lorentzian peak; 1=trapezoidal integration of the oversampled modeled lorentzian peak; 2=two halves ; 3=integrated area of the NOT oversampled baseline-corrected SPECTRA
NORMALIZAR=0; % divides each spectrum to the summed intensities, relative to a total averaged summed intensity, so the overall intensities do not change much
RETAIN_FIRST_FIT = 1; % if 1, these parameters are calculated for first spatial point (idx1=idx2=1) and used for all remaining points:  pixleft, pixright, pixleftbkg, pixrightbkg, 

%% DEBUGGING OPTIONS
%Nparam1=20;  % to speed things up for debugging, only this spatial point is processed
DEBUG_DISCARDING=0; % shows information of the spectra-discarding step for averaging, with specific statistics of idx1=DEBUG_DISCARDING
DEBUG_PEAKS=0; % plot peak fitting, any number not zero and will stop at that idx and show the plots
DEBUG_BP=0; % plot Boltzmann-plots any number not zero is the idx1 value to display
DEBUG_SBP=0; % plot Saha-boltzmann-plots
DEBUG_NE_HISTOGRAM=0; % plot Ne from Stark histogram
DEBUG_STABILITY=0; % OJO!! not really implemented yet, info about stability of calcs point by point (such as number of valid peaks in each BP...)


%% CALIBRATION-FREE DATABASE 
CFNlines = 600; % peaks to process, max number just for matrices allocation
Halfa=1;CaI=2;CaII=3;MgI=4;MgII=5;SrI=6;SrII=7;ZnI=8;ZnII=9;CuI=10;CuII=11;AlI=12;AlII=13;SnI=14;SnII=15;VI=16;VII=17;FeI=18;FeII=19;NiI=20;NiII=21;MnI=22;MnII=23;TiI=24;TiII=25;HgI=26;HgII=27;ArI=28;ArII=29;NeI=30;NeII=31;CI=32;CII=33;
species = [ Halfa CaI CaII MgI MgII SrI SrII ZnI ZnII CuI CuII AlI AlII SnI SnII VI VII FeI FeII NiI NiII MnI MnII TiI TiII HgI HgII ArI ArII NeI NeII CI CII];
speciesStr = { 'Halfa',  'CaI ', 'CaII', 'MgI ', 'MgII', 'SrI ', 'SrII', 'ZnI','ZnII','CuI','CuII',   'AlI','AlII',   'SnI','SnII',     'VI','VII',       'FeI','FeII', 'NiI','NiII'       ,'MnI','MnII',  'TiI','TiII' , 'HgI','HgII'  , 'ArI','ArII' , 'NeI','NeII' , 'CI','CII'    }; % each specie
atomicWeights = [1.008 , 40.078, 40.078 , 24.305, 24.305, 87.62 , 87.62, 65.38, 65.38, 63.546, 63.546, 26.982, 26.982, 118.71 , 118.71 , 50.9415, 50.9415, 55.845,55.845,58.6934, 58.69334 , 54.938,54.938, 47.867,47.867 , 0,0 ,             0,0,         0,0,    12.0107, 12.0107  ];
Nspecies = max(species); % number of species
EACH_BP=1;SINGLE_TE_BP=2;SAHA_BP=3;SINGLE_TE_SAHA_BP=4; % enumeration of options for conc calc
% other species not in composition[] , pi_usar()=0, are not proccesed
%  'composition' should be: [neutral, ion, neutral, ion...] (neutrals first, all ions included)
% AND the species with reliable Te shoulde be the first one

%limpets
composition = [CaI CaII MgI MgII ]; % Nordic gold 89% copper, 5% aluminium, 5% zinc, and 1% tin
whichSpecieForTe0 = CaI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot

%Zn+Cu
composition = [CuI CuII ZnI ZnII ]; % laton
realComposition = [ 0.643 0.357 ]; % allow us to calculate a5n error metric
whichSpecieForTe0 = CuI; %  the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot

%20cents coin nordic gold
composition = [CuI CuII AlI AlII ZnI ZnII SnI SnII]; % Nordic gold 89% copper, 5% aluminium, 5% zinc, and 1% tin
realComposition = [ 0.89 0.05 0.05 0.01]; % allow us to calculate an error metric

%20cents coin nordic gold WITHOUT Sn (not enough peaks?)
composition = [CuI CuII AlI AlII ZnI ZnII ]; % Nordic gold 89% copper, 5% aluminium, 5% zinc, and 1% tin
realComposition = [ 0.90 0.05 0.05 ]; % allow us to calculate an error metric

whichSpecieForTe0 = ZnI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot AND the specie to iterate the Te in SahaBP, 0=use the one from singleBP and single_SahaBP
whichSpeciesForSingleBPTe = [CuI ZnI AlI ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CuI ZnI AlI ]; % only these species are used for the simultaneus fitting of all SAHA BPs 

%limpets
composition = [CaI CaII MgI MgII SrI SrII]; % Nordic gold 89% copper, 5% aluminium, 5% zinc, and 1% tin
realComposition = [ 0.97 0.02 0.01 ]; % allow us to calculate an error metric
whichSpecieForTe0 = 0; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [MgI MgII CaI CaII SrI SrII]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [MgI CaI]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SAHA_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP



%Hg-Ar Calibration Lamp lines - FOR CALIBRATION 
%composition = [HgI HgII ArI ArII ]; 
%realComposition = [ 0.5 0,.5 ]; 
%whichSpecieForTe0 = 0; 
%whichSpeciesForSingleBPTe = [];
%whichSpeciesForSingleSahaBPTe = []; 
%getConcentrationsFrom = SINGLE_TE_BP; 

%Ne calibration Lamp lines - FOR CALIBRATION 
%composition = [NeI NeII]; 
%realComposition = [ 1.0  ]; 
%whichSpecieForTe0 = 0; 
%whichSpeciesForSingleBPTe = [];
%whichSpeciesForSingleSahaBPTe = []; 
%getConcentrationsFrom = SINGLE_TE_BP; 



composition = [CI CII FeI FeII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = 0; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII FeI FeII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII FeI FeII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP

composition = [CI CII CuI CuII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = 0; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII CuI CuII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII CuI CuII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP

%limpets ONLY Ca & Mg
composition = [CaI CaII MgI MgII ]; 
realComposition = [ 0.9847 0.0153 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = 0; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [MgI MgII CaI CaII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [MgI MgII CaI CaII ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP


composition = [CI CII FeI FeII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = CI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII FeI FeII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII FeI FeII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP




composition = [CI CII FeI FeII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = CI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII FeI FeII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII FeI FeII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP

composition = [CI CII NiI NiII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = NiI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII NiI NiII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII NiI NiII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP


composition = [CI CII CuI CuII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = 0; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII CuI CuII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII CuI CuII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP



composition = [CI CII NiI NiII ]; 
realComposition = [ 0.9 0.1 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight
whichSpecieForTe0 = NiI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [CI CII NiI NiII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [CI CII NiI NiII  ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SINGLE_TE_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP

%limpets ONLY Ca & Mg
composition = [CaI CaII MgI MgII ]; 
realComposition = [ 0.9939 0.0061 ]; % 25 mmol/mol of Mg, equals 1.53% in weight, 20 mmol/mol are 1.22% in weight, 10mmol/mol 0.61%
whichSpecieForTe0 = CaI; % the Te from BP of this specie is used as Te0 for Saha-Boltzmann plot
whichSpeciesForSingleBPTe = [MgI MgII CaI CaII ]; % only these species are used for the simultaneus fitting of all BPs 
whichSpeciesForSingleSahaBPTe = [MgI MgII CaI CaII ]; % only these species are used for the simultaneus fitting of all SAHA BPs 
getConcentrationsFrom = SAHA_BP;  % EACH_BP or SINGLE_TE_BP or SAHA_BP or SINGLE_TE_SAHA_BP


%% end of configuration section  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
piquees=zeros(CFNlines,1); %identificador de especie
piUse=zeros(CFNlines,1); %si se usan o no, por defecto es 0
pila=zeros(CFNlines,1); %lambdas
pilaobs=zeros(CFNlines,1); %observed lambda for a line (miscalibration...), if not zero, pila() is substituted by pilaobs
pipx=zeros(CFNlines,1); %pixeles
piGA=zeros(CFNlines,1); %producto gxA en el NIST
piEm=zeros(CFNlines,1); %energía del nivel superior
piStark=zeros(CFNlines,1); %cociente Nreferencia/omega de ensanchamiento HWHM stark de la línea, en nm
pigk=zeros(CFNlines,1); % integer degenerate value g_k (for partition function)
piAcc=cell(CFNlines,''); % accuracy of lines from NIST, strings using their format: A+, A, B+ ...
pi_lngAI=zeros(CFNlines,1); %valores de Ln(g.A/I) para el boltzmann-plot as for 13dic2018, WITHOUT LAMBDA
corr_pi_lngAI=zeros(CFNlines,1); %valores de Ln(g.A/I) corregidos para el saha-boltzmann-plot
corr_piEm=zeros(CFNlines,1); %energía del nivel superior
pi_z=zeros(CFNlines,1); % valor de ionización para Saha-Boltzmann, z=0 para I, z=1 para II
NeStarkNe0=1E17; % initival Ne for iterative Ne 

%quick and dirty definition of parameteres for the multiple-line BP fitting
mdl1 = @(beta,x) beta(1)*x+beta(2);
mdl2 = @(beta,x) beta(1)*x+beta(3);
mdl3 = @(beta,x) beta(1)*x+beta(4);
mdl4 = @(beta,x) beta(1)*x+beta(5);
mdl5 = @(beta,x) beta(1)*x+beta(6);
mdl6 = @(beta,x) beta(1)*x+beta(7);
mdl7 = @(beta,x) beta(1)*x+beta(8);
mdl8 = @(beta,x) beta(1)*x+beta(9);
mdl9 = @(beta,x) beta(1)*x+beta(10);
mdl10 = @(beta,x) beta(1)*x+beta(11);


%% Table of emission lines
% i is just an arbitrary index
% piusar=1 if the peak should be used for CF calculations 
% piquees is the ID of the species ( speciesStr(piquees()) gives the name)
% pila is the wavelength
% piEi is Ei lower level from ASD in eV units (default units are cm-1)
% piEm is Ek upper level from ASD in eV units (default units are cm-1)
% piGA is the product g*A, should be chosen from the configuration screen
% pi_Z is the ionization state, 0 for "I" , 1 for "II" ...
% pi_gk is the statistical weight of the upper leve, from ASD (check "g" in configuration screen, lower-right corner) 
% pi_Stark is the linear Stark parameter as Ne_REF/HWHM (nm) units . LIBS++
%    database has a 1E16 ref value , half width values, and wavelength in AMSTRONGS
% IRSAC_reference(specie) = i; 
% Reasons to remove
DISCARD_Ek0 = 0;
DISCARD_TEST = 0;
DISCARD_NOTSEEN = 0;
DISCARD_VAR = 0; % sometimes appears, introduces variability in concentrations

%FIRST peak is special: Halfa for Ne estimation. Its stark parameter is larger but not linearly dependent on Ne
i=1;   piUse(i)=1;piquees(i)=Halfa;pila(i)=656.279; piEm(i)=12.088; piGA(i)=7.94E8; pi_z(i)=0;pigk(i)=18;piStark(i)=1; % Stark parameter not valid, H uses another formula

% 10ene2021: lines marked as 0*0*1 are peaks that not always are valid across sequence for PD.541, trying to obtain a more stable sequence
%CaI
i=2;   piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=299.496; PiEi(i)=1.879;piEm(i)=6.01789; piGA(i)=1.1E8; pi_z(i)=0;pigk(i)=3;  % son dos líneas muy juntas
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=299.731;  PiEi(i)=1.886;piEm(i)=6.0211; piGA(i)=1.2E8; pi_z(i)=0;pigk(i)=5;  % dos muy juntas
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=299.964;  PiEi(i)=1.886;piEm(i)=6.01789; piGA(i)=8.37E7; pi_z(i)=0;pigk(i)=3;  % dos muy juntas
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=300.086;  PiEi(i)=1.886;piEm(i)=6.01622; piGA(i)=1.58e8; pi_z(i)=0;pigk(i)=1;  % dos muy juntas
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=300.686; PiEi(i)=1.899;piEm(i)=6.0211; piGA(i)=3.8E8; pi_z(i)=0;pigk(i)=5; % pico Ca ventana primigenia, realmente son tres mezclados de intensidad parecida
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=300.921; PiEi(i)=1.899;piEm(i)=6.017898; piGA(i)=1.29E8; pi_z(i)=0;pigk(i)=3;  
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=335.03;  piEm(i)=5.585; piGA(i)=4.335E7; pi_z(i)=0;pigk(i)=3;  % NO USAR: dos solapadas OJO! PRUEBA CON EL VALOR MEDIO DE Em, no, sale muy mal el BP
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=336.19;  piEm(i)=0; piGA(i)=0E8; pi_z(i)=0;pigk(i)=7;  % NO USAR: dos solapadas
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=362.41;  piEi(i)=1.879;piEm(i)=5.29946; piGA(i)=6.36E7; pi_z(i)=0;pigk(i)=3;  % 
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=363.075; piEi(i)=1.885;piEm(i)=5.300; piGA(i)=11.25E7; pi_z(i)=0;pigk(i)=5;  % NO USAR: dos solapadas OJO!! PRUEBA CON EL VALOR MEDIO DE Em, no, sale muy mal el BP
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=364.44;  piEm(i)=5.30; piGA(i)=0E8; pi_z(i)=0;pigk(i)=7;  % NO USAR: tres solapadas
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=422.67;  piEi(i)=0.0;piEm(i)=2.93; piGA(i)=6.54E8; pi_z(i)=0; piStark(i)=1E16/0.0063; pigk(i)=3;% CaI ojo! Ek=000000 NO USAR REF: https://griem.obspm.fr/index.php?page=pages/result.php&element=Ca&base=1
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=428.30;  piEi(i)=1.886;piEm(i)=4.7798; piGA(i)=2.17E8; pi_z(i)=0;pigk(i)=5;% 
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=428.93;  piEi(i)=1.879;piEm(i)=4.769; piGA(i)=1.8E8; pi_z(i)=0;pigk(i)=3;%  
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=430.253; piEm(i)=4.779; piGA(i)=6.8E8; pi_z(i)=0;pigk(i)=5;% CaI dos subpicos lados, seems reversed in alternateDelay of apex
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=430.77;  piEm(i)=4.7631; piGA(i)=1.99E8; pi_z(i)=0;pigk(i)=1;% CaI 
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=431.86;  piEi(i)=1.8989;piEm(i)=4.769; piGA(i)=2.2E8; pi_z(i)=0;pigk(i)=3;piStark(i)=1E16/0.00077; % CaI REF: BD LIBS++ (Amstr-1E16-halfwidth)
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=442.54;  piEm(i)=4.68; piGA(i)=1.49E8; pi_z(i)=0;pigk(i)=3;piStark(i)=1E16/0.00145; % CaI % CaI del artículo de pandhija  442.5nm
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=443.569; piEm(i)=0.0; piGA(i)=0; pi_z(i)=0;pigk(i)=3;% CaI  NO USAR dos solapadas OJO!! 
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=445.47;  piEm(i)=4.681; piGA(i)=6.1E8; pi_z(i)=0;pigk(i)=5;% CaI% CaI del artículo de pandhija  445.5nm
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=504.16;  piEm(i)=5.167; piGA(i)=9.9E7; pi_z(i)=0;pigk(i)=3;% CaI aislado *
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=518.88;  piEm(i)=5.321; piGA(i)=2E8;pi_z(i)=0; pigk(i)=5;% CaI aislado , too overlapped in apex
i=i+1; piUse(i)=0;piquees(i)=CaI;pila(i)=526.56;  piEm(i)=4.877; piGA(i)=1.3E8; pi_z(i)=0;pigk(i)=3;% CaI aislado *
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=551.298; piEm(i)=5.181; piGA(i)=1.1E8; pi_z(i)=0;pigk(i)=1;
i=i+1; piUse(i)=1;piquees(i)=CaI;pila(i)=558.197; piEm(i)=4.7435; piGA(i)=0.42E8; pi_z(i)=0;pigk(i)=7; % selected CaI lines in [Praher2010] below self-absorption threshold, an outlier in BP apex
IRSAC_reference(CaI) = i; % the above line has, according to Praher2010, the lowest selfabsorption coefficient.
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=560.129; piEm(i)=4.7386; piGA(i)=0.43E8; pi_z(i)=0;pigk(i)=5; % selected CaI lines in [Praher2010] below self-absorption threshold  reversed in alternateDelay in Apex
i=i+1; piUse(i)=0*0*1;piquees(i)=CaI;pila(i)=560.285; piEm(i)=4.7353; piGA(i)=0.42E8; pi_z(i)=0;pigk(i)=3; % selected CaI lines in [Praher2010] below self-absorption threshold  reversed in alternateDelay in Apex


% CaII
i=i+1; piUse(i)=1;piquees(i)=CaII;pila(i)=315.89; piEm(i)=7.047; piGA(i)=1.2E9; pi_z(i)=1;pigk(i)=4;piStark(i)=1E16/0.00292; % CaII flojo REF: LIBS++
i=i+1; piUse(i)=1;piquees(i)=CaII;pila(i)=317.93 ; piEm(i)=7.050; piGA(i)=2.2E9; pi_z(i)=1;pigk(i)=6;piStark(i)=1E16/0.00292;% CaII  REF:LIBS++ reversed in alternateDelay
i=i+1; piUse(i)=0;piquees(i)=CaII;pila(i)=318.127; piEm(i)=7.047; piGA(i)=2.3E8; pi_z(i)=1;pigk(i)=4;% CaII , selected CaII line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=CaII;pila(i)=370.603; piEi(i)=3.123;piEm(i)=6.4679; piGA(i)=1.8E8; pi_z(i)=1;pigk(i)=2;piStark(i)=1E16/0.0035; % CaII aislado *
i=i+1; piUse(i)=1;piquees(i)=CaII;pila(i)=373.69 ; piEi(i)=3.151;piEm(i)=6.4679; piGA(i)=3.4E8; pi_z(i)=1;pigk(i)=2;% CaII aislado * , selected CaII line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=CaII;pila(i)=393.36 ; piEi(i)=0.0;  piEm(i)=3.151; piGA(i)=5.88E8; pi_z(i)=1;pigk(i)=4;piStark(i)=1E16/0.00065; % CaII  OJO!! Ek=0 pero hay pocos, reversed in alternateDelay of apex
i=i+1; piUse(i)=0;piquees(i)=CaII;pila(i)=396.84 ; piEi(i)=0.0;  piEm(i)=3.123; piGA(i)=2.8E8; pi_z(i)=1;pigk(i)=2;piStark(i)=1E16/0.00065;% CaII Ek=0 totally reversed in alternateDelay apex
i=i+1; piUse(i)=1;piquees(i)=CaII;pila(i)=849.802 ; piEi(i)=1.692;  piEm(i)=3.151; piGA(i)=4.44E6; pi_z(i)=1;pigk(i)=4;piStark(i)=0;% CaII from ASD-LIBS
%IRSAC_reference(CaII) = i; 


%MgI
i=i+1; piUse(i)=0;piquees(i)=MgI;pila(i)=277.983;  piEi(i)=2.715; piEm(i)=7.173; piGA(i)=8E8; piStark(i)=0; pi_z(i)=0; pigk(i)=4;  %MgI DON'T USE: two overlapped lines, parameters are the mean of the two
i=i+1; piUse(i)=0*1;piquees(i)=MgI;pila(i)=285.2127; piEi(i)=0;piEm(i)=4.3458; piGA(i)=1.47E9; piStark(i)=1.28E17/0.00025; pi_z(i)=0; pigk(i)=3;%MgI REF: ASD
i=i+1; piUse(i)=0;piquees(i)=MgI;pila(i)=333.2146; piEi(i)=2.7116;piEm(i)=6.4314; piGA(i)=3.06E7; piStark(i)=0; pi_z(i)=0; pigk(i)=3;%MgI REF: LIBS++ , selected MgI line in [Praher2010] below self-absorption threshold
IRSAC_reference(MgI) = i; % the 333.21nm line has, according to Praher2010, the lowest selfabsorption coefficient, but is not seen in the limpet shell.
i=i+1; piUse(i)=0;piquees(i)=MgI;pila(i)=333.6674; piEi(i)=2.7166;piEm(i)=6.4314; piGA(i)=5.10E7; piStark(i)=0; pi_z(i)=0; pigk(i)=3;%MgI , selected MgI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=382.9355; piEi(i)=2.7091;piEm(i)=5.9459; piGA(i)=2.70E8; piStark(i)=0; pi_z(i)=0; pigk(i)=3;%MgI , selected MgI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=383.23; piEi(i)=2.711;piEm(i)=5.9459; piGA(i)=6.05e8; pi_z(i)=0; piStark(i)=1E16/0.0107;pigk(i)=5; %MgI
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=383.83; piEi(i)=2.716;piEm(i)=5.9459; piGA(i)=1.13E9; pi_z(i)=0; piStark(i)=1E16/0.0107;pigk(i)=7; %MgI*  , selected MgI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=516.7322; piEi(i)=2.709;piEm(i)=5.1078; piGA(i)=3.39E7; piStark(i)=1E16/0.00328; pi_z(i)=0; pigk(i)=3; %MgI , selected MgI line in [Praher2010] below self-absorption threshold, stark from LIBS++
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=517.2684; piEi(i)=2.7115;piEm(i)=5.1078; piGA(i)=1.019E8; piStark(i)=1E16/0.00328; pi_z(i)=0; pigk(i)=3; %MgI , from ASD-LIBS, STark from LIBS++
i=i+1; piUse(i)=1;piquees(i)=MgI;pila(i)=518.360; piEi(i)=2.7166;piEm(i)=5.1078; piGA(i)=1.68E8; piStark(i)=0; pi_z(i)=0; pigk(i)=3; %MgI , from ASD-LIBS


%MgII
i=i+1; piUse(i)=DISCARD_VAR;piquees(i)=MgII;pila(i)=279.078;piEi(i)=4.422;piEm(i)=8.8637; piGA(i)=1.6e9;piStark(i)=1E16/0.00089; pi_z(i)=1;pigk(i)=4;%MgII  from ASD-LIBS
i=i+1; piUse(i)=0;piquees(i)=MgII;pila(i)=279.55;piEi(i)=0.0;piEm(i)=4.4338; piGA(i)=1.04e9;piStark(i)=1E16/0.000218; pi_z(i)=1;pigk(i)=4;%MgII    Ek=0!!  REF: LIBS++ should check with https://griem.obspm.fr/index.php?page=pages/result.php&element=Mg&base=2   gives too large Ne values (>1E18, surely selfabsorbed)
i=i+1; piUse(i)=DISCARD_VAR;piquees(i)=MgII;pila(i)=279.7998;piEi(i)=4.4338;piEm(i)=8.8637; piGA(i)=2.87e9;piStark(i)=1E16/0.00089; pi_z(i)=1;pigk(i)=6;%MgII  from ASD-LIBS
i=i+1; piUse(i)=0;piquees(i)=MgII;pila(i)=280.2704;piEi(i)=0.0;piEm(i)=4.4224; piGA(i)=5.14E8;piStark(i)=1E16/0.00019; pi_z(i)=1;pigk(i)=2;%Mg Ek=00000!!!!!!  gives too large Ne values (>1E18, surely selfabsorbed)
i=i+1; piUse(i)=DISCARD_NOTSEEN;piquees(i)=MgII;pila(i)=292.8634;piEi(i)=4.4224;piEm(i)=8.6547; piGA(i)=2.3e8;piStark(i)=0; pi_z(i)=1;pigk(i)=2;%MgII  from ASD-LIBS
i=i+1; piUse(i)=DISCARD_NOTSEEN;piquees(i)=MgII;pila(i)=293.6509;piEi(i)=4.4337;piEm(i)=8.6547; piGA(i)=4.6e8;piStark(i)=0; pi_z(i)=1;pigk(i)=2;%MgII  from ASD-LIBS
%IRSAC_reference(MgII) = i;


%SrI
i=i+1; piUse(i)=1;piquees(i)=SrI;pila(i)=330.1734; piEi(i)=1.775;piEm(i)=5.529; piGA(i)=1.8E8; pi_z(i)=0;pigk(i)=3;  % from ASD
i=i+1; piUse(i)=1;piquees(i)=SrI;pila(i)=460.73; piEi(i)=0.0;piEm(i)=2.69; piGA(i)=6.03E8; pi_z(i)=0;pigk(i)=3;  % 460.6nm SrI del artículo de pandhija easily reversed
i=i+1; piUse(i)=1;piquees(i)=SrI;pila(i)=481.18; piEm(i)=4.423; piGA(i)=4.5E8; pi_z(i)=0;pigk(i)=5;  % 481.1nm SrI del artículo de pandhija

%SrII
i=i+1; piUse(i)=0;piquees(i)=SrII;pila(i)=338.0711; piEi(i)=2.940;piEm(i)=6.6066; piGA(i)=0.0000;piStark(i)=0; pi_z(i)=1;pigk(i)=4; %SrII  from ASD, gAk not in ASD, not seen in ASD-LIBS
i=i+1; piUse(i)=1;piquees(i)=SrII;pila(i)=346.4457; piEi(i)=3.04;piEm(i)=6.6174; piGA(i)=1.9E9;piStark(i)=0; pi_z(i)=1;pigk(i)=6; %SrII  from ASD-LIBS
i=i+1; piUse(i)=1;piquees(i)=SrII;pila(i)=407.7; piEi(i)=0.0;piEm(i)=3.04; piGA(i)=5.64E8;piStark(i)=1E16/0.0016; pi_z(i)=1;pigk(i)=4; %SrII  REF: LIBS++ EASILY REVERSED
i=i+1; piUse(i)=1;piquees(i)=SrII;pila(i)=416.1796; piEi(i)=2.940;piEm(i)=5.9186; piGA(i)=1.3E8;piStark(i)=1E16/0.0016; pi_z(i)=1;pigk(i)=2; %SrII  from ASD-LIBS
i=i+1; piUse(i)=0;piquees(i)=SrII;pila(i)=421.5524; piEi(i)=0.0;piEm(i)=2.9403; piGA(i)=2.55E8;piStark(i)=1E16/0.0016; pi_z(i)=1;pigk(i)=2; %SrII  from ASD-LIBS Ek=0


%ZnI 
i=i+1; piUse(i)=0;piquees(i)=ZnI;pila(i)= 213.857; piEi(i)=0.0; piEm(i)=5.796; piGA(i)=2.14E9; piStark(i)=0; pi_z(i)=0;pigk(i)=3; piResonant(i)=1;  % highest in ASD-LIBS, interference with left peak,discarded, too much interference, Ei=0, but using it reduces the error!
%i=i+1; piUse(i)=0;piquees(i)=ZnI;pila(i)= 280.086; piEm(i)=8.503; piGA(i)=1E9; piStark(i)=0; pi_z(i)=0;pigk(i)=7;  % DO NOT USE: no Ak in ASD or elsewhere, arbitrary value just to try
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 307.59; piEi(i)=0.0;piEm(i)=4.03; piGA(i)=1.1E5; piStark(i)=0; pi_z(i)=0;pigk(i)=3;  % ASD, strong in brass spectra, Ei=0 but low transition probability
IRSAC_reference(ZnI) = i; % chosen due to low gAk
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 328.23; piEi(i)=4.006; piEm(i)=7.782; piGA(i)=2.7E8; piStark(i)=0; pi_z(i)=0;pigk(i)=3;  % from colao 2004 (too much background in coin spectra, though) removed due to interference of a peak on the left (coin)
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 330.26; piEi(i)=4.029; piEm(i)=7.7827; piGA(i)=6.0E8; piStark(i)=1E16/0.007; pi_z(i)=0;pigk(i)=5;  % from LIBS++ % three overlapped lines? - according to ASD-LIBS, there is no interferences, rechecked: if used, everything burn into flames; zhao2018 uses these two lines, gAk=5.35E8
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 334.50; piEi(i)=4.077; piEm(i)=7.7833; piGA(i)=1.19E9; piStark(i)=1E16/0.0087; pi_z(i)=0;pigk(i)=7;  % from LIBS++ - according to ASD-LIBS, there is no interferences, Zhao2018 uses these two lines, gAk=1.05E9
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 468.014; piEi(i)=4.006; piEm(i)=6.655; piGA(i)=3*1.4E7; piStark(i)=1E16/0.0014; pi_z(i)=0;pigk(i)=3;  % Ak from GENIE, added again because gA value was wrong; not seen in ASD-LIBS nor nordic gold
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 472.21; piEi(i)=4.029; piEm(i)=6.6545; piGA(i)=3*4.2E7; piStark(i)=1E16/0.001335; pi_z(i)=0;pigk(i)=3;  % no Ak in ASD, source for alternate Aki: https://www-amdis.iaea.org/GENIE/ ; not seen in ASD-LIBS nor nordic gold
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 481.053; piEi(i)=4.077; piEm(i)=6.6545; piGA(i)=3*7.0E7;piStark(i)=1E16/0.001022; pi_z(i)=0;pigk(i)=3;  %   no Ak in ASD, source for alternate Aki: https://www-amdis.iaea.org/GENIE/; not seen in ASD-LIBS nor nordic gold
i=i+1; piUse(i)=1;piquees(i)=ZnI;pila(i)= 636.24; piEi(i)=5.796; piEm(i)=7.7438; piGA(i)=2.4E8; piStark(i)=0; pi_z(i)=0;pigk(i)=5;  % from ASD larger intensities, it shows up at 636.15nm but there is no other possible line in this wavelength range.; not seen in ASD-LIBS but seen clearly in nordic gold

%ZnII
i=i+1; piUse(i)=0;piquees(i)=ZnII;pila(i)= 206.2; piEm(i)=6.0108; piGA(i)=1.32E9; pi_z(i)=1;pigk(i)=4;  %  discarded <230nm too noisy
i=i+1; piUse(i)=0;piquees(i)=ZnII;pila(i)= 202.548; piEm(i)=6.119; piGA(i)=1.63E9; pi_z(i)=1;pigk(i)=4;  % ASD larger relative intensities  , discarded <230nm too noisy
i=i+1; piUse(i)=1;piquees(i)=ZnII;pila(i)= 250.20; piEm(i)=10.965; piGA(i)=3.94E8; pi_z(i)=1;pigk(i)=2;  % ASD larger relative intensities   removed because fitting is wrong, could be used if another window/baseline processing fix that.
IRSAC_reference(ZnII) = i; % the above line has been selected due to low Ak and high Ek,
i=i+1; piUse(i)=0;piquees(i)=ZnII;pila(i)= 255.795; piEm(i)=10.965; piGA(i)=7.82E8; pi_z(i)=1;pigk(i)=2;  % ASD larger relative intensities  
i=i+1; piUse(i)=1;piquees(i)=ZnII;pila(i)= 491.16; piEm(i)=14.539; piGA(i)=1.09E9; pi_z(i)=1;pigk(i)=6;  %  wide weird peak
i=i+1; piUse(i)=0;piquees(i)=ZnII;pila(i)= 492.401; piEm(i)=14.538; piGA(i)=2.18E9; pi_z(i)=1;pigk(i)=8;  %  seen in nordic gold
i=i+1; piUse(i)=0;piquees(i)=ZnII;pila(i)= 589.43; piEm(i)=8.114; piGA(i)=0; pi_z(i)=1;pigk(i)=4;  %  Ak is not in ASD, not found elsewhere

%CuI  i=64
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)=216.509; piEi(i)=0.0;piEm(i)=5.724; piGA(i)=2.2E8; piStark(i)=0; pi_z(i)=0;pigk(i)=4; piResonant(i)=1;  % from ADS-LIBS, gives a bad vertical value,but i do not know why... wait! Ek=0!!!
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)=217.894; piEi(i)=0.0 ;piEm(i)=5.68; piGA(i)=2.0E8; piStark(i)=0; pi_z(i)=0;pigk(i)=4; piResonant(i)=1;  % from ADS-LIBS , discarded: overlapped with CuII
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)=219.975; piEi(i)=1.6422;piEm(i)=7.2768; piGA(i)=1.2564E9; piStark(i)=0; pi_z(i)=0;pigk(i)=4;piResonant(i)=1;  % reference line in [sun2009] Ak from LIBS++ (units of 1E8). overlapped with 219.959
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)=261.837; piEi(i)=1.389;piEm(i)=6.123; piGA(i)=1.23E8; piStark(i)=0; pi_z(i)=0;pigk(i)=4;  % ASD, strong in brass spectra
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)=276.64; piEi(i)=1.642;piEm(i)=6.123; piGA(i)=3.8E7; piStark(i)=0; pi_z(i)=0;pigk(i)=4;piResonant(i)=1;  % ASD, strong in brass spectra
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 306.341; piEi(i)=1.642 ;piEm(i)=5.688; piGA(i)=6.20E6; piStark(i)=0; pi_z(i)=0;pigk(i)=4;  % ASD
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 319.410; piEi(i)=1.642 ;piEm(i)=5.5228; piGA(i)=6.2E6; piStark(i)=0; pi_z(i)=0;pigk(i)=4;  % ASD
i=i+1; piUse(i)=0*1;piquees(i)=CuI;pila(i)= 324.7540;piEi(i)=0.0000000;piEm(i)=3.8166920;piGA(i)=5.580e+08;piAcc{i}='AA';pi_z(i)=0;pigk(i)=4;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1; piUse(i)=0*1;piquees(i)=CuI;pila(i)= 327.395; piEi(i)=1.642 ;piEm(i)=3.7859; piGA(i)=2.752E8; piStark(i)=0; pi_z(i)=0;pigk(i)=2;piResonant(i)=1;   % from Tognoni 2007 Ek=0!!! discarded ad-hoc becouse got worse Te
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 333.785; piEi(i)=1.389 ;piEm(i)=5.102; piGA(i)=3.0E6; piStark(i)=0; pi_z(i)=0;pigk(i)=8;  % ASD
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 402.26; piEi(i)=3.786 ;piEm(i)=6.867; piGA(i)=7.60E7; piStark(i)=0; pi_z(i)=0;pigk(i)=4;  % ASD
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 406.264; piEi(i)=3.817;piEm(i)=6.868; piGA(i)=1.26E8; piStark(i)=0; pi_z(i)=0;pigk(i)=6;  % ASD
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 424.896; piEi(i)=5.076;piEm(i)=7.994; piGA(i)=3.9E7; piStark(i)=0; pi_z(i)=0;pigk(i)=2;  % ASD, strong in brass spectra
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 427.51; piEi(i)=4.838;piEm(i)=7.73; piGA(i)=2.76E8; piStark(i)=1E16/0.004; pi_z(i)=0;pigk(i)=8;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 448.03; piEi(i)=3.786;piEm(i)=6.55; piGA(i)=6.0E6; piStark(i)=1E16/0.012; pi_z(i)=0;pigk(i)=2;  % from LIBS++%
%IRSAC_reference(CuI) = i; % the above line has been selected due to low Ak and high Ek, the reference in [sun2009] @ 219.98nm is not seen 
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 453.07; piEi(i)=3.817;piEm(i)=6.55; piGA(i)=1.7E7; piStark(i)=1E16/0.011; pi_z(i)=0;pigk(i)=2;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 453.96; piEi(i)=5.153;piEm(i)=7.88; piGA(i)=8.48E7; piStark(i)=1E16/0.09; pi_z(i)=0;pigk(i)=4;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 458.697; piEi(i)=5.102;piEm(i)=7.805; piGA(i)=1.92E8; piStark(i)=1E16/0.07; pi_z(i)=0;pigk(i)=6;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 465.11; piEi(i)=5.072;piEm(i)=7.737; piGA(i)=3.04E8; piStark(i)=1E16/0.00435; pi_z(i)=0;pigk(i)=8;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 470.46; piEi(i)=5.102;piEm(i)=7.737; piGA(i)=4.04E7; piStark(i)=0; pi_z(i)=0;pigk(i)=8;  % ASD, seen in nordic gold

i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 501.66; piEi(i)=5.522;piEm(i)=7.993; piGA(i)=2*0.1618E8; piStark(i)=0; pi_z(i)=0;pigk(i)=2;  % from LIBS++
i=i+1; piUse(i)=0*1;piquees(i)=CuI;pila(i)= 510.55; piEi(i)=1.389;piEm(i)=3.8167; piGA(i)=8.0E6; piStark(i)=1E16/0.00215; pi_z(i)=0;pigk(i)=4;  % from LIBS++
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 511.191; piEi(i)=5.569;piEm(i)=7.993; piGA(i)=2*0.1049E8; piStark(i)=0; pi_z(i)=0;pigk(i)=2;  % from ASD, gAk from LIBS++
i=i+1; piUse(i)=0*1;piquees(i)=CuI;pila(i)= 515.32; piEi(i)=3.786;piEm(i)=6.1911; piGA(i)=2.46E8; piStark(i)=1E16/0.0095; pi_z(i)=0;pigk(i)=4;  % 
i=i+1; piUse(i)=1;piquees(i)=CuI;pila(i)= 521.8202; piEi(i)=3.816692; piEm(i)=6.1920251; piGA(i)=4.5e+08; piAcc{i}='C+'; pi_z(i)=0; pigk(i)=6; piIntens(i)=2500;piStark(i)=0; % auto imported from NIST ASD Intensity>400
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 529.25; piEi(i)=5.395;piEm(i)=7.737; piGA(i)=8.72E7; piStark(i)=0; pi_z(i)=0;pigk(i)=8;  % ASD, seen in nordic gold coin with DP
IRSAC_reference(CuI) = i;
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 555.494; piEi(i)=5.506;piEm(i)=7.737; piGA(i)=8*0.006201E8; piStark(i)=0; pi_z(i)=0;pigk(i)=8;  % ASD, seen in nordic gold coin with DP
i=i+1; piUse(i)=0;piquees(i)=CuI;pila(i)= 570.024; piEi(i)=1.642;piEm(i)=3.8167; piGA(i)=9.6E5; piStark(i)=0; pi_z(i)=0;pigk(i)=4;  % ASD ; REMOVED DUE TO MISCALIBRATION OF THAT CHANNEL, PEAK FITTING IS WRONG
i=i+1; piUse(i)=0*1;piquees(i)=CuI;pila(i)= 578.21; piEi(i)=1.642;piEm(i)=3.7859; piGA(i)=3.3E6; piStark(i)=1E16/0.0036; pi_z(i)=0;pigk(i)=2;  % from LIBS++;  REMOVED DUE TO MISCALIBRATION OF THAT CHANNEL, PEAK FITTING IS WRONG
%i=i+1; piUse(i)=1;piquees(i)=CuI;pila(i)=793.31;%piEm(i)=5.3483;piGA(i)=???; piStark(i)=1E16/0.016;pi_z(i)=0;pigk(i)=2;  %%from LIBS++, gk not known; 


%CuII we need lines above 250nm. according to ASD-LIBS for Cu, no CuII
%above 248nm
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 204.38; piEm(i)=8.783; piGA(i)=9.94E8; pi_z(i)=1;pigk(i)=7;  % from ASD-LIBS, not too much interference, 
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 205.498; piEm(i)=8.864; piGA(i)=8.25E8; pi_z(i)=1;pigk(i)=5;  % from ASD-LIBS, not too much interference, 
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 207.86; piEm(i)=14.20; piGA(i)=1.9E9;  pi_z(i)=1;pigk(i)=3;  % from ASD larger intensities, discarded <230nm too noisy
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 211.21; piEm(i)=9.125; piGA(i)=1.0E9;  pi_z(i)=1;pigk(i)=3;  % from ASD larger intensities, discarded <230nm too noisy
i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 213.598; piEm(i)=8.52; piGA(i)=4.13E9; pi_z(i)=1;pigk(i)=9;  % from ASD-LIBS, not too much interference, overlapped with ZnI
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 214.898; piEm(i)=8.486; piGA(i)=6.2E8; pi_z(i)=1;pigk(i)=7;  % from ASD-LIBS, not too much interference, 
i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 219.22; piEm(i)=8.486; piGA(i)=2.0E9;  pi_z(i)=1;pigk(i)=7;  % from ASD-LIBS, not too much interference, but adding this peak got worse results
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 221.81; piEm(i)=8.42; piGA(i)=1.0E9;   pi_z(i)=1;pigk(i)=3;  % from ASD-LIBS, not too much interference, but adding this peak got worse results
i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 224.70; piEm(i)=8.235; piGA(i)=1.6E9;  pi_z(i)=1;pigk(i)=5;  % from ASD larger intensities, removed due to interference with peak to the left, could be used if another window/baseline processing fix that.
i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 236.99; piEi(i)=3.256;piEm(i)=8.486; piGA(i)=3.7E8;  pi_z(i)=1;pigk(i)=7;  % ASD, in coind with DP, I think the peak is 236.93 of AlI
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 240.01; piEi(i)=3.256;piEm(i)=8.420; piGA(i)=2.1E7;  pi_z(i)=1;pigk(i)=3;  % ASD, seen in coin spectra with DP, no interference verified
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 240.33; piEi(i)=8.2349;piEm(i)=13.392; piGA(i)=7.4E8;  pi_z(i)=1;pigk(i)=7;  % ASD, seen in coin spectra with DP, no interference verified
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 248.58; piEi(i)=8.6625;piEm(i)=13.649; piGA(i)=4.26E8;  pi_z(i)=1;pigk(i)=3;  % ASD, seen in coin spectra with DP, no interference verified
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 248.97; piEi(i)=3.2564;piEm(i)=8.235; piGA(i)=5.0E6;  pi_z(i)=1;pigk(i)=5;  % ASD, seen in coin spectra with DP, no interference verified
IRSAC_reference(CuII) = i; % chosen randomly, check
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 250.63; piEm(i)=13.43; piGA(i)=1.0E9;  pi_z(i)=1;pigk(i)=5;  % from ASD larger intensities, seen in coin spectra, verified no interferences
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 254.48; piEm(i)=13.39; piGA(i)=1.36E9; pi_z(i)=1;pigk(i)=7;  % from ASD larger intensities, checked no interferences with ASD-LIBS and coin spectra

i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=273.97662;piEi(i)=9.12471797;piEm(i)=13.64873556;piGA(i)=1.2e+07;piAcc{i}='D';pi_z(i)=1;pigk(i)=3;piIntens(i)=86000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=274.52710;piEi(i)=8.91695539;piEm(i)=13.43190338;piGA(i)=7.00e+07;piAcc{i}='C+';pi_z(i)=1;pigk(i)=5;piIntens(i)=140000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=276.96690;piEi(i)=8.91695539;piEm(i)=13.39213290;piGA(i)=4.7e+08;piAcc{i}='C+';pi_z(i)=1;pigk(i)=7;piIntens(i)=310000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=283.73682;piEi(i)=9.06349786;piEm(i)=13.43190338;piGA(i)=1.2e+08;piAcc{i}='C+';pi_z(i)=1;pigk(i)=5;piIntens(i)=150000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=287.76996;piEi(i)=9.12471797;piEm(i)=13.43190338;piGA(i)=1.2e+08;piAcc{i}='C+';pi_z(i)=1;pigk(i)=5;piIntens(i)=170000;piStark(i)=0;% auto imported from NIST ASD Intensity>400

i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 329.04; piEi(i)=14.328;piEm(i)=18.096; piGA(i)=7.7E8; pi_z(i)=1;pigk(i)=13;  % ASD
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 336.56; piEi(i)=14.423;piEm(i)=18.105; piGA(i)=2.6E8; pi_z(i)=1;pigk(i)=9;  % ASD,
i=i+1; piUse(i)=0;piquees(i)=CuII;pila(i)= 338.071; piEi(i)=14.696;piEm(i)=18.362; piGA(i)=2.4E8; pi_z(i)=1;pigk(i)=7;  % ASD, seen in nordic gold
i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 422.7; piEi(i)=14.423;piEm(i)=18.105; piGA(i)=2.6E8; pi_z(i)=1;pigk(i)=9;  % ASD, not seen, removed

i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=490.973351;piEi(i)=14.32872908;piEm(i)=16.85329769;piGA(i)=2.65e+09;piAcc{i}='B+';pi_z(i)=1;pigk(i)=13;piIntens(i)=160000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=491.83779;piEi(i)=14.59881073;piEm(i)=17.1189423;piGA(i)=2.6e+09;piAcc{i}='C';pi_z(i)=1;pigk(i)=9;piIntens(i)=54000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=493.16982;piEi(i)=14.34032994;piEm(i)=16.85365487;piGA(i)=2.1e+09;piAcc{i}='C';pi_z(i)=1;pigk(i)=11;piIntens(i)=140000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=495.37244;piEi(i)=14.61564135;piEm(i)=17.1177911;piGA(i)=3.4e+09;piAcc{i}='C';pi_z(i)=1;pigk(i)=11;piIntens(i)=82000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=498.550499;piEi(i)=14.39211337;piEm(i)=16.87831320;piGA(i)=8.7e+08;piAcc{i}='C+';pi_z(i)=1;pigk(i)=9;piIntens(i)=70000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=505.179210;piEi(i)=14.42818740;piEm(i)=16.88176497;piGA(i)=1.70e+09;piAcc{i}='C+';pi_z(i)=1;pigk(i)=11;piIntens(i)=120000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=506.545858;piEi(i)=14.69012052;piEm(i)=17.13707847;piGA(i)=1.45e+09;piAcc{i}='C+';pi_z(i)=1;pigk(i)=9;piIntens(i)=70000;piStark(i)=0;% auto imported from NIST ASD Intensity>400

i=i+1; piUse(i)=0*1;piquees(i)=CuII;pila(i)= 715.776; piEi(i)=13.392;piEm(i)=15.124; piGA(i)=2.6E6; pi_z(i)=1;pigk(i)=5;  % ASD, seen in nordic gold

i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=740.43561;piEi(i)=14.88955769;piEm(i)=16.56357367;piGA(i)=1.4e+08;piAcc{i}='D+';pi_z(i)=1;pigk(i)=7;piIntens(i)=55000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0*1;piquees(i)=CuII;pila(i)=766.46451;piEi(i)=14.96299725;piEm(i)=16.58016355;piGA(i)=1.8e+08;piAcc{i}='D+';pi_z(i)=1;pigk(i)=5;piIntens(i)=52000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=780.76534;piEi(i)=14.97602748;piEm(i)=16.56357367;piGA(i)=2.5e+08;piAcc{i}='C+';pi_z(i)=1;pigk(i)=7;piIntens(i)=82000;piStark(i)=0;% auto imported from NIST ASD Intensity>400
i=i+1;piUse(i)=0;piquees(i)=CuII;pila(i)=782.56528;piEi(i)=13.39213290;piEm(i)=14.97602748;piGA(i)=4.7e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=9;piIntens(i)=59000;piStark(i)=0;% auto imported from NIST ASD Intensity>400



%AlI
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 237.31;piEi(i)=0.014; piEm(i)=5.237; piGA(i)=5.44E8;pi_z(i)=0;pigk(i)=6; piStark(i)=0;%  from ASD-LIBS
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 256.798; piEi(i)=0.0;piEm(i)=4.826; piGA(i)=7.68E7;pi_z(i)=0;pigk(i)=4; piStark(i)=0;%  from ASD-LIBS larger intensities of Al-I
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 257.509; piEi(i)=0.014;piEm(i)=4.827; piGA(i)=2.16E8;pi_z(i)=0;pigk(i)=6; piStark(i)=0;%  from ASD-LIBS larger intensities of Al-I
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 265.24; piEi(i)=0;piEm(i)=4.673; piGA(i)=2.84E7;pi_z(i)=0;pigk(i)=2; piStark(i)=0;%  from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 266.04; piEi(i)=0.014;piEm(i)=4.673; piGA(i)=5.68E7;pi_z(i)=0;pigk(i)=2; piStark(i)=0;%  from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 305.47; piEm(i)=7.6557; piGA(i)=7.08E5;pi_z(i)=0;pigk(i)=4; piStark(i)=0;%  reference line for IRSAC [sun2009], gkAk value from GENIE, not sure it is Ok, discarded due to strong peak to the right 
IRSAC_reference(AlI) = i; % it is the line with highest coeff in the first trial, it could be the one with lowest SA;  lowest gA, reference line in [sun2009] is 305.47 which dos not have gAk value
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 308.21;piEi(i)=0; piEm(i)=4.022; piGA(i)=2.35E8;pi_z(i)=0;pigk(i)=4; piStark(i)=1E16/0.00225; %  from ASD larger intensities, stark from LIBS++
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 309.27;piEi(i)=0.014; piEm(i)=4.022; piGA(i)=4.37E8;pi_z(i)=0;pigk(i)=6; piStark(i)=0;%  from ASD larger intensities discarded, are two peaks, simmilar intensity,the second one with Ei aprox. =0
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 394.35; piEi(i)=0;piEm(i)=3.1427; piGA(i)=9.98E7;pi_z(i)=0;pigk(i)=2; piStark(i)=1E16/0.00165;%  LIBS++, fitting fails, changed wavelength to 394.35, warning! Ne from stark of this peak is an outlier, check
i=i+1; piUse(i)=1;piquees(i)=AlI;pila(i)= 396.15; piEi(i)=0.014;piEm(i)=3.1427; piGA(i)=1.97E8;pi_z(i)=0;pigk(i)=2; piStark(i)=1E16/0.00165;%  LIBS++  warning! Ne from stark of this peak is an outlier, check
%IRSAC_reference(AlI) = 0; % the only possible line is 305.47 but gives way too low values of SA coefficient


%AlII
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=266.916; piEi(i)=0.0;piEm(i)=4.644; piGA(i)=9.84E3; pi_z(i)=1;pigk(i)=3; PiStark(i)=0;  % from ASD-Libs larger intensities, not seen in nordic gold
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=281.62; piEi(i)=7.42;piEm(i)=11.822; piGA(i)=3.57E8; pi_z(i)=1;pigk(i)=1; PiStark(i)=1E16/0.00212;  % El Sherbini 2006, stark from LIBS++
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=358.707; piEi(i)=11.84;piEm(i)=15.302; piGA(i)=1.46E9; pi_z(i)=1;pigk(i)=7; PiStark(i)=0;  % from ASD-Libs larger intensities, not seen in nordic gold, mix of several lines with simmilar parameters
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=466.31; piEi(i)=10.60;piEm(i)=13.256; piGA(i)=1.74E8; pi_z(i)=1;pigk(i)=3; PiStark(i)=1E16/0.003015; %  El Sherbini 2006, stark from LIBS++
IRSAC_reference(AlII) = i; % not many to choose from
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=704.206; piEi(i)=11.317;piEm(i)=13.077; piGA(i)=2.89E8; pi_z(i)=1;pigk(i)=5; PiStark(i)=0; %  from ASD-LIBS larger intensities
i=i+1; piUse(i)=1;piquees(i)=AlII;pila(i)=705.66; piEi(i)=11.317;piEm(i)=13.073; piGA(i)=1.72E8; pi_z(i)=1;pigk(i)=3; PiStark(i)=0; %  from ASD-LIBS larger intensities

%SnI
i=i+1; piUse(i)=0;piquees(i)=SnI;pila(i)=235.48;piEm(i)=5.47; piGA(i)=8.5E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  piResonant(i)=1; %  from ASD larger intensities, checked no interferences with ASD_LIBS, totally overlapped
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=242.16;piEm(i)=6.186; piGA(i)=1.8E9;pi_z(i)=0;pigk(i)=7; piStark(i)=0; piResonant(i)=1;  %  from ASD larger intensities, checked no interferences with ASD_LIBS
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=242.95;piEm(i)=5.53; piGA(i)=1E9;pi_z(i)=0;pigk(i)=7; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, checked no interferences with ASD_LIBS
i=i+1; piUse(i)=1*1;piquees(i)=SnI;pila(i)=249.57;piEm(i)=6.03; piGA(i)=3.1E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;%  from ASD larger intensities, checked with ASD-LIBS: interference with CuI, discarded
i=i+1; piUse(i)=1*1;piquees(i)=SnI;pila(i)=266.12;piEm(i)=4.87; piGA(i)=3.3E7;pi_z(i)=0;pigk(i)=3; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, checked no interferences with ASD-LIBS, but was an outlier, removing this peak overall error decrease a lot
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=283.998;piEm(i)=4.789;piGA(i)=8.58E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, peaks not confirmed, it seems isolated , checked with ASD-LIBS; gAk is 1.7E8 in GENIE and 1.57 in Griffoni 2016 ¿?
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=286.33;piEm(i)=4.329;piGA(i)=1.6E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, peaks not confirmed, it seems isolated , checked
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=317.50;piEm(i)=4.329;piGA(i)=3.0E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, peaks not confirmed, it seems isolated , checked
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=326.23;piEm(i)=4.87;piGA(i)=8.1E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0; piResonant(i)=1; %  from ASD larger intensities, peaks not confirmed, it seems isolated , does not appears in ASD-LIBS
i=i+1; piUse(i)=1*0;piquees(i)=SnI;pila(i)=380.10;piEm(i)=4.329;piGA(i)=8.4E7;pi_z(i)=0;pigk(i)=3; piStark(i)=0; piResonant(i)=1;%  from ASD larger intensities, peaks not confirmed, it seems isolated , checked
i=i+1; piUse(i)=1;piquees(i)=SnI;pila(i)=452.47;piEm(i)=4.867;piGA(i)=7.8E7;pi_z(i)=0;pigk(i)=3; piStark(i)=0;%  from ASD larger intensities, peaks not confirmed, it seems isolated , does not appear in ASD-LIBS
IRSAC_reference(SnI) = i; % most are resonant


%SnII - none with enough intensity in the coin spectrum
i=i+1; piUse(i)=0;piquees(i)=SnII;pila(i)=328.31;piEm(i)=11.07; piGA(i)=1.02E9;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from ASD larger intensities, could be CuII @328.32, too much background and interferences, IT'S ZnI
IRSAC_reference(SnII) = 0; % not many to choose from

%VI - Vanadium
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=292.362;piEi(i)=0.0686;piEm(i)=4.308; piGA(i)=6.6E8;pi_z(i)=0;pigk(i)=8;  % from Cremers book p.264
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=318.396;piEi(i)=0.04;piEm(i)=3.933; piGA(i)=2.8E9;pi_z(i)=0;pigk(i)=10;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=370.357;piEi(i)=0.3;piEm(i)=3.647; piGA(i)=8.96E8;pi_z(i)=0;pigk(i)=8;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=411.178;piEi(i)=0.3;piEm(i)=3.315; piGA(i)=1.00E9;pi_z(i)=0;pigk(i)=10;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=438.47;piEi(i)=0.287;piEm(i)=3.113; piGA(i)=9.2E8;pi_z(i)=0;pigk(i)=10;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=438.998;piEi(i)=0.275;piEm(i)=3.098; piGA(i)=9.2E8;pi_z(i)=0;pigk(i)=10;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=437.92;piEi(i)=0.3;piEm(i)=3.131; piGA(i)=1.38E9;pi_z(i)=0;pigk(i)=12;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=446.03;piEi(i)=0.3;piEm(i)=3.079; piGA(i)=2.09E8;pi_z(i)=0;pigk(i)=8;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=459.41;piEi(i)=0.06;piEm(i)=2.766; piGA(i)=6.8E7;pi_z(i)=0;pigk(i)=12;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=488.156;piEi(i)=0.068;piEm(i)=2.607; piGA(i)=6.2E7;pi_z(i)=0;pigk(i)=8;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=572.705;piEi(i)=1.08;piEm(i)=3.244; piGA(i)=1.94E8;pi_z(i)=0;pigk(i)=10;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VI;pila(i)=624.31;piEi(i)=0.3;piEm(i)=2.286; piGA(i)=1.94E7;pi_z(i)=0;pigk(i)=10;  % from Cremers book


%VII 
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=289.332;piEi(i)=0.368;piEm(i)=4.652; piGA(i)=8.61E8;pi_z(i)=1;pigk(i)=7;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=309.31;piEi(i)=0.392;piEm(i)=4.399; piGA(i)=2.6E9;pi_z(i)=1;pigk(i)=13;  % from Cremers book
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=310.23;piEi(i)=0.368;piEm(i)=4.363; piGA(i)=1.96E9;pi_z(i)=1;pigk(i)=11;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=311.07;piEi(i)=0.348;piEm(i)=4.332; piGA(i)=1.42E9;pi_z(i)=1;pigk(i)=9;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=311.838;piEi(i)=0.333;piEm(i)=4.307; piGA(i)=1.03E9;pi_z(i)=1;pigk(i)=7;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=312.528;piEi(i)=0.323;piEm(i)=4.289; piGA(i)=7.50E8;pi_z(i)=1;pigk(i)=5;  % from ASD larger intensities
i=i+1; piUse(i)=1;piquees(i)=VII;pila(i)=327.61;piEi(i)=1.128;piEm(i)=4.911; piGA(i)=1.87E9;pi_z(i)=1;pigk(i)=11;  % from Cremers book

%FeI -- Iron
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=346.59;piEm(i)=3.686; piGA(i)=0.119E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0; % Shah 2012, piEi=0.110,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=355.49;piEm(i)=6.319; piGA(i)=1.40E8;pi_z(i)=0;pigk(i)=13; piStark(i)=0; % Shah 2012, piEi=2.832,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=361.88;piEm(i)=4.415; piGA(i)=0.722E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % Shah 2012, piEi=0.990,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=367.9913;piEm(i)=3.36826; piGA(i)=1.24E7;pi_z(i)=0;pigk(i)=9; piStark(i)=0; % , selected FeI line in [Praher2010] below self-absorption threshold, Ei=0!!!
IRSAC_reference(FeI) = i; % 367.9913 lowest selfabsorption in Praher2010
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=368.7456;piEm(i)=4.2203; piGA(i)=7.2E7;pi_z(i)=0;pigk(i)=9; piStark(i)=0; % , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=370.556;piEm(i)=3.39; piGA(i)=0.0322E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.052, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=370.925;piEm(i)=4.25; piGA(i)=0.156E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.91, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels  , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=371.994;piEm(i)=3.33; piGA(i)=0.162E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.00, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=372.762;piEm(i)=4.28; piGA(i)=0.225E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.96, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels, , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=373.486;piEm(i)=4.18; piGA(i)=0.901E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0;  %  from Aragon 2014, piEi=0.86
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=373.713;piEm(i)=3.36; piGA(i)=0.141E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.052, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=374.3362;piEm(i)=4.301278; piGA(i)=7.80E7;pi_z(i)=0;pigk(i)=3; piStark(i)=0; % , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=374.556;piEm(i)=3.39; piGA(i)=0.115E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.087, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=374.589;piEm(i)=3.43; piGA(i)=0.0733E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.12, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=374.826;piEm(i)=3.41; piGA(i)=0.0915E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  S.M. Pershin 2012, piEi=0.11, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=374.948;piEm(i)=4.22; piGA(i)=0.763E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0;  %  from Aragon 2014, piEi=0.91
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=375.823;piEm(i)=4.26; piGA(i)=0.634E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  from Aragon 2014, piEi=0.96
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=376.3789;piEm(i)=4.28; piGA(i)=0.544E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  from Aragon 2014, piEi=0.99
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=376.5539;piEm(i)=6.53; piGA(i)=0.951E8;pi_z(i)=0;pigk(i)=15; piStark(i)=0;  %  from Aragon 2014, piEi=3.24 ,pigi=13
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=376.7192;piEm(i)=4.30; piGA(i)=0.639E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0;  %  from Aragon 2014, piEi=1.01 ,
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=379.5002;piEm(i)=4.26; piGA(i)=0.115E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  from Aragon 2014, piEi=0.99 ,pigi=5
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=379.7515;piEm(i)=6.5; piGA(i)=0.457E8;pi_z(i)=0;pigk(i)=13; piStark(i)=0;  %  from Aragon 2014, piEi=3.24 ,
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=380.5343;piEm(i)=6.56; piGA(i)=0.860E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0;  %  from Aragon 2014, piEi=3.30 ,pigi=9
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=381.584;piEm(i)=4.73; piGA(i)=1.12E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0;  %  from Aragon 2014, piEi=1.48 ,pigi=9
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=382.0425;piEm(i)=4.10; piGA(i)=0.667E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0;  %  from Aragon 2014, piEi=0.86 ,pigi=11
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=382.1178;piEm(i)=6.51; piGA(i)=0.554E8;pi_z(i)=0;pigk(i)=13; piStark(i)=0;  %  from Aragon 2014, piEi=3.27 ,pigi=11
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=382.4444;piEm(i)=3.24097; piGA(i)=1.98E7;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=382.7823;piEm(i)=4.80; piGA(i)=1.05E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  from Aragon 2014, piEi=1.56 ,pigi=7
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=383.4222;piEm(i)=4.19; piGA(i)=0.452E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  from Aragon 2014, piEi=0.96 ,pigi=7
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=384.1048;piEm(i)=4.83; piGA(i)=1.36E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0;  %  from Aragon 2014, piEi=1.61 ,pigi=5
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=385.6372;piEm(i)=3.27; piGA(i)=0.0464E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  from Aragon 2014, piEi=0.05 ,pigi=7
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=385.9212;piEm(i)=5.62; piGA(i)=0.0725E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0;  %  from Aragon 2014, piEi=2.40 ,pigi=13
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=385.9911;piEm(i)=3.21; piGA(i)=0.0969E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0;  %  from Aragon 2014, piEi=0.00 
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=386.5523;piEm(i)=4.22; piGA(i)=0.155E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0;  %  from Aragon 2014, piEi=1.01 
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=387.2501;piEm(i)=4.19; piGA(i)=0.105E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0;  %  from Aragon 2014, piEi=0.99
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=387.3760;piEm(i)=5.63; piGA(i)=0.0657E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0;  %  from Aragon 2014, piEi=2.43 , pigi=11
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=389.9707;piEm(i)=3.2657; piGA(i)=1.29E7;pi_z(i)=0;pigk(i)=5; piStark(i)=0; % , selected FeI line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=1;piquees(i)=FeI; pila(i)=396.93;piEm(i)=4.608; piGA(i)=0.226E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % Shah 2012, piEi=1.485,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=425.08;piEm(i)=4.473; piGA(i)=0.102E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % Shah 2012, piEi=1.557,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=426.05;piEm(i)=5.308; piGA(i)=0.399E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0; % Shah 2012, piEi=2.399,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=432.58;piEm(i)=4.473; piGA(i)=0.516E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % Shah 2012, piEi=1.608,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeI; pila(i)=492.05;piEm(i)=5.532; piGA(i)=0.358E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0; % Shah 2012, piEi=2.832,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy

%FeII
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=233.28;piEm(i)=5.361; piGA(i)=1.31E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0; % Shah 2012, piEi=0.048,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=233.80;piEm(i)=5.408; piGA(i)=1.13E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0; % Shah 2012, piEi=0.107,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=236.8596;piEm(i)=5.58; piGA(i)=0.606E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0;  %  from Aragon 2014, piEi=0.35, pigi=6
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=237.3736;piEm(i)=5.22; piGA(i)=0.425E8;pi_z(i)=1;pigk(i)=10; piStark(i)=0;  %  from Aragon 2014, piEi=0.00, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=237.5193;piEm(i)=5.60; piGA(i)=0.981E8;pi_z(i)=1;pigk(i)=2; piStark(i)=0;  %  from Aragon 2014, piEi=0.39, pigi=4
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=238.2093;piEm(i)=5.20; piGA(i)=3.13E8;pi_z(i)=1;pigk(i)=12; piStark(i)=0;  %  from Aragon 2014, piEi=0.00, pigi=10
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=239.1478;piEm(i)=5.48; piGA(i)=0.0377E8;pi_z(i)=1;pigk(i)=10; piStark(i)=0;  %  from Aragon 2014, piEi=0.30, pigi=8
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=239.5626;piEm(i)=5.22; piGA(i)=1.926E8;pi_z(i)=1;pigk(i)=14; piStark(i)=0;  %  from Aragon 2014, piEi=0.05, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=240.4887;piEm(i)=5.24; piGA(i)=1.96E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=0.08, pigi=6
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=241.0520;piEm(i)=5.25; piGA(i)=1.55E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=0.11, pigi=4
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=241.1069;piEm(i)=5.26; piGA(i)=2.37E8;pi_z(i)=1;pigk(i)=2; piStark(i)=0;  %  from Aragon 2014, piEi=0.12, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=241.33;piEm(i)=5.257; piGA(i)=1.02E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0; % Shah 2012, piEi=0.121,  ARTICULO: Quantitative elemental analysis of steel using calibration-free laser-induced breakdown spectroscopy
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=243.0079;piEm(i)=7.93; piGA(i)=1.91E8;pi_z(i)=1;pigk(i)=10; piStark(i)=0;  %  from Aragon 2014, piEi=2.83, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=243.2262;piEm(i)=7.94; piGA(i)=1.57E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=2.84, pigi=6
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=243.9302;piEm(i)=8.23; piGA(i)=2.25E8;pi_z(i)=1;pigk(i)=14; piStark(i)=0;  %  from Aragon 2014, piEi=3.15, pigi=12
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=244.4516;piEm(i)=7.65; piGA(i)=2.78E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=2.58, pigi=6
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=258.5876;piEm(i)=4.79; piGA(i)=0.894E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=0.00, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=259.1543;piEm(i)=5.82; piGA(i)=0.5725E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=1.04,
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=259.2785;piEm(i)=8.86; piGA(i)=2.74E8;pi_z(i)=1;pigk(i)=16; piStark(i)=0;  %  from Aragon 2014, piEi=4.08, pigi=14
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=259.8370;piEm(i)=4.82; piGA(i)=1.43E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=0.05, pigi=8
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=259.9396;piEm(i)=4.77; piGA(i)=2.35E8;pi_z(i)=1;pigk(i)=10; piStark(i)=0;  %  from Aragon 2014, piEi=0.00, 
i=i+1; piUse(i)=1;piquees(i)=FeII; pila(i)=261.1874;piEm(i)=4.79; piGA(i)=1.20E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=0.05, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=261.7618;piEm(i)=4.82; piGA(i)=0.488E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=0.08, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=272.754;piEm(i)=5.58; piGA(i)=0.857E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0;  %  S.M. Pershin 2012, piEi=1.04, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=273.0734;piEm(i)=5.62; piGA(i)=0.279E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0;  %  from Aragon 2014, piEi=1.08, 
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=273.6966;piEm(i)=5.60; piGA(i)=1.22E8;pi_z(i)=1;pigk(i)=2; piStark(i)=0;  %  from Aragon 2014, piEi=1.08, pigi=4
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=273.9548;piEm(i)=5.51; piGA(i)=2.21E8;pi_z(i)=1;pigk(i)=8; piStark(i)=0;  %  from Aragon 2014, piEi=0.99,
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=274.3197;piEm(i)=5.62; piGA(i)=1.97E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0;  %  from Aragon 2014, piEi=1.10,pigi=2
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=274.6484;piEm(i)=5.59; piGA(i)=2.05E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=1.08,pigi=4
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=274.6982;piEm(i)=5.55; piGA(i)=1.69E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0;  %  from Aragon 2014, piEi=1.04,
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=275.3288;piEm(i)=7.77; piGA(i)=1.89E8;pi_z(i)=1;pigk(i)=12; piStark(i)=0;  %  from Aragon 2014, piEi=3.27,pigi=10
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=275.5737;piEm(i)=5.48; piGA(i)=2.15E8;pi_z(i)=1;pigk(i)=10; piStark(i)=0;  %  from Aragon 2014, piEi=0.99,pigi=8
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=276.181;piEm(i)=5.58; piGA(i)=0.11E8;pi_z(i)=1;pigk(i)=4; piStark(i)=0;  %  S.M. Pershin 2012, piEi=1.09, ARTICULO: Physics of selective evaporation of components during laser ablation of stainless steels
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=298.482;piEm(i)=5.823; piGA(i)=2.6E8;pi_z(i)=1;pigk(i)=6; piStark(i)=0; % , selected FeII line in [Praher2010] below self-absorption threshold
i=i+1; piUse(i)=0;piquees(i)=FeII; pila(i)=298.554;piEm(i)=5.8756; piGA(i)=9.56E7;pi_z(i)=1;pigk(i)=4; piStark(i)=0; % , selected FeII line in [Praher2010] below self-absorption threshold
IRSAC_reference(FeII) = i; % 298.554 lowest selfabsorption in Praher2010

%NiI Nickel
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=310.1557;piEi(i)=0.10908;piEm(i)=4.105; piGA(i)=4.4E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % from ? there is another interfering line at 310.1878nm
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=341.476;piEi(i)=0.0254;piEm(i)=3.6552; piGA(i)=5.0E8;pi_z(i)=0;pigk(i)=9; piStark(i)=1E16/0.000381; % from Cremer's book p. 264 & LIBS++ Stark & Reinhard !!! Ej too small
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=344.626;piEi(i)=0.109;piEm(i)=3.706; piGA(i)=2.2E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=345.289;piEi(i)=0.109;piEm(i)=3.699; piGA(i)=6.9E7;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=345.847;piEi(i)=0.212;piEm(i)=3.796; piGA(i)=3.0E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=347.254;piEi(i)=0.109;piEm(i)=3.679; piGA(i)=8.4E7;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=349.296;piEi(i)=0.109;piEm(i)=3.657; piGA(i)=2.9E8;pi_z(i)=0;pigk(i)=3; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=351.034;piEi(i)=0.212;piEm(i)=3.743; piGA(i)=1.2E8;pi_z(i)=0;pigk(i)=1; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=351.505;piEi(i)=0.109;piEm(i)=3.635; piGA(i)=2.9E8;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=352.454;piEi(i)=0.02539;piEm(i)=3.542; piGA(i)=5.0E8;pi_z(i)=0;pigk(i)=5; piStark(i)=0; % from Cremer's book p. 264
i=i+1; piUse(i)=1;piquees(i)=NiI; pila(i)=361.939;piEi(i)=0.423;piEm(i)=3.847; piGA(i)=4.6E8;pi_z(i)=0;pigk(i)=7; piStark(i)=1E16/0.001; % from Cremer's & LIBS++ Stark
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=385.8297;piEi(i)=0.423;piEm(i)=3.635; piGA(i)=4.8E7;pi_z(i)=0;pigk(i)=7; piStark(i)=0; % from NIST higher intensities
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=388.97;piEi(i)=0.2748;piEm(i)=3.46455; piGA(i)=1.2E4;pi_z(i)=0;pigk(i)=3; piStark(i)=0; % reference line in [sun2009] IRSAC
IRSAC_reference(NiI) = i; % the above line has been selected as reference in [sun2009] IRSAC
i=i+1; piUse(i)=0;piquees(i)=NiI; pila(i)=471.441;piEi(i)=3.380;piEm(i)=6.09; piGA(i)=5.1E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0; % from a table of a book?

%NiII 
i=i+1; piUse(i)=0;piquees(i)=NiII; pila(i)=220.672;piEi(i)=1.254;piEm(i)=6.871; piGA(i)=1.33E9;pi_z(i)=1;pigk(i)=8; piStark(i)=0; % NIST higher intensities
IRSAC_reference(NiII) = i; % the reference line in [sun2009] IRSAC 221.32 does not exist. the one above has the lowest gA and higer Em
i=i+1; piUse(i)=0;piquees(i)=NiII; pila(i)=221.648;piEi(i)=1.0407;piEm(i)=6.632; piGA(i)=4.1E9;pi_z(i)=1;pigk(i)=12; piStark(i)=0; % NIST higher intensities
i=i+1; piUse(i)=0;piquees(i)=NiII; pila(i)=231.604;piEi(i)=1.0407;piEm(i)=6.538734; piGA(i)=2.3E9;pi_z(i)=1;pigk(i)=8; piStark(i)=0; % Reinhard book  p. 529

%MnI Manganese 
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=279.482;piEi(i)=0.000;piEm(i)=4.435; piGA(i)=3.0E9;pi_z(i)=0;pigk(i)=8; piStark(i)=0; % Cremers !!! Ei=0 
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=280.106;piEi(i)=0.000;piEm(i)=4.429; piGA(i)=1.5E9;pi_z(i)=0;pigk(i)=4; piStark(i)=0; % NIST higher intensities !!! Ei=0 
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=380.6711;piEi(i)=2.114;piEm(i)=5.370; piGA(i)=7.1E8;pi_z(i)=0;pigk(i)=12; piStark(i)=0; % NIST higher intensities 
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=404.136;piEi(i)=2.114;piEm(i)=5.181; piGA(i)=7.87E8;pi_z(i)=0;pigk(i)=10; piStark(i)=0; % Cremers ,  there is another interfering line at 310.1878nm
IRSAC_reference(MnI) = i; % the above line has been selected as reference in [sun2009] IRSAC
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=403.076;piEi(i)=0.000;piEm(i)=3.075; piGA(i)=1.4E8;pi_z(i)=0;pigk(i)=8; piStark(i)=0'; % cremers book apendix c.2 !!! Ei=0
i=i+1; piUse(i)=1;piquees(i)=MnI; pila(i)=482.352;piEi(i)=2.319;piEm(i)=4.889; piGA(i)=3.99E8;pi_z(i)=0;pigk(i)=8; piStark(i)=0'; % cremers book apendix c.2


%MnII
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=257.61;piEi(i)=0;piEm(i)=4.811; piGA(i)=2.52E9;pi_z(i)=1;pigk(i)=9; piStark(i)=0; % Reinhard book !!! Ei=0
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=261.020;piEi(i)=3.415;piEm(i)=8.164; piGA(i)=4.5E9;pi_z(i)=1;pigk(i)=15; piStark(i)=0; % NIST higher intensities
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=261.815;piEi(i)=3.418;piEm(i)=8.152; piGA(i)=3.8E9;pi_z(i)=1;pigk(i)=13; piStark(i)=0; % NIST higher intensities
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=263.817;piEi(i)=3.421;piEm(i)=8.119; piGA(i)=1.9E9;pi_z(i)=1;pigk(i)=7; piStark(i)=0; % Reinhard book p. 516
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=270.84;piEi(i)=3.4199;piEm(i)=7.9963; piGA(i)=1.5E9;pi_z(i)=1;pigk(i)=9; piStark(i)=0; % reference line in [sun2009] IRSAC
IRSAC_reference(MnII) = i; % the above line has been selected as reference in [sun2009] IRSAC
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=293.305;piEi(i)=1.1745;piEm(i)=5.400; piGA(i)=6.12E8;pi_z(i)=1;pigk(i)=3; piStark(i)=0; % Reinhard book p. 270
i=i+1; piUse(i)=1;piquees(i)=MnII; pila(i)=294.92;piEi(i)=1.1745;piEm(i)=5.377; piGA(i)=1.37E9;pi_z(i)=1;pigk(i)=7; piStark(i)=0; % Reinhard book p. 517 

%TiI Titanium
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=334.853;piEi(i)=0.000;piEm(i)=3.701; piGA(i)=2.7E6;pi_z(i)=0;pigk(i)=3; piStark(i)=0; % Cremers !!! Ei=0 
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=365.349;piEi(i)=0.048;piEm(i)=3.440; piGA(i)=9.5E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0; % Cremers !!! Ei=0 
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=498.17;piEi(i)=0.848;piEm(i)=3.336; piGA(i)=8.58E8;pi_z(i)=0;pigk(i)=13; piStark(i)=0; % Cremers C.2 
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=499.11;piEi(i)=0.836;piEm(i)=3.319; piGA(i)=6.42E8;pi_z(i)=0;pigk(i)=11; piStark(i)=0; % paper ¿?
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=503.646;piEi(i)=1.443;piEm(i)=3.905; piGA(i)=3.55E8;pi_z(i)=0;pigk(i)=9; piStark(i)=0; % NIST higher intensities, not Ei=0
i=i+1; piUse(i)=1;piquees(i)=TiI; pila(i)=521.038;piEi(i)=0.048;piEm(i)=2.427; piGA(i)=3.5E7;pi_z(i)=0;pigk(i)=9; piStark(i)=0; % NIST higher intensities
% Cremers table C.2: 351.9, 398.92nm cant't find in NIST

%TiII
i=i+1; piUse(i)=1;piquees(i)=TiII; pila(i)=308.803;piEi(i)=0.04878;piEm(i)=4.062; piGA(i)=1.2E9;pi_z(i)=1;pigk(i)=8; piStark(i)=0; % paper ¿? !!! Ei=0
i=i+1; piUse(i)=1;piquees(i)=TiII; pila(i)=323.452;piEi(i)=0.0487;piEm(i)=3.881; piGA(i)=1.71E9;pi_z(i)=1;pigk(i)=10; piStark(i)=0; % Cremers !!! Ei=0
i=i+1; piUse(i)=1;piquees(i)=TiII; pila(i)=336.1;piEi(i)=0.0280;piEm(i)=3.716; piGA(i)=1.58E9;pi_z(i)=1;pigk(i)=10; piStark(i)=0; % paper ¿? !!! Ei=0
i=i+1; piUse(i)=1;piquees(i)=TiII; pila(i)=439.503;piEi(i)=1.084;piEm(i)=3.904; piGA(i)=7.5E7;pi_z(i)=1;pigk(i)=8; piStark(i)=0; % paper ¿? 

%Hg lines FOR CALIBRATION PORPOUSES mainly
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=184.94994;piEi(i)=0.0000000;piEm(i)=6.7036623;piGA(i)=2.24e+09;piAcc{i}='A';pi_z(i)=0;pigk(i)=3;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=234.543;piEi(i)=4.6673829;piEm(i)=9.951960;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=237.8324;piEi(i)=4.6673829;piEm(i)=9.878884;piGA(i)=1.1e+07;piAcc{i}='D';pi_z(i)=0;pigk(i)=3;piIntens(i)=4000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=248.2001;piEi(i)=4.8864946;piEm(i)=9.880321;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=300;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=248.3815;piEi(i)=4.8864946;piEm(i)=9.876651;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=170;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=253.4772;piEi(i)=4.6673829;piEm(i)=9.5572499;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=253.65210;piEi(i)=0.0000000;piEm(i)=4.8864946;piGA(i)=2.52e+07;piAcc{i}='A+';pi_z(i)=0;pigk(i)=3;piIntens(i)=900000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=257.6285;piEi(i)=4.8864946;piEm(i)=9.697563;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=265.2039;piEi(i)=4.8864946;piEm(i)=9.5601496;piGA(i)=2.0e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=1600;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=265.3690;piEi(i)=4.8864946;piEm(i)=9.5572499;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=6000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=265.5134;piEi(i)=4.8864946;piEm(i)=9.5547145;piGA(i)=5.50e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=269.8828;piEi(i)=5.4606247;piEm(i)=10.053254;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=7;piIntens(i)=200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=269.9376;piEi(i)=5.4606247;piEm(i)=10.052323;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=275.2777;piEi(i)=4.6673829;piEm(i)=9.1700120;piGA(i)=1.8e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=280.3466;piEi(i)=5.4606247;piEm(i)=9.881854;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=7;piIntens(i)=180;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=289.36010;piEi(i)=4.8864946;piEm(i)=9.1700120;piGA(i)=4.71e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=296.72830;piEi(i)=4.6673829;piEm(i)=8.8445373;piGA(i)=1.4e+08;piAcc{i}='D';pi_z(i)=0;pigk(i)=3;piIntens(i)=3000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=302.15040;piEi(i)=5.4606247;piEm(i)=9.5628233;piGA(i)=3.6e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=7;piIntens(i)=1200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=302.3471;piEi(i)=5.4606247;piEm(i)=9.5601496;piGA(i)=4.7e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=300;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=312.56740;piEi(i)=4.8864946;piEm(i)=8.8519848;piGA(i)=3.3e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=4000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=0;piquees(i)=HgI;pila(i)=313.15550;piEi(i)=4.8864946;piEm(i)=8.8445373;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=3000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=0;piquees(i)=HgI;pila(i)=313.18440;piEi(i)=4.8864946;piEm(i)=8.8441713;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=4000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=334.14840;piEi(i)=5.4606247;piEm(i)=9.1700120;piGA(i)=5.04e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=700;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=365.01580;piEi(i)=5.4606247;piEm(i)=8.8563375;piGA(i)=9.03e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=7;piIntens(i)=9000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=365.48420;piEi(i)=5.4606247;piEm(i)=8.8519848;piGA(i)=9.20e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=3000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=366.28870;piEi(i)=5.4606247;piEm(i)=8.8445373;piGA(i)=1.1e+07;piAcc{i}='C';pi_z(i)=0;pigk(i)=3;piIntens(i)=500;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=366.32840;piEi(i)=5.4606247;piEm(i)=8.8441713;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=404.65650;piEi(i)=4.6673829;piEm(i)=7.7304550;piGA(i)=6.21e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=12000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=407.78370;piEi(i)=4.8864946;piEm(i)=7.9260766;piGA(i)=4.0e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=434.74945;piEi(i)=6.7036623;piEm(i)=9.5547145;piGA(i)=4.2e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=150;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=435.83350;piEi(i)=4.8864946;piEm(i)=7.7304550;piGA(i)=1.7e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=12000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=496.010;piEi(i)=9.539998;piEm(i)=12.03889;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=7;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=535.4034;piEi(i)=7.7304550;piEm(i)=10.045526;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=5;piIntens(i)=130;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=546.07500;piEi(i)=5.4606247;piEm(i)=7.7304550;piGA(i)=1.5e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=6000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=567.581;piEi(i)=7.7304550;piEm(i)=9.914254;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=600;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=576.96100;piEi(i)=6.7036623;piEm(i)=8.8519848;piGA(i)=1.18e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=579.06700;piEi(i)=6.7036623;piEm(i)=8.8441713;piGA(i)=1.6e+08;piAcc{i}='D';pi_z(i)=0;pigk(i)=5;piIntens(i)=900;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=580.3782;piEi(i)=7.9260766;piEm(i)=10.0617497;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=585.9254;piEi(i)=7.7304550;piEm(i)=9.8459093;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=130;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=671.634;piEi(i)=7.9260766;piEm(i)=9.771569;piGA(i)=7.2e+05;piAcc{i}='D';pi_z(i)=0;pigk(i)=3;piIntens(i)=600;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=690.746;piEi(i)=7.7304550;piEm(i)=9.524891;piGA(i)=8.0e+06;piAcc{i}='D';pi_z(i)=0;pigk(i)=5;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=708.1901;piEi(i)=7.7304550;piEm(i)=9.4806917;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=3;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=HgI;pila(i)=709.1860;piEi(i)=7.7304550;piEm(i)=9.4782338;piGA(i)=0;piAcc{i}='';pi_z(i)=0;pigk(i)=1;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ADD Intensity>100

%Argon lines for calibration
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=415.8590;piEi(i)=11.54835433;piEm(i)=14.52891337;piGA(i)=7.00e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=419.0713;piEi(i)=11.54835433;piEm(i)=14.50606752;piGA(i)=1.40e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=419.8317;piEi(i)=11.62359262;piEm(i)=14.57594866;piGA(i)=2.57e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=420.0674;piEi(i)=11.54835433;piEm(i)=14.49905352;piGA(i)=6.8e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=7;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=425.9362;piEi(i)=11.82807106;piEm(i)=14.73811524;piGA(i)=3.98e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=1;piIntens(i)=200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=426.6286;piEi(i)=11.62359262;piEm(i)=14.52891337;piGA(i)=1.6e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=427.2169;piEi(i)=11.62359262;piEm(i)=14.52491318;piGA(i)=2.4e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=150;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=430.0101;piEi(i)=11.62359262;piEm(i)=14.50606752;piGA(i)=1.88e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=433.3561;piEi(i)=11.82807106;piEm(i)=14.68829018;piGA(i)=2.8e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=451.0733;piEi(i)=11.82807106;piEm(i)=14.57594866;piGA(i)=1.18e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=1;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=667.7282;piEi(i)=11.62359262;piEm(i)=13.47988670;piGA(i)=2.36e+05;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=675.2834;piEi(i)=12.90701519;piEm(i)=14.74254073;piGA(i)=9.65e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=5;piIntens(i)=150;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=687.1289;piEi(i)=12.90701519;piEm(i)=14.71089798;piGA(i)=8.34e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=3;piIntens(i)=150;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=696.5431;piEi(i)=11.54835433;piEm(i)=13.32785693;piGA(i)=1.9e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=703.0251;piEi(i)=13.07571560;piEm(i)=14.83881088;piGA(i)=1.34e+07;piAcc{i}='C';pi_z(i)=0;pigk(i)=5;piIntens(i)=150;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=706.7218;piEi(i)=11.54835433;piEm(i)=13.30222736;piGA(i)=1.9e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=706.8736;piEi(i)=13.09487245;piEm(i)=14.84836887;piGA(i)=6.0e+06;piAcc{i}='D+';pi_z(i)=0;pigk(i)=3;piIntens(i)=100;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=714.7042;piEi(i)=11.54835433;piEm(i)=13.28263891;piGA(i)=1.9e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=727.2936;piEi(i)=11.62359262;piEm(i)=13.32785693;piGA(i)=5.49e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=737.2118;piEi(i)=13.07571560;piEm(i)=14.7570515;piGA(i)=1.7e+07;piAcc{i}='D+';pi_z(i)=0;pigk(i)=9;piIntens(i)=200;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=738.3980;piEi(i)=11.62359262;piEm(i)=13.30222736;piGA(i)=4.2e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=750.3869;piEi(i)=11.82807106;piEm(i)=13.47988670;piGA(i)=4.5e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=751.4652;piEi(i)=11.62359262;piEm(i)=13.27303799;piGA(i)=4.0e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=15000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=763.5106;piEi(i)=11.54835433;piEm(i)=13.17177759;piGA(i)=1.22e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=25000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=772.3761;piEi(i)=11.54835433;piEm(i)=13.15314376;piGA(i)=1.6e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=15000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=772.4207;piEi(i)=11.72316029;piEm(i)=13.32785693;piGA(i)=3.51e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=794.8176;piEi(i)=11.72316029;piEm(i)=13.28263891;piGA(i)=5.58e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=800.6157;piEi(i)=11.62359262;piEm(i)=13.17177759;piGA(i)=2.4e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=801.4786;piEi(i)=11.54835433;piEm(i)=13.09487245;piGA(i)=4.6e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=25000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=810.3693;piEi(i)=11.62359262;piEm(i)=13.15314376;piGA(i)=7.50e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=811.5311;piEi(i)=11.54835433;piEm(i)=13.07571560;piGA(i)=2.3e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=7;piIntens(i)=35000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=826.4522;piEi(i)=11.82807106;piEm(i)=13.32785693;piGA(i)=4.59e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=840.8210;piEi(i)=11.82807106;piEm(i)=13.30222736;piGA(i)=1.12e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=15000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=842.4648;piEi(i)=11.62359262;piEm(i)=13.09487245;piGA(i)=1.08e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=852.1442;piEi(i)=11.82807106;piEm(i)=13.28263891;piGA(i)=4.17e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=15000;piStark(i)=0;% auto imported from NIST ADD Intensity>100
i=i+1;piUse(i)=1;piquees(i)=ArI;pila(i)=866.7944;piEi(i)=11.72316029;piEm(i)=13.15314376;piGA(i)=7.29e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=4500;piStark(i)=0;% auto imported from NIST ADD Intensity>100


% Neon lines for calibration
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=341.79031;piEi(i)=16.67082679;piEm(i)=20.29727974;piGA(i)=4.4e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=5;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=344.77022;piEi(i)=16.61906995;piEm(i)=20.21417951;piGA(i)=1.06e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=346.65781;piEi(i)=16.71538094;piEm(i)=20.29091559;piGA(i)=3.75e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=3;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=347.25706;piEi(i)=16.61906995;piEm(i)=20.18843456;piGA(i)=1.18e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=350.12154;piEi(i)=16.67082679;piEm(i)=20.21098944;piGA(i)=3.60e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=3;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=351.51900;piEi(i)=16.67082679;piEm(i)=20.19691626;piGA(i)=3.4e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=5;piIntens(i)=2000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=352.04714;piEi(i)=16.84805355;piEm(i)=20.36885387;piGA(i)=7.58e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=1;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=359.35263;piEi(i)=16.84805355;piEm(i)=20.29727974;piGA(i)=4.44e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=359.36385;piEi(i)=16.84805355;piEm(i)=20.29717103;piGA(i)=2.2e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=3;piIntens(i)=3000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=503.13484;piEi(i)=18.55510789;piEm(i)=21.01865501;piGA(i)=9.38e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=2500;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=503.77512;piEi(i)=18.55510789;piEm(i)=21.01552421;piGA(i)=4.00e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=9;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=508.03830;piEi(i)=18.57583585;piEm(i)=21.01560530;piGA(i)=2.37e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=1500;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=511.65032;piEi(i)=18.38162308;piEm(i)=20.80416986;piGA(i)=1.30e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=1500;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=512.22565;piEi(i)=18.69335927;piEm(i)=21.11318520;piGA(i)=9.80e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=5;piIntens(i)=1500;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=540.05616;piEi(i)=16.67082679;piEm(i)=18.96595353;piGA(i)=9.00e+05;piAcc{i}='B+';pi_z(i)=0;pigk(i)=1;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=576.44188;piEi(i)=18.55510789;piEm(i)=20.70536489;piGA(i)=1.01e+08;piAcc{i}='B';pi_z(i)=0;pigk(i)=9;piIntens(i)=7000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=582.01558;piEi(i)=18.57583585;piEm(i)=20.70550140;piGA(i)=6.02e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=585.24878;piEi(i)=16.84805355;piEm(i)=18.96595353;piGA(i)=6.15e+07;piAcc{i}='AA';pi_z(i)=0;pigk(i)=1;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=588.18950;piEi(i)=16.61906995;piEm(i)=18.72638130;piGA(i)=3.45e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=594.48340;piEi(i)=16.61906995;piEm(i)=18.70407102;piGA(i)=5.65e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=597.55343;piEi(i)=16.61906995;piEm(i)=18.69335927;piGA(i)=1.05e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=6000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=602.99968;piEi(i)=16.67082679;piEm(i)=18.72638130;piGA(i)=1.68e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=607.43376;piEi(i)=16.67082679;piEm(i)=18.71137652;piGA(i)=6.03e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=1;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=609.61630;piEi(i)=16.67082679;piEm(i)=18.70407102;piGA(i)=9.05e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=3000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=614.30627;piEi(i)=16.61906995;piEm(i)=18.63679141;piGA(i)=1.41e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=616.35937;piEi(i)=16.71538094;piEm(i)=18.72638130;piGA(i)=4.38e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=621.72812;piEi(i)=16.61906995;piEm(i)=18.61270512;piGA(i)=1.91e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=626.64952;piEi(i)=16.71538094;piEm(i)=18.69335927;piGA(i)=7.47e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=633.44276;piEi(i)=16.61906995;piEm(i)=18.57583585;piGA(i)=8.05e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=638.29914;piEi(i)=16.67082679;piEm(i)=18.61270512;piGA(i)=9.63e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=640.22480;piEi(i)=16.61906995;piEm(i)=18.55510789;piGA(i)=3.604e+08;piAcc{i}='AAA';pi_z(i)=0;pigk(i)=7;piIntens(i)=20000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=650.65277;piEi(i)=16.67082679;piEm(i)=18.57583585;piGA(i)=1.50e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=15000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=659.89528;piEi(i)=16.84805355;piEm(i)=18.72638130;piGA(i)=6.96e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=667.82766;piEi(i)=16.84805355;piEm(i)=18.70407102;piGA(i)=1.16e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=5000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=692.94672;piEi(i)=16.84805355;piEm(i)=18.63679141;piGA(i)=8.70e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=100000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=702.40500;piEi(i)=16.84805355;piEm(i)=18.61270512;piGA(i)=5.67e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=34000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=703.24128;piEi(i)=16.61906995;piEm(i)=18.38162308;piGA(i)=7.98e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=85000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=705.12922;piEi(i)=18.38162308;piEm(i)=20.13945716;piGA(i)=7.14e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=3;piIntens(i)=2200;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=705.91079;piEi(i)=18.38162308;piEm(i)=20.13751108;piGA(i)=3.73e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=10000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=717.39380;piEi(i)=16.84805355;piEm(i)=18.57583585;piGA(i)=1.44e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=77000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=724.51665;piEi(i)=16.67082679;piEm(i)=18.38162308;piGA(i)=3.03e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=77000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=743.88981;piEi(i)=16.71538094;piEm(i)=18.38162308;piGA(i)=7.41e+06;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=60000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=747.24383;piEi(i)=18.38162308;piEm(i)=20.04038629;piGA(i)=8.88e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=3;piIntens(i)=3100;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=748.88712;piEi(i)=18.38162308;piEm(i)=20.03674654;piGA(i)=1.16e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=32000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=753.57739;piEi(i)=18.38162308;piEm(i)=20.02644506;piGA(i)=9.18e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=3;piIntens(i)=28000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=754.40439;piEi(i)=18.38162308;piEm(i)=20.02464191;piGA(i)=3.87e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=13000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=811.85495;piEi(i)=18.61270512;piEm(i)=20.13945716;piGA(i)=1.16e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=3;piIntens(i)=3800;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=826.60769;piEi(i)=18.63679141;piEm(i)=20.13629502;piGA(i)=2.33e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=7200;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=837.76070;piEi(i)=18.55510789;piEm(i)=20.03464876;piGA(i)=4.49e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=9;piIntens(i)=76000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=841.71597;piEi(i)=18.57583585;piEm(i)=20.04842432;piGA(i)=7.49e+06;piAcc{i}='C+';pi_z(i)=0;pigk(i)=7;piIntens(i)=2700;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=849.53591;piEi(i)=18.57583585;piEm(i)=20.03486930;piGA(i)=2.75e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=7;piIntens(i)=69000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=859.12583;piEi(i)=18.69335927;piEm(i)=20.13610655;piGA(i)=1.56e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=41000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=863.46472;piEi(i)=18.61270512;piEm(i)=20.04820273;piGA(i)=1.14e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=35000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=878.37539;piEi(i)=18.72638130;piEm(i)=20.13751108;piGA(i)=1.62e+08;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=43000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500
i=i+1;piUse(i)=1;piquees(i)=NeI;pila(i)=885.38669;piEi(i)=18.63679141;piEm(i)=20.03674654;piGA(i)=8.20e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=5;piIntens(i)=27000;piStark(i)=0;% auto imported from NIST ADD Intensity>1500

%CI
% carbon lines
i=i+1;piUse(i)=0*1;piquees(i)=CI;pila(i)=193.090540;piEi(i)=1.2637284;piEm(i)=7.68476771;piGA(i)=1.02e+09;piAcc{i}='A';pi_z(i)=0;pigk(i)=3;piIntens(i)=20000000000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=199.3627;piEi(i)=1.2637284;piEm(i)=7.48277591;piGA(i)=2.5e+05;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=13000000000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=1;piquees(i)=CI;pila(i)=247.85612;piEi(i)=2.6840136;piEm(i)=7.68476771;piGA(i)=8.4e+07;piAcc{i}='C+';pi_z(i)=0;pigk(i)=3;piIntens(i)=640000000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=258.29010;piEi(i)=2.6840136;piEm(i)=7.48277591;piGA(i)=1.7e+04;piAcc{i}='C+';pi_z(i)=0;pigk(i)=3;piIntens(i)=530000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=296.48460;piEi(i)=0.0020354130;piEm(i)=4.1826219;piGA(i)=4.e+01;piAcc{i}='D';pi_z(i)=0;pigk(i)=5;piIntens(i)=6900000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=296.72240;piEi(i)=0.0053825826;piEm(i)=4.1826219;piGA(i)=1.0e+02;piAcc{i}='D';pi_z(i)=0;pigk(i)=5;piIntens(i)=13000000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=477.173346;piEi(i)=7.48779862;piEm(i)=10.08537735;piGA(i)=4.0e+06;piAcc{i}='C';pi_z(i)=0;pigk(i)=5;piIntens(i)=69000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=493.202553;piEi(i)=7.68476771;piEm(i)=10.19792594;piGA(i)=6.0e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=73000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=505.214919;piEi(i)=7.68476771;piEm(i)=10.13817181;piGA(i)=1.3e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=160000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=538.03308;piEi(i)=7.68476771;piEm(i)=9.98852470;piGA(i)=5.58e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=3;piIntens(i)=120000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=601.322;piEi(i)=8.64716033;piEm(i)=10.708468;piGA(i)=9.0e+06;piAcc{i}='D';pi_z(i)=0;pigk(i)=5;piIntens(i)=79000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=601.322;piEi(i)=8.64716033;piEm(i)=10.7084389;piGA(i)=4.0e+06;piAcc{i}='D';pi_z(i)=0;pigk(i)=9;piIntens(i)=79000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=711.3180;piEi(i)=8.64716033;piEm(i)=10.38970413;piGA(i)=2.22e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=9;piIntens(i)=110000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=711.519;piEi(i)=8.64039597;piEm(i)=10.38244864;piGA(i)=4.4e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=1;piIntens(i)=110000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=711.519;piEi(i)=8.64302361;piEm(i)=10.38507435;piGA(i)=1.53e+07;piAcc{i}='B';pi_z(i)=0;pigk(i)=7;piIntens(i)=110000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=786.088;piEi(i)=8.85066275;piEm(i)=10.42745830;piGA(i)=7.65e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=110000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0;piquees(i)=CI;pila(i)=805.859;piEi(i)=8.85066275;piEm(i)=10.38877079;piGA(i)=5.45e+06;piAcc{i}='B';pi_z(i)=0;pigk(i)=5;piIntens(i)=59000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000
i=i+1;piUse(i)=0*1;piquees(i)=CI;pila(i)=833.51443;piEi(i)=7.68476771;piEm(i)=9.17184604;piGA(i)=3.51e+07;piAcc{i}='B+';pi_z(i)=0;pigk(i)=1;piIntens(i)=5900000;piStark(i)=0;% auto imported from NIST ASD Intensity>50000

% Carbon II

i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=250.912;piEi(i)=13.715649;piEm(i)=18.655492;piGA(i)=1.88e+08;piAcc{i}='A';pi_z(i)=1;pigk(i)=4;piIntens(i)=250;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=251.206;piEi(i)=13.720781;piEm(i)=18.654858;piGA(i)=3.37e+08;piAcc{i}='A';pi_z(i)=1;pigk(i)=6;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=274.649;piEi(i)=16.331740;piEm(i)=20.844687;piGA(i)=1.74e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=250;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=283.671;piEi(i)=11.963699;piEm(i)=16.333123;piGA(i)=1.32e+08;piAcc{i}='B+';pi_z(i)=1;pigk(i)=4;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=283.760;piEi(i)=11.963699;piEm(i)=16.331740;piGA(i)=6.58e+07;piAcc{i}='B+';pi_z(i)=1;pigk(i)=2;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=391.898;piEi(i)=16.331740;piEm(i)=19.494540;piGA(i)=1.27e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=2;piIntens(i)=570;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=392.069;piEi(i)=16.333123;piEm(i)=19.494540;piGA(i)=2.54e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=2;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=426.700;piEi(i)=18.045809;piEm(i)=20.950644;piGA(i)=1.34e+09;piAcc{i}='C+';pi_z(i)=1;pigk(i)=6;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=426.726;piEi(i)=18.045987;piEm(i)=20.950644;piGA(i)=1.90e+09;piAcc{i}='C+';pi_z(i)=1;pigk(i)=8;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=426.726;piEi(i)=18.045987;piEm(i)=20.950644;piGA(i)=9.54e+07;piAcc{i}='C+';pi_z(i)=1;pigk(i)=6;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=513.294;piEi(i)=20.701286;piEm(i)=23.116071;piGA(i)=1.56e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=513.328;piEi(i)=20.704212;piEm(i)=23.118840;piGA(i)=1.68e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=6;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=514.349;piEi(i)=20.704212;piEm(i)=23.114045;piGA(i)=1.55e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=2;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=514.516;piEi(i)=20.709788;piEm(i)=23.118840;piGA(i)=3.89e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=6;piIntens(i)=570;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=515.109;piEi(i)=20.709788;piEm(i)=23.116071;piGA(i)=1.66e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=400;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=564.807;piEi(i)=20.704212;piEm(i)=22.898763;piGA(i)=7.88e+07;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=250;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=566.247;piEi(i)=20.709788;piEm(i)=22.898763;piGA(i)=1.17e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=588.977;piEi(i)=18.045987;piEm(i)=20.150478;piGA(i)=1.26e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=4;piIntens(i)=570;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=589.159;piEi(i)=18.045809;piEm(i)=20.149650;piGA(i)=6.98e+07;piAcc{i}='B';pi_z(i)=1;pigk(i)=2;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=657.805;piEi(i)=14.448827;piEm(i)=16.333123;piGA(i)=1.47e+08;piAcc{i}='A';pi_z(i)=1;pigk(i)=4;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=658.288;piEi(i)=14.448827;piEm(i)=16.331740;piGA(i)=7.32e+07;piAcc{i}='A';pi_z(i)=1;pigk(i)=2;piIntens(i)=570;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=678.390;piEi(i)=20.709788;piEm(i)=22.536906;piGA(i)=2.92e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=8;piIntens(i)=250;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=711.563;piEi(i)=22.532398;piEm(i)=24.274338;piGA(i)=2.88e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=8;piIntens(i)=250;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=711.990;piEi(i)=22.536906;piEm(i)=24.277799;piGA(i)=4.19e+08;piAcc{i}='B';pi_z(i)=1;pigk(i)=10;piIntens(i)=350;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=723.132;piEi(i)=16.331740;piEm(i)=18.045809;piGA(i)=1.40e+08;piAcc{i}='A';pi_z(i)=1;pigk(i)=4;piIntens(i)=800;piStark(i)=0;% auto imported from NIST ASD Intensity>1
i=i+1;piUse(i)=0;piquees(i)=CII;pila(i)=723.642;piEi(i)=16.333123;piEm(i)=18.045987;piGA(i)=2.51e+08;piAcc{i}='A';pi_z(i)=1;pigk(i)=6;piIntens(i)=1000;piStark(i)=0;% auto imported from NIST ASD Intensity>1


%%%%%%%%%%%
%disp('*********************************************using external piUses');
%piUses_optim_apice_31ago2019_1248_for_nparam1;

%partitions functions, interpolated using NIST data from
%https://physics.nist.gov/PhysRefData/ASD/levels_form.html Te=6000..20000K
poly_partfunc=zeros(Nspecies,4,1); % 3 degree polynomial interpolation of partition functions, from 6000K to 20000K
poly_partfunc(CaI,:) = [+1.51811E-11, -3.61057E-08, -1.37503E-03, +7.85876E+00];
poly_partfunc(CaII,:)= [+1.21705E-13, +5.11919E-09, +1.88296E-04, +1.04049E+00];
poly_partfunc(MgI,:) = [+5.84800E-12, -1.02634E-07, +6.38554E-04, -3.37430E-01];
poly_partfunc(MgII,:)= [+2.33535E-13, -4.25505E-09, +3.23152E-05, +1.90703E+00];
poly_partfunc(SrI,:) = [+6.11127E-12, +1.50741E-07, -0.0024826489,+9.8401229316];  
poly_partfunc(SrII,:)= [+3.54626E-13, +2.27577E-09, +1.68447E-04, + 1.13143];
poly_partfunc(ZnI,:) = [+2.83044E-12, -6.62522E-08, +5.34237E-04, -4.38897E-01];
poly_partfunc(ZnII,:)= [+1.64423E-13, -3.58172E-09, +2.67447E-05, +1.93217E+00];
poly_partfunc(CuI,:) = [+6.56177E-12, -1.25898E-07, +1.13081E-03, -1.09017E+00];
poly_partfunc(CuII,:)= [-2.95249E-13, +2.58956E-08, -2.00221E-04, +1.40651E+00];
poly_partfunc(AlI,:) = [+3.77444E-12, +1.10894E-08, -6.48764E-04, +8.65612E+00];
poly_partfunc(AlII,:)= [+1.00996E-13, +6.24378E-10, -1.95496E-05, +1.07319E+00];
poly_partfunc(SnI,:) = [+6.71496E-12, -1.40628E-07, +1.55177E-03, +1.85903E-01];
poly_partfunc(SnII,:)= [+7.55230E-13, -2.91627E-08, +5.03192E-04, +1.31003E+00];
poly_partfunc(VI,:) =  [-2.73579E-11, +2.10201E-06, -1.37134E-02, +6.88887E+01];
poly_partfunc(VII,:)=  [-1.46240E-12, +2.31560E-07, +4.38807E-03, +1.53628E+01];
poly_partfunc(FeI,:) = [+1.70038E-11, +4.28415E-07, -3.27306E-03, +3.24825E+01];
poly_partfunc(FeII,:)= [+2.43755E-12, +1.36340E-07, +2.17327E-03, +2.90773E+01];
poly_partfunc(NiI,:) = [+6.14267E-12, +1.01750E-07, -6.24182E-04, +3.15010E+01];
poly_partfunc(NiII,:)= [+1.45490E-12, -2.71607E-08, +1.90850E-03, +1.56002E+00];
poly_partfunc(MnI,:) = [+1.02761E-11, +1.49377E-07, -2.35897E-03, +1.36964E+01];
poly_partfunc(MnII,:)= [+2.04381E-12, +1.65890E-07, -1.52921E-03, +1.12160E+01];
poly_partfunc(TiI,:) =  [1];
poly_partfunc(TiII,:)=  [1];
poly_partfunc(CI,:) =  [+3.66793E-12, -1.03972E-07, +1.1629E-03, +5.31929 ];
poly_partfunc(CII,:)=  [+1.07323E-13, -4.65368E-10, -1.87951E-06, +5.94667];
%poly_partfunc(I,:) =  [];
%poly_partfunc(II,:)=  [];


%Einfinite for saha-bolztmann plot from http://www.webelements.com/calcium/atoms.html
% better yet: https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H-DS+i&units=1&at_num_out=on&el_name_out=on&shells_out=on&level_out=on&e_out=0&unc_out=on&biblio=on
Einfinite = zeros(Nspecies,1);
Einfinite(MgI)=7.65;% magnesium
Einfinite(MgII)=7.65; 
Einfinite(CaI)=6.11;% calcium
Einfinite(CaII)=6.11;
Einfinite(ZnI)= 9.394;%zinc
Einfinite(ZnII)= 9.394;
Einfinite(CuI)=7.726;%copper
Einfinite(CuII)=7.726;
Einfinite(AlI)=5.986; %aluminium
Einfinite(AlII)=5.986;
Einfinite(SnI)=7.35; % Tin
Einfinite(SnII)=7.35; %
Einfinite(VI) = 6.746187; % Vanadium
Einfinite(VII) = 6.746187;
Einfinite(FeI) =  7.9025; % Iron
Einfinite(FeII) = 7.9025; 
Einfinite(NiI) =  7.639878; % Nickel
Einfinite(NiII) = 7.639878; 
Einfinite(MnI) =  7.4340379; % Manganese
Einfinite(MnII) = 7.4340379; 
Einfinite(TiI) = 6.828120 ; % Titanium
Einfinite(TiII) = 6.828120 ; 
Einfinite(CI) = 11.2602880 ; % Carbon
Einfinite(CII) = 11.2602880 ; 

%Einfinite(I) =  ; % 
%Einfinite(II) = ; 
 

% discard lines from species not included in 'composition'
for s=1:Nspecies
    if ~any(composition(:)==s) && s~=Halfa % Halfa is always used for Ne estimation
        indexes = find(piquees==s);
        piUse(indexes)=0; % mark as not to be used
        %disp(['discarding ' strjoin(speciesStr(s))]);
    end
end
% discard resonant lines 
if (REMOVE_RESONANT_LINES)
    piUse(piResonant(:)==1)=0;
end
% discard  lines with lower energy level < 0.2
if (REMOVE_LOWER_LEVEL_ZERO)
    disp(['  Warning!! removal of lines with zero Ei, be sure ALL the lines has the value declared, otherwise it is zero by default!']);
    piUse(piEi(:)<0.2)=0;
end

if (USAR_HR2000)
  ValSaturation = 16000;   % expected maximum value without saturation
  ValMinHeightAllPeaks = 1000; % if no peak is above this value, the spectrum is considered not-valid. noise values are typ. <1000
  ValMinHeightRelative = 300; % altura RELATIVA mínima de un pico para que sea válido
end
if (USAR_PIMAX3)
  ValSaturation = 65000;   % expected maximum value without saturation
  ValMinHeightAllPeaks = 5000; % if no peak is above this value, the spectrum is considered not-valid. noise values are typ. <1000
  ValMinHeightSinglePeak = 1200; % new v20 19mar2020: if DISCARD_ALSO_WEAK_PEAKS ==1 and the raw height of any single peak is below this value, the entire spectrum is marked as not valid  
  ValMinHeightRelative = 200; % altura RELATIVA mínima de un pico para que sea válido
end
if (USAR_AVANTES)
  ValSaturation = 65000;   % expected maximum value without saturation
  ValMinHeightAllPeaks = 2000; % if no peak is above this value, the spectrum is considered not-valid. noise values are typ. <1000
  ValMinHeightSinglePeak = 1800; % new v20 19mar2020: if DISCARD_ALSO_WEAK_PEAKS ==1 and the raw height of any single peak is below this value, the entire spectrum is marked as not valid
  ValMinHeightRelative = 100; % altura RELATIVA mínima de un pico para que sea válido
end
ValSaturationCorrected = ValSaturation*ones(Npixeles,1); % si no hay corrección de irradiancia, se queda así 

if ( strcmp(param1Name,'deltaT') == 1)   modo2D = 1; else     modo2D = 0; end

if (USAR_HR2000)  % spectra from HR2000 is in another variable
 Npixeles=2048;
 spectra = hr2000_spectra;
 lambdas=hr2000_lambdas;
 ancho = lambdas(2048)-lambdas(1);
 central = ceil(lambdas(1) + ancho/2);
 nm2pixels=2048/ancho;   
 grating=['hr2000']; 
 lampara_FWHM = 1.0;  % según datasheet
end

if (USAR_PIMAX3)  % spectra from HR2000 is in another variable
   
    if exist('REMOTE_PIMAX','var')  % spectra from PIMAX3 is from remote machine
      Npixeles=1024;
      %spectra = pimax3_spectra;
      %lambdas=pimax3_lambdas;
      spectra = spectra_remote_pimax;
      lambdas=lambdas_remote_pimax;
      
    end
   
   ancho = lambdas(1024)-lambdas(1);
   central = ceil(lambdas(1) + ancho/2);
   nm2pixels=1024/ancho;   
   if (ancho > 250) % grating 150lpp
      grating=['pimax3 150lpp@' num2str(central) 'nm']; 
      lampara_FWHM = 0.9;  % pico de 435 nm con grating 150lpp @ 400nm, otras anchuras de otros picos son 0.86nm, 0.94nm
   else
      if (ancho > 30) % graging 1200lpp
        grating=['pimax3 1200lpp@' num2str(central) 'nm'];
        lampara_FWHM = 0.133;
      
      else
        grating=['pimax3 1800lpp@' num2str(central) 'nm'];
        lampara_FWHM = 0.1;
      end
   end

end

 % instrumental width FWHM (measured with the HgAr calibration lamp)
 % Avantes:  178.5-255.3	0.097;  251.3-316.9	0.065 ; 312.4-365.6	0.0661;
 % 360.7-403.5	0.053 ; 398.8-492.4	0.085 ; 488.9-564.9	0.0819 ;
 % 560.2-680.2	0.072; 681.7-889.3	0.142

if (USAR_AVANTES)  % spectra from Avantes is already the main spectra() variable 
 ancho = lambdas(Npixeles)-lambdas(1);
 central = ceil(lambdas(1) + ancho/2);
 nm2pixels=Npixeles/ancho;   % this value is not per channel and it is an aproximation
 grating=['avantes']; 
 lampara_FWHM = [ 251.3 0.097 312.4 0.065 360.7 0.0661 398.8 0.053 488.9 0.085 560.2 0.0819 681.7 0.072 889.3 0.142 ]; % hasta 251.3 0.097nm, hasta 312.4 0.065nm, etc. 
 % these are the wavelength windows for each channel: h5 = 178.63 to 251.32
 % h6 = 251.34 to 312.34 %h7 = 312.35 to 360.73 %h8 = 360.74 to 398.74
 %h1 = 398.75 to 488.81 % h2 = 488.85 to 560.17 %h3 = 560.17 to 680.17
 %h4 = 681.57 to 889.26nm.
end

 %% spectra stitching for multiwindows measurements  *new* 12mar2017
 if (MULTI_WINDOW==1)
   NpixelesEachWindow = Npixeles;
   Npixeles = Nwindows*NpixelesEachWindow;
   spectra = zeros(Npixeles,Npulses, Nparam1, Nparam2,'single'); % nueva variable con los espectros pegados
   lambdas = zeros(Npixeles,1);
   lastLambda = 0.0;
   lastIndex = 0;
   for nWindow =1:Nwindows % se ejecuta una sola vez en el modo no MULTI_WINDOW
     ix = find(MW_lambdas(:,nWindow)>lastLambda); % primera lambda mayor que la última
     lambdas(lastIndex+1:lastIndex+1+NpixelesEachWindow-ix(1)) = MW_lambdas(ix(1):NpixelesEachWindow,nWindow);
     spectra(lastIndex+1:lastIndex+1+NpixelesEachWindow-ix(1),:,:,:) = MW_spectra(ix(1):NpixelesEachWindow,:,:,:,nWindow);
     lastIndex = lastIndex+1+NpixelesEachWindow-ix(1);
     lastLambda = MW_lambdas(NpixelesEachWindow,nWindow);
   end
   lambdas = lambdas(1:lastIndex);
   spectra = spectra(1:lastIndex,:,:,:);
   Npixeles = lastIndex;
   grating=['multiwindow@' num2str(central) 'nm'];
   lampara_FWHM = 0.133;  % se supone que es el de 1200lpp, medido experimentalmente, modelando los picos, y valor medio
 end

%pixex=ceil(AvantesFWHMpixeles(lampara_FWHM)/2); % número de píxeles hacia
%los lados que ocupa el ancho FWHM, para integrar, en el Avantes es
%variable ; esto ya no se puede hacer porque es variable. La anchura
%instrumental sale siempre de 2 pixeles, es muy pequeña.

FEATURES='';
FEATURES = [ FEATURES 'CalibrationFree,'];
FEATURES = [ FEATURES 'REMOVE_BACKGROUND=' num2str(REMOVE_BACKGROUND) ','];
FEATURES = [ FEATURES 'REMOVE_BASELINE=' num2str(REMOVE_BASELINE) ','];
if (USAR_HR2000~=0)
      FEATURES = [ FEATURES 'HR2000,'];
end
FEATURES = [ FEATURES grating ','];

% fix the spectral response of the setup, should be done for CFLIBS, but in case there is no calibration data, it can be skipped
if (CORREGIR_IRRADIANCIA==1 ) %%% Cambié la condicion para poder usar carpetas antiguas
   %24august2019 polynomials have large error around 250nm (emission max),
   %there is a new DH2000.MAT file with dh2000_irradiance_both, onlyUV,
   %onlyIR, split and lambdas
   poly_deuterio = fliplr([ -318.1404847, 4.735948769, -0.02636499, 6.98226E-05, -8.81643E-08, 4.2219E-11 ]); % el primero es p0, al revés que polyfit y polyval 
   poly_halog = fliplr([ 1.275454122, -0.002887867, - 2.17019E-05, 6.94333E-08, -3.7258E-11 ]);
   poly_emision_combinada = [ -8.08165E-09, 1.89917E-05, - 1.26612E-02, 2.71727 ]; % p0 is the last value, 6jun2018 este polinomio no modela bien en el borde izquierdo (200nm)
   load ('dh2000.mat'); % 24august2019 contains the irradiance values from calibration sheet. values below  200nm are extrapolated

   if (strcmp(grating,'pimax3 150lpp@390nm')) % august2019, G=10, gate=20ms
       load('respuesta_espectrometro_pimax_390nm_azul_agosto2019.mat');
       emision = timeseries(dh2000_irradiance_onlyUV,dh2000_lambdas);
       emision2= resample(emision,lambdas);
       respuesta_espectrometro_raw = measured_irradiance_dh2000spectrum_noBackground ./ emision2.Data ;
       respuesta_espectrometro = smooth(respuesta_espectrometro_raw,100,'lowess');
       dividir_respuesta_espectrometro  = respuesta_espectrometro ./ max(respuesta_espectrometro); % this is the vector to divide to
   elseif (strcmp(grating,'multiwindow@290nm')) % medidas multiwindow 1200lpp 6 ventana
      load('respuesta_espectrometro_1200lpp_6ventanas_backgroundquitado_mar2017.mat');
      i = find(ganancias_respuesta_espectrometro==gainIntensifier);
      if (i==0)
       error('There is not an spectral correction file for the gain apparently used for the measurements');
      end
      dividir_respuesta_espectrometro = respuesta_espectrometro(:,i); % hay que dividir entre este vector
      disp('Warning!!! spectral response correction for Avantes is not correct below 200nm (not enough light)');
   elseif (strcmp(grating,'avantes'))
       % 24ago2019 new calibration with SPLIT spectrum: irradiance captured
       % with only deuterium, only halogen, and then stitched, we should
       % use the var dh2000_irradiance_split (not both on at once)
       load('respuesta_espectrometro_avantes_roja_1mmSR_con_lentes_split_22ene2021'); %updated 22ene2021
       emision = timeseries(dh2000_irradiance_split,dh2000_lambdas);
       emision2= resample(emision,lambdas);
       respuesta_espectrometro_raw = measured_irradiance_dh2000spectrum_noBackground ./ emision2.Data ;
       % smothing create artifacts in the stitches of the channels, I think
       % is better with no smoothing, it is already done in a per-channel basis
       %respuesta_espectrometro =
       %smooth(respuesta_espectrometro_raw,100,'lowess'); % already done
       respuesta_espectrometro = respuesta_espectrometro_raw;
       dividir_respuesta_espectrometro  = respuesta_espectrometro ./ max(respuesta_espectrometro); % this is the vector to divide to
   else    
     error 'There is no spectral correction file for this particular spectrometer/spectral window configuration'
   end
   % al corregir en irradiancia el valor de saturación ya no es aplicable,
   % voy a hacer una matriz de valores de saturación "corregidos"
   %25ago2019 there are values in "dividir_..." less than 0, subtitute to
   %0.005 (max is 1)
   dividir_respuesta_espectrometro(dividir_respuesta_espectrometro<0.005) = 0.005;
   which_lambdas_have_too_low_sensitivity = dividir_respuesta_espectrometro<0.005;
   ValSaturationCorrected = ValSaturation ./ dividir_respuesta_espectrometro;  
end

PNbarras=strrep(PN,'\','/'); % gira un poco las barras :)
para_el_titulo = [' ' grating ' ' num2str(floor(LAMBDA_MAGNESIO)) '_' num2str(floor(LAMBDA_CALCIO)) 'nm' ' I=' num2str(INTEGRAR) ' B=' num2str(REMOVE_BACKGROUND) ' ' que_ratio_es_texto];

%cálculo de los polinomios lambda<->pixel, permite usar decimales de pixel
%y recalibrar el espectro
pixeles = (1:size(lambdas))';
%poly_lambda2pixel = polyfit(lambdas,pixeles,2);  % esto con el de avantes
%ya no vale
%poly_pixel2lambda = polyfit(pixeles,lambdas,2);

% SPLIT_PARAM2 guarda la matriz spectra() en bloques para cada iteración param2
if (SPLIT_PARAM2==1)
   spectra = zeros(Npixeles,Npulses, Nparam1, Nparam2,'single');
   for param2IndexReal= 1:Nparam2
      param2 = param2From + (param2IndexReal-1)*param2Step; 
      fn_split = strcat('spectra_idx2=',sprintf('%04i',param2IndexReal),'_', param2Name,'=',sprintf('%0.2f',param2));
      fn_split = strrep(fn_split,'.',','); %replaces dots with commas, just in case 
      if (exist(strcat(PN,'\',fn_split,'.mat'), 'file')==2)
        tmp_spectra = load(strcat(PN,'\',fn_split),'spectra'); % ead the partial spectra for this param2 iteration
        spectra(:,:,:,param2IndexReal) = tmp_spectra.spectra;
%      else if(param2IndexReal==4)
%          spectra(:,:,:,param2IndexReal) = backup_spectra;   
%          end
      end
      
   end
end

% SPLIT la dimensión 1 en trozos más pequeños si REARRANGE!=0, vale para secuencias de pulsos muy largas en profundidad, promediar cada pocos espectros
if (REARRANGE > 0)
  ra_factor = Npulses / REARRANGE; % factor en el que se incrementa dimensión 1 y se reduce Npulses
  ra_spectra = zeros(Npixeles,REARRANGE, Nparam1 * ra_factor, Nparam2,'single');
  k=1;
  for j=1:Nparam1
   for i=1:ra_factor
    ra_spectra(:,:,k,:) = spectra(:,1+(i-1)*REARRANGE:i*REARRANGE,j,:);
    k=k+1;
   end
  end
  Npulses = REARRANGE;
  Nparam1 = Nparam1*ra_factor;
  param1Name = [ param1Name '_rearranged!'];
  param1Step = param1Step / ra_factor;
  spectra = ra_spectra;
  spectra_is_valid = ones(REARRANGE, Nparam1 * ra_factor, Nparam2,'int8');
end

%ojo!!hay que revisar esto porque algo falla en la definición del eje
if (REARRANGE > 0)
     ejeParam1 = param1From:param1Step:Nparam1*param1Step+param1From-1; % he centralizado aquí la definición del eje horizontal en v04 para cuadrarlo con lo de REARRANGE
else
     ejeParam1 = param1From:param1Step:Nparam1*param1Step+param1From-1; % 1setp2019 gave error, I don't know if changing this have colateral damages
end
ejeParam1 = param1From:param1Step:param1To;
ejeParam2 = param2From:param2Step:param2To;

ratio_MgCa_todos = NaN(Npulses,Nparam1,Nparam2); % se calcula el ratio espectro a espectro
intensidad_Mg = NaN(Npulses,Nparam1,Nparam2); % esto es para el cálculo del SNR
intensidad_Ca = NaN(Npulses,Nparam1,Nparam2);  
intensidad_backg_Mg = NaN(Npulses,Nparam1,Nparam2); 
intensidad_backg_Ca = NaN(Npulses,Nparam1,Nparam2); 
SNR_Mg = NaN(Nparam1,Nparam2); % signal to noise ratio of both lines
SNR_Ca = NaN(Nparam1,Nparam2); 
SNB_Mg = NaN(Nparam1,Nparam2); % signal to BACKGROUND ratio, as the all-shots average at the middle WL in the background4xx range
SNB_Ca = NaN(Nparam1,Nparam2); 
ratio_MgCa_promediados = zeros(Nparam1,Nparam2); % primero todos los ratios con espectros individuales, luego se promedian los ratios
ratio_MgCa_del_espectro_promedio = zeros(Nparam1,Nparam2); % ratio de los espectros promediados previamente
ratio_MgCa_del_promedio_de_ratios = zeros(Nparam1,Nparam2); %ratio obtenido de todos los ratios individuales
ratio_resina_del_promedio_de_ratios = zeros(Nparam1,Nparam2); % primero todos los ratios con espectros individuales, luego se promedian los ratios
rsd_ratio2l_bruto = NaN(Nparam1,Nparam2);
rsd_linea_bruto = NaN(Nparam1,Nparam2);
sem_ratio2l_bruto = NaN(Nparam1,Nparam2);
std_errorbars = NaN(Nparam1,Nparam2); % guarda el STD para las barras de error
sem_errorbars = NaN(Nparam1,Nparam2); % guarda el Standard Error of the Mean para las barras de error
intensidad_todos = zeros(Npulses,Nparam1,Nparam2); % se guarda la intensidad integrada del espectro para descartar outliers
intensidad_outliers = zeros(Npulses,Nparam1,Nparam2); % diferencia entre la intensidad del espectro y el promedio de todos
mlineaMg = zeros(Nparam1,Nparam2); % valor medio de la intensidad de la línea en cada punto de medida
mlineaCa = zeros(Nparam1,Nparam2); % valor medio de la intensidad de la línea en cada punto de medida
pixelRSD = buscaPixel(lambdas,LAMBDA_RSD(1),mean(mean(spectra,2),3),BUSCAR_PIXEL); % devuelve el píxel más cercano a ese valor, luego comprueba que sea un máximo
tmp_spectra = mean(spectra,2);
pixelResina = buscaPixel(lambdas,LAMBDA_RESINA(1),tmp_spectra(:,1),2); % igual hay que poner flag=2 para que no se mueva porque el pico al principio puede ser inexistente si se mide desde el ápice y se acaba en la resina. Se hace con el primer punto suponiendo está en la resina
mean_spectra_only_valids=zeros(Npixeles,Nparam1,Nparam2);   % mean spectra using only valid spectra

% specific for CF-LIBS
 kB_eVK = 8.6173324E-5;
 kB_JK = 1.38064852E-23;
 h_Js = 6.62607E-34;
 h_eVs = 4.135667E-15;
 me_MeVc2 = 0.511E6/9E16;
 me_gr = 9.1093836E-31 * 1E3;
 peaks_isvalid=zeros(CFNlines,Npulses,Nparam1,Nparam2); % if a peak model is valid or not
 peaks_heights=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the calculated peaks height
 peaks_fwhm=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the calculated FWHM of peaks
 peaks_area=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the calculated integrated intensity of peaks
 peaks_stark=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the stark broadenings
 peaks_central=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the fitted central lambdas
 peaks_selfabscoeff=zeros(CFNlines,Npulses,Nparam1,Nparam2); % selfabsorption coefficients (Sun2009IRSAC or PRaher2010)
 peaks_Ne_fromStark=zeros(CFNlines,Npulses,Nparam1,Nparam2); % all the Ne from stark broadening
 % calculado para el espectro promedio en Npulses
 peaks_heights_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the calculated peaks height
 peaks_heights_frommean_SAcorrected=zeros(CFNlines,Nparam1,Nparam2); % all the calculated peaks height, after SA correction
 peaks_fwhm_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the calculated FWHM of peaks
 peaks_intensity_frommean=zeros(CFNlines,Nparam1,Nparam2); % the intensity, depending on LINE_INTENSITY, it could be the height, area...  used for BP, SBP, ratios ...
 peaks_twohalfareas_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the calculated integrated instensity of peaks, two times the right half area (to overcom asymmetry problem)
 peaks_area_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the calculated integrated instensity of peaks
 peaks_stark_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the stark broadenings of mean spectra
 peaks_central_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the fitted central lambdas
 peaks_fit_error_frommean=zeros(CFNlines,Nparam1,Nparam2); % the error (relative squared sum of residuals) of the fitting
 peaks_fit_pixleftbkg1=zeros(CFNlines,Nparam1,Nparam2); % the pixel to the left used for background substraction (area of the line)
 peaks_fit_pixrightbkg1=zeros(CFNlines,Nparam1,Nparam2); % the pixel to the right used for for background substraction (area of the line)
 peaks_fit_pixleft=zeros(CFNlines,Nparam1,Nparam2); % the pixel to the left used for peak fitting
 peaks_fit_pixright=zeros(CFNlines,Nparam1,Nparam2); % the pixel to the right used for peak fitting
 peaks_Ne_fromStark_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the Ne from stark broadening
 peaks_Ne_frommean=zeros(Nparam1,Nparam2); % could be from Halfa, other peaks' stark, or anything
 peaks_Ne_fromLineRatio_frommean=zeros(CFNlines,Nparam1,Nparam2); % all the Ne from ????
 peaks_selfabscoeff_frommean=zeros(CFNlines,Nparam1,Nparam2); % selfabsorption coefficients (Sun2009IRSAC or Praher2010) <1, DIVIDE
 % results from individual spectra, 
 slope_BP=zeros(Nspecies,Npulses,Nparam1,Nparam2); % bolztmann plot slope
 intercept_BP=zeros(Nspecies,Npulses,Nparam1,Nparam2); % intercept qs of bolztmann plot 
 NumPeaks_BP=zeros(Nspecies,Npulses,Nparam1,Nparam2); %number of peaks used in each BP
 Te_BP=zeros(Nspecies,Npulses,Nparam1,Nparam2); % calculated Te with boltzmann-plot 
 Te_SahaBP=zeros(Nspecies,Npulses,Nparam1,Nparam2); % calculated Te with saha-boltzmann-plot 
 %results from mean spectrum
 peaks_isvalid_frommean=zeros(CFNlines,Nparam1,Nparam2); % if a peak model is valid or not
 slope_BP_frommean=zeros(Nspecies,Nparam1,Nparam2); % bolztmann plot slope
 slope_single_BP_frommean=zeros(Nparam1,Nparam2); % singles slopes from BP of all the species at once
 intercept_single_BP_frommean=zeros(Nspecies,Nparam1,Nparam2); % bolztmann plot slope for each species BUT with common slope
 intercept_BP_frommean=zeros(Nspecies,Nparam1,Nparam2); % intercept qs of bolztmann plot 
 Te_BP_frommean=zeros(Nspecies,Nparam1,Nparam2); % calculated Te with boltzmann-plot 
 Te_single_BP_frommean=zeros(Nparam1,Nparam2); % calculated Te with a single boltzmann-plot of all the species
 NumPeaks_BP_frommean=zeros(Nspecies,Nparam1,Nparam2); %number of peaks used in each BP
 slope_SahaBP_frommean=zeros(Nspecies,Nparam1,Nparam2); % bolztmann plot slope
 slope_single_SahaBP_frommean=zeros(Nparam1,Nparam2); % saha-bolztmann plot slope with single Te
 intercept_SahaBP_frommean=zeros(Nspecies,Nparam1,Nparam2); % intercept qs of saha bolztmann plot 
 intercept_single_SahaBP_frommean=zeros(Nspecies,Nparam1,Nparam2); % intercept qs of saha bolztmann plot  single Te
 Te_SahaBP_frommean=zeros(Nspecies,Nparam1,Nparam2); % calculated Te with saha-boltzmann-plot 
 Te_single_SahaBP_frommean=zeros(Nparam1,Nparam2); % calculated Te with saha-boltzmann-plot single Te
 NumPeaks_SahaBP_frommean=zeros(Nspecies,Nparam1,Nparam2); %number of peaks used in each SahaBP
 peaks_area_frommean_fitted_for_reference_line=zeros(Nspecies,Nparam1,Nparam2); % fitted values of IRSAC reference lines, after first BP
 
 %concentrations and ratios
 concentration=zeros(Nspecies,Nparam1,Nparam2); %
 concentration_z1=zeros(Nspecies,Nparam1,Nparam2); % with z=1 substitution
 concentration_weight  =zeros(Nspecies,Nparam1,Nparam2); % in weight
 concentration_error = zeros(Nparam1,Nparam2); % in weight
 ratio_MgCa_CFLibs = zeros(Nparam1,Nparam2); % ratio from CF calculations (mmols/mol)
 opts = optimset('Display','off');

% se calculan los píxeles de cada lambda
    pixelesCa = zeros(size(LAMBDA_CALCIO,2),1);pixelesMg = zeros(size(LAMBDA_MAGNESIO,2),1);
    for n=1:size(LAMBDA_CALCIO,2)
      pixelesCa(n) = buscaPixel(lambdas,LAMBDA_CALCIO(n),mean(mean(spectra,2),3),BUSCAR_PIXEL);
    end
    pixelesMg = zeros(size(LAMBDA_MAGNESIO,2),1);
    for n=1:size(LAMBDA_MAGNESIO,2)
      pixelesMg(n) = buscaPixel(lambdas,LAMBDA_MAGNESIO(n),mean(mean(spectra,2),3),BUSCAR_PIXEL);
    end

    % estas lambdas se usan con INTEGRAR=3 para calcular la recta de
    % regresión del background
    pixeles_background4lambdas_mg = zeros(size(background4lambdas_mg));
    pixeles_background4lambdas_ca = zeros(size(background4lambdas_ca));
    pixeles_background4lambdas_resina = zeros(size(background4lambdas_resina));
    for lin=1:size(background4lambdas_mg,1)
        for j=1:size(background4lambdas_mg,2)
           pixeles_background4lambdas_mg(lin,j) =    buscaPixel(lambdas,background4lambdas_mg(lin,j),mean(mean(spectra,2),3),2); % no se mueve, es background
        end
    end
    for lin=1:size(background4lambdas_ca,1)
        for j=1:size(background4lambdas_ca,2)
           pixeles_background4lambdas_ca(lin,j) =    buscaPixel(lambdas,background4lambdas_ca(lin,j),mean(mean(spectra,2),3),2); % no se mueve, es background
        end
    end
    for lin=1:size(background4lambdas_resina,1)
        for j=1:size(background4lambdas_resina,2)
           pixeles_background4lambdas_resina(lin,j) =    buscaPixel(lambdas,background4lambdas_resina(lin,j),mean(mean(spectra,2),3),2); % no se mueve, es background
        end
    end
    
spectra_original = spectra; % before changing spectra with background removal and irradiance correction, save the original one just in case to verify things

%% CF-Step 1   calc firstLaserpulse for each spatial point, to remove background with more precision
%antes de preprocesar voy a calcular el valor de ncleaningshots con precisión porque el ruido crece linealmente hasta el disparo 6, para no fallar porque creo que quitar el ruido aditivo es importante
%Ncleaningshots = NpulsesClean ; % más o menos el laser empieza a disparar 2 pulsos antes de lo especificado, por si acaso
%modificado v13: esto es importante en los espectros de la PIMAX3, pero el
%primero llega a valores de 11000.
 disp(['CF-Step 1: calc of firstLaserPulse for every spatial point...']);    
 maxspectravalues=zeros(Npulses,Nparam1,Nparam2);
 firstLaserpulse=zeros(Nparam1,Nparam2);
 for idx2= 1:Nparam2
  for idx1= 1:Nparam1
    for npulse=2:Npulses
      maxspectravalues(npulse,idx1,idx2)=max(spectra(:,npulse,idx1,idx2));
    end
    f = find(maxspectravalues(:,idx1,idx2)>ValMinHeightAllPeaks,1);  % no intensity above this value means noise spectrum
    if (~isempty(f))
        firstLaserpulse(idx1,idx2) = f;
    else
        firstLaserpulse(idx1,idx2) = 1; % actually there is no spectra
    end
  end
 end
 
%% CF-Step 2   remove background with the electronic (additive) noise of the first spectra, METHOD #1
% preprocesado para quitar el background AQUI, el método 1 usa los primeros pulsos sin láser (NpulsesClean)
% la amplitud del ruido va creciendo asintóticamente, con Npulsesclean=7
% (5) no se estabiliza todavía
% v12 14dic2017 el ruido crece linealmente, no es estable hasta el 6, y no
% se garantiza que empiece a disparar exactamente en NpulsesClean+1
if (REMOVE_BACKGROUND== 1) 
 if (USAR_PIMAX3==1)
    first_noise = 6;
 else
    first_noise = 1; % only pimax3 has this effect of high-noise at the first pulses
 end
 disp(['CF-Step 2: Removing background using electronic noise from first spectra starting with ' num2str(first_noise) '...']);    
 
 % calculamos un background medio con todos los espectros de ruido válidos
 specback = zeros (size(lambdas,1),1);
 num_specbacks = 0;
 for idx2= 1:Nparam2
  for idx1= 1:Nparam1
   n_espectros = firstLaserpulse(idx1,idx2) - first_noise;
   if (n_espectros > (NpulsesClean-first_noise-2) && n_espectros < (NpulsesClean-first_noise+3)) % el número de espectros de ruido es razonable
     tmp_mean_noise =   sum(spectra(:,first_noise:firstLaserpulse(idx1,idx2)-1,idx1,idx2),2);
      if max(tmp_mean_noise) > 3000 * n_espectros
          disp('Warning: background with peaks, something is wrong with the estimation of the first few spectra with no LIBS emission');
      end
     specback = specback + tmp_mean_noise;
     num_specbacks = num_specbacks + n_espectros;
   end
  end
 end
 
 specback = specback / num_specbacks;
 Ncleaningshots = mode(firstLaserpulse(find(firstLaserpulse>0))); %% warning! this single value is the MOST frequent, it could be different in every spatial point due to timming issues in the q-switch controller
 for idx2= 1:Nparam2
  for idx1= 1:Nparam1
      specback = nanmeanLocal(spectra(:,first_noise:firstLaserpulse(idx1,idx2)-1,idx1,idx2),2); % electronic noise averaged
      for npulse=firstLaserpulse(idx1,idx2):Npulses
        if (spectra_is_valid(npulse,idx1,idx2)) % en principio no está descartado
         %dos opciones, restar el espectro de ruido promedio global o el local para este disparo ¿?¿
         %tmp_spectra = spectra(:,npulse,idx1,idx2) - specback;
         tmp_spectra = spectra(:,npulse,idx1,idx2) - specback;
         %plot(lambdas,specback,lambdas,spectra(:,npulse,idx1,idx2),lambdas,tmp_spectra);axis([ lambdas(1) lambdas(end) 0 6000]);
         spectra(:,npulse,idx1,idx2) = tmp_spectra;
         
        else
          %spectra(:,npulse,idx1,idx2) = zeros(Npixeles,1); verificado
          %29mayo2015 que no se usan luego para procesar
        end   
      end
  end
 end
 disp(['Removed background substracting ' num2str(num_specbacks) ' noise spectra']);
 
 B=spectra;
 
end

 %% Divides each spectrum by its total intensity. Warning! this change the valSaturation value and makes it different for each spectrum 
if (NORMALIZAR==1)
 disp('Total area normalization for each spectrum...');
 

 SUMATOTAL=(sum(spectra,1));
 meanIntegratedIntensity = nanmeanLocal(nanmeanLocal(SUMATOTAL(1,10:end,:),2),3); % quick and dirty average of valid spectra
 spectra=spectra./ SUMATOTAL .* meanIntegratedIntensity; % this way de overal vertical scale is not changed so much, yet it is normalized
 
 elseif(NORMALIZAR==2)
    
end

%% CF-Step 3   Irradiance correction 
% irradiance correction, the actual correction is perfomed now because substraction of additive noise should be performed with the raw spectra 

if (CORREGIR_IRRADIANCIA == 1 || DO_CF==1) 
 disp('CF-Step 3: Irradiance correction...');    
 for idx2= 1:Nparam2
  for idx1= 1:Nparam1
   for npulse=firstLaserpulse(idx1,idx2):Npulses
    %irradiance response correction using DH2000 data 
    spectra(:,npulse,idx1,idx2) = spectra(:,npulse,idx1,idx2) ./ dividir_respuesta_espectrometro;
   end
  end
 end
end

% wavelength sorted, too close lines discarded, regardless of being used  
  whichLambdas=[];whichPeaks=[];
  for s=composition
    whichLambdas = vertcat(whichLambdas,  pila(piquees==s));
    whichPeaks = vertcat(whichPeaks,  find(piquees==s));
  end
  ordered = sortrows([whichLambdas whichPeaks]);
  for i=1:size(ordered,1)-1
      if (abs(ordered(i,1)-ordered(i+1,1)) < 3*AvantesFWHM(lampara_FWHM,pila(ordered(i,2))) && piUse(ordered(i,2))==1 && piUse(ordered(i+1,2))==1 ) % three times the FHWM is the threshold, it is about 0.2.0.3nm apart for the Avantes. Risk of overlapping or misclasification
         if (REMOVE_OVERLAPPED_LINES==0)
           disp(['  Warning!! two overlapped lines IN USE: ', strjoin(speciesStr(piquees(ordered(i,2)))), '@',num2str(pila(ordered(i,2))), ' and ' ,  strjoin(speciesStr(piquees(ordered(i+1,2)))), '@',num2str(pila(ordered(i+1,2))), ' not discarded (REMOVE_OVERLAPPED=0)']);
         else 
           disp(['  Warning!! two overlapped lines IN USE: ', strjoin(speciesStr(piquees(ordered(i,2)))), '@',num2str(pila(ordered(i,2))), ' and ' ,  strjoin(speciesStr(piquees(ordered(i+1,2)))), '@',num2str(pila(ordered(i+1,2))) , ' DISCARDED!!']);
           piUse(ordered(i,2)) =0;
           piUse(ordered(i+1,2)) =0;
         end
      end
  end
  
% also for peaks in lambdas with sensitivity very low ("dividir_..." < 0.005)
% a warning is shown
if (CORREGIR_IRRADIANCIA>0) % only with irradiance correction the sensitivy of any lambda is known
 for i=1:size(ordered,1)
      px = lambda2pixel(lambdas,pila(ordered(i,2)));
      if (which_lambdas_have_too_low_sensitivity(px)==1)
         disp(['  Warning!! peak ', strjoin(speciesStr(piquees(ordered(i,2)))), '@',num2str(pila(ordered(i,2))), ' is at an channel edge with poor sensitivity - irradiance correction could not be Ok']);
      end
 end
end

% pila() is substituted by pilaobs() for shifted peaks
for p=1:CFNlines
    if(pilaobs(p)~=0)
      pila(p) = pilaobs(p);
    end
end


%% CF-Step 4.0   Mark as non-valid for weak or saturated spectra and mean spectrum
% we don't know in advance which species (peaks) are going to analyze. 
% One possibility is to look for ALL peaks in the CF list. if at least one
% is saturated, mark as not valid. if the height of the highest is under a
% threshold, the entire spectrum is week, mark as not valid.
% new v16! the mean is calculate from FIRST_SHOT_FOR_AVERAGING to LAST_SHOT_FOR_AVERAGING
% new v16! we indeed know ('composition') which lines are we going to use.
% for spectra with peaks always saturated, you didn't get any.
whichLambdas=pila(piUse>0); % only peaks in composition, other peaks has been set to piUse(i)=0
whichPeaks = find(piUse>0);
doNotUse = zeros(size(piUse));
doNotUse(:)=1; % let's start with all peaks Ok
pixels2look_orig = lambda2pixel(lambdas,whichLambdas'); % transpose, if not, non-singleton error  
discardedSpectraSat=zeros(Nparam1,Nparam2);
discardedSpectraAllLow=zeros(Nparam1,Nparam2);
discardedSpectraSinglePeakLow=zeros(Nparam1,Nparam2);
totalValidSpectra=zeros(Nparam1,Nparam2);
totalSaturatedPeaks=zeros(Nparam1,Nparam2,size(whichPeaks,1)); % count of accumulated discarded peaks due to at least one saturated
totalAllLowPeaks=zeros(Nparam1,Nparam2,size(whichPeaks,1)); % stores the npulses of spectra that has all peaks below a threshold 
totalSinglePeakLow=zeros(Nparam1,Nparam2,size(whichPeaks,1)); % stores the npulses of spectra that has at least one peak below a threshold

% creates strings to name emission lines
lineName = {};
for i=1:size(whichPeaks)
     lineName{i} = [ '#' num2str(whichPeaks(i)) ' ', strjoin(speciesStr(piquees(whichPeaks(i)))), '@',num2str(pila(whichPeaks(i)))];
end

if (DISCARD_NONVALID_SPECTRA > 0) 
 fprintf('CF-Step 4: Marking as non-valid weak or saturated spectra (only CF PEAKS in composition checked):\n');
       if (FIRST_SHOT_FOR_AVERAGING>1 || LAST_SHOT_FOR_AVERAGING < Npulses)
           disp(['  Warning!! averaging spectrum only for ', num2str(FIRST_SHOT_FOR_AVERAGING), ':' ,num2str(min(LAST_SHOT_FOR_AVERAGING,Npulses)),' pulses (out of ' ,num2str(Npulses), ' pulses).']);
      end
 for idx2= 1:Nparam2  
  for idx1= 1:Nparam1
      which_spectrum_has_saturated_peaks = {}; 
      for npulse=1:Npulses
        if (spectra_is_valid(npulse,idx1,idx2)) % not previously marked, should be all valid
          %pixels2look are "tuned" looking for higher neighbour pixels, some pixels are on the peak wings well below the max
          pixels2look = pixels2look_orig;
          for i=1:size(pixels2look_orig,2)
            [maxl,maxi] = max(mean_spectra_only_valids(pixels2look_orig(i)-2:pixels2look_orig(i)+2,idx1,idx2));
            pixels2look(i) = pixels2look_orig(i) + maxi - 3; % the pixel around central with max value
          end
          peaksSat =  spectra_original(pixels2look,npulse,idx1,idx2) > ValSaturation;
          numPeaksSat = sum(peaksSat==1);
          allpeaksLow =  spectra_original(pixels2look,npulse,idx1,idx2) < ValMinHeightAllPeaks;
          numPeaksNotLow = sum(allpeaksLow==0);
          singlepeakLow =  spectra_original(pixels2look,npulse,idx1,idx2) < ValMinHeightSinglePeak; % added 19mar2020 v20 
          numPeaksLow = sum(singlepeakLow==1); % number of peaks below the thresdhold
          if (numPeaksSat > 0) % at least one peak in the ORIGINAL non-processed spectrum is saturated
              totalSaturatedPeaks(idx1,idx2,:) = totalSaturatedPeaks(idx1,idx2,:) + reshape(peaksSat,[1 1 size(peaksSat,1)]);
              discardedSpectraSat(idx1,idx2) = discardedSpectraSat(idx1,idx2)  +1 ; % not really discarded, only marked
              which_spectrum_has_saturated_peaks {discardedSpectraSat(idx1,idx2)} = npulse; 
          end
          if (numPeaksNotLow == 0) % ALL of the peaks in the ORIGINAL non-processed spectrum are low
              totalAllLowPeaks(idx1,idx2,:) = totalAllLowPeaks(idx1,idx2,:) + reshape(allpeaksLow,[1 1 size(allpeaksLow,1)]);
              discardedSpectraAllLow(idx1,idx2) = discardedSpectraAllLow(idx1,idx2) +1 ;
              spectra_is_valid(npulse,idx1,idx2)=0; 
          end
          if (numPeaksLow > 0 && DISCARD_ALSO_WEAK_PEAKS && spectra_is_valid(npulse,idx1,idx2)~=0) % we want to discard the entire spectrum because at least one peak is too low, ONLY COUNTS SPECTRA NOT DISCARDED AS NOISE
              totalSinglePeakLow(idx1,idx2,:) = totalSinglePeakLow(idx1,idx2,:) + reshape(singlepeakLow,[1 1 size(singlepeakLow,1)]);
              discardedSpectraSinglePeakLow(idx1,idx2) = discardedSpectraSinglePeakLow(idx1,idx2) +1 ;
              spectra_is_valid(npulse,idx1,idx2)=0; 
          end
          if (spectra_is_valid(npulse,idx1,idx2)~=0) % everything's ok
              totalValidSpectra(idx1,idx2) = totalValidSpectra(idx1,idx2) + 1;
          end
          if (spectra_is_valid(npulse,idx1,idx2) ==0) % se ha descartado ahora
              %disp(['  Not valid spectrum(' , num2str(npulse),',',num2str(idx1),',',num2str(idx2),')'  ]);
              %fprintf([num2str(npulse) '-' num2str(idx1) '-' num2str(idx2) ' ']);
          end
        end   
      end % npulse
      % all the Npulses evaluated, are most of them discarded (LOW)?
      if sum(spectra_is_valid(:,idx1,idx2)) < Npulses/5 % 80% of spectra in this spatial points have been discarded
          disp(['  Warning!! only ' , num2str(sum(spectra_is_valid(:,idx1,idx2))) ,' valid spectra in spatial point (' , num2str(idx1),',',num2str(idx2),')'  ]);
      end
      
      % new v20 19mar2020: if AVERAGE_ONLY_FIRST > 0 only the N first valid spectra are marked as valid
      if AVERAGE_ONLY_FIRST > 0 
        for n=1:Npulses
          if spectra_is_valid(n,idx1,idx2) 
            for  i=n+AVERAGE_ONLY_FIRST:Npulses
              spectra_is_valid(i,idx1,idx2)=0; % mark as nonvalid everything else
            end
            break
          end
        end
      end

      % what about individual peaks with some of them saturated?
      for i=1:size(totalSaturatedPeaks,3)
          if (totalSaturatedPeaks(idx1,idx2,i)>0) % at least one peak was saturated
              if (discardedSpectraSat(idx1,idx2) > totalValidSpectra(idx1,idx2)) % the number of spectra with saturated peaks is more than 50%, it is better to avoid using this peak...
                disp(['  Discarding line #' num2str(whichPeaks(i)) ' of ', strjoin(speciesStr(piquees(whichPeaks(i)))), '@',num2str(pila(whichPeaks(i))), ' (',num2str(idx1),',',num2str(idx2),') with ' , num2str(discardedSpectraSat(idx1,idx2)), ' saturated spectra out of ', num2str(discardedSpectraSat(idx1,idx2)+totalValidSpectra(idx1,idx2))]);
                doNotUse(whichPeaks(i))=0; 
              else % ... otherwise, it is better to unv-alid spectra
                str_spectra = '';
                for n=1:totalSaturatedPeaks(idx1,idx2,i)
                    spectra_is_valid(which_spectrum_has_saturated_peaks{n},idx1,idx2)=0;
                    totalValidSpectra(idx1,idx2) = totalValidSpectra(idx1,idx2) - 1; % was not marked as invalid before
                    str_spectra = strcat(str_spectra,num2str(which_spectrum_has_saturated_peaks{n}));
                end
                disp(['  Discarding spectra (', str_spectra, ') of spatial point (' ,num2str(idx1),',',num2str(idx2),') because ' strjoin(speciesStr(piquees(whichPeaks(i)))), '@',num2str(pila(whichPeaks(i))), ' (',num2str(idx1),',',num2str(idx2),') has ' , num2str(discardedSpectraSat(idx1,idx2)), ' saturated spectra out of ', num2str(discardedSpectraSat(idx1,idx2)+totalValidSpectra(idx1,idx2))]);
              end
              %disp('not discarding actually');
              
          end
      end
      piUse = piUse & doNotUse;
      spectra_is_valid(1:FIRST_SHOT_FOR_AVERAGING-1,idx1,idx2)  = 0;
      spectra_is_valid(LAST_SHOT_FOR_AVERAGING+1:Npulses,idx1,idx2)  = 0;
      mean_spectra_only_valids(:,idx1,idx2) = nanmeanLocal(spectra(:,spectra_is_valid(:,idx1,idx2)==1,idx1,idx2),2);  
       % 17mar2020 debugging information
    if (DEBUG_DISCARDING==idx1) % shows statistics of spectra and peaks for this particular spatial point
     figure(41);
     bar([ squeeze(totalAllLowPeaks(idx1,idx2,:)) squeeze(totalSinglePeakLow(idx1,idx2,:)) squeeze(totalSaturatedPeaks(idx1,idx2,:))] ,'stacked');
     title({['(idx1=', num2str(idx1) , ') Number of times a peak is < ValMinHeightAllPeaks (blue)'],'< ValMinHeightSinglePeak (red), saturated (other)'});
     set(gca,'xtick',[1:size(whichPeaks)],'xticklabel',lineName)
     xtickangle(90);
     if (size(find(spectra_is_valid(:,idx1,idx2)==1),1)>0)
       figure(42);
       plot(lambdas,spectra_original(:,spectra_is_valid(:,idx1,idx2)==1,idx1,idx2));
       title(['Valid spectra (',num2str(size(find(spectra_is_valid(:,idx1,idx2)==1),1)),') to average for idx1=' num2str(idx1)]);
     end
     if (size(find(spectra_is_valid(:,idx1,idx2)==0),1)>0)
       figure(43);
       plot(lambdas,spectra_original(:,spectra_is_valid(:,idx1,idx2)==0,idx1,idx2));
       title(['NON valid discarded spectra (', num2str(size(find(spectra_is_valid(:,idx1,idx2)==0),1)),') for idx1=' num2str(idx1)]);
     end
     figure(44);
     imagesc(1:Npulses,1:size(whichPeaks),spectra_original(pixels2look,:,idx1,idx2));
     set(gca,'ytick',[1:size(whichPeaks)],'yticklabel',lineName)
     xlabel('number of pulse');
     title( ['Raw height value of peaks for every pulse (idx1=' , num2str(idx1) ,')']);
     colormap(pink);
     caxis( [ 1000 3000]);
   end

  end % idx1
 end % idx2
 
 if (DEBUG_DISCARDING>0) % shows plot of number of valid spectra for all idx1, idx2=0
     figure(40);
     clf;
     hold on;
     plot(ejeParam1,discardedSpectraSat(:,1),'LineWidth',2,'Color','r');    
     plot(ejeParam1,discardedSpectraAllLow(:,1),'LineWidth',2,'Color','b');    
     plot(ejeParam1,discardedSpectraSinglePeakLow(:,1),'LineWidth',2,'Color','k');    
     plot(ejeParam1,totalValidSpectra(:,1),'LineWidth',2,'Color','g');    
     title({'Number of spectra with:','at least one saturated peak from composition (red)','all peaks zero (blue)','single peak below threshold (black)','valid spectra (green)'});   
     xlabel(param1Name);
     ylabel('number of spectra out of ');
     hold off;
 end
end

%% CF-Step 5   substract spectra if ALTERNATE_DELAY_WHICH_ONE!=0
% ALTERNATE_DELAY_WHICH_ONE=1 every alternate shot is at different delay. 
% after this step, each spatial point has individual difference spectra
% every two shots, and mean_only_valids() has the MEAN difference spectrum
if ALTERNATE_DELAY_WHICH_ONE>0 && USAR_AVANTES==1
 fprintf('CF-Step 5: subtracting MEAN spectra for ALTERNATE_DELAY_WHICH_ONE=');
 if ALTERNATE_DELAY_WHICH_ONE==1
      fprintf('1 (Substract hotter-colder)\n');
 elseif ALTERNATE_DELAY_WHICH_ONE==2
      fprintf('2 (only hotter)\n');     
 elseif ALTERNATE_DELAY_WHICH_ONE==3
      fprintf('3 (only colder)\n');     
 else
      fprintf('??? unknown\n');     
 end
 for idx2= 1:Nparam2  
  for idx1= 1:Nparam1
      hotter_spectra = zeros(size(lambdas,1),1);
      colder_spectra = zeros(size(lambdas,1),1);
      npfirst = max(firstLaserpulse(idx1,idx2),FIRST_SHOT_FOR_AVERAGING);
%      hotter_intensity = sum(spectra(:,npfirst,idx1,idx2),1);
%      colder_intensity = sum(spectra(:,npfirst+1,idx1,idx2),1);
%      if (hotter_intensity > colder_intensity)
%        ; % the first pulse is the hotter (shorter delay) 
%      else
%        npfirst = npfirst + 1; % start always with hotter spectra (shorter delay, larger & hotter spectrum)
%      end
      num_spectra = 0;
      for npulse=npfirst:2:min(LAST_SHOT_FOR_AVERAGING,Npulses)-1
          if (spectra_is_valid(npulse,idx1,idx2)==1 && spectra_is_valid(npulse+1,idx1,idx2)==1)
            diff_spectrum = spectra(:,npulse,idx1,idx2)-spectra(:,npulse+1,idx1,idx2);
            hotter_spectra = hotter_spectra + spectra(:,npulse,idx1,idx2);
            colder_spectra = colder_spectra + spectra(:,npulse+1,idx1,idx2);
            num_spectra = num_spectra+1;
          else
            disp(['  Spectrum (' , num2str(npulse),' or ', num2str(npulse+1),',',num2str(idx1),',',num2str(idx2),') is not valid within AlternateDelay averaging and differencing range, not used'  ]);
          end
      end
      % calc difference of MEAN hotter and colder spectra. It is safer to
      % get the substract right with the mean spectra, not with the two
      % first
      hotter_intensity = sum(hotter_spectra,1);
      colder_intensity = sum(colder_spectra,1);
      if (ALTERNATE_DELAY_WHICH_ONE==1)
        if (hotter_intensity > colder_intensity)
         mean_spectra_only_valids(:,idx1,idx2) = (hotter_spectra - colder_spectra) ./ num_spectra;
        else
         mean_spectra_only_valids(:,idx1,idx2) = (colder_spectra - hotter_spectra) ./ num_spectra; 
         spectra_is_valid(npfirst,idx1,idx2)=0; % this cannot be substracted
         npfirst = npfirst+1; % the hotter one is the next one
        end
      elseif (ALTERNATE_DELAY_WHICH_ONE==2)
           mean_spectra_only_valids(:,idx1,idx2) = hotter_spectra ./ num_spectra;
      elseif (ALTERNATE_DELAY_WHICH_ONE==3)
           mean_spectra_only_valids(:,idx1,idx2) = colder_spectra ./ num_spectra;
      end
      if (0)
          figure;
          subplot(3,1,1);   plot(lambdas,mean_spectra_only_valids(:,idx1,idx2));title('short delay spectrum - hotter');
          subplot(3,1,2);   plot(lambdas,mean_spectra_only_valids(:,idx1+1,idx2));title('long delay spectrum - colder');
          subplot(3,1,3);   plot(lambdas,mean_spectra_only_valids(:,idx1,idx2)- mean_spectra_only_valids(:,idx1+1,idx2));title('subtracted spectrum');
      end
      % now that we know which spectrum is the hotter and the colder, let's
      % substract per-pair. needed to calculate intensity
      if (ALTERNATE_DELAY_WHICH_ONE==1)
       for npulse=npfirst:2:min(LAST_SHOT_FOR_AVERAGING,Npulses)-1
         diff_spectrum = spectra(:,npulse,idx1,idx2)-spectra(:,npulse+1,idx1,idx2);
         spectra(:,npulse,idx1,idx2) = diff_spectrum;
         spectra(:,npulse+1,idx1,idx2) = diff_spectrum;
       end
       elseif (ALTERNATE_DELAY_WHICH_ONE==2)
        for npulse=npfirst:2:min(LAST_SHOT_FOR_AVERAGING,Npulses)-1
         spectra(:,npulse+1,idx1,idx2) = spectra(:,npulse,idx1,idx2); % duplicates hotter
        end
       elseif (ALTERNATE_DELAY_WHICH_ONE==3)
        for npulse=npfirst:2:min(LAST_SHOT_FOR_AVERAGING,Npulses)-1
         spectra(:,npulse,idx1,idx2) = spectra(:,npulse+1,idx1,idx2); % duplicates colder
        end
      end
   end
 end
 
 
end


%% Figure 1 with all the spectra in 2D, to check there are useful spectra
if (0)
me = spectra;
 for idx2= 1:Nparam2
  for idx1= 1:Nparam1
   for npulse=1:Npulses   
     me(:,npulse,idx1,idx2) = me(:,npulse,idx1,idx2)/max(me(:,npulse,idx1,idx2));
   end
  end
 end
figure(1);
me_promedio = squeeze(mean(me,2));
h=pcolor(me_promedio);set(h, 'EdgeColor', 'none');
title([PNbarras(size(PNbarras,2)-25:size(PNbarras,2)) para_el_titulo]);   
end


%% CF-Step 6   remove baseline (emission background) with different algorithms 

 if (REMOVE_BASELINE== 1)  % method Yu2009 
  disp(['CF-Step 6: Removing baseline (emission background) with Yu2009 method...']);      
  for idx2= 1:Nparam2
   for idx1= 1:Nparam1
     if ONLY_CF_MEAN==0 % should be applied to ALL individual spectra
       for npulse=1:Npulses
        if (spectra_is_valid(npulse,idx1,idx2)) % it is not discarded yet
          tmp_spectra = remove_background_Yu2009(lambdas, spectra(:,npulse,idx1,idx2),0); % last argument could be 1 to plot the corrected spectrum 
          spectra(:,npulse,idx1,idx2) = tmp_spectra;
        end
       end
     else % mean spectra has been already obtained
       tmp_spectra = remove_background_Yu2009(lambdas, mean_spectra_only_valids(:,idx1,idx2),0); % last argument could be 1 to plot the corrected spectrum 
       mean_spectra_only_valids(:,idx1,idx2) = tmp_spectra;
     end   
   end
  end
 end

%% CF-LIBS calculations

if ONLY_CF_MEAN==0
%% CF-Step 7  Calcs everything for every valid single spectrum. Output variables are peaks_????
  disp('the per-spectrum algorithm was so outdated that it has been removed 1setp2019, should be copied and adapted from mean algorithm');
end
  
%% CF-Step 7   Calcs everything AGAIN for each MEAN spectrum at every spatial point. Output variables are peaks_??_frommean
% 
% Now, we know wich peaks are valid, Everything again but using MEAN spectra of ONLY VALID peaks
% wait! which peaks are valid has been determined from every single spectrum, we should start from scratch using the mean spectrum (only valid ones) ?
% 27nov2018 rewriten using the real MEAN spectrum (only for spectra_is_valid)
  for idx2= 1:Nparam2
    for idx1= 1:Nparam1
        if (~isnan(mean_spectra_only_valids(1,idx1,idx2))) % mean spectrum from non-data gives NaN
          disp(['CF-Step 7: Doing peak fitting for (', num2str(idx1),',',num2str(idx2),')']);
          for p=1:CFNlines
           if (piUse(p)==1 && pila(p)>= lambdas(1) && pila(p)<=lambdas(end)  ) % marked to be used, AND in composition, P==1 (Halways fitted)
              if (piquees(p)==Halfa)
                 widthCriteriaForDiscard=15; % for discarding, Halfa peak stark broadening is much larger
                 HalfaLargerWidthCorrection=7; % factor to use a larger wavelength range for (much wider) Halfa line
                 fit_excursion=20; % FWHM times, how far the peak is considered to extend, max, tested with Halfa 28feb2020, UPDATE: 9ene2021 for PD.541 experiments was too large, changed from 30 to 20
              else
                 widthCriteriaForDiscard=4;
                 HalfaLargerWidthCorrection=1; % for "normal" peaks
                 fit_excursion=8;
              end     
              % FIT-Step1 first calc an estimated center, and the extend of
              % the peak. Central pixel is searched up to 2 pixels left and
              % right to compensate for miscalibration. piobs() could be used to compensate larger shifts 
              pixcentral=lambda2pixel(lambdas,pila(p));
              [maxl,maxi] = max(mean_spectra_only_valids(pixcentral-2:pixcentral+2,idx1,idx2));
              pixcentral = pixcentral + maxi - 3; % the pixel around central with max value
              heightcentral = mean_spectra_only_valids(pixcentral,idx1,idx2); 
              
              pixleft=lambda2pixel(lambdas,lambdas(pixcentral) - 0.5*AvantesFWHM(lampara_FWHM,lambdas(pixcentral))*HalfaLargerWidthCorrection); % starting search point is FWHM/2
              max_excursion = pixleft - abs(lambda2pixel(lambdas,lambdas(pixcentral)+fit_excursion*AvantesFWHM(lampara_FWHM,lambdas(pixcentral)))-lambda2pixel(lambdas,lambdas(pixcentral))); % number of pixels max for peak extend
              while(pixleft>1 && pixleft>max_excursion && mean_spectra_only_valids(pixleft,idx1,idx2)>heightcentral/40) % stop if we are too far or intensity is <2.5%
                if(mean_spectra_only_valids(pixleft,idx1,idx2)<mean_spectra_only_valids(pixleft-1,idx1,idx2)) %once it starts to increase again...
                    break;
                end
                pixleft=pixleft-1;
              end
              pixright=lambda2pixel(lambdas,lambdas(pixcentral)+ 0.5*AvantesFWHM(lampara_FWHM,lambdas(pixcentral))*HalfaLargerWidthCorrection); % starting search point is FWHM/2
              max_excursion = pixright + abs(lambda2pixel(lambdas,lambdas(pixcentral)+fit_excursion*AvantesFWHM(lampara_FWHM,lambdas(pixcentral)))-lambda2pixel(lambdas,lambdas(pixcentral))); % number of pixels max for peak extend
              while(pixright<Npixeles-1 && pixright<max_excursion && mean_spectra_only_valids(pixright,idx1,idx2)>heightcentral/40) % stop if we are too far or intensity is <5%
                if(mean_spectra_only_valids(pixright,idx1,idx2)<mean_spectra_only_valids(pixright+1,idx1,idx2)) %once it starts to increase again...
                    break;
                end
                pixright=pixright+1;
              end
              % FIT-Step2 we have a [pixleft,pixright] window for peak FITTING, now
              % we search a background zone , stop when intensity builds up
              % it is repited twice, after first iteration, peak extend is
              % reevaluated if a point of the peak is below the baseline
              for bf=1:2
               pixleftbkg1=pixleft;
               pixleftbkg2=pixleft;
               max_excursion = pixleftbkg1 - abs(lambda2pixel(lambdas,lambdas(pixcentral)+fit_excursion*AvantesFWHM(lampara_FWHM,lambdas(pixcentral)))-lambda2pixel(lambdas,lambdas(pixcentral))); % max pixel to the left of background 
               while(pixleftbkg2>1 && pixleftbkg2>max_excursion && (mean_spectra_only_valids(pixleftbkg2,idx1,idx2) - mean_spectra_only_valids(pixleftbkg1,idx1,idx2)) < heightcentral/10) % stop when intensity builds up 
                pixleftbkg2=pixleftbkg2-1;
               end
               pixrightbkg1=pixright;
               pixrightbkg2=pixright;
               max_excursion = pixrightbkg1 + abs(lambda2pixel(lambdas,lambdas(pixcentral)+fit_excursion*AvantesFWHM(lampara_FWHM,lambdas(pixcentral)))-lambda2pixel(lambdas,lambdas(pixcentral))); % max pixel to the right of background 
               while(pixrightbkg2<Npixeles-1 && pixrightbkg2<max_excursion && (mean_spectra_only_valids(pixrightbkg2,idx1,idx2) - mean_spectra_only_valids(pixrightbkg1,idx1,idx2)) < heightcentral/10) % stop when intensity builds up 
                 pixrightbkg2=pixrightbkg2+1;
               end
               % v22: if RETAIN_FIRST_FIT=1, the first calc for pixleftbkg1, pixrightbkg1 is used
               if RETAIN_FIRST_FIT==1 && (idx1~=1 || idx2~=1)
                   pixrightbkg1 = peaks_fit_pixrightbkg1(p,1,1);
                   pixleftbkg1 = peaks_fit_pixleftbkg1(p,1,1);
               end
               % FIT-Step3 linear of quadratic fit of the baseline
               % quadratic with background points (over pixels space)
               %pBkg = fit([lambdas(pixleftbkg2:pixleftbkg1); lambdas(pixrightbkg1:pixrightbkg2)], [mean_spectra_only_valids(pixleftbkg2:pixleftbkg1,idx1,idx2); mean_spectra_only_valids(pixrightbkg1:pixrightbkg2,idx1,idx2)],'poly2' );
               % linear with two points
               %pBkg = fit([lambdas(pixleftbkg2:pixleftbkg1); lambdas(pixrightbkg1:pixrightbkg2)], [mean_spectra_only_valids(pixleftbkg2:pixleftbkg1,idx1,idx2); mean_spectra_only_valids(pixrightbkg1:pixrightbkg2,idx1,idx2)],'poly2' );
               pBkgPx = fit([pixleftbkg1; pixrightbkg1], [mean_spectra_only_valids(pixleftbkg1,idx1,idx2); mean_spectra_only_valids(pixrightbkg1,idx1,idx2)],'poly1' );
               % the baseline over lambda is not really needed
               %pBkgLa = fit([lambdas(pixleftbkg1); lambdas(pixrightbkg1)], [mean_spectra_only_valids(pixleftbkg1,idx1,idx2); mean_spectra_only_valids(pixrightbkg1,idx1,idx2)],'poly1' ); % baseline in lambda space, needed to interpolate at the central wavelenth
               % FIT-Step4 corrects peak extend IF peak wings points are UNDER the baseline
               %DEBUG:  plot(pixleft:pixright,mean_spectra_only_valids(pixleft:pixright,idx1,idx2),pixleft:pixright,feval(pBkg,pixleft:pixright)); title([strjoin(speciesStr(piquees(p))),'@' num2str(pila(p))]);
               baselineOk=true;
               for i=pixcentral:-1:pixleft
                  if(mean_spectra_only_valids(i,idx1,idx2) < feval(pBkgPx,i)*0.99) % if the wing of the peak is below the interpolated baseline
                      pixleft=i;
                      baselineOk=false;
                      break; % look no further
                  end
               end
               for i=pixcentral:1:pixright
                  if(mean_spectra_only_valids(i,idx1,idx2) < feval(pBkgPx,i)*0.99) % if the wing of the peak is below the interpolated baseline
                      pixright=i;
                      baselineOk=false;
                      break; % look no further
                  end
               end
               if baselineOk
                     break; % no need to reevaluate
               end
              end % for 1:2 reevaluate if peak extend is reduced
              if RETAIN_FIRST_FIT==1 && (idx1~=1 || idx2~=1)
                   pixright = peaks_fit_pixright(p,1,1);
                   pixleft = peaks_fit_pixleft(p,1,1);
               end
              %v22 we keep pixels used for fitting, we could check for changes that could explain variability of some parameters
              peaks_fit_pixleftbkg1(p,idx1,idx2)=pixleftbkg1;
              peaks_fit_pixrightbkg1(p,idx1,idx2)=pixrightbkg1;
              peaks_fit_pixleft(p,idx1,idx2)=pixleft;
              peaks_fit_pixright(p,idx1,idx2)=pixright;
              % lorentzian fitting
              peak_lambdas = lambdas(pixleft:pixright);
              if (size(peak_lambdas,1)>4) % 28feb2020 if the area has no peak and just a smooth valley, pixleft=pixright=pixcentral
               peak_intensities =  mean_spectra_only_valids(pixleft:pixright,idx1,idx2); % the peak
               peak_intwobaseline = mean_spectra_only_valids(pixleft:pixright,idx1,idx2) - feval(pBkgPx,pixleft:pixright); % the peak minus interpolated baseline
               ws = warning('off','all');  % Turn off warning messages of lorentzfit
               [yprime1 loparams resnorm1 residual1]=lorentzfit(double(peak_lambdas),double(peak_intwobaseline),[],[],[],opts);warning(ws)  % Turn it back on.
               % store fitting error (relative, %)
               peaks_fit_error_frommean(p,idx1,idx2) = 100*norm(residual1) / sqrt(size(peak_intensities,1)) / max(peak_intensities); % in percentage, same metric as peakfit() from Tom O'Haver
               if (peaks_fit_error_frommean(p,idx1,idx2)>20) % the fitting error is too large
                   disp(['      Warning!! fitting error of line #' num2str(p) ' of ', strjoin(speciesStr(piquees(p))), '@',num2str(pila(p)), ' (',num2str(idx1),',',num2str(idx2),') is ' , num2str(peaks_fit_error_frommean(p,idx1,idx2)), '%']);
               end
               % store height
               peaks_heights_frommean(p,idx1,idx2)=loparams(1)/loparams(3)+loparams(4); % total height
               % store central lambda
               peaks_central_frommean(p,idx1,idx2)=loparams(2);
               error_lambda_central = peaks_central_frommean(p,idx1,idx2)-pila(p);
               if (abs(error_lambda_central) > 2.0*AvantesFWHM(lampara_FWHM,pila(p))) % the error in lambda_central is too large               
                   disp(['      Warning!! central wavelength error for fitting of line #' num2str(p) ' of ', strjoin(speciesStr(piquees(p))), '@',num2str(pila(p)), ' (',num2str(idx1),',',num2str(idx2),') is ' , num2str(error_lambda_central), 'nm']);
               end
               % store total fwhm width and calculate Stark width
               % 14dic2018 there are two possible methods: cuadratic and linear,according to LIBS-tutorial-partII only cuadratic is valid
               %peaks_stark_frommean(selectedPeaks,idx1,idx2) = -(AvantesFWHM(lampara_FWHM,pila(selectedPeaks)).^2 - peaks_fwhm_frommean(selectedPeaks,idx1,idx2).^2)./peaks_fwhm_frommean(selectedPeaks,idx1,idx2);
               peaks_fwhm_frommean(p,idx1,idx2)=2*sqrt(loparams(3)); 
              % 14dic2018 we should test different approaches: linear subs, cuadratic, torres2007 formula...
              % instrumental deconvolution should be done here
              peaks_stark_frommean(p,idx1,idx2) = -(AvantesFWHM(lampara_FWHM,pila(p)).^2 - peaks_fwhm_frommean(p,idx1,idx2).^2)./peaks_fwhm_frommean(p,idx1,idx2);
              %peaks_stark_frommean(p,idx1,idx2) = sqrt(peaks_fwhm_frommean(p,idx1,idx2).^2 - AvantesFWHM(lampara_FWHM,pila(p)).^2 ); 
              if ( peaks_stark_frommean(p,idx1,idx2) < 0 || ~isreal(peaks_stark_frommean(p,idx1,idx2)) )                    
                  peaks_stark_frommean(p,idx1,idx2) = 0; % we can get negative values, remove them
                  disp(['      Warning!! stark width is less than zero for fitting of line #' num2str(p) ' of ', strjoin(speciesStr(piquees(p))), '@',num2str(pila(p)), ' (',num2str(idx1),',',num2str(idx2),')']);
              end
              % calculate and store AREA of the peak
              % new v21 5ago2020 yo can choose wich method to use with the variable AREA_METHOD =
              oversampled_lambdas = lambdas(pixleft):0.001:lambdas(pixright);
              oversampled_lorentzPeak = loparams(1)./((oversampled_lambdas - loparams(2)).^2 + loparams(3)) + loparams(4); % oversampled lorentzian shape
              if AREA_METHOD==0 % analytical formula of modeled lorentzian peak
                peaks_area_frommean(p,idx1,idx2)= pi*peaks_heights_frommean(p,idx1,idx2)*peaks_fwhm_frommean(p,idx1,idx2)/2; % area of the lorentzian function, according to Demtröder, ec. 3.10d, p. 64
              elseif AREA_METHOD==1 % trapezoidal integration of the oversampled modeled lorentzian peak
                peaks_area_frommean(p,idx1,idx2)= trapz(oversampled_lambdas,oversampled_lorentzPeak);
              elseif AREA_METHOD==2 % two halves area from modelled lorentzian peak  
                %the area of the right half could be a better estimation for heavily asymmetric peaks (some channels of our spectrometer have this problem)
                %5ago2020: last thought is that, as coma aberration does not change the total area (energy conservation), either sides of the modeled peak are wrong
                idx_central = find(oversampled_lorentzPeak== max(oversampled_lorentzPeak));
                half_osl = oversampled_lambdas(idx_central:end);
                half_osi = oversampled_lorentzPeak(idx_central:end);
                peaks_twohalfareas_frommean(p,idx1,idx2)= 2*trapz(half_osl,half_osi);
                peaks_area_frommean(p,idx1,idx2)= peaks_twohalfareas_frommean(p,idx1,idx2);
              elseif AREA_METHOD==3 % trapezoidal integration of the NOT oversampled baseline-corrected SPECTRA
                peaks_area_frommean(p,idx1,idx2)= trapz(peak_lambdas,peak_intwobaseline ); 
              else
                error('AREA_METHOD not yet implemented');
              end
             
              %27nov2018 sometimes the fitting produces complex numbers, when there is not clear peaks
              if (~isreal(peaks_fwhm_frommean(p,idx1,idx2)))
                   %figure(33);    plot(trozo_lambdas,trozo_pico,trozo_lambdas,yprime1)
                   %disp('imaginario!');
                   peaks_fwhm_frommean(p,idx1,idx2) = 0;
              end
              else % no enough points for fitting
                  peaks_heights_frommean(p,idx1,idx2)=0; % will be rendered invalid
              end
              % criteria for discard: height < minimum height, saturated,width < instrumental fwhm/1.5 (to include a margin of error, witdh> instrumental fwhm*WidthCriteria (to include Halfa peak width)
              if(peaks_heights_frommean(p,idx1,idx2) < ValMinHeightRelative || peaks_heights_frommean(p,idx1,idx2) > ValSaturationCorrected(lambda2pixel(lambdas,peaks_central_frommean(p,idx1,idx2))) || peaks_fwhm_frommean(p,idx1,idx2) < AvantesFWHM(lampara_FWHM,pila(p))/1.5 || peaks_fwhm_frommean(p,idx1,idx2) > AvantesFWHM(lampara_FWHM,peaks_central_frommean(p,idx1,idx2))*widthCriteriaForDiscard || ~isreal(peaks_fwhm_frommean(p,idx1,idx2))) %discarded
                   disp(['      Discarded peak #', num2str(p) , ' of ', strjoin(speciesStr(piquees(p))), '@',num2str(pila(p)), ' (', num2str(idx1),',',num2str(idx2),') ', '> height=', sprintf('%5.0f',peaks_heights_frommean(p,idx1,idx2)), '   width=', sprintf('%1.2f',peaks_fwhm_frommean(p,idx1,idx2))  ]);
                   peaks_isvalid_frommean(p,idx1,idx2)=0; % tag as not valid
               else   
                   peaks_isvalid_frommean(p,idx1,idx2)=1; % valid
               end
                             
               if(DEBUG_PEAKS==idx1) % && peaks_isvalid_frommean(p,idx1,idx2)==1) 
               %if p==44    
                    figure(33);clf;hold on;
                    plot(lambdas(pixleftbkg2:pixrightbkg2),mean_spectra_only_valids(pixleftbkg2:pixrightbkg2,idx1,idx2),'b'); title([strjoin(speciesStr(piquees(p))),'@' num2str(pila(p))]); % the original peak
                    plot(lambdas(pixleft),mean_spectra_only_valids(pixleft,idx1,idx2),'r*'); % extend of peak for fitting
                    plot(lambdas(pixright),mean_spectra_only_valids(pixright,idx1,idx2),'r*');
                    plot(lambdas(pixleftbkg2),mean_spectra_only_valids(pixleftbkg2,idx1,idx2),'g+'); % baseline extend
                    plot(lambdas(pixrightbkg2),mean_spectra_only_valids(pixrightbkg2,idx1,idx2),'g+');
                    plot(lambdas(pixleftbkg2:pixrightbkg2),feval(pBkgPx,pixleftbkg2:pixrightbkg2),'r--'); % interpolated baseline
                    if(size(peak_lambdas,1)>4) % fitting is not performed for less than 5 points, next line would give am error
                      plot(peak_lambdas,yprime1,'g');  % the fitting
                      plot(oversampled_lambdas,oversampled_lorentzPeak,'k'); % the oversampled lorentzian peak
                      plot(peak_lambdas,peak_intwobaseline,'m.');  % the peaks minus baseline
                      disp(['      ',strjoin(speciesStr(piquees(p))), sprintf('@%3.2f',pila(p)),'> height=', sprintf('%5.0f',peaks_heights_frommean(p,idx1,idx2)), '   measured width=', sprintf('%1.4f',peaks_fwhm_frommean(p,idx1,idx2)),'   stark width=', sprintf('%1.4f',peaks_stark_frommean(p,idx1,idx2)),'   area=',sprintf('%5.1f',peaks_area_frommean(p,idx1,idx2)),'   area(2xhalf)=',sprintf('%5.1f',peaks_twohalfareas_frommean(p,idx1,idx2)),'   error=',sprintf('%2.2f',peaks_fit_error_frommean(p,idx1,idx2)),'%'  ]); 
                      hold off;
                      figure(34);clf;hold on;
                      peakfit([peak_lambdas peak_intwobaseline],lambdas(pixcentral),2.0,1,2,0,99,0,0,0,1); % peak fitting by the Tom O'haver
                    end
                    pause; 
               end
            end % if piusar (we want to use that peak)
          end % p peaks

       
          disp(['CF-Step 8: Calculating Te from Boltzmann-plot for point (', num2str(idx1),',',num2str(idx2),')']);
          % Calculate Te using Boltzmann-plot fron the MEAN spectrum
          % if self-absorption correction, Te is iterated as correction
          % depends on Te
          if(DEBUG_BP==idx1) % clear all figures
            for s=composition
                figure(30+s);clf;
            end
            figure(50);clf;
          end
          % decides the line intensity to use from calcs: height, area ...
          if (LINE_INTENSITY==0) % height
             peaks_intensity_frommean = peaks_heights_frommean; 
          elseif (LINE_INTENSITY==1) % area
             peaks_intensity_frommean = peaks_area_frommean;
          elseif (LINE_INTENSITY==3) % twice the half right area, from v21, this is done with AREA_METHOD=2 
              error('LINE_INTENSITY value deprecated, use AREA_METHOD instead'); % %peaks_intensity_frommean = peaks_twohalfareas_frommean;
          else
             error('LINE_INTENSITY value not implemented.');
          end
          new_B_plot=1;BP_iteration = 1;
          peaks_selfabscoeff_frommean(:,idx1,idx2) = 1; % by default, coefficient should be 1 if the specie has no reference line
          lastTe = 0;
          while(1) % Te iteration with self-absorption correction
            cellcount=1; x_allspecies={}; y_allspecies={}; s_allspecies={};
            for s=composition
             selectedPeaks_specie= find(piquees==s);%<<<<<<< specie to be used
             selectedPeaks_douse = find(piUse==1);
             selectedPeaks_valid = find(peaks_isvalid_frommean(:,idx1,idx2)==1);
             selectedPeaks = intersect(intersect(selectedPeaks_specie,selectedPeaks_douse),selectedPeaks_valid);
             %pi_lngAI = log(pila(selectedPeaks).*peaks__frommean(selectedPeaks,idx1,idx2)./piGA(selectedPeaks)); % for boltzmann plot
             pi_lngAI = log(peaks_intensity_frommean(selectedPeaks,idx1,idx2) ./ peaks_selfabscoeff_frommean(selectedPeaks,idx1,idx2) ./ piGA(selectedPeaks)); % for boltzmann plot, area instead of height, lambda removed
             % WITH LAMBDA 
             %pi_lngAI = log(pila(selectedPeaks).*peaks_intensity_frommean(selectedPeaks,idx1,idx2) ./ peaks_selfabscoeff_frommean(selectedPeaks,idx1,idx2) ./ piGA(selectedPeaks)); % for boltzmann plot, area instead of height, lambda removed
             x = piEm(selectedPeaks); % valores de Ek para el eje horizontal
             if (size(selectedPeaks,1) > 1 && (max(x)-min(x) > 0.1 )) % only if there is at least two peaks AND an enough range of Ek 
              NumPeaks_BP_frommean(s,idx1,idx2)=size(selectedPeaks,1);
              y = pi_lngAI; 
              p = polyfit(x,y,1); % linear regression 
              yfit = polyval(p,x); % to get R^2
              yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
              r2 = 1 - SSresid/SStotal;
              slope_BP_frommean(s,idx1,idx2) = p(1);
              intercept_BP_frommean(s,idx1,idx2) = p(2);
              Te_BP_frommean(s,idx1,idx2)= - 1 / (8.617E-5 * slope_BP_frommean(s,idx1,idx2));
              disp([ '              ' , strjoin(speciesStr(s)), ' Te=' sprintf('%5.0f',Te_BP_frommean(s,idx1,idx2)),'K']);
              % store x,y in a cell for "single-BP" Te estimation ONLY if Te for this specie is reasonable: 2000K-30000K
              % AND it is in the whichSpeciesForSingleBPTe list
              if (ismember(s,whichSpeciesForSingleBPTe) && Te_BP_frommean(s,idx1,idx2) > 2000 && Te_BP_frommean(s,idx1,idx2)<40000)
                s_allspecies{cellcount} = s; % store which species is each x,y in order to store the intercept later on
                x_allspecies{cellcount} = x;
                y_allspecies{cellcount} = y;cellcount=cellcount+1;
              end
              if(DEBUG_BP==idx1)  % optional plot for debugging
                figure(30+s);which_subplot=1; %BP_iteration; % clf removed
                if which_subplot>9 which_subplot=1; end
                subplot(1,1,which_subplot);clf;hold on;plot(x,y,'*');plot(x,yfit,'--g');title(['BP ',' #',num2str(BP_iteration),' ',strjoin(speciesStr(s)),' Te=' sprintf('%5.0f',Te_BP_frommean(s,idx1,idx2)) 'K R^2=' sprintf('%1.2f',r2)]);hold off;
              end
             end %size(selpicos)
            end %species Bolztmann-plot from the mean spectra
            % now a common BP for all the species, so we get a single Te
            % value, stored in single_slope_BP_frommean(idx1,idx2)
            mdl_cell = {mdl1,mdl2,mdl3,mdl4,mdl5,mdl6,mdl7,mdl8,mdl9,mdl10};
            beta0 = [-1, 0.5 , 0.5, 0.5 , 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
            % function nlinmultifit from https://es.mathworks.com/matlabcentral/fileexchange/40613-multiple-curve-fitting-with-common-parameters-using-nlinfit
            try
              [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(x_allspecies, y_allspecies, mdl_cell(1:length(x_allspecies)), beta0(1:1+length(x_allspecies))); % if error, probably the number of species is more than 10
            catch
              beta(:) = -1E20 ; % 0K slope
              disp([ '               > Warning! no fit for single-Te BP']);
            end
            slope_single_BP_frommean(idx1,idx2) = beta(1);
            for i=1:length(x_allspecies)
              intercept_single_BP_frommean(s_allspecies{i},idx1,idx2) = beta(i+1); % this is the intercept of each species BUT with common slope (Te) fitting
            end
            Te_single_BP_frommean(idx1,idx2)= - 1 / (8.617E-5 * slope_single_BP_frommean(idx1,idx2));
            disp([ '              mean      Te=' sprintf('%5.0f',nanmeanLocal(nonzeros(Te_BP_frommean(whichSpeciesForSingleBPTe,idx1,idx2)))),'K']); % leaving zeros out
            disp([ '              Single-BP Te=' sprintf('%5.0f',Te_single_BP_frommean(idx1,idx2)),'K']);
            if(DEBUG_BP==idx1)  % optional plot for debugging
                figure(50);
                r2=zeros(length(x_allspecies),1);
                for i=1:length(x_allspecies)
                 p = [beta(1) beta(i+1)];
                 x=x_allspecies{i};
                 y=y_allspecies{i};
                 yfit = polyval(p,x); % to get R^2
                 yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
                 r2(i) = 1 - SSresid/SStotal;
                 if(r2(i)<0)
                      disp(['            Warning!! ' strjoin(speciesStr(s_allspecies{i})) ' in single-Te BP (mean spectra) has a negative r^2']); % meaning the fitting is soooo bad
                      %stop % the fitting is very bad, something is wrong,
                      %but ... life must go on
                 end
                 hold on; plot(x,y,'*');plot(x,yfit,'--g');
                end
                title(['BP ',' #',num2str(BP_iteration),' all species Te=' sprintf('%5.0f',Te_single_BP_frommean(idx1,idx2)) 'K, mean R^2=' sprintf('%1.2f',mean(r2))]);hold off;
            end
            if(DEBUG_BP==idx1)
              pause;
            end
            % now there is a Te estimation, if self-absorption correction
            % is enabled, it should be iterated
            if (CORRECT_SELFABSORPTION>0)
              if (CORRECT_SELFABSORPTION==2)             % self-absorption correction from sun2009 IRSAC
                disp(['           Selfabsorption correction using Sun2009 IRSAC for iteration #',num2str(BP_iteration)]); 
                for s=composition
                  selectedPeaks_valid=find(peaks_isvalid_frommean(:,idx1,idx2)==1 & piquees(:)==s); % only valid peaks
                  if (BP_iteration==1 && IRSAC_reference(s)>0 ) % the first time, we need to calculate the "fitted value of the reference line", which is the intensity from the BP at the Em value
                      %some BP has no points, left intensity as is,
                      %otherwise SAcoeff goes to zero
                      if (slope_BP_frommean(s,idx1,idx2)==0 && intercept_BP_frommean(s,idx1,idx2)==0) % there is no BP
                          peaks_intensity_frommean_fitted_for_reference_line(s,idx1,idx2) =  peaks_intensity_frommean(IRSAC_reference(s),idx1,idx2);
                      else   
                          peaks_intensity_frommean_fitted_for_reference_line(s,idx1,idx2) = exp(slope_BP_frommean(s,idx1,idx2) * piEm(IRSAC_reference(s)) + intercept_BP_frommean(s,idx1,idx2))*piGA(IRSAC_reference(s));
                      end
                  end
                  % a single value of Te is used from the simultaneous
                  % fitting of all species in the BP
                  if (IRSAC_reference(s)==0 || Te_BP_frommean(s,idx1,idx2)==0 ) % there is no reliable reference line, for example, Al in nordic gold; OR there is no BP for this specie
                    peaks_selfabscoeff_frommean(selectedPeaks_valid,idx1,idx2) = 1;
                  else
                    peaks_selfabscoeff_frommean(selectedPeaks_valid,idx1,idx2) =  peaks_intensity_frommean(selectedPeaks_valid,idx1,idx2) ./ peaks_intensity_frommean_fitted_for_reference_line(s,idx1,idx2) .* piGA(IRSAC_reference(s)) ./ piGA(selectedPeaks_valid) .* exp(-(piEm(IRSAC_reference(s))-piEm(selectedPeaks_valid)) ./ (kB_eVK .* Te_BP_frommean(s,idx1,idx2) ));
                  end
                  if(1)
                    for n=1:size(selectedPeaks_valid,1)
                         if (peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2) > 1)
                           disp([ '               ', strjoin(speciesStr(piquees(selectedPeaks_valid(n)))) , '@' , num2str(pila(selectedPeaks_valid(n))) , ' SAcoeff=', num2str(peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)) ' >> 1.0']);
                           peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)=1.0;
                         else
                           disp([ '               ', strjoin(speciesStr(piquees(selectedPeaks_valid(n)))) , '@' , num2str(pila(selectedPeaks_valid(n))) , ' SAcoeff=', num2str(peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)) ]);
                         end
                    end
                  end
                end
              end
              %Te = Te_single_BP_frommean(idx1,idx2);
              Te = nanmeanLocal(nonzeros(Te_BP_frommean(whichSpeciesForSingleBPTe,idx1,idx2))); % Te from mean of each Te, so rsd<5% can be performed
              if (abs(Te-lastTe)/Te < 0.1 ) % 10% of difference
                  % ok, Te has converged, but, according to sun2009, the
                  % RSD of temperatures should be within 5%, otherwise, the
                  % highest Te should be substitud by mean Te and iterate
                  % again
                actual_temps = Te_BP_frommean(whichSpeciesForSingleBPTe,idx1,idx2);
                TeRsd = std(actual_temps) / nanmeanLocal(actual_temps);
                if (TeRsd > 0.05) % more than 5% (sun2009)
                  [Temax,cualTemax] = max(nonzeros(Te_BP_frommean(whichSpeciesForSingleBPTe,idx1,idx2)));
                  disp([ '               rsd of Te>5% substituing Te=', num2str(Te_BP_frommean(whichSpeciesForSingleBPTe(cualTemax),idx1,idx2)),'K of ', strjoin(speciesStr(whichSpeciesForSingleBPTe(cualTemax))), ' with mean value Te=',num2str(nanmeanLocal(actual_temps)),'K']);
                  Te_BP_frommean(whichSpeciesForSingleBPTe(cualTemax),idx1,idx2) = nanmeanLocal(actual_temps);
                  % the next step will be calc BPs again, so we perform the
                  % SA correction now, the code is duplicated but ...
                  % I miss the GOTO instruction...
                  % CODE DUPLICATED FROM ABOVE, should be refactored 
                  if (CORRECT_SELFABSORPTION==2)             % self-absorption correction from sun2009 IRSAC
                   disp(['           Selfabsorption correction using Sun2009 IRSAC for iteration #',num2str(BP_iteration)]); 
                   for s=composition
                   selectedPeaks_valid=find(peaks_isvalid_frommean(:,idx1,idx2)==1 & piquees(:)==s); % only valid peaks
                   if (BP_iteration==1 && IRSAC_reference(s)>0 ) % the first time, we need to calculate the "fitted value of the reference line", which is the intensity from the BP at the Em value
                      peaks_intensity_frommean_fitted_for_reference_line(s,idx1,idx2) = exp(slope_BP_frommean(s,idx1,idx2) * piEm(IRSAC_reference(s)) + intercept_BP_frommean(s,idx1,idx2))*piGA(IRSAC_reference(s));
                   end
                   % a single value of Te is used from the simultaneous
                   % fitting of all species in the BP
                   if (IRSAC_reference(s)==0 ||  Te_BP_frommean(s,idx1,idx2)==0) % there is no reliable reference line, for example, Al in nordic gold; or there is no BP for this specie
                     peaks_selfabscoeff_frommean(selectedPeaks_valid,idx1,idx2) = 1;
                   else
                     peaks_selfabscoeff_frommean(selectedPeaks_valid,idx1,idx2) =  peaks_intensity_frommean(selectedPeaks_valid,idx1,idx2) ./ peaks_intensity_frommean_fitted_for_reference_line(s,idx1,idx2) .* piGA(IRSAC_reference(s)) ./ piGA(selectedPeaks_valid) .* exp(-(piEm(IRSAC_reference(s))-piEm(selectedPeaks_valid)) ./ (kB_eVK .* Te_BP_frommean(s,idx1,idx2) ));
                   end
                   if(1)
                    for n=1:size(selectedPeaks_valid,1)
                         if (peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2) > 1)
                           disp([ '               ', strjoin(speciesStr(piquees(selectedPeaks_valid(n)))) , '@' , num2str(pila(selectedPeaks_valid(n))) , ' SAcoeff=', num2str(peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)) ' >> 1.0']);
                           peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)=1.0;
                         else
                           disp([ '               ', strjoin(speciesStr(piquees(selectedPeaks_valid(n)))) , '@' , num2str(pila(selectedPeaks_valid(n))) , ' SAcoeff=', num2str(peaks_selfabscoeff_frommean(selectedPeaks_valid(n),idx1,idx2)) ]);
                         end
                    end
                   end
                  end
                 end
                 % END OF CODE DUPLICATED FROM ABOVE, should be refactored 
                else
                 break; % we are done
                end
              end
              lastTe = Te;
            else
               break; % No SA correction, no iterations
            end  % if(CORRECT_...
            BP_iteration = BP_iteration +1;  
          end % while(1)
          % at this point, we have:
          % -Te_BP_frommean: Te for each specie
          % -Te_single_BP_frommean: single Te for all species
          % -intercept_BP_frommean: intercept of the BP of each specie
          % -intercept_single_BP_frommean: intercept of EACH specie from the
          %      single (common) BP
          
          
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % 1setp2019 removed, I am not sure it works
          % Ne from the intensity ratio of preselected lines
          % With Calcicum: peak 33 CaII 470nm and all available CaI
          % With Copper:  peak 71 CuII 254.48nm and all available CuI
          %xxxxxx  DUPLICATE IN STEP 5 ONCE IT WORKS!!!
          %ojo!! sale 1E14 pero estable, ¿?¿?¿ tendrá algún error?
          %whichPeakIon_forNeLineRatio = find(pila==254.48); % CuII , should be one peak of a II ion, CaII 470nm or CuII  or ...
          %whichSpecieNeutral_forNeLineRatio = CuI; % neutral species          
          %disp(['           Calculating Ne from line ratios (not corrected!!) of ' strjoin(speciesStr(whichSpecieNeutral_forNeLineRatio))]);
          %selectedPeaks = intersect(intersect(find(piquees==whichSpecieNeutral_forNeLineRatio),find(piUse==1)),find(peaks_isvalid_frommean(:,idx1,idx2)==1));
          %Te =  Te_BP_frommean(whichSpecieNeutral_forNeLineRatio,idx1,idx2);
          %peaks_Ne_fromLineRatio_frommean(selectedPeaks,idx1,idx2) = (peaks_heights_frommean(selectedPeaks,idx1,idx2)./peaks_heights_frommean(whichPeakIon_forNeLineRatio,idx1,idx2)).*((pila(selectedPeaks)./pila(whichPeakIon_forNeLineRatio)).*(piGA(whichPeakIon_forNeLineRatio)./piGA(selectedPeaks)))*(2.*3.1415927.*me_MeVc2.*kB_eVK./h_eVs.^2.*Te).^1.5.*exp((-piEm(whichPeakIon_forNeLineRatio)+piEm(selectedPeaks)-Einfinite(s))./kB_eVK./Te)./1E6; % gives m-1
          
          % initial value of Te0 to create the SahaBP and to calculate Ne from STark
          if (whichSpecieForTe0==0) % use the Te from Te_single_BP_frommean()
                  Te =  Te_single_BP_frommean(idx1,idx2); % the one from the single BP
          else % use the Te from the chosen specie
                  Te =  Te_BP_frommean(whichSpecieForTe0,idx1,idx2); % the chosen one 
          end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Ne from stark width using the mean spectra of every peak 
          % we know now a reliable Te
          disp('           Ne from Stark broadening');  
          selectedPeaks=intersect(find(piStark~=0),find(piUse~=0)); % peaks with stark parameters AND piusar=1 (some with large selfabsorption and piusar=0 produced large Ne values
          selectedPeaks=intersect(selectedPeaks,find(peaks_isvalid_frommean(:,1,1)==1));  % only valid peaks should be used to calculate Ne
          peaks_Ne_fromStark_frommean(selectedPeaks,idx1,idx2)=peaks_stark_frommean(selectedPeaks,idx1,idx2)./2.*piStark(selectedPeaks); %linear Ne-Stark broadening relationship
          % different formula for Halfa, % from Pardini 2013 & El Sherbini 2006
          %if(peaks_isvalid_frommean(1,idx1,idx2))
          %1sept2019 this is tricky, the formula below should be used
          %regardless if the peak is valid or not, but fitting could be
          %wrong enough to give an error 
            a = sqrt(4033.8-24.45*(log(Te))^2 ); 
            b = 1.029E8+174576.3*Te;
            alfa = 1 /(a+b/sqrt(NeStarkNe0));
            peaks_Ne_fromStark_frommean(1,idx1,idx2) = 8.02E12 *(peaks_stark_frommean(1,idx1,idx2)*10/alfa)^(3/2); % p=1 is Halfa
            %disp(['             Ne_Halfa=', num2str(peaks_Ne_fromStark_frommean(1,idx1,idx2))]);
          %end
          idx_nz = find(peaks_Ne_fromStark_frommean(:,idx1,idx2)>0); 
          idx_nz = idx_nz(2:end);% Halfa (always idx=1) should not be here
          % we get an alternative Ne value from the MEDIAN of start width
          % of peaks with stark parameters, removing outliers
          temp_Ne =  median(peaks_Ne_fromStark_frommean(idx_nz,idx1,idx2),1);
          idx_nosolarge = find(peaks_Ne_fromStark_frommean(:,idx1,idx2)<temp_Ne*7);
          idx_nososmall = find(peaks_Ne_fromStark_frommean(:,idx1,idx2)>temp_Ne/7);
          idx_nz_noOutliers = intersect(idx_nosolarge,idx_nososmall); % the range is x7 and /7 to be considered not an outlier
          idx_nz_noOutliers = idx_nz_noOutliers(2:end);% Halfa (always idx=1) should not be here
          temp_Ne = median(peaks_Ne_fromStark_frommean(idx_nz_noOutliers,idx1,idx2),1);
          disp(['             Ne_Halfa=' , sprintf('%1.2e',peaks_Ne_fromStark_frommean(1,idx1,idx2)),'     Ne_median_no_outliers=', sprintf('%1.2e',temp_Ne), '(',num2str(size(idx_nz_noOutliers,1)), ' peaks)']);
          %1setp2019 Halfa is very unreliable in the experiment (apice
          %2014-13), maybe with another one...
          peaks_Ne_frommean(idx1,idx2) = (100-HALFA_PERCENTAGE)/100*temp_Ne + (HALFA_PERCENTAGE)/100*peaks_Ne_fromStark_frommean(1,idx1,idx2);
          
          %peaks_Ne_frommean(idx1,idx2) = 1.61E16;disp('OJO!!! valor Ne fijo'); % OJO!! pruebo a dejarlo fijo
          
          if (DEBUG_NE_HISTOGRAM==idx1)
             figure(54);
             subplot(2,1,1);
             histogram(peaks_Ne_fromStark_frommean(idx_nz,idx1,idx2),50);
             subplot(2,1,2);
             histogram(peaks_Ne_fromStark_frommean(idx_nz_noOutliers,idx1,idx2),50);
             pause;
          end

          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % SAHA-bolztman plot of elements. code imported from normalization algorithm
          % 
          % 
          new_SH_plot=1; % allow to insert multiples plots in one
          if (whichSpecieForTe0==0)
              FEATURES = [FEATURES 'Te0 for Saha-BP from singleBP, ' ];
              disp(['           Doing Saha-Boltzmann plot using Te0 from singleBP and Ne from Halfa/median_others_no_outliers']);
          else
              FEATURES = [FEATURES 'Te0 for Saha-BP from ' strjoin(speciesStr(whichSpecieForTe0)) ', '];
              disp(['           Doing Saha-Boltzmann plot using Te0 from BP of ', strjoin(speciesStr(whichSpecieForTe0)), ' and Ne from Halfa/median_others_no_outliers']);
          end

           Ne = peaks_Ne_frommean(idx1,idx2); % and Ne already selected from Starks measurements (Halfa or whatever)
           %or
           %Ne = nanmeanLocal( peaks_Ne_fromLineRatio_frommean(peaks_Ne_fromLineRatio_frommean(:,idx1,idx2)~=0,idx1,idx2),1); % NE from line intensity ratio of CaI and CaII <<<< CHANGE?
           %disp('Using Ne from line ratios');
           if (Ne==0 || isnan(Ne) || Te > 50000 || Te < 2000 || isnan(Te)) % some background-only spectra could end here 
             disp(['Cannot perform SahaBP: Te0='  num2str(Te),' Ne0=' , num2str(Ne), '  (' , num2str(idx1),',',num2str(idx2),')'  ]);
             continue;
           end
           
           SahaBP_iteration=1;
           while(1) % iterative process to converge to a new Te
           % new V16 the for s=composition loop is INSIDE the while(1), as
           % in the BP
           cellcount=1; x_allspecies={}; y_allspecies={}; s_allspecies={}; % to accumlate data for simultaneous fitting of all the SahaBP
           for s = composition(mod(composition,2)==0) % we will use each specie AND the next index, CaI+CaII , MgI+MgII... 'composition' should have neutral, ion, neutral, ion...
              selectedPeaks_specie = find(piquees==s);%<<<<<<< specie to be used
              selectedPeaks_specie_ion = find(piquees==s+1);%the ion II
              selectedPeaks_douse = find(piUse==1);
              selectedPeaks_valid = find(peaks_isvalid_frommean(:,idx1,idx2)==1);
              selectedPeaks = intersect(intersect([selectedPeaks_specie ; selectedPeaks_specie_ion],selectedPeaks_douse),selectedPeaks_valid);
              selectedPeaks_z0 = intersect(intersect(selectedPeaks_specie,selectedPeaks_douse),selectedPeaks_valid);
              selectedPeaks_z1 = intersect(intersect(selectedPeaks_specie_ion,selectedPeaks_douse),selectedPeaks_valid);
              %pi_lngAI_z0 = log(pila(selectedPeaks_z0).* (selectedPeaks_z0,idx1,idx2)./ peaks_selfabscoeff_frommean(selectedPeaks_z0,idx1,idx2) ./piGA(selectedPeaks_z0)); % for boltzmann plot
              %pi_lngAI_z1 = log(pila(selectedPeaks_z1).* peaks_area_frommean(selectedPeaks_z1,idx1,idx2)./ peaks_selfabscoeff_frommean(selectedPeaks_z1,idx1,idx2) ./piGA(selectedPeaks_z1)) - log(2*(2*3.1415927*0.511E6/9E16*8.6173324E-5/(4.135668E-15)^2)^1.5*Te^1.5 / (Ne*1E6)); %eq. 2 from Aguilera2007
              pi_lngAI_z0 = log( peaks_intensity_frommean(selectedPeaks_z0,idx1,idx2)./ peaks_selfabscoeff_frommean(selectedPeaks_z0,idx1,idx2) ./piGA(selectedPeaks_z0)); % for boltzmann plot
              pi_lngAI_z1 = log( peaks_intensity_frommean(selectedPeaks_z1,idx1,idx2)./ peaks_selfabscoeff_frommean(selectedPeaks_z1,idx1,idx2) ./piGA(selectedPeaks_z1)) - log(2*(2*3.1415927*0.511E6/9E16*8.6173324E-5/(4.135668E-15)^2)^1.5*Te^1.5 / (Ne*1E6)); %eq. 2 from Aguilera2007
              pi_Em_corr = piEm(selectedPeaks_z1)+ pi_z(selectedPeaks_z1)*Einfinite(s+1); % corrected Em for ion species (II)
              x = [ piEm(selectedPeaks_z0) ; pi_Em_corr ] ; % Ek values for X axis
              y = [ pi_lngAI_z0 ; pi_lngAI_z1]; % y axis
              if (size(selectedPeaks,1) > 1 && (max(x)-min(x) > 0.1 )) % only if there is at least two peaks AND an enough range of Ek 
                % if not set, Te is calculated again for every specie 
               p = polyfit(x,y,1); % linear regression 
               yfit = polyval(p,x); % to get R^2
               yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
               r2 = 1 - SSresid/SStotal;
               slope_SahaBP_frommean(s,idx1,idx2) = p(1);
               intercept_SahaBP_frommean(s,idx1,idx2) = p(2);
               Te_SahaBP_frommean(s,idx1,idx2)= - 1 / (8.617E-5 * slope_SahaBP_frommean(s,idx1,idx2));              
               % OJO!! check those two lines
               %slope_SahaBP_frommean(s,idx1,idx2)= - 1 / (8.617E-5 * Te_SahaBP_frommean(whichSpecieForTe0,idx1,idx2));  % slope from the first specie, which should be whichSpecieForTe
               %intercept_SahaBP_frommean(s,idx1,idx2) = +(mean(y) - slope_SahaBP_frommean(s,idx1,idx2)*mean(x)); % linear fit with fixed slope https://lists.gnu.org/archive/html/help-octave/2010-05/msg00308.html
               yfit = polyval([slope_SahaBP_frommean(s,idx1,idx2) intercept_SahaBP_frommean(s,idx1,idx2)],x); % this is only to get R^2
               yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
               r2 = 1 - SSresid/SStotal;
               % store x,y in a cell for "single-sahaBP" Te estimation ONLY if Te for this specie is reasonable: 2000K-30000K
               % AND it is in the whichSpeciesForSingleSahaBPTe list
               if (ismember(s,whichSpeciesForSingleSahaBPTe) && Te_SahaBP_frommean(s,idx1,idx2) > 2000 && Te_SahaBP_frommean(s,idx1,idx2)<30000)
                s_allspecies{cellcount} = s; % store which species is each x,y in order to store the intercept later on
                x_allspecies{cellcount} = x;
                y_allspecies{cellcount} = y;cellcount=cellcount+1;
               end
               disp([ '              Saha-Boltzmann of ' ,strjoin(speciesStr(s)), '+' ,strjoin(speciesStr(s+1)),' newTe=', sprintf('%5.0f',Te_SahaBP_frommean(s,idx1,idx2)) , '  Ne=', sprintf('%1.2e',Ne)  ]);                 
               if(DEBUG_SBP==idx1)
                  figure(15);
                  if (new_SH_plot==1)
                     clf;which_subplot=1;
                  end    
                     subplot(2,2,which_subplot);hold on;plot(x,y,'*');plot(x,yfit,'--g');title(['SahaBP ' , strjoin(speciesStr(s)),'+',strjoin(speciesStr(s+1)),' Te=' sprintf('%5.0f',Te_SahaBP_frommean(s,idx1,idx2)) 'K R^2=' sprintf('%1.2f',r2)]);hold off;pause;
                     which_subplot = which_subplot +1;
                     if which_subplot>4 which_subplot=1; end
                     if (new_SH_plot==1) new_SH_plot=0; end
               end
             else
               disp([ '               > Not enough peaks or Ek-Ei range for ' ,strjoin(speciesStr(s)), '+' ,strjoin(speciesStr(s+1))]);
               %break;% there are not enough number of peaks
             end % if there are peaks for the SahaBP
           end % for s=coomposition loop
            
            % all the species has been processed, let's do a singleSahaBP
            mdl_cell = {mdl1,mdl2,mdl3,mdl4,mdl5,mdl6,mdl7,mdl8,mdl9,mdl10};
            beta0 = [-1, 0.5 , 0.5, 0.5 , 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
            % function nlinmultifit from https://es.mathworks.com/matlabcentral/fileexchange/40613-multiple-curve-fitting-with-common-parameters-using-nlinfit
            try
             [beta,r,J,Sigma,mse,errorparam,robustw] = nlinmultifit(x_allspecies, y_allspecies, mdl_cell(1:length(x_allspecies)), beta0(1:1+length(x_allspecies))); % if error, probably the number of species is more than 10
            catch
              beta(:) = -1E20; % 0K slope
              disp([ '               > Warning! no fit for single-Te SahaBP']);
            end
            slope_single_SahaBP_frommean(idx1,idx2) = beta(1);
            for i=1:length(x_allspecies)
              intercept_single_SahaBP_frommean(s_allspecies{i},idx1,idx2) = beta(i+1); % this is the intercept of each species BUT with common slope (Te) fitting
            end
            Te_single_SahaBP_frommean(idx1,idx2)= - 1 / (8.617E-5 * slope_single_SahaBP_frommean(idx1,idx2));
            disp([ '              All species Te=' sprintf('%5.0f',Te_single_SahaBP_frommean(idx1,idx2)),'K']);
            if(DEBUG_SBP==idx1)  % optional plot for debugging
                figure(51);
                r2=zeros(length(x_allspecies),1);
                for i=1:length(x_allspecies)
                 p = [beta(1) beta(i+1)];
                 x=x_allspecies{i};
                 y=y_allspecies{i};
                 yfit = polyval(p,x); % to get R^2
                 yresid = y - yfit; SSresid = sum(yresid.^2); SStotal = (length(y)-1) * var(y);
                 r2(i) = 1 - SSresid/SStotal;
                 if(r2(i)<0)
                      disp(['  Warning!! ' strjoin(speciesStr(s_allspecies{i})) ' in single-Te SahaBP (mean spectra) has a negative r^2']); % meaning the fitting is soooo bad
                      %stop % the fitting is very bad, something is wrong,
                      %but ... life must go on
                 end
                 hold on; plot(x,y,'*');plot(x,yfit,'--g');
                end
                title(['SahaBP ',' #',num2str(SahaBP_iteration),' all species Te=' sprintf('%5.0f',Te_single_SahaBP_frommean(idx1,idx2)) 'K, mean R^2=' sprintf('%1.2f',mean(r2))]);hold off;
            end
            if(DEBUG_SBP==idx1)
              pause;
            end
            %update Te
            lastTe=Te;
            if (whichSpecieForTe0==0) % we iterate the single Te 
                  Te =  Te_single_SahaBP_frommean(idx1,idx2); % the one from the single SahaBP Te
            else % we trust in this specie for the Te iteration
                  Te =  Te_SahaBP_frommean(whichSpecieForTe0,idx1,idx2); % the chosen one 
            end
            % we should update also Ne  ¿?
            %{
                    CaIIpeakToBeUsedForIon = 33; % CaII 470nm
                    selectedPeaks_specie= find(piquees==s);%<<<<<<< specie to be used
                    selectedPeaks_douse = find(piusar==1);
                    selectedPeaks_valid = find(peaks_isvalid_frommean(:,idx1,idx2)==1);
                    selectedPeaks = intersect(intersect(selectedPeaks_specie,selectedPeaks_douse),selectedPeaks_valid);
                    tmp_Ne = (peaks_heights_frommean(selectedPeaks,idx1,idx2)./peaks_heights_frommean(CaIIpeakToBeUsedForIon,idx1,idx2)).*((pila(selectedPeaks)./pila(CaIIpeakToBeUsedForIon)).*(piGA(CaIIpeakToBeUsedForIon)./piGA(selectedPeaks)))*(2.*3.1415927.*me_MeVc2.*kB_eVK./h_eVs.^2.*Te).^1.5.*exp((-piEm(CaIIpeakToBeUsedForIon)+piEm(selectedPeaks)-Einfinite(s))./kB_eVK./Te)./1E6; % gives m-1
                    Ne = nanmeanLocal( tmp_Ne(tmp_Ne~=0),1); % NE from line intensity ratio of CaI and CaII <<<< CHANGE?
            %}
            if (abs(Te-lastTe)/Te < 0.10 ) % if Te changes less than 10%, let finish 
                break;
            end
            SahaBP_iteration = SahaBP_iteration +1;  
            end  %while(1)
        else
          disp(['           No valid mean spectrum for (', num2str(idx1),',',num2str(idx2),')']);
        end %all CF calcs for a spatial points 
        %disp(['After idx1=' num2str(idx1) ' the number of  valid lines is ' num2str(size(find(peaks_isvalid_frommean(:,idx1,idx2)==1),1))]);
    end %idx1
  end %idx2 of the mean spectra loop

% plot with the evolution of valid peaks discarded by the peak fitting
% algorithm
if (1)
 figure(52); % peaks that are valid in each spatial points, could be used again to restrict BP for only those peaks, in all spatial points
 which_peaks_start = find(peaks_isvalid_frommean(:,1,1));
 common_peaks_isvalid_frommean = peaks_isvalid_frommean(:,1,1);
 for idx2= 1:Nparam2
  for idx1= 2:Nparam1
    if (idx1~=43)
      common_peaks_isvalid_frommean = and(common_peaks_isvalid_frommean,peaks_isvalid_frommean(:,idx1,1));
    end
  end
 end
 subplot(1,2,2);
 pivfm_only = peaks_isvalid_frommean(which_peaks_start,:,:);
 imagesc(pivfm_only);
 ylabel(num2str(which_peaks_start));
 axis on;

 subplot(1,2,1);
 imagesc(common_peaks_isvalid_frommean(which_peaks_start));
 %ylabel(num2str(which_peaks_start));
 %title(num2str(which_peaks_start))
 % var 'common_peaks_isvalid_frommean' has the mask of peaks ALWAYS valid
 % in all the Nparam iterations
end


  %FEATURES = [FEATURES 'Concentrations from Te_BP_mean & Ne_fromStark_mean, WITH z=1 substitution,'];
  FEATURES = [FEATURES 'Concentrations from Te_??BP_mean & Ne_fromStark_Halfa_mean, WITH z=1 substitution,'];
  % calculate concentrations
  % Te and intercepts could be from:
  %    Intercept of individual BP of each species, at
  %         intercept_BP_frommean() and Te_BP_frommean()
  %    single Te after selfabsorption correction with single slope fitting,
  %         at intercept_single_BP_frommean() and  Te_single_BP_frommean()
  %    Te for every element from SahaBP, at  
  %         intercept_SahaBP_frommean() and Te_SahaBP_frommean()
  % Ne can be from:
  %    H-alpha line   peaks_Ne_fromStark_frommean(1,idx1,idx2)
  %    a specific line peaks_Ne_fromStark_frommean(i,idx1,idx2)
  partition_func = zeros(max(composition),1);
  CsF = zeros(max(composition),1);
  for idx2= 1:Nparam2
    for idx1= 1:Nparam1
      if(getConcentrationsFrom==SAHA_BP)           Te = Te_SahaBP_frommean(whichSpecieForTe0,idx1,idx2); end % Te from Saha-Boltzmann plot <:<:<:<:<:<:<
      if(getConcentrationsFrom==SINGLE_TE_SAHA_BP) Te = Te_single_SahaBP_frommean(idx1,idx2); end % Te from Saha-Boltzmann plot <:<:<:<:<:<:<
      if(getConcentrationsFrom==EACH_BP)           Te = Te_BP_frommean(whichSpecieForTe0,idx1,idx2);  end % Te from individual BP, with SA correction Te is the same for all species
      if(getConcentrationsFrom==SINGLE_TE_BP)      Te = Te_single_BP_frommean(idx1,idx2); end % Te from commen BP, with SA correction Te is a single one 
      Ne = peaks_Ne_frommean(idx1,idx2); % and Ne from peak 1 of H_alfa or median(other peaks)
      for s=composition
        partition_func(s) = polyval(poly_partfunc(s,:),Te);
        if (getConcentrationsFrom== SAHA_BP)           intercept = intercept_SahaBP_frommean(s,idx1,idx2); end
        if (getConcentrationsFrom== SINGLE_TE_SAHA_BP) intercept = intercept_single_SahaBP_frommean(s,idx1,idx2); end
        if (getConcentrationsFrom== EACH_BP)           intercept = intercept_BP_frommean(s,idx1,idx2); end
        if (getConcentrationsFrom== SINGLE_TE_BP)      intercept = intercept_single_BP_frommean(s,idx1,idx2); end
        if (intercept==0) % actually, there is no intercept
            disp(['  Warning!! No intercept for ' strjoin(speciesStr(s)) ]);
            intercept = -1E20; % should give 0%
        end
        CsF(s) = partition_func(s) * exp (intercept); % intercept from SahaBP (it seems its for neutrals only)
      end
      F=sum(CsF);
      concentration(composition,idx1,idx2) = CsF(composition) ./ F;
      disp(['concentrations for (', num2str(idx1),',',num2str(idx2),')'])
      for s = composition(mod(composition,2)==0) % these are the neutrals z=0
        disp(['  ',strjoin(speciesStr(s)),'+',strjoin(speciesStr(s+1)),'=', num2str(concentration(s,idx1,idx2)+concentration(s+1,idx1,idx2))]);
      end
      
      %% "by hand" substitution of unreliable values using the Saha eq. for ions (z=1)
      for s = composition(mod(composition,2)==0) % these are the neutrals z=0
        CsF(s+1) = CsF(s)*2*partition_func(s+1)/partition_func(s)*(2*3.1415927*me_MeVc2*kB_eVK/h_eVs^2*Te)^1.5*exp(-Einfinite(s)/kB_eVK/Te)/(Ne*1E6); %MeV
      end      
      F=sum(CsF);
      concentration_z1(composition,idx1,idx2) = CsF(composition) ./ F;
      disp(['concentrations (z=1 substitution) for (', num2str(idx1),',',num2str(idx2),')'])
      for s = composition(mod(composition,2)==0) % these are the neutrals z=0
        disp(['  ',strjoin(speciesStr(s)),'+',strjoin(speciesStr(s+1)),'=', num2str(concentration_z1(s,idx1,idx2)+concentration_z1(s+1,idx1,idx2))]);
      end

      %% composition in weight from composition_z1
      for s = composition % neutrals and ions
        concentration_weight(s,idx1,idx2) = concentration_z1(s,idx1,idx2) * atomicWeights(s);
      end      
      FC = sum(concentration_weight(:,idx1,idx2));
      concentration_weight(composition,idx1,idx2) = concentration_weight(composition,idx1,idx2) ./ FC;
      disp(['concentrations in weight (z=1 substitution) for (', num2str(idx1),',',num2str(idx2),')'])
      for s = composition(mod(composition,2)==0) % these are the neutrals z=0
        disp(['  ',strjoin(speciesStr(s)),'+',strjoin(speciesStr(s+1)),'=', num2str(concentration_weight(s,idx1,idx2)+concentration_weight(s+1,idx1,idx2))]);
      end
      
      % error , cuadratic difference between real composition and weight
      % composition with z=1 substitution
      err=0; idx=1;
      for s = composition(mod(composition,2)==0) % these are the neutrals z=0
        err = err + (concentration_weight(s,idx1,idx2)+concentration_weight(s+1,idx1,idx2) - realComposition(idx))^2;
        idx = idx+1;
      end
      concentration_error(idx1,idx2) = sqrt(err); % error in composition (real-calculated)
      disp(['   Error: ',num2str(sqrt(err)*100),'%']);      
      
      %C_Ca = (CsF(CaI) + CsF(CaII)) / F;
      %C_Mg = (CsF(MgI) + CsF(MgII)) / F;
      ratio_MgCa_CFLibs(idx1,idx2) = (concentration(MgI,idx1,idx2)+concentration(MgII,idx1,idx2)) / (concentration(CaI,idx1,idx2)+concentration(CaII,idx1,idx2)) * 1000; % milimols/mol
      ratio_MgCa_CFLibs_z1(idx1,idx2) = (concentration_z1(MgI,idx1,idx2)+concentration_z1(MgII,idx1,idx2)) / (concentration_z1(CaI,idx1,idx2)+concentration_z1(CaII,idx1,idx2)) * 1000; % milimols/mol
      
      %C_Sr = (CsF(SrI) + CsF(SrII)) / F;
      %ratio_SrCa_CFLibs(idx1,idx2) = C_Sr / C_Ca * 1000; % milimols/mol
      %C_Cu = (CsF(CuI) + CsF(CuII)) / F;
      %C_Zn = (CsF(ZnI) + CsF(ZnII)) / F;
      %ratio_CuZn_CFLibs(idx1,idx2) = C_Zn / C_Cu * 1000; % milimols/mol
    end %idx1
   end %idx2

 %statistics about the stability of the results point-by-point
 if(DEBUG_STABILITY)
    % number of valid peaks in each spatial point for every specie peak
    figure(23);surf(1:Nparam1,1:CFNlines,squeeze(sum(peaks_isvalid(:,:,:,1),2)));title('for each emission line, number of valid peaks in each spatial point');
    % number of points used for the SahaBP
    figure(24);surf(1:Nparam1,1:Nspecies,NumPeaks_SahaBP_frommean);title('for each specie, number of points used in SahaBP_frommean');
 end


% new v21 8ene2021 the variable "ratio_to_export"  will have the choosen ratio values (from WHICH_RATIO_TO_EXPORT)
%  'lir_from_mean' = line intensity ratio from mean spectrum; 'mean_lir' = mean of ratios from individual spectra; 'molar_cf' = molar concentration from cf libs; 'molar_cf_z1' with z1 substitution ; 'weight_cf_z1' weight concentration (not molar)

if (strcmp(WHICH_RATIO_TO_EXPORT,'lir_from_mean') || strcmp(WHICH_RATIO_TO_EXPORT , 'mean_lir')) && DO_CF==1
    disp('WARNING!!, ratio_to_export not available for CF=1, changing to molar_cf_z1');
    WHICH_RATIO_TO_EXPORT = 'molar_cf_z1';
end
if (strcmp(WHICH_RATIO_TO_EXPORT , 'molar_cf') || strcmp(WHICH_RATIO_TO_EXPORT , 'molar_cf_z1') || strcmp(WHICH_RATIO_TO_EXPORT , 'weight_cf_z1'))  && DO_CF==0
    disp('WARNING!!, ratio_to_export not available for CF=0, changing to mean_lir');
    WHICH_RATIO_TO_EXPORT = 'mean_lir';
end
if strcmp(WHICH_RATIO_TO_EXPORT , 'lir_from_mean')
    ratio_to_export = ratio_MgCa_del_espectro_promedio;
elseif strcmp(WHICH_RATIO_TO_EXPORT , 'mean_lir')
    ratio_to_export = ratio_MgCa_del_promedio_de_ratios;
elseif strcmp(WHICH_RATIO_TO_EXPORT , 'molar_cf')
    ratio_to_export = ratio_MgCa_CFLibs;
elseif strcmp(WHICH_RATIO_TO_EXPORT , 'molar_cf_z1')
    ratio_to_export = ratio_MgCa_CFLibs_z1;
elseif strcmp(WHICH_RATIO_TO_EXPORT , 'weight_cf_z1')
    ratio_to_export = ratio_MgCa_CFLibs_weight_z1;
else
    error('ratio to export not implemented');
end

%v22 RSD of the ratio sequence to export (ratio_to_export)

rsd_of_seq = nanstdLocal(ratio_to_export)/nanmeanLocal(ratio_to_export);
disp(['STD of ratio_to_export:  ' num2str(rsd_of_seq*100) '%']);

%% esto está copiado de calidadSecuencia el 23/7/2015 para meter en el excel los resultados
% genera una fila excel con datos del experimento y alguna estimación de su
% "calidad"
hay_ssa=0;
if (DO_SSA==1)
    
if (modo2D==1)
    xnan = nanmeanLocal (ratio_to_export ,1);
else
    xnan = ratio_to_export(:,1);
end
x=xnan(~isnan(xnan));
sn_ssa = 0; senal_ssa=0;


if (size(x,1) >  16 )  % es posible que no haya espectros válidos suficientes, 
 disp('Doing FFT & SSA to smooth sequence...');
 fs = 1 / deltaR;
 if (size(x,1) < 32)
  NFFT=16; % para que no se pare con secuencias pequeñas, aunque el resultado no será comparable
 else
  NFFT=32;  
 end

 xfft = x(1:NFFT) - mean(x(1:NFFT)); % para que no tenga continua
 L=length(xfft);	 	 
 X=fft(xfft,NFFT);	 	 
 Px=X.*conj(X)/(NFFT*L); %Power of each freq components	 	 
 fVals=fs*(0:NFFT/2-1)/NFFT;	 	 
 figure(20);
 plot(fVals,Px(1:NFFT/2),'b','LineSmoothing','on','LineWidth',1);	 	 
 title('One Sided Power Spectral Density');	 	 
 xlabel('Frequency (Hz)')	 	 
 ylabel('PSD');

  % SSA

 xssa = x - mean(x);
 M=8; % embedding dimension
 [E,V,A,R]=ssa_matlab(xssa, M); 
 senal_ssa = sum(R(:,1:SSA_SUAVIZADO),2);
 ruido_ssa = sum(R(:,SSA_SUAVIZADO+1:8),2);

 figure(21);
 eje = ejeParam1(1:length(xssa));
 plot(eje,xssa,eje,senal_ssa);
 title('SSA reconstruida sum(1:2) señal');	 	 
 figure(22);
 plot(eje,ruido_ssa);
 title('SSA reconstruida ruido');	 	 

% vamos a probar un estimador
% para la FFT, dejarla fija siempre a 16, esto deja el 2º componente en
% 0.3125 mm-1 , es decir, 3 mm por ciclo; el 3º en 0.625 mm-1, 1.6 mm por ciclo, sigue siendo
% un valor razonable para un cambio brusco. voy a dejar 2 primeras
% frecuencias como señal, del 3 al 8 como ruido.; cambio a 32 para mirar un
% parte significativa de más secuencias, 1:3 y 4:15
% para el SSA, voy a empezar haciendo el sum(rms) de los componentes 1:2 de
% 6 como señal, y el resto como ruido, o quizá 8 dimensiones para que sea
% comparable a la FFT???, si venga, 8 dimensiones,
%
 sn_fft = sum(Px(1:3)) / sum(Px(4:NFFT/2-1)); % 32 -> 16, 16->8 etc.

 sum_senal_ssa = sum(senal_ssa .^ 2);
 sum_ruido_ssa = sum(ruido_ssa .^ 2);
 sn_ssa = sum_senal_ssa / sum_ruido_ssa;
 if (size(senal_ssa) > 1 )
      hay_ssa = 1;
 else
     hay_ssa = 0;
 end
else
    hay_ssa = 0;
end % no se hace FFT of SSA si hay pocos puntos
end

%columnas en la excel empezando en C:   muestra, path, fecha, arq/moderna
%lamina/sección	isotopos previos	número de medida superpuesta	puntos
%deltaR ventana lineas lente/fibra ganancia t_retraso t_captura Npulses, integrado background RSD minimo	RSD medio	S/N FFT
%S/N FFT RSD de los espectros integrados , rsd_secuencia, experiment.txt
%
%
puntos = size(x,1); % número de puntos = ratios disponibles realmente 
lineas =  [ num2str(floor(LAMBDA_MAGNESIO)) '/' num2str(floor(LAMBDA_CALCIO)) 'nm'];
amplitud = (max(ratio_to_export(:,1)) - min(ratio_to_export(:,1))) / nanmeanLocal(ratio_to_export(:,1));

rsd_secuencia = nanstdLocal(x) /  nanmeanLocal(x);

rsd_minimo = min(rsd_ratio2l_bruto(rsd_ratio2l_bruto>0.01)); % a veces el rsd sale cero , pero queremos el valor mínimo, también hay espúreos, así que he puesto un threshold del 1%


%exportación a Excel
% columna1 : eje param1
% columna2 : ratio Mg/Ca
% columna3 : RSD del ratio (%)
% columna4 : intensidad media de la línea de magnesio en cada punto
% columna5 :  "" de calcio
% columna6 : ratio resina (v09)

if (DO_EXCEL==1)
 disp('Writing Excel file...');
 ahoraDatePath = strcat(datestr(now,'ddmmyyyy-HHMMSS'),SUBEXPERIMENT_NAME);

 %guarda el excel
 xlsFileName = strcat(PN, '\ratio_' , que_ratio_es_texto, '_', ahoraDatePath, '.xlsx');

  xlswrite(xlsFileName,[{'Carpeta:'} {PN}] , 'Hoja1', 'A1');
  xlswrite(xlsFileName,[{'Procesado:'} {para_el_titulo} ], 'Hoja1', 'A2');
  xlswrite(xlsFileName,[{'Descripción:'} {EXPERIMENT_DESCRIPTION} ], 'Hoja1', 'A3');
  if (DO_SSA==1)
    xlswrite(xlsFileName,[{'SN_SSA:'} {sn_ssa} ], 'Hoja1', 'A4');
  end
  xlswrite(xlsFileName,[{'RSD_secuencia:'} {rsd_secuencia} ], 'Hoja1', 'A5');  

  if (modo2D==1)
   cadena = cell(size(ejeParam1));
   for idx1= 1:Nparam1+1
     if (idx1==Nparam1+1)
       cadena(idx1)= { 'promedio' };
     else
       cadena(idx1)= { strcat(param1Name,'=', num2str(ejeParam1(idx1))) };
     end
   end
   xlswrite(xlsFileName, cadena , 'Hoja1', 'C10');
   for idx1= 1:Nparam1+1
     celda = strcat('C',num2str(idx1+10));
     if (idx1==Nparam1+1)
       xlswrite(xlsFileName,[ejeParam2' nanmeanLocal(ratio_to_export,1)'  ], 'Hoja1',celda);    
     else
       xlswrite(xlsFileName,[ejeParam2' ratio_to_export(idx1,:)' ], 'Hoja1',celda);  
   end
  end
  
 
 else
  xlswrite(xlsFileName,[ {param1Name} {['ratio-' WHICH_RATIO_TO_EXPORT]} {'rsd_ratio [%]'} {'media línea Mg'} {'media línea Ca'} {'ratio suavizado'} {'resina (ratio Si/Ca)'} ] , 'Hoja1', 'C10');
  xlswrite(xlsFileName,[ejeParam1' ratio_to_export(:,1) 100*rsd_ratio2l_bruto(:,1) mlineaMg(:,1) mlineaCa(:,1) ratio_resina_del_promedio_de_ratios(:,1) ], 'Hoja1','C11');
  if (DO_SSA==1)
    xlswrite(xlsFileName,senal_ssa + mean(x), 'Hoja1','H11');   
  end
 end 


 %crea el gráfico
 %http://www.mathworks.com/matlabcentral/newsreader/view_thread/298622 
 % esto está sacado de http://stackoverflow.com/questions/18165593/plot-excel-bar-chart-from-matlab
 e = actxserver('excel.application');
 eWs = e.Workbooks.Open(xlsFileName); % Open workbook
 eS = eWs.ActiveSheet;
 e.Visible = 0;

 %añade el gráfico
 eCO = eS.ChartObjects.Add(500, 30, 400, 250);
 eC = eCO.Chart;
 eC.SeriesCollection.NewSeries;
 howManyData = length(ratio_to_export(ratio_to_export>0)); % número real de datos
 datosRatio = ratio_to_export(ratio_to_export>0);

 if (hay_ssa==1)
  datosSuave = senal_ssa + mean(x);
  datosSuave = datosSuave(datosSuave>0);
  RangoArriba = round(((max(max([datosRatio datosSuave]))+0.05)*100))/100;
  RangoAbajo = round(((min(min([datosRatio datosSuave]))-0.05)*100))/100;
 else
  RangoArriba = round(((max(max(datosRatio ))+0.05)*100))/100;
  RangoAbajo = round(((min(min(datosRatio ))-0.05)*100))/100;
 end

 %eC.SeriesCollection.NewSeries;
 eC.SeriesCollection(1).Values = eS.Range(['D11:D' num2str(howManyData+10)]);
 eC.SeriesCollection(1).XValue = eS.Range(['C11:C' num2str(howManyData+10)]);
 eC.SeriesCollection(1).Name  = '=Hoja1!D10';
 eC.SeriesCollection.NewSeries;
 eC.SeriesCollection(2).Values = eS.Range(['H11:H'  num2str(howManyData+10)]);
 eC.SeriesCollection(2).XValue = eS.Range(['C11:C'  num2str(howManyData+10)]);
 eC.SeriesCollection(2).Name  = '=Hoja1!H10';

 eCO.Chart.ChartType = 74; %scatter with lines http://it.toolbox.com/wiki/index.php/EXCEL_Chart_Type_Enumeration 

 %ejes
 ChartAxes = invoke(eCO.Chart,'Axes',1);
 set(ChartAxes,'HasTitle',1);
 set(ChartAxes.AxisTitle,'Caption',param1Name');
 ChartAxes = invoke(eCO.Chart,'Axes',2);
 set(ChartAxes,'HasTitle',2);
 set(ChartAxes.AxisTitle,'Caption','ratio');
 eCO.Chart.Axes(2).MinimumScale = RangoAbajo;
 eCO.Chart.Axes(2).MaximumScale = RangoArriba;
 eCO.Chart.HasTitle = true;
 eCO.Chart.ChartTitle.Text = ['Ratio '  que_ratio_es_texto  ' '  EXPERIMENT_NAME]; 
 eWs.Save;
 eWs.Close(false);
 e.Quit;
 delete(e);
end % DO_EXCEL

%makes a copy of this script, matched to the XLS with now() timestamp
disp('saving a copy of this script with a timestamp just in case...');
copyfile([ mfilename('fullpath') '.m'] ,[ PN '\\' mfilename '_' datestr(now,'ddmmyyyy-HHMMSS') '.m']); %copy the current program to the experiment folder, easier to repeat

disp('done!');


