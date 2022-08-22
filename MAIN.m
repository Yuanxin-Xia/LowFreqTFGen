close all
clc
clear
fid = fopen('TF.txt','w');


% BASIC PARAMETER SETTING
% Wall impedance setting
Room = [
3 4.5 2.7 
3.5 10 2.7 
2.4 15 2.7 
4 6.75 2.7 
3.8 5.15 2.7 
4.2 5 2.7 
4 8 3 
5 5 3 
3.3 3.3 3.3 
2.4 20 2.4
];

% Set impedamce for each surface, if not set, the default value is absorption coefficient = 0.01
S1 = 'R';
S2 = 'R';
S3 = 'R';
S4 = 'R';
S5 = 'R';
S6 = 'R';

% Study frequency
UFreq = 354; % Upper frequency
LFreq = 1;  % Lower frequency
IFreq = 1;  % Frequency interval

% Source and receiver number
nums = 6; 
numr = 216;

% Initial transfer functions contaner
room = size(Room,1);
TF = zeros(room*nums*numr,3+floor((UFreq-LFreq)/IFreq+1)); % Transfer functions between one source and multiple receivers


import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress('COMSOL_Progress.txt');

%% LOOP THE ROOM
for r = 1:1:size(Room,1)
    Lx = Room(r,1);
    Ly = Room(r,2);
    Lz = Room(r,3);

[Source,Receiver] = Point_generation(Lx,Ly,Lz,0);


%% LOOP SOURCE POINTs

for u = 1:size(Source,1) % Loop every source
    disp(strcat('This is room: ',num2str(r),' ,and this is source: ',num2str(u)))
    tic % set a timer
    
    % Create a new model and show progress
    model = ModelUtil.create(strcat('Model',num2str(r)));
    model.modelPath(pwd);

    % Frequency vector
    Freq = strcat('range(',num2str(LFreq),',',num2str(IFreq),',',num2str(UFreq),')');

    % Source 1
    model.param.set('S1_X', strcat(num2str(Source(u,1)),'[m]'));
    model.param.set('S1_Y', strcat(num2str(Source(u,2)),'[m]'));
    model.param.set('S1_Z', strcat(num2str(Source(u,3)),'[m]'));

    % Set receiver 1 to N
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:size(Receiver,1)
        model.param.set(strcat('M',num2str(i),'_X'), strcat(num2str(Receiver(i,1)),'[m]'));
        model.param.set(strcat('M',num2str(i),'_Y'), strcat(num2str(Receiver(i,2)),'[m]'));
        model.param.set(strcat('M',num2str(i),'_Z'), strcat(num2str(Receiver(i,3)),'[m]'));
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Room geometry
    model.param.set('LX', strcat(num2str(Lx),'[m]'));
    model.param.set('LY', strcat(num2str(Ly),'[m]'));
    model.param.set('LZ', strcat(num2str(Lz),'[m]'));

    % Ambient parameter
    model.param.set('rm', '10.9', 'Specific flow resistivity');
    model.param.set('rho', '1.2[kg*m^-3]', 'Air density');
    model.param.set('c', '343[m/s]', 'sound speed');
    model.param.label('Parameters 1');

    model.component.create('comp1', true);

    model.component('comp1').geom.create('geom1', 3);

    model.component('comp1').label('Component 1');

    tbl1 = model.result.table.create('tbl1', 'Table');

    % Impedance interpolation
    % Impedance A
    model.func.create('int1', 'Interpolation');
    model.func('int1').set('source', 'file');
    model.func('int1').label('ZA Interpolation');
    model.func('int1').set('funcs', {'ZAreal' '1'; 'ZAimag' '2'});
    model.func('int1').set('filename', strcat(pwd,'/ImpedanceFile/A.txt'));
    model.func('int1').set('fununit', {'Pa*s/m' 'Pa*s/m'});
    model.func('int1').set('argunit', {'Hz'});
    
    % Impedance B
    model.func.create('int2', 'Interpolation');
    model.func('int2').set('source', 'file');
    model.func('int2').label('ZB Interpolation');
    model.func('int2').set('funcs', {'ZBreal' '1'; 'ZBimag' '2'});
    model.func('int2').set('filename', strcat(pwd,'/ImpedanceFile/B.txt'));
    model.func('int2').set('fununit', {'Pa*s/m' 'Pa*s/m'});
    model.func('int2').set('argunit', {'Hz'});

    % Impedance C
    model.func.create('int3', 'Interpolation');
    model.func('int3').set('source', 'file');
    model.func('int3').label('ZR Interpolation');
    model.func('int3').set('funcs', {'ZRreal' '1'; 'ZRimag' '2'});
    model.func('int3').set('filename', strcat(pwd,'/ImpedanceFile/R.txt'));
    model.func('int3').set('fununit', {'Pa*s/m' 'Pa*s/m'});
    model.func('int3').set('argunit', {'Hz'});
 
    %% CREATE GEOMETRY
    model.component('comp1').mesh.create('mesh1');
	
    % Create a block

    model.component('comp1').geom('geom1').geomRep('comsol');
    model.component('comp1').geom('geom1').create('blk1', 'Block');
    model.component('comp1').geom('geom1').feature('blk1').set('size', {'LX' 'LY' 'LZ'});
    
    
    % Set the room geometry to the block
    model.component('comp1').geom('geom1').feature('blk1').set('size', {'LX' 'LY' 'LZ'});

    % Create points: source and receiver
    
    % Source
    model.component('comp1').geom('geom1').create('pt1', 'Point');
    model.component('comp1').geom('geom1').feature('pt1').label('Source1');
    model.component('comp1').geom('geom1').feature('pt1').set('p', {'S1_X' 'S1_Y' 'S1_Z'});

    % Receivers 1 to N
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:size(Receiver,1)
        model.component('comp1').geom('geom1').create(strcat('pt',num2str(i+1)), 'Point');
        model.component('comp1').geom('geom1').feature(strcat('pt',num2str(i+1))).label(strcat('Mic',num2str(i)));
        model.component('comp1').geom('geom1').feature(strcat('pt',num2str(i+1))).set('p', {strcat('M',num2str(i),'_X') strcat('M',num2str(i),'_Y') strcat('M',num2str(i),'_Z')});
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Form Union
    model.component('comp1').geom('geom1').feature('fin').label('Form Union');
    model.component('comp1').geom('geom1').run;
    model.component('comp1').geom('geom1').run('fin');

    %% SELECTION POINTS
    
    % Rerank source and receiver. Reason? Because the receiver and
    % source are merged together in COMSOL modelling process, so we need to
    % find the source again and assign the rigth number to it
    
    order = find(Receiver(:,1)>Source(u,1),1); % Rerank source and receiver
    
    
    % Select source
    model.component('comp1').selection.create('sel1', 'Explicit');
    model.component('comp1').selection('sel1').geom('geom1', 0);
    model.component('comp1').selection('sel1').set([4+order]);

    % Select receivers and jump source number
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:size(Receiver,1)
        if i < order 
        model.component('comp1').selection.create(strcat('sel',num2str(i+1)), 'Explicit');
        model.component('comp1').selection(strcat('sel',num2str(i+1))).geom('geom1', 0);
        model.component('comp1').selection(strcat('sel',num2str(i+1))).set([4+i]);
        elseif i>=order
        model.component('comp1').selection.create(strcat('sel',num2str(i+1)), 'Explicit');
        model.component('comp1').selection(strcat('sel',num2str(i+1))).geom('geom1', 0);
        model.component('comp1').selection(strcat('sel',num2str(i+1))).set([5+i]);
        end
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Assign a name to source
    model.component('comp1').selection('sel1').label('Source 1');

    % Assign names to receivers
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:size(Receiver,1)
        model.component('comp1').selection(strcat('sel',num2str(i+1))).label(strcat('Mic',num2str(i)));
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    %% MATERIAL SETTING
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
    model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
    model.component('comp1').material('mat1').propertyGroup.create('NonlinearModel', 'Nonlinear model');
    model.component('comp1').material('mat1').propertyGroup.create('idealGas', 'Ideal gas');
    model.component('comp1').material('mat1').propertyGroup('idealGas').func.create('Cp', 'Piecewise');

    %% Create boundary condition, probe and material
    model.component('comp1').physics.create('acpr', 'PressureAcoustics', 'geom1');

    % Create boundary condition (Initialization)

    % Surface impedance 1 imp1
    model.component('comp1').physics('acpr').create('imp1', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp1').selection.set([1 2 3 4 5 6]);

    % Surface impedance 2 imp2
    model.component('comp1').physics('acpr').create('imp2', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp2').selection.set([1]);

    % Surface impedance 3 imp3
    model.component('comp1').physics('acpr').create('imp3', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp3').selection.set([2]);

    % Surface impedance 4 imp4
    model.component('comp1').physics('acpr').create('imp4', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp4').selection.set([3]);

    % Surface impedance 5 imp5
    model.component('comp1').physics('acpr').create('imp5', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp5').selection.set([4]);

    % Surface impedance 6 imp6
    model.component('comp1').physics('acpr').create('imp6', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp6').selection.set([5]);

    % Surface impedance 7 imp7
    model.component('comp1').physics('acpr').create('imp7', 'Impedance', 2);
    model.component('comp1').physics('acpr').feature('imp7').selection.set([6]);

    % Point source
    model.component('comp1').physics('acpr').create('pp1', 'PressurePoint', 0);
    model.component('comp1').physics('acpr').feature('pp1').selection.named('sel1');

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Point probe
    for i = 1:size(Receiver,1)
        model.component('comp1').probe.create(strcat('point',num2str(i)), 'Point');
        model.component('comp1').probe(strcat('point',num2str(i))).selection.named(strcat('sel',num2str(i+1)));
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Material
    model.thermodynamics.label('Thermodynamics');

    model.material.label('Materials');
    model.component('comp1').material('mat1').label('Air');
    model.component('comp1').material('mat1').set('family', 'air');
    model.component('comp1').material('mat1').propertyGroup('def').label('Basic');
    model.component('comp1').material('mat1').propertyGroup('def').func('eta').label('Piecewise');
    model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
    model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
    model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
    model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
    model.component('comp1').material('mat1').propertyGroup('def').func('Cp').label('Piecewise 2');
    model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
    model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
    model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
    model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').label('Analytic');
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa' 'K'});
    model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '293.15'});
    model.component('comp1').material('mat1').propertyGroup('def').func('k').label('Piecewise 3');
    model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
    model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
    model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
    model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').label('Analytic 2');
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
    model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').label('Analytic 1');
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa' 'K'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').label('Analytic 2a');
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
    model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
    model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
    model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
    model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
    model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
    model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
    model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
    model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
    model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
    model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
    model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
    model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').label('Refractive index');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').info('category').label('Information');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
    model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
    model.component('comp1').material('mat1').propertyGroup('NonlinearModel').label('Nonlinear model');
    model.component('comp1').material('mat1').propertyGroup('NonlinearModel').info('category').label('Information');
    model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
    model.component('comp1').material('mat1').propertyGroup('idealGas').label('Ideal gas');
    model.component('comp1').material('mat1').propertyGroup('idealGas').info('category').label('Information');
    model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
    model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('arg', 'T');
    model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
    model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
    model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
    model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
    model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
    model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
    model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897');
    model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
    model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
    model.component('comp1').material('mat1').materialType('nonSolid');

    model.component('comp1').coordSystem('sys1').label('Boundary System 1');

    model.common('cminpt').label('Default Model Inputs');

    % label
    model.component('comp1').physics('acpr').label('Pressure Acoustics, Frequency Domain');

    %% Set boundary condition
    % 1. Basic: Pressure acoustic
    model.component('comp1').physics('acpr').feature('fpam1').label('Pressure Acoustics 1');
    model.component('comp1').physics('acpr').feature('fpam1').featureInfo('info').label('Equation View');

    % 2. Basic: Sound Hard Boundary (Wall)
    model.component('comp1').physics('acpr').feature('shb1').label('Sound Hard Boundary (Wall) 1');
    model.component('comp1').physics('acpr').feature('shb1').featureInfo('info').label('Equation View');

    % 3. Basic: Initial Values
    model.component('comp1').physics('acpr').feature('init1').label('Initial Values 1');
    model.component('comp1').physics('acpr').feature('init1').featureInfo('info').label('Equation View');

    % 4. Basic: Continuity
    model.component('comp1').physics('acpr').feature('dcont1').label('Continuity 1');
    model.component('comp1').physics('acpr').feature('dcont1').featureInfo('info').label('Equation View');

    % 5. Rigid wall absorption coefficient
    model.component('comp1').physics('acpr').feature('imp1').set('ImpedanceModel', 'AbsorptionCoefficient');
    model.component('comp1').physics('acpr').feature('imp1').label('Rigid wall');
    model.component('comp1').physics('acpr').feature('imp1').featureInfo('info').label('Equation View');

    % 6. Wall impedance
    
    % wall 1
    if exist('S1','var')
        model.component('comp1').physics('acpr').feature('imp2').set('Zi', strcat('Z',S1,'real(freq)+i*Z',S1,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp2').label('Impedance 1');
        model.component('comp1').physics('acpr').feature('imp2').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp2').active(false);
    end
    
    % wall 2
    if exist('S2','var')
        model.component('comp1').physics('acpr').feature('imp3').set('Zi', strcat('Z',S2,'real(freq)+i*Z',S2,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp3').label('Impedance 2');
        model.component('comp1').physics('acpr').feature('imp3').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp3').active(false);
    end
    
    % wall 3
    if exist('S3','var')
        model.component('comp1').physics('acpr').feature('imp4').set('Zi', strcat('Z',S3,'real(freq)+i*Z',S3,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp4').label('Impedance 3');
        model.component('comp1').physics('acpr').feature('imp4').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp4').active(false);
    end
    
    % wall 4
    if exist('S4','var')
        model.component('comp1').physics('acpr').feature('imp5').set('Zi', strcat('Z',S4,'real(freq)+i*Z',S4,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp5').label('Impedance 4');
        model.component('comp1').physics('acpr').feature('imp5').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp5').active(false);
    end
    
    % wall 5
    if exist('S5','var')
        model.component('comp1').physics('acpr').feature('imp6').set('Zi', strcat('Z',S5,'real(freq)+i*Z',S5,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp6').label('Impedance 5');
        model.component('comp1').physics('acpr').feature('imp6').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp6').active(false);
    end
    
    % wall 6
    if exist('S6','var')
        model.component('comp1').physics('acpr').feature('imp7').set('Zi', strcat('Z',S6,'real(freq)+i*Z',S6,'imag(freq)'));
        model.component('comp1').physics('acpr').feature('imp7').label('Impedance 6');
        model.component('comp1').physics('acpr').feature('imp7').featureInfo('info').label('Equation View');
    else
        model.component('comp1').physics('acpr').feature('imp7').active(false);
    end
    
    % 7. Source term
    model.component('comp1').physics('acpr').feature('pp1').set('p0', 1);
    model.component('comp1').physics('acpr').feature('pp1').label('Source');
    model.component('comp1').physics('acpr').feature('pp1').featureInfo('info').label('Equation View');

    disp('Begin mesh!')
    disp( datestr(now, 'yy/mm/dd-HH:MM'))
    %% MESH
    model.component('comp1').mesh('mesh1').label('Mesh 1');

    disp('End mesh!')
    disp( datestr(now, 'yy/mm/dd-HH:MM'))
    %% PROBE SETTING

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i = 1:size(Receiver,1)
        model.component('comp1').probe(strcat('point',num2str(i))).label(strcat('Point Probe ',num2str(i)));
        model.component('comp1').probe(strcat('point',num2str(i))).set('expr', 'acpr.Lp_t');
        model.component('comp1').probe(strcat('point',num2str(i))).set('unit', 'dB');
        model.component('comp1').probe(strcat('point',num2str(i))).set('descr', 'Total sound pressure level');
        model.component('comp1').probe(strcat('point',num2str(i))).set('table', 'tbl1');
        model.component('comp1').probe(strcat('point',num2str(i))).set('window', 'window2');
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
 
    %% STUDY STEETING
    %
    model.study.create('std1');
    model.study('std1').create('freq', 'Frequency');

    % solver configuration
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('p1', 'Parametric');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').create('d1', 'Direct');
    model.sol('sol1').feature('s1').create('i1', 'Iterative');
    model.sol('sol1').feature('s1').create('i2', 'Iterative');
    model.sol('sol1').feature('s1').create('i3', 'Iterative');
    model.sol('sol1').feature('s1').create('i4', 'Iterative');
    model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i3').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature('i4').create('dd1', 'DomainDecomposition');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').create('mg1', 'Multigrid');
    model.sol('sol1').feature('s1').feature.remove('fcDef');

    %% RESULT SETTING
    % dataset 
    model.result.dataset.create('dset2', 'Solution');
    model.result.dataset('dset2').set('probetag', 'point2');

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       for i = 1:size(Receiver,1)
            if i< order
                model.result.dataset.create(strcat('avh',num2str(i)), 'Average');
                % average value
                model.result.dataset(strcat('avh',num2str(i))).set('probetag', strcat('point',num2str(i)));
                model.result.dataset(strcat('avh',num2str(i))).set('data', 'dset2'); % keep dset2
                model.result.dataset(strcat('avh',num2str(i))).selection.geom('geom1', 0); % keep geom1
                model.result.dataset(strcat('avh',num2str(i))).selection.set([4+i]);
                model.result.numerical.create(strcat('pev',num2str(i)), 'EvalPoint');
                model.result.numerical(strcat('pev',num2str(i))).set('probetag', strcat('point',num2str(i)));
            elseif i>=order
                model.result.dataset.create(strcat('avh',num2str(i)), 'Average');
                % average value
                model.result.dataset(strcat('avh',num2str(i))).set('probetag', strcat('point',num2str(i)));
                model.result.dataset(strcat('avh',num2str(i))).set('data', 'dset2'); % keep dset2
                model.result.dataset(strcat('avh',num2str(i))).selection.geom('geom1', 0); % keep geom1
                model.result.dataset(strcat('avh',num2str(i))).selection.set([5+i]);
                model.result.numerical.create(strcat('pev',num2str(i)), 'EvalPoint');
                model.result.numerical(strcat('pev',num2str(i))).set('probetag', strcat('point',num2str(i)));
            end
        end

    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % PlotGroup1D
    model.result.create('pg4', 'PlotGroup1D');
    model.result('pg4').set('probetag', 'window2_default');
    model.result('pg4').create('tblp1', 'Table');
    
    
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    string = '';
    for i = 1:size(Receiver,1)
        string = strcat('point',num2str(i),',',string);
    end
    string = string(1:end-1);

    model.result('pg4').feature('tblp1').set('probetag', string);

    for i = 1:size(Receiver,1)
        model.component('comp1').probe(strcat('point',num2str(i))).genResult([]);
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    %% STUDY
    model.study('std1').label('Study 1');
    model.study('std1').feature('freq').label('Frequency Domain');
    model.study('std1').feature('freq').set('plist', Freq);

    model.sol('sol1').attach('std1');
    model.sol('sol1').label('Solution 1');
    model.sol('sol1').feature('st1').label([native2unicode(hex2dec({'7f' '16'}), 'unicode')  native2unicode(hex2dec({'8b' 'd1'}), 'unicode')  native2unicode(hex2dec({'65' 'b9'}), 'unicode')  native2unicode(hex2dec({'7a' '0b'}), 'unicode') ': Frequency Domain']);
    model.sol('sol1').feature('v1').label('Dependent Variables 1');
    model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
    model.sol('sol1').feature('v1').set('cname', {'freq'});
    model.sol('sol1').feature('v1').set('clist', {strcat(Freq,'Hz')});



    model.sol('sol1').feature('s1').label('Stationary Solver 1');
    model.sol('sol1').feature('s1').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('aDef').label('Advanced');
    model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
    model.sol('sol1').feature('s1').feature('pDef').label([native2unicode(hex2dec({'53' 'c2'}), 'unicode')  native2unicode(hex2dec({'65' '70'}), 'unicode')  native2unicode(hex2dec({'53' '16'}), 'unicode') ' 2']);
    model.sol('sol1').feature('s1').feature('pDef').set('pname', {'freq'});
    model.sol('sol1').feature('s1').feature('pDef').set('plistarr', {Freq});

    model.sol('sol1').feature('s1').feature('pDef').set('punit', {'Hz'});
    model.sol('sol1').feature('s1').feature('pDef').set('pcontinuationmode', 'no');
    model.sol('sol1').feature('s1').feature('pDef').set('preusesol', 'auto');
    model.sol('sol1').feature('s1').feature('pDef').set('uselsqdata', false);
    model.sol('sol1').feature('s1').feature('p1').label('Parametric 1');
    model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
    model.sol('sol1').feature('s1').feature('p1').set('plistarr', {Freq});

    model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
    model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
    model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
    model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
    model.sol('sol1').feature('s1').feature('d1').label('Suggested Direct Solver (acpr)');
    model.sol('sol1').feature('s1').feature('i1').label('Suggested Iterative Solver (GMRES with GMG) (acpr)');
    model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver');
    model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('i2').label('Suggested Iterative Solver (FGMRES with GMG) (acpr)');
    model.sol('sol1').feature('s1').feature('i2').set('linsolver', 'fgmres');
    model.sol('sol1').feature('s1').feature('i2').feature('ilDef').label('Incomplete LU');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').label('Multigrid 1');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').label('Presmoother');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').label('Postsmoother');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').label('Coarse Solver');
    model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('i3').label('Suggested Iterative Solver (Shifted Laplace) (acpr)');
    model.sol('sol1').feature('s1').feature('i3').feature('ilDef').label('Incomplete LU');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').label('Multigrid 1');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('mcasegen', 'coarse');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('scale', 3);
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('slaplacemain', true);
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('epsslaplacemain', {'acpr' '0.4'});
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('slaplacemg', true);
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').set('epsslaplacemg', {'acpr' '0.4'});
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('pr').label('Presmoother');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('po').label('Postsmoother');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('po').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('cs').label('Coarse Solver');
    model.sol('sol1').feature('s1').feature('i3').feature('mg1').feature('cs').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('i4').label('Suggested Iterative Solver (Domain Decomposition) (acpr)');
    model.sol('sol1').feature('s1').feature('i4').feature('ilDef').label('Incomplete LU');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').label('Domain Decomposition (Schwarz) 1');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('prefun', 'ddhyb');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('dompernodemaxactive', true);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('ddolhandling', 'ddrestricted');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('userac', false);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('usecoarse', false);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('ddboundary', 'absorbing');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('alphaabsorbing', {'acpr' '1'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('sndorderabsorbing', {'acpr' 'on'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('betaabsorbing', {'acpr' '0.1'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('slaplacemain', true);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('epsslaplacemain', {'acpr' '0.4'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('slaplacemg', true);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').set('epsslaplacemg', {'acpr' '0.4'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('cs').label('Coarse Solver');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('cs').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').label('Domain Solver');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('dDef').label('Direct');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').label('Multigrid 1');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').set('slaplacemg', true);
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').set('epsslaplacemg', {'acpr' '0.4'});
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('pr').label('Presmoother');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('po').label('Postsmoother');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('po').feature('soDef').label('SOR 1');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('cs').label('Coarse Solver');
    model.sol('sol1').feature('s1').feature('i4').feature('dd1').feature('ds').feature('mg1').feature('cs').feature('dDef').label('Direct');
    disp('Begin study!')
    disp( datestr(now, 'yy/mm/dd-HH:MM'))
    model.sol('sol1').runAll;


    str = mphtable(model,'tbl1'); 
    tbl_data = str.data;
    tb1_data1 = tbl_data(:,2:end)';

    toc

    % The first three colomn of transfer functions is room size
        for i = 1:size(Receiver,1)
            TF((r-1)*nums*numr+(u-1)*numr+i,:) = [Lx,Ly,Lz,tb1_data1(i,:)];
            fprintf(fid,'%d\t',[Lx,Ly,Lz,tb1_data1(i,:)]);
            fprintf(fid,'\n');
        end


disp('End study!')
disp( datestr(now, 'yy/mm/dd-HH:MM'))   
disp('---------------------------------------------------------------------')
end
    disp('*******************************************')
    disp(strcat('*********Room ',num2str(r),'finished***********'))
    disp('*******************************************')
end
fclose(fid);
save('TF.mat','TF')
writematrix(TF,'TF.xls','Sheet',1)
