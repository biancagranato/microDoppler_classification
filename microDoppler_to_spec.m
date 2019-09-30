%% Introduction

% This script was used to simulate micro-Doppler effect signatures of three
%  objects: pedestrian, helicopter, parked car.
% Most of the simulation functions (beginning in line 145) were taken from 
%  MATLAB's implementation of Introduction to Micro Doppler effects:
%  https://www.mathworks.com/help/phased/examples/introduction-to-micro-doppler-effects.html
%
% The rest of the code was written by me and is under MIT License. 
% If you have any questions, don't hesitate to get in touch. 

%% Initial conditions
%Parameters of interest
tic
windows = {'flattopwin', 'hamming'}; %Window function
wsize = [64 256]; %Window size
olap = [0.2 0.8]; %Overlap (proportion)
clrs = [4 256]; %Number of colours in colour map

%% Make directories for each group of images
mkdir('Images') %Make folder where images are saved
for w = 1:numel(windows) %Window type
    for z = 1:numel(wsize) %Window size
        for o = 1:numel(olap) %Overlap (proportion)
            for c = 1:numel(clrs) %Number of colours in colour map
                %Person
                dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/person',...
                    windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c)); %Get folder name
                [~,msg] = mkdir(['Images/' dirnam]); %Make directory
                if ~isempty(msg) %Display message if any
                    fprintf('Person - W: %s, Z: %u, O = %u, C = %u\n%s\n',...
                        windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c),msg)
                end
                %Parked car
                dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/car',...
                    windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c)); %Get folder name
                [~,msg] = mkdir(['Images/' dirnam]); %Make directory
                if ~isempty(msg) %Display message if any
                    fprintf('W: %s, Z: %u, O = %u, C = %u\n%s\n',...
                        windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c),msg)
                end
                
                %Helicopter
                dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/helicopter',...
                    windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c)); %Get folder name
                [~,msg] = mkdir(['Images/' dirnam]); %Make directory
                if ~isempty(msg) %Display message if any
                    fprintf('W: %s, Z: %u, O = %u, C = %u\n%s\n',...
                        windows{w},wsize(z),fix(olap(o)*wsize(z)),clrs(c),msg)
                end
            end
        end
    end
end

%% The epic loop
prf = 2e4; %pulse repeat frequency
for x = [5 10] %vary starting x position
    for y = [2 8] %vary starting y position
        for z = [0 1] %vary starting z position
            for vx = [1 2]
                for vy = [0 1]
                    for vz = [0 1]
                        %Get micro doppler return
                        parkedcar_pos = [x+2;y;z-1]; %parked car position
                        ped_pos = [y-1;z;x]; %person's position
                        ped_vel = [vx, vy, vz]; %person's speed
                        tgtvel = [vx;vy;vz]*10; %helicopter is 10x as fast
                        tgtinitpos = [x;y;z]; %helicopter position
                        [mped, mcar] = getMdop(prf,parkedcar_pos,ped_pos,ped_vel); %call function
                        mdop = getHelicopter(tgtinitpos,tgtvel,prf);
                        
                        parfor wdw = 1:numel(windows)
                            for sz = wsize
                                for o = olap
                                    for c = clrs
                                        
                                        w = str2func(windows{wdw}); %convert window name to function
                                        
                                        %%%Person%%%
                                        %Get spectrogram
                                        S = spectrogram(mped,w(sz),fix(o*sz),512,prf,...
                                            'yaxis','centered'); %get spectrogram data
                                        
                                        %Convert spectrogram data to image
                                        imData = ind2rgb(im2uint8(abs(S)), parula(c)); %convert to image
                                        imData = imresize(imData, [227 227]); %resize
                                        
                                        %Save image
                                        dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/person/',...
                                            windows{wdw},sz,fix(o*sz),c); %Get directory
                                        sv = sprintf('pers_%s%usz%uolp%uclrs_%ux%uy%uz_%uvx%uvy%uvz.png',...
                                            windows{wdw},sz,fix(o*sz),c,x,y,z,vx,vy,vz); %Name to save image
                                        imwrite(imData,['Images/' dirnam sv]) %Save
                                        
                                        %%%Car%%%
                                        %Get spectrogram
                                        S = spectrogram(mcar,w(sz),fix(o*sz),512,prf,...
                                            'yaxis','centered'); %get spectrogram data
                                        
                                        %Convert spectrogram data to image
                                        imData = ind2rgb(im2uint8(abs(S)), parula(c)); %convert to image
                                        imData = imresize(imData, [227 227]); %resize
                                        
                                        %Save image
                                        dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/car/',...
                                            windows{wdw},sz,fix(o*sz),c); %Get directory
                                        sv = sprintf('car_%s%usz%uolp%uclrs_%ux%uy%uz_%uvx%uvy%uvz.png',...
                                            windows{wdw},sz,fix(o*sz),c,x,y,z,vx,vy,vz); %Name to save image
                                        imwrite(imData,['Images/' dirnam sv]) %Save
                                        
                                        
                                        %%%Helicopter%%%
                                        %Get spectrogram
                                        S = spectrogram(mdop,w(sz),fix(o*sz),512,prf,...
                                            'yaxis','centered'); %get spectrogram data
                                        
                                        %Convert spectrogram data to image
                                        imData = ind2rgb(im2uint8(abs(S)), parula(c)); %convert to image
                                        imData = imresize(imData, [227 227]); %resize
                                        
                                        %Save image
                                        dirnam = sprintf('%s/%u_size/%u_olap/%u_clrs/helicopter/',...
                                            windows{wdw},sz,fix(o*sz),c); %Get directory
                                        sv = sprintf('helicopter_%s%usz%uolp%uclrs_%ux%uy%uz_%uvx%uvy%uvz.png',...
                                            windows{wdw},sz,fix(o*sz),c,x,y,z,vx,vy,vz); %Name to save image
                                        imwrite(imData,['Images/' dirnam sv]) %Save
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
toc

%% micro-Doppler simulation
%The micro-Doppler simulation functions below were mostly taken from the MATLAB 
% introduction to the micro Doppler effect, which can be found here:
% https://www.mathworks.com/help/phased/examples/introduction-to-micro-doppler-effects.html
% The code above, however, is mostly written by me.
function [mped, mcar] = getMdop(prf,parkedcar_pos,ped_pos,ped_vel)
c = 2.998e+8;
fc = 5e9;
fs = 1e6;
pWdth = 2e-6;
ped_heading = 90;
ped_height = 1.8;

wav = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',pWdth,'PRF',prf);
tx  = phased.Transmitter;
rx  = phased.ReceiverPreamp;

egocar_pos = [0;0;0];
egocar_vel = [0;0;0];
egocar = phased.Platform('InitialPosition',egocar_pos,'Velocity',egocar_vel,...
    'OrientationAxesOutputPort',true);

parkedcar_vel = [0;0;0];
parkedcar = phased.Platform('InitialPosition',parkedcar_pos,'Velocity',parkedcar_vel,...
    'OrientationAxesOutputPort',true);
parkedcar_tgt = phased.RadarTarget('PropagationSpeed',c,'OperatingFrequency',fc,'MeanRCS',5);

ped = phased.BackscatterPedestrian('InitialPosition',ped_pos,'InitialHeading',ped_heading,...
    'PropagationSpeed',c,'OperatingFrequency',fc,'Height',ped_height,...
    'WalkingSpeed',sqrt(sum(ped_vel.^2)));

chan_ped = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,...
    'TwoWayPropagation',true,'SampleRate',fs);
chan_pcar = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,...
    'TwoWayPropagation',true,'SampleRate',fs);

Tsamp = 0.001;
Niter = 1e4/2;
npulse = round(fs/prf);
xr_ped = complex(zeros(npulse,Niter));
xr_pcar = complex(zeros(npulse,Niter));

for m = 1:Niter
    [pos_ego,vel_ego,~] = egocar(Tsamp);
    [pos_pcar,vel_pcar,~] = parkedcar(Tsamp);
    [pos_ped,vel_ped,ax_ped] = move(ped,Tsamp,0);
    [~,angrt_ped] = rangeangle(pos_ego,pos_ped,ax_ped);
    
    x = tx(wav());
    xt_ped = chan_ped(repmat(x,1,size(pos_ped,2)),pos_ego,pos_ped,vel_ego,vel_ped);
    xt_pcar = chan_pcar(x,pos_ego,pos_pcar,vel_ego,vel_pcar);
    xt_ped = reflect(ped,xt_ped,angrt_ped);
    xt_pcar = parkedcar_tgt(xt_pcar);
    xr_ped(:,m) = rx(xt_ped);
    xr_pcar(:,m) = rx(xt_pcar);
end

xd_ped = conj(dechirp(xr_ped,x));
xd_pcar = conj(dechirp(xr_pcar,x));

mped = sum(xd_ped);
mcar = sum(xd_pcar);

end

function out = getHelicopter(tgtinitpos,tgtvel,prf)
%Get variables
reflectivity = [10 .1 .1 .1 .1];
radarpos = [0;0;0];
radarvel = [0;0;0];
Nblades   = 4;
bladeang  = (0:Nblades-1)*2*pi/Nblades;
bladelen  = 6.5;
bladerate = deg2rad(4*360);  % rps -> rad/sec
fc = 5e9;
fs     = 1e6;
tx  = phased.Transmitter;
rx  = phased.ReceiverPreamp;
c = 2.998e+8;
lambda = c/fc;
NSampPerPulse = round(fs/prf);
Niter = 1e4;
y     = complex(zeros(NSampPerPulse,Niter));
pWdth = 2e-6;


tgtmotion  = phased.Platform('InitialPosition',tgtinitpos,'Velocity',tgtvel);
helicop = phased.RadarTarget('MeanRCS',reflectivity,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
wav = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',pWdth,'PRF',prf);
ura = phased.URA('Size',4,'ElementSpacing',lambda/2);
env = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,...
    'TwoWayPropagation',true,'SampleRate',fs);
txant = phased.Radiator('Sensor',ura,'PropagationSpeed',c,'OperatingFrequency',fc);
rxant = phased.Collector('Sensor',ura,'PropagationSpeed',c,'OperatingFrequency',fc);

rng(2018);
for m = 1:Niter
    % update helicopter motion
    t = (m-1)/prf;
    [scatterpos,scattervel,scatterang] = helicopmotion(t,tgtmotion,bladeang,...
        bladelen,bladerate,prf);
    
    % simulate echo
    x  = txant(tx(wav()),scatterang);                    % transmit
    xt = env(x,radarpos,scatterpos,radarvel,scattervel); % propagates to/from scatterers
    xt = helicop(xt);                                    % reflect
    xr = rx(rxant(xt,scatterang));                       % receive
    y(:,m) = sum(xr,2);                                  % beamform
end

mfcoeff = getMatchedFilter(wav);
mf  = phased.MatchedFilter('Coefficients',mfcoeff);
ymf = mf(y);
[~,ridx] = max(sum(abs(ymf),2)); % detection via peak finding along range

out = ymf(ridx,:);
end

function [scatterpos,scattervel,scatterang] = helicopmotion(...
    t,tgtmotion,BladeAng,ArmLength,BladeRate,prf)


radarpos = [0;0;0];
Nblades  = size(BladeAng,2);

[tgtpos,tgtvel] = tgtmotion(1/prf);

RotAng     = BladeRate*t;
scatterpos = [0 ArmLength*cos(RotAng+BladeAng);0 ArmLength*sin(RotAng+BladeAng);zeros(1,Nblades+1)]+tgtpos;
scattervel = [0 -BladeRate*ArmLength*sin(RotAng+BladeAng);...
    0 BladeRate*ArmLength*cos(RotAng+BladeAng);zeros(1,Nblades+1)]+tgtvel;

[~,scatterang] = rangeangle(scatterpos,radarpos);

end