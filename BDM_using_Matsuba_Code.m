%% Replicate Matsuba analysis for signature 1000 ADCP 
% 
% This code has been developed replicate the code that was used in:
% 
% Matsuba, Y., Roelvink, D., Reniers, A. J. H. M., Rijnsdorp, D. P., 
% & Shimozono, T. (2022). Reconstruction of directional spectra of 
% infragravity waves. Journal of Geophysical Research: Oceans, 127, 
% e2021JC018273. https://doi.org/10.1029/2021JC018273
%
% but to use the beam velocities from a Nortek Signature 1000, rather than
% the ADCP used in their study. The main motivation for this revision is
% the additional beams in the signature 1000, as well as the non-level
% deployment which makes the beam geometry slightly more complicated
%
% The data used in that paper as well as the MATLAB codes are available at:
% https://figshare.com/articles/dataset/Dataset_for_Reconstruction_of_Directional_Spectra_of_Infragravity_Waves/17157902/1
%
% 
%% Main differences between their code and this code
%
% The main differences between their code and this version are:
%
% * Pitch, roll and heading angles are fed to the code (as opposed to hard
% coded)
% * that inbuild function CPSD is used to calculate the cross-spectral
% matrix rather than Matsuba's cross spectrum estimates
% * this code was developed to take signature data which has been processed
% through the Nortek data conversion tool, so some conventions that are in
% the sample data vary from the signature data (beam velocity direction is
% away from the instrument for a signature, also the bins are labelled
% nearest to farthest from the instrument... to name a few specific
% differences)
% 
%% Set-up for the analysis. This is the only section requiring human
% This section included some parameters which are likely to change, such as
% sampling rate, frequency limits, how many and which velocity bins to use
% from the signature. Where practicable, data are drawn from the config
% files that are created in the data conversion process, so those need to
% be available too. 
%
% This _should_ be the only section of the code that need to be modified
% when loading data from different instruments

clear

% Constants
g = 9.81; % gravity
rho = 1025; % water density 
deviation = 12.89; % the deviation, positive for East. Used in the magnetic to true conversion

% Time series / spectral estimation parametes

dt = 0.5; %the sampling rate of the downsampled data
fmin = 0.05; % the minimum frequency used in calculations +/- df 
fmax = 0.20; % the maximum frequency used in calculations +/- df 
window_seconds = 200; % window length in seconds
nfft = window_seconds / dt ; %the number of points in each window 
noverlap = nfft/2; % the overlap used 
window = hanning(nfft); % the window for cpsd estimates 
dtheta = 5; % the width of the angular bins for the final directional spectrum 

% Instrument specific parameters
alpha = 25; % the angle that the non-vertical beams are tilted from vertical 
first_bin = 1; % the first depth bin to use (for signature this is the deepest)
bin_skip = 3; % thin out the number of bins that are used for the analysis
min_bin_depth = -3; % the minimum acceptable depth (relative record mean depth) to use
Z0 = 0.3; % the distance above the bottom that the instrument was mounted  
number_of_beams = 4;
record_has_pressure = true;
record_has_altimeter = false;
velocity_save_format = 2;
    % There are a few different ways that velocity data can be saved
    % 1 == saved as 3 dimensional matrix: Burst_Velocity_Beam: [578972×nbeams×nbins single]
    % 2 == saved as 1 dimensional matrix: Burst_VelBeam1: [543072×20 single]
    %                                     Burst_VelBeam2: [543072×20 single]
    %                                     Burst_VelBeam3: [543072×20 single]
    %                                     Burst_VelBeam4: [543072×20 single]

% setup the way you want to ingest the datafiles. This is setup for the
% case where all the ad2cp files are in the same folder, and the site is
% identified by a unique code. In addition, they are chronologically
% indexed. 

% Sample Data 

data_directory = pwd;
instrument_codes = "Data_subset_CR10*small.mat";
save_file_code =   "CR10_dspec_BDM";
data_start = datetime(2021,08,30,00,00,00); % start time of good data
data_stop = datetime(2021,08,30,11,59,59); % end time of good data
dt_inst = 0.25; % sampling rate of the instrument

time_skip = dt/dt_inst;

%% Check that the Matsuba functions are available:

if exist("det_overshoot.m","file") ~= 2
    fprintf("The functions from Matsuba et al must be downloaded and in \n" + ...
        "the working directory. They are available from: \n" + ...
        "https://figshare.com/articles/dataset/Dataset_for_Reconstruction_of_Directional_Spectra_of_Infragravity_Waves/17157902/1\n")
    error("Matsuba functions are not in the path (see note above)")
end

%% Use the configuration above to setup for the data ingest 
%
% A lot of additional variables are created here which are dependent on how
% the user has configured above. While spelling them out here is somewhat
% redundant, it aids readability

% Setup some parameters based on the configuration above:
LT = 360/dtheta; %the number of bins, based on dtheta
theta = linspace(-pi,pi,LT+1); %theta in radians
    theta = theta(2:end); %this removes the double up at +/- pi
fs = 1/dt; % sampling frequency
fN = 0.5 / dt; % nyquist 
F = 0:1/nfft/dt:fN; % frequency vector based on nfft and dt
df = F(2)-F(1); 
index_f = (F > fmin - df & F < fmax + df); % frequencies within the nominal range
% index_FIG = F > fmin - df & F <= fIG;
% index_FSS = F > fIG & F < fmax + df;
f = F(index_f);
% f_SS = F(index_FSS);
% f_ig = F(index_FIG);

LF = length(f); % the length of the frequecny vector, it is called often
% LF_ig = length(f_ig);
% LF_SS = length(f_SS);

% the instrument relative unit vector of the beams 
if number_of_beams == 4 % 4 beams
    beams_inst_relative = [sind(alpha), 0,            -sind(alpha), 0;
                       0          , -sind(alpha), 0,            sind(alpha);
                       cosd(alpha), cosd(alpha),  cosd(alpha), cosd(alpha)];
elseif number_of_beams == 5
    beams_inst_relative = [sind(alpha), 0,            -sind(alpha), 0, 0;
                       0          , -sind(alpha), 0,            sind(alpha),0;
                       cosd(alpha), cosd(alpha),  cosd(alpha), cosd(alpha),1];
else; error("The number of beams is not a valid entry. Only 4 and 5 beams are supported")
end
    
% create an index to load the files for each location 
working_directory = pwd; 
addpath(data_directory)
cd(data_directory)
for ii = 1 : length(instrument_codes)
    data_file{ii} = dir(instrument_codes(ii));
end
cd(working_directory)

for site_code = 1 : length(instrument_codes)

% Now pull out data from the config file 
load(data_file{site_code}(1).name, "Config")
blanking_distance = Config.Burst_BlankingDistance;
cell_size = Config.Burst_CellSize;
n_cells = single(Config.Burst_NCells);
cell_range = double(blanking_distance + (1:n_cells)*cell_size);
slant_range = cell_range ./ cosd(alpha);
Burst_NSample = Config.Burst_NSample;

% remove the skipped bins (first_bin and bin_skip defined above):
slant_range = slant_range(first_bin:bin_skip:end);
number_of_bins = length(slant_range);
% a time vector which covers all the valid hours of data 
dtime = data_start(site_code) : hours(1) : data_stop(site_code);
LTIME = length(dtime); 

% preallocate for speed
D_f_theta = nan(LF,LT,LTIME); % directional spread
S_f_theta = nan(LF,LT,LTIME); % directional spectrum
S_f = nan(LF,LTIME); % surface energy density spectrum
% S_f_theta_bound_ig = nan(LF_ig, LT, LTIME); % bound IG waves
% S_f_theta_free_ig  = nan(LF_ig, LT, LTIME); % free IG waves
p =       nan(Burst_NSample,1);
eta =     nan(Burst_NSample,1);
pitch =   nan(size(p));
roll =    nan(size(p));
heading = nan(size(p));
u       = nan(Burst_NSample, number_of_beams, length(slant_range));

% now the fun bit, reading in the data from the nortek .mat output files 
load(data_file{site_code}(1).name, "Data")
file_number = 2;

file_dtime = datetime(Data.Burst_Time, "ConvertFrom", "datenum");

for record_time = 1 : LTIME
fprintf("file %s, record %d of %d \n", instrument_codes(site_code), record_time, LTIME)
good = file_dtime >= dtime(record_time) & file_dtime <= dtime(record_time) + hours(1);
n_good = sum(good);
if record_has_pressure; p(1:n_good) = Data.Burst_Pressure(good); end
if record_has_altimeter; eta(1:n_good) = Data.Burst_AltimeterDistanceAST(good); end %#ok<UNRCH> 
pitch(1:n_good) =     Data.Burst_Pitch(good);
roll(1:n_good) =      Data.Burst_Roll(good);
heading(1:n_good) =   Data.Burst_Heading(good);
if velocity_save_format == 2
    for ii = 1 : number_of_beams
        eval("u(1:n_good,ii,:) = Data.Burst_VelBeam"+ii+"(good,first_bin:bin_skip:end);");
    end
else; u(1:n_good,:,:) = Data.Burst_Velocity_Beam(1:n_good,:,first_bin:bin_skip:end);
end

    % this activates if there is not a complete record, and reads in the
    % next file 
    if sum(good) < Burst_NSample %& record_time ~= 424
        load(data_file{site_code}(file_number).name, "Data")
        file_number = file_number + 1;
        file_dtime = datetime(Data.Burst_Time, "ConvertFrom", "datenum");
        good = file_dtime >= dtime(record_time) & file_dtime <= dtime(record_time) + hours(1);
        loop_number = 1; % an escape variable for the next loop
        while sum(good) == 0 % in case the data isn't yet found. 
            load(data_file{site_code}(file_number).name, "Data")
            file_number = file_number + 1;
            file_dtime = datetime(Data.Burst_Time, "ConvertFrom", "datenum");
            good = file_dtime >= dtime(record_time) & file_dtime <= dtime(record_time) + hours(1);
            loop_number = loop_number + 1;
            if loop_number > length(data_file{site_code})
                error("Data not found in files")
            end
        end

        if record_has_pressure; p(n_good+1:end) = Data.Burst_Pressure(good); end
        if record_has_altimeter; eta(n_good+1:end) = Data.Burst_AltimeterDistanceAST(good); end %#ok<UNRCH> 

        pitch(n_good+1:end) =     Data.Burst_Pitch(good);
        roll(n_good+1:end) =      Data.Burst_Roll(good);
        heading(n_good+1:end) =   Data.Burst_Heading(good);
        if velocity_save_format == 2
        for ii = 1 : number_of_beams
        eval("u(n_good+1:end,ii,:) = Data.Burst_VelBeam"+ii+"(good,first_bin:bin_skip:end);");
        end
        else; u(n_good+1:end,:,:) = Data.Burst_Velocity_Beam(1:n_good,:,first_bin:bin_skip:end);
        end

    end

% down sample all the data
P = p(1:time_skip:end);
U = u(1:time_skip:end, :, :);
PITCH = pitch(1:time_skip:end);
ROLL = roll(1:time_skip:end);
HEADING = heading(1:time_skip:end);

% Ok so at this stage we now have one hour of data. From here the next
% thing to do is use the pitch, roll and heading data to calculate a normal
% vector for each beam so that we can later do the transformation
% 

h = mean(P) + Z0;
pitch_mean = mean(PITCH);
roll_mean = mean(ROLL);
heading_mean = mean(HEADING);

rot_head = [sind(heading_mean + deviation), -cosd(heading_mean + deviation), 0;
            cosd(heading_mean + deviation), sind(heading_mean + deviation),  0;
            0            , 0               1];
rot_pitch = [cosd(pitch_mean), 0, -sind(pitch_mean);
             0          , 1,  0;
             sind(pitch_mean), 0, cosd(pitch_mean) ];
rot_roll = [1, 0         , 0;
            0, cosd(roll_mean), -sind(roll_mean);
            0, sind(roll_mean), cosd(roll_mean) ];
T = rot_head * rot_pitch * rot_roll;
beam_ENU_unit_vectors = T * beams_inst_relative;

%
% calculate the location (x,y,z) relative to the instrument, with x
% increasing east, y increasing north, and the origin at the adcp
%

bin_loc = [0;0;0];
for ii = 1 : length(slant_range)
    bin_loc = [bin_loc beam_ENU_unit_vectors * slant_range(ii)];
end

% now convert the z - coordinate from instrument relative to an actual
% depth value 
instrument_depth = -h +Z0;
bin_loc(3,:) = instrument_depth + bin_loc(3,:);
good_bin = bin_loc(3,:) <= min_bin_depth;

% if altimeter is being used, then change the location back to the surface.
if record_has_altimeter; bin_loc(:,1) = [0;0;0]; end %#ok<UNRCH>

% Now generate a matrix of unit vectors with the same length as bin_loc 
bin_norm = nan(size(bin_loc));
bin_norm(:,1) = [0;0;1]; % normal vector for eta/pressure measurement 
beam = 1;
for ii = 2 : number_of_beams * number_of_bins + 1
    bin_norm(:,ii) = beam_ENU_unit_vectors(:,beam);
    beam = beam + 1;
    if beam == number_of_beams + 1
        beam = 1;
    end
end

% conduct cpsd on the data in the with the same parameters as the paper
M = sum(good_bin); %this differs because pressure is now used, not eta

N = M * (M+1) / 2; % this is the number elements needed to describe the
    %cross spectral matrix (given the symmetry due to conjugates)

% prepare the data into one array for simpler indexing
% the negative value for the velocity is used because that is how Matsuba
% does it, shouldn't matter, but seems odd. 
signal = nan(length(P), number_of_beams * number_of_bins + 2);

if record_has_altimeter; signal(:,1) = eta; else; signal(:,1) = P; end %#ok<UNRCH> 

for ii_beam = 1 : number_of_beams
    for ii_bin = 1 : number_of_bins
        ind = 1 + ii_beam + 4 * (ii_bin - 1 );
        signal(:,ind) = U(:,ii_beam,ii_bin);
    end
end

% now clear out the data from the bins which are above the threshold
% depth. At the same time, I will clear off the associated normal vectors
% and bin locations. 
bin_loc = bin_loc(:,good_bin);
bin_norm = bin_norm(:,good_bin);
signal = signal(:,good_bin);
signal = detrend(signal,2);

% The elevation spectrum, no need to do a transformation, I will need to
% come back to this and use pressure to ensure that there are no ghosts in
% the pressure transformation 
if record_has_pressure
Spp = cpsd(signal(:,1), signal(:,1), window, noverlap, nfft, fs);
Spp = Spp(index_f);
% need an Seta to use to scale the beam velocity, 
kw = wavek(1./f, h);
Kp = cosh(kw.*(h+bin_loc(end,1))) ./ cosh(kw .* h);
Seta = Spp'./(Kp.^2);
Seta = Seta(:);
S_f(:,record_time) = Seta;
else
Seta = cpsd(signal(:,1), signal(:,1), window, noverlap, nfft, fs); %#ok<UNRCH> 
Seta = Seta(index_f);
S_f(:,record_time) = Seta;
end

%calculate the power spectra of each signal

Su = nan(LF,M-1);
for ii = 1 : M-1
    X = cpsd(signal(:,ii+1),signal(:,ii+1),window,noverlap,nfft,fs);
    Su(:,ii) = X(index_f);
end

% the power spectra of all signals, ee is Matsuba terminology
if record_has_altimeter;  ee = [Spp, Su]; %#ok<UNRCH> 
else; ee = [Seta, Su];
end

% calculate the cross-spectra - phi

% note that they do not save their cross spectra as a square, but rather as
% a 1*(M*(M+1)/2) array. This is only a cosmetic difference, but I will use
% a square 

phi = nan(2*N, LF);
DW = nan(N, LF);

for ii = 1 : M
    for jj = 1 : ii

        counter = ii*(ii-1)/2+jj;

        X = cpsd(signal(:,ii),signal(:,jj),window,noverlap,nfft,fs);

        phi(counter,:) = real(X(index_f))./Seta;
        phi(counter + N, :) = -imag(X(index_f))./Seta;

        DW(counter,:) = sqrt(ee(:,ii).*ee(:,jj))./Seta;

    end
end

[S_f_theta(:,:,record_time), D_f_theta(:,:,record_time)] = dspec_BDM_pU(Seta, ...
    phi, DW, bin_loc, bin_norm, f, theta, h);


% % % Sometimes this is useful, especially when bug chasing
% % save each step:
% S_f_theta_temp = S_f_theta(:,:,record_time);
% D_f_theta_temp = D_f_theta(:,:,record_time);
% S_f_temp = S_f(:,record_time);
% cd(temp_folder)
% fname = "temp_file_"+save_file_code(site_code)+"_timestep_" + num2str(record_time,'%.3d');
% save(fname, "S_f_theta_temp", "D_f_theta_temp", "S_f_temp")
% cd(working_directory)

end

dir_naut_unsorted= mod(270 - rad2deg(theta),360);
[dir_naut, ind] = sort(dir_naut_unsorted);

dspec.S_f_theta_BDM = S_f_theta(:,ind,:);
dspec.D_f_theta_BDM = D_f_theta(:,ind,:);
dspec.S_f = S_f;
dspec.f = f(:);
dspec.direction_nautical = dir_naut(:);
dspec.dtime = dtime(:);
dspec.nfft = nfft;
dspec.analysis_period_hours = 1;

save(save_file_code(site_code),"dspec")

end

function [S_f_theta, D_f_theta] = dspec_BDM_pU(S_f, ...
    phi, DW, bin_loc, bin_norm, f, theta, h)

M = length(bin_loc);
N = M * (M+1) / 2;
dtheta = theta(2)-theta(1);

LT = length(theta);
LF = length(f);

S = nan(LF, LT);

D = -2.*eye(LT,LT);
D(1,2)=1;
D(1,end)=1;
D(end,end-1)=1;
D(end,1)=1;
for lt=2:(LT-1)
    D(lt, lt-1)=1;
    D(lt, lt+1)=1;
end

for lf = 1 : LF
% lf
ff = f(lf); % pulls out one frequency at a time
k = wavek(1./ff, h);
omega = 2 * pi * ff;

% now calculate a transfer function for this one frequency
H = zeros(M,LT);
% H(1,:) = 1; If eta is measure, otherwise use pressure 
H(1,:) = cosh(k*(h+bin_loc(3,1)))/cosh(k*h);
for jj = 2 : M
        
% extract the bin location 
    x = bin_loc(1,jj); y = bin_loc(2,jj); z = bin_loc(3,jj);

% extract the bin normal vector
    nx = bin_norm(1,jj); ny = bin_norm(2,jj); nz = bin_norm(3,jj);

% Transfer function
    H(jj,:) = omega .* cosh(k.*(h+z)) ./ sinh(k.*h) ...
        .*(nx.*cos(theta) + ny.*sin(theta) - 1i .* nz .*tanh(k.*(h+z))) ...
        .*exp(1i .* (k .* cos(theta) .* x + k .* sin(theta) .* y));

end

% use only the BDM model, so I have not included the if statement that is
% containd in Matsuba's code. From here out, it is just the Matsuba Code 

A0 = complex(zeros(2*N, LT));
B0 = complex(zeros(2*N, 1));
    for m=1:M
    for n=1:m
        counter=m*(m-1)/2+n;
        A0(counter,:)=real(H(m,:).*conj(H(n,:))).*dtheta./DW(counter,lf);
        A0(N+counter,:)=imag(H(m,:).*conj(H(n,:))).*dtheta./DW(counter,lf);
        B0(counter)=phi(counter,lf)./DW(counter,lf);
        B0(counter+N)=phi(counter+N,lf)./DW(counter,lf);
    end
    end

indexNaN=find(~isnan(B0));
A=A0(indexNaN, :);
B=B0(indexNaN, :);
N2=length(indexNaN);

d=[B;zeros(LT,1)];
ABIC=Inf;
i0=-Inf;

for i=-10:20
    u0=10^0.*(0.5.^(i-1));
    C=[A;abs(u0).*D];
    x0 = lsqnonneg(C, d);
    CI=transpose(A)*A+u0.^2*transpose(D)*D;
     [~, det0, r0]=det_overshoot(CI);
    logdetC=log(det0)+sum(log(r0));
    R=(sum((A*x0-B).^2)+u0.^2.*sum((D*x0).^2))/N2;
    ABIC0=N2*log(R)-LT*log(u0^2)+logdetC;
    if ABIC0<ABIC && ABIC0>-Inf
        ABIC=ABIC0;
        x=x0;
        i0=i;
    end
end
if i0>-Inf
    x=x./sum(x.*dtheta);
    S(lf,:)=S_f(lf).*x;
end
end

% rescale this S to extract the directional spread function and the
% directional spectra 

D_f_theta = nan(size(S));
S_f_theta = nan(size(S));

for lf = 1 : LF

%need to wrap the spread function around so that the integral is correct
X = [S(lf,end),S(lf,:)];
THETA = [-theta(end), theta];

area = trapz(THETA,X);

D_f_theta(lf,:) = S(lf,:)./area;
S_f_theta(lf,:) = D_f_theta(lf,:) .* S_f(lf);

end

end






