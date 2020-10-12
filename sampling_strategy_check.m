% THe goal of this script is to kind fo act as a second version to
% config_check.m. In config_check.m I just made a million plots trying to
% figure out which were the ones that actually showed what we wanted to
% see. So I am going to pick pieces from the code and consolidate it here
% as a more surefire sampling strategy check.

% Notes:
% Right now we are trying to check a few particular time periods which si
% why I wil choose some files and not others

% What we don't know:
% - What beams actually belong in our solution?
% - Actual depth range of the instrument due to tilt of mount and tilt of
% glider

%%
% Primers
working_dir = '/Users/coakley/Documents/MATLAB/ADCP_Wave/';
deployment = 'ru29_20190906T1535';
data_dir = [working_dir 'Data/' deployment '/deployment_1/'];

%% Convert all pd0s to .mats

if ~isfolder([data_dir 'mat/'])
% Make sure the file names do not have any spaces in them.

    cd([data_dir 'PD0/'])

    % Get a list of file names
    pd0 = dir(pwd);
    pd0 = {pd0(3:end).name};

    % Process all pd0s
    for ii = 1:length(pd0)
        read_ADCP_PD0_v3(pd0{ii}, [pd0{ii}(1:end-3) 'mat']);
        mkdir ../mat
        movefile('*.mat', '../mat');
    end

else
    disp('You already processed your pd0s my guy')
end

%% Check particular time stamps
% Our desired profiles to check
desired_date = datenum(2019,9,13,3,0,0);  %1
%desired_date = datenum(2019,9,21,18,0,0); %2
%desired_date = datenum(2019,9,23,21,0,0); %3

% Get a list of the data files
cd([data_dir 'mat/']);
mat = dir(pwd);
mat = {mat(3:end).name};

% This loop finds the first file that contains the desired_time
for ii = 1:length(mat)
    % Load in a file
    load(mat{ii});

    % Get start and end times for the segment
    ts = data.time(1);
    te = data.time(end);
    
    % Find the differences between the start and end times from our desired
    % time
    dif_start = ts - desired_date;
    dif_end   = te - desired_date;
    
    if dif_start <= 0 && dif_end >= 0
        desired_file = mat{ii};
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking Ruth's pd0%%%%%%%
% read_ADCP_PD0_v3('ruths.PD0', ['ruths.mat']);
% load ruths.mat

%% Check the data

% DVL does automatic quality checks so for now I will stick with just those
% unless instructed otherwise
% figure(1)
% pcolor(data.u1+data.u2+data.u3+data.u4)
% shading flat
% set(gca,'YDir','reverse')

% So for know, I will combine all beams together when checking the validity
% of a least square solution
%all_beams = data.u1+data.u2+data.u3+data.u4;
all_beams = data.u1+data.u3+data.u4;

% From looking at the data, ther are actual fill values of -32768 and where
% the instruments qc did not catch weak returns looks like 0's in all beams
all_beams_nan = all_beams;
all_beams_nan(all_beams_nan == 0 | all_beams_nan <= -10000) = NaN;

% I am also going to assume that all data below a bad return is also bad
for ii = 1:size(all_beams_nan,2)
    ind = find(isnan(all_beams_nan(:,ii)), 1, 'first');
    all_beams_nan(ind:end,ii) = NaN;
end

% This tells us how many bins to use in checking the condition of the G
% matrix
bins_by_ping = sum(~isnan(all_beams_nan),1);

figure(2)
plot(bins_by_ping)
xlabel('Ping number')
ylabel('# of good bins')

avg_good_bins = floor(mean(bins_by_ping));

%% Check condition of the G matrix
% bs = double(cfg.binsize)/100; %Bin size in m
% nbin = double(cfg.nbins);
% t = floor(nanmean(diff(data.time*(24*3600)))); % Convert time stamps from days to seconds
% w = 0.1;
% dived = floor(max(abs(data.transdep)));
% ping_num = 1;
% flip = 0;
% 
% dz = 75;
% [~, sol] = settings2Gmatrix(bs, nbin, t, w, dived, dz, ping_num, flip);
% dz

% For desired time 1, singular matrix at dz = 
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% For desired time 2, singular matrix at dz = 
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% For desired time 3, singular matrix at dz = 
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% Keep testing different dz's for accuracy and for different times

%% Rotate the beam velocities // remap depth bins // check cond of G matrix


% Removing fill values from the data
data.u1(data.u1 < -10000) = NaN;
data.u2(data.u2 < -10000) = NaN;
data.u3(data.u3 < -10000) = NaN;
data.u4(data.u4 < -10000) = NaN;

% Throw out any data with extreme behavior
% Pitch/Roll
ind = abs(data.pitch) > 30 | abs(data.pitch) < 20 | abs(data.roll) > 10;
data.u1(ind) = NaN;
data.u2(ind) = NaN;
data.u3(ind) = NaN;
data.u4(ind) = NaN;

% data.ei1(isnan(data.u1)) = NaN;
% data.ei2(isnan(data.u2)) = NaN;
% data.ei3(isnan(data.u3)) = NaN;
% data.ei4(isnan(data.u4)) = NaN;

% This function uses a linear extrapolation to place the velocities into
% the proper depth bin. This function as 
binmap_data = mapdepthcell_EH(data,cfg);

% QC should be done before conversion to beam coordinates
% Depth cells should also be remapped before this step because now we have
% ENU velocities in these bins
binmap_data = beamtoearth_EH(binmap_data,cfg);

% Create a grid of the depths of each bin
binmap_data.z = repmat(binmap_data.bins,[size(binmap_data.u1,2),1 ])';
for ii= 1:size(binmap_data.z,2)
    binmap_data.z(:,ii) = binmap_data.z(:,ii) + binmap_data.transdep(ii);
end

% Remove full profiles of NaNs
nans1 = sum(isnan(binmap_data.VE),1);
nans2 = sum(isnan(binmap_data.VN),1);

ind = nans1 == 15 | nans2 == 15;

binmap_data.VE = binmap_data.VE(:,~ind);
binmap_data.VN = binmap_data.VN(:,~ind);
binmap_data.z = binmap_data.z(:,~ind);

% Now we can run the actual inversion for 0 surface current
dz = 200;
[O_ls,G_ls,bnew,C] = inversion_leastSquare_sparse_2019(...
   binmap_data.VE',binmap_data.VN',binmap_data.z',dz,max(abs(binmap_data.transdep)),[0 0]);

% For desired time 1, singular matrix at dz = 60
%                   , poorly conditioned at dz = 
%                   , good for dz =
% For desired time 2, singular matrix at dz = 60, 100
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% For desired time 3, singular matrix at dz = 100
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% Keep testing different dz's for accuracy and for different times

% Having removed beam 2:
% For desired time 1, singular matrix at dz = 20, 50, 200
%                   , poorly conditioned at dz = 
%                   , good for dz =
% For desired time 2, singular matrix at dz = 200
%                   , poorly conditioned at dz = 
%                   , good for dz = 
% For desired time 3, singular matrix at dz = 200
%                   , poorly conditioned at dz = 
%                   , good for dz = 

%% Plots of echo intensity

% Make a scatter plot of ei1,2,3,4 in subplots
% Plot the mean value for each bin as well as the spread across each
% profile
xl = [40 140];

figure(3)
bin_num_rep = repmat([1:cfg.nbins],[size(binmap_data.ei1,2),1])';

subplot(151)
scatter(binmap_data.ei1(:),bin_num_rep(:),'bx','MarkerEdgeAlpha',0.1);
hold on
plot(nanmean(binmap_data.ei1,2), 1:cfg.nbins,'r-o')
xlim(xl)
title('Echo intensity beam 1')
set(gca,'YDir','reverse')

subplot(152)
scatter(binmap_data.ei2(:),bin_num_rep(:),'bx','MarkerEdgeAlpha',0.1);
hold on
plot(nanmean(binmap_data.ei2,2), 1:cfg.nbins,'r-o')
xlim(xl)
title('Echo intensity beam 2')
set(gca,'YDir','reverse')

subplot(153)
scatter(binmap_data.ei3(:),bin_num_rep(:),'bx','MarkerEdgeAlpha',0.1);
hold on
plot(nanmean(binmap_data.ei3,2), 1:cfg.nbins,'r-o')
xlim(xl)
title('Echo intensity beam 3')
set(gca,'YDir','reverse')

subplot(154)
scatter(binmap_data.ei4(:),bin_num_rep(:),'bx','MarkerEdgeAlpha',0.1);
hold on
plot(nanmean(binmap_data.ei4,2), 1:cfg.nbins,'r-o')
xlim(xl)
title('Echo intensity beam 4')
set(gca,'YDir','reverse')

subplot(155)
title('Mean EI by Beam')
plot(nanmean(binmap_data.ei1,2), 1:cfg.nbins,'r-o')
hold on
plot(nanmean(binmap_data.ei2,2), 1:cfg.nbins,'b-o')
plot(nanmean(binmap_data.ei3,2), 1:cfg.nbins,'g-o')
plot(nanmean(binmap_data.ei4,2), 1:cfg.nbins,'k-o')
legend('Mean Beam 1', 'Mean Beam 2', 'Mean Beam 3', 'Mean Beam 4','Location','southoutside') 
set(gca,'YDir','reverse')
