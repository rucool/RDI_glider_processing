function RDI_batch_pd02nc(in, out, deployment)
% This function is used to process pd0 files from RDI-ADCP instruments and
% convert them to one netcdf.
% in = path to the directory containing the individual pd0s
% out= desired location of the .mat
% Sam Coakley
% 5/10/19


%Testing settings (this data is bad)
% in = '/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/ooi_cp_glider_389/dvl/';
% out = '/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/ooi_cp_glider_389/dvl/';
% deployment ='ooi_cp_glider_389';

%% Make sure the inputs are good
if ~strcmp(in(end),'/')
    in = [in '/'];
end
if ~strcmp(out(end),'/')
    out = [out '/'];
end
if ~isfolder(out)
    mkdir(out)
end

%% Read in and process the pd0s to mat files
raw_files=dir([in '*.pd0']);

% Prevent re-processing of files in the directory
for ii=1:length(raw_files) 
    if ~isfile([in raw_files(ii).name(1:end-3) 'mat'])
        temp=read_ADCP_PD0_v3([in raw_files(ii).name], [in raw_files(ii).name(1:end-4) '.mat']);
    end
end

%% Load the mats
proc_files=dir([in '*.mat']);
test_data=[];
for ii=1:length(proc_files)
    test_data=[test_data; load([proc_files(ii).name],'cfg','data')];
end

%% Write the mats to a netcdf

% Concat all the data from each segment
time=[];heading=[];pitch=[];roll=[];pressure=[];u1=[];
u2=[];u3=[];u4=[];c1=[];c2=[];c3=[];c4=[];pg1=[];pg2=[];pg3=[];pg4=[];
for ii=1:length(test_data)
    time=[time test_data(ii).data.time];
    heading=[heading test_data(ii).data.heading];
    pitch=[pitch test_data(ii).data.pitch];
    roll=[roll test_data(ii).data.roll];
    pressure=[pressure test_data(ii).data.pressure];
    u1=[u1 test_data(ii).data.u1];
    u2=[u2 test_data(ii).data.u2];
    u3=[u3 test_data(ii).data.u3];
    u4=[u4 test_data(ii).data.u4];
    c1=[c1 test_data(ii).data.c1];
    c2=[c2 test_data(ii).data.c2];
    c3=[c3 test_data(ii).data.c3];
    c4=[c4 test_data(ii).data.c4];
    pg1=[pg1 test_data(ii).data.pg1];
    pg2=[pg2 test_data(ii).data.pg2];
    pg3=[pg3 test_data(ii).data.pg3];
    pg4=[pg4 test_data(ii).data.pg4];
end

% Change filled values to NaN
u1(u1 == -32768) = NaN;
u2(u2 == -32768) = NaN;
u3(u3 == -32768) = NaN;
u4(u4 == -32768) = NaN;

% This would the place to do QA/QC for the OOI data
%%% Beams with correlation values below 64 in any given bin were removed
%%% from averaging. But just to make sure we will remove all data where
%%% there was any beam with too low of a correlation value
% corr_lim = 64;
% c1_low = c1 < corr_lim;
% c2_low = c2 < corr_lim;
% c3_low = c3 < corr_lim;
% c4_low = c4 < corr_lim;
% clow = c1_low + c2_low + c3_low + c4_low;
% u1(clow > 0) = NaN;
% u2(clow > 0) = NaN;
% u3(clow > 0) = NaN;
% u4(clow > 0) = NaN;
% clear c*_low clow

% Behavior exclusion
% Not sure if I should use NaNs or just chuck the data
% Since this is step 1 in the RDI processing I will leave it in as NaNs and
% keep the time and pressure stamps
% behav_ind = abs(pitch)<=15 | abs(roll)>=10;
% u1(:,behav_ind) = NaN;
% u2(:,behav_ind) = NaN;
% u3(:,behav_ind) = NaN;
% u4(:,behav_ind) = NaN;
% pitch(behav_ind) = NaN;
% roll(behav_ind) = NaN;
% heading(behav_ind) = NaN;

prof_num=size(u1,2);
bin_num = size(u1,1);
fname=[out deployment '_adcp.nc'];

nccreate(fname, 'time', 'Dimensions',{'Profile', prof_num}, 'Datatype', 'double');
ncwrite(fname, 'time', time);
nccreate(fname, 'heading', 'Dimensions',{'Profile', prof_num}, 'Datatype', 'double');
ncwrite(fname, 'heading', heading);
nccreate(fname, 'pitch', 'Dimensions',{'Profile', prof_num}, 'Datatype', 'double');
ncwrite(fname, 'pitch', pitch);
nccreate(fname, 'roll', 'Dimensions',{'Profile', prof_num}, 'Datatype', 'double');
ncwrite(fname, 'roll', roll);
nccreate(fname, 'pressure', 'Dimensions',{'Profile', prof_num}, 'Datatype', 'double');
ncwrite(fname, 'pressure', pressure);
nccreate(fname, 'u1', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'u1', u1');
nccreate(fname, 'u2', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'u2', u2');
nccreate(fname, 'u3', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'u3', u3');
nccreate(fname, 'u4', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'u4', u4');
nccreate(fname, 'c1', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'c1', c1');
nccreate(fname, 'c2', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'c2', c2');
nccreate(fname, 'c3', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'c3', c3');
nccreate(fname, 'c4', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'c4', c4');
nccreate(fname, 'pg1', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'pg1', pg1');
nccreate(fname, 'pg2', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'pg2', pg2');
nccreate(fname, 'pg3', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'pg3', pg3');
nccreate(fname, 'pg4', 'Dimensions',{'Profile', prof_num, 'Bin', bin_num}, 'Datatype', 'double');
ncwrite(fname, 'pg4', pg4');

% Write the cfg file as a global attribute
cfields=fieldnames(test_data(1).cfg);
cstruct=struct2cell(test_data(1).cfg);
for ii=1:length(cfields)
    if isa(cstruct{ii},'uint8') || isa(cstruct{ii},'uint16') ||...
            isa(cstruct{ii},'unit32') || isa(cstruct{ii},'unit64')
        ncwriteatt(fname,'/',cfields{ii},single(cstruct{ii}))
    else
        ncwriteatt(fname,'/',cfields{ii},single(cstruct{ii}))
    end
end

%% Delete the mat files
for ii=1:length(proc_files)
    delete([proc_files(ii).folder '/' proc_files(ii).name])
end
end