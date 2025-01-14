






%% Set environment


close all
clearvars
clc

addpath(genpath('/Users/jpohaku/Documents/MATLAB/mTRF-Toolbox-master/'));
addpath  '/Users/jpohaku/Documents/MATLAB/speechenv_toolbox/';
addpath  '/Users/jpohaku/Documents/MATLAB/mocaptoolbox-1.6.6/';


homefolder='/Users/jpohaku/Documents/gocsg/';

dsh=homefolder(end);
cd(homefolder)
addpath(genpath(homefolder))
format long g

%% Get kinemative data summary

load([ homefolder 'setup_files/jointsegmentkeys.mat']);

fileidx=dir([homefolder 'TRFoutput/Study1/']);
fileidx = fileidx(~[fileidx.isdir]);
clearvars TRF_study
for f=1:size(fileidx,1)
    TRF_study.(strrep(fileidx(f).name,'.mat',''))=load([fileidx(f).folder '/' fileidx(f).name]);
end
filenames=fieldnames(TRF_study);
outputnames=fieldnames(TRF_study.(filenames{1}).TRF_master);
outputnames(end)=[];
outputnames(4:end)=[];


TRF_study_vel=TRF_study;
for t=1:length(outputnames)
        counter=0;
        for f=1:length(filenames)
            fns=fieldnames(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest);
            for fn=1:length(fns)
                    for bp=1:20
TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)=gradient(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp));
                    end
            end
        end
end
TRF_study_accel=TRF_study;


TRF_study=TRF_study_vel;
clear meanholdall stdholdall pholdall nholdall varholdall_vel
for t=1:length(outputnames)
    clear mh_p_all mh_n_all
    for bp=1:20
        clear mh_p mh_n 
        counter=0;
        for f=1:length(filenames)
            fns=fieldnames(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).test);
            for fn=1:length(fns)
                counter=counter+1;
                if counter==1
                    mh_p=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)>0,bp);
                    mh_n=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)<0,bp);
                else
                    mh_p=vertcat(mh_p,TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)>0,bp));
                    mh_n=vertcat(mh_n,TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)<0,bp));
                end
            end
        end
        mh_p_all(bp,1)= mean(mh_p);
        mh_p_all(bp,2)=std(mh_p)/sqrt(counter); 
        mh_p_all(bp,3)=var(mh_p); 

        mh_n_all(bp,1)= mean(mh_n);
        mh_n_all(bp,2)=std(mh_n)/sqrt(counter);
        mh_n_all(bp,3)= var(mh_n);

        var_bp(bp,1)=var(vertcat(mh_p,mh_n));



    end
    pholdall.(outputnames{t})=mh_p_all;
    nholdall.(outputnames{t})= mh_n_all;
    varholdall.(outputnames{t})= var_bp;

end
pholdall_vel=pholdall;
nholdall_vel=nholdall;
varholdall_vel=varholdall;

TRF_study=TRF_study_accel;
clear meanholdall stdholdall pholdall nholdall varholdall_accel
for t=1:length(outputnames)
    clear mh_p_all mh_n_all
    for bp=1:20
        clear mh_p mh_n 
        counter=0;
        for f=1:length(filenames)
            fns=fieldnames(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).test);
            for fn=1:length(fns)
                counter=counter+1;
                if counter==1
                    mh_p=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)>0,bp);
                    mh_n=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)<0,bp);
                else
                    mh_p=vertcat(mh_p,TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)>0,bp));
                    mh_n=vertcat(mh_n,TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).rtest.(fns{fn})(:,bp)<0,bp));
                end
            end
        end
        mh_p_all(bp,1)= mean(mh_p);
        mh_p_all(bp,2)=std(mh_p)/sqrt(counter); 
        mh_p_all(bp,3)=var(mh_p); 

        mh_n_all(bp,1)= mean(mh_n);
        mh_n_all(bp,2)=std(mh_n)/sqrt(counter);
        mh_n_all(bp,3)= var(mh_n);
                var_bp(bp,1)=var(vertcat(mh_p,mh_n));

    end
    pholdall.(outputnames{t})=mh_p_all;
    nholdall.(outputnames{t})= mh_n_all;
        varholdall.(outputnames{t})= var_bp;

end
pholdall_accel=pholdall;
nholdall_accel=nholdall;
varholdall_accel=varholdall;


outputlabels={'Horizontal (X)','Depth (Z)','Vertical (Y)'}; 
close all
figure;
for t=1:length(outputnames)
    subplot(2,3,t)
    b= bar(1:20,varholdall_vel.(outputnames{t}),'FaceColor', [0 0 0]);
    hold on
    ylim([0 6])
    set(gca,'xtick',1:20,'xticklabel',{1:20}),  grid on
    title(outputlabels{t})
    set(gca, 'FontSize', 8);
    if t==1
        ylabel('Movement \sigma^2 (Velocity)', 'Interpreter', 'tex')
    end
    hold off
end
for t=1:length(outputnames)
    subplot(2,3,t+3)
    b= bar(1:20,varholdall_accel.(outputnames{t}),'FaceColor', [0 0 0]);
    hold on
   ylim([0 0.75])
    set(gca,'xtick',1:20,'xticklabel',{[1:20]'}),  grid on
     if t==1
        ylabel('Movement \sigma^2 (Acceleration)', 'Interpreter', 'tex')
    end
    set(gca, 'FontSize', 8);
    hold off
end
sgtitle('Kinematic Data Summary', 'Interpreter', 'tex')
set(gcf, 'Position', [126,300,1266,404]);




%% Setup TRF loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Use to recreate Study 1 analysis %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name TRF outputs
readme.saveTRFname='Study1';
readme.Afolder=[homefolder 'audiofiles/'];
readme.Mfolder=[homefolder 'skelfiles/'];
readme.outputfolder=[homefolder 'TRFoutput/'];
readme.setupfolder=[homefolder 'setup_files/'];

% Model hyperparameters
readme.Dir = -1; % direction of causality
readme.tmin = -1000; % minimum time lag (ms)
readme.tmax = 200; % maximum time lag (ms)
readme.lambda = 10.^(-6:2:6);
readme.fs=30;
readme.zeropad=0;
% readme.audiorep='env1'; % F0 env1

% Analysis options
readme.foldlength=15; % in seconds
% readme.bmotavatar='mc_avatar'; % 'db_avatar' or 'mc_avatar'
%  readme.bmotagg='separate' ;  % 'sum' or 'separate'
% readme.bmotdirtype='split'; % 'concat' or 'split'
% readme.bmottype='all'; % 'kinetics' or 'kinematics' or'all'
% readme.speechtype='env'; % env# or dv
% readme.axistype='concat'; % split or concat
readme.bmotlist={'dempXvel1','dempYvel2','dempZvel3',...
    'dempXspeed1','dempYspeed2','dempZspeed3'};
readme.flippred='no';

% Run TRF loop
TRF_master=gocsg_TRFloop(readme);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Use to recreate Study 2 analysis %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name TRF outputs
readme.saveTRFname='Study2';
readme.Afolder=[homefolder 'audiofiles/'];
readme.Mfolder=[homefolder 'skelfiles/'];
readme.outputfolder=[homefolder 'TRFoutput/'];
readme.setupfolder=[homefolder 'setup_files/'];

% Model hyperparameters
readme.Dir = -1; % direction of causality
readme.tmin = -1000; % minimum time lag (ms)
readme.tmax = 200; % maximum time lag (ms)
readme.lambda = 10.^(-6:2:6);
readme.fs=30;
readme.zeropad=0;
% readme.audiorep='env1'; % F0 env1

% Analysis options
readme.foldlength=15; % in seconds
% readme.bmotavatar='mc_avatar'; % 'db_avatar' or 'mc_avatar'
% readme.bmotagg='separate' ;  % 'sum' or 'separate'
% readme.bmotdirtype='split'; % 'concat' or 'split'
% readme.bmottype='all'; % 'kinetics' or 'kinematics' or'all'
% readme.scaleup='no'; % 'DBI' 'DBII' 'all'
% readme.speechtype='env'; % env# or dv
% readme.axistype='concat'; % split or concat
readme.bmotlist={'dempXYZspeed','dempXYZvel','dempXYZaccel'};
readme.flippred='yes';

% Run TRF loop
TRF_master=gocsg_TRFloop(readme);



%% Study 1 analysis


fileidx=dir([homefolder 'TRFoutput/Study1/']);
fileidx = fileidx(~[fileidx.isdir]);
clearvars TRF_study
for f=1:size(fileidx,1)
    TRF_study.(strrep(fileidx(f).name,'.mat',''))=load([fileidx(f).folder '/' fileidx(f).name]);
end
filenames=fieldnames(TRF_study);
outputnames=fieldnames(TRF_study.(filenames{1}).TRF_master);
outputnames(end)=[];

clearvars tbl
current_row=0;
for f=1:length(filenames)
    for t=1:length(outputnames)
        fns=fieldnames(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).test);
        for fn=1:length(fns)
            current_row=current_row+1;
            tbl{current_row,1}=filenames{f};
            tbl{current_row,2}=outputnames{t};
            tbl{current_row,3}=fns{fn};
            tbl{current_row,4}=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).test.(fns{fn}).r;
            [~,tbl{current_row,5}]=max(squeeze(sum(abs(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).forwardmodel.(fns{fn}).w),2)));
            %   squeeze(sum(abs(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).forwardmodel.(fns{fn}).w),2));
            %   TRF_study.(filenames{f}).TRF_master.readme
            if contains(outputnames{t},'vel')
                tbl{current_row,6}='vector';
            end
            if contains(outputnames{t},'speed')
                tbl{current_row,6}='scalar';
            end
            if contains(outputnames{t},'X')
                tbl{current_row,7}='X';
            end
            if contains(outputnames{t},'Y')
                tbl{current_row,7}='Y';
            end
            if contains(outputnames{t},'Z')
                tbl{current_row,7}='Z';
            end

        end
    end
end
tbl=cell2table(tbl,'VariableNames',{'File','Bmottype','Fold','test_r','maxbp','Encoding','Axis'});

clearvars lme lme_output modelnames modelequations results lmeresults
modelnames={'Intercept','Additive','Interaction'};
modelequations={'test_r ~ 1  + (1|File)',...
    'test_r ~ 1  + Axis + Encoding + (1|File)',...
    'test_r ~ 1  + Axis*Encoding + (1|File)'}; 

% Loop to run lmers
for m=1:length(modelnames)
    lme.(modelnames{m}) = fitlme(tbl,modelequations{m}, 'FitMethod','REML'); % or +(1|itemID_A)
    lme_output.(modelnames{m})=dataset2cell(lme.(modelnames{m}).Coefficients);
end
cmpnames={'m1m2','m2m3'};
for ii=1:length(cmpnames)
    lmeresults.(cmpnames{ii})=dataset2cell(compare(lme.(modelnames{ii}),lme.(modelnames{ii+1})));
    lmeresults.(cmpnames{ii}){4,3}=lmeresults.(cmpnames{ii}){3,3}-lmeresults.(cmpnames{ii}){2,3};
    lmeresults.(cmpnames{ii}){4,4}=lmeresults.(cmpnames{ii}){3,4}-lmeresults.(cmpnames{ii}){2,4};
    lmeresults.(cmpnames{ii}){4,5}=lmeresults.(cmpnames{ii}){3,5}-lmeresults.(cmpnames{ii}){2,5};
end
lmeresults.all=vertcat(lmeresults.(cmpnames{1}),lmeresults.(cmpnames{2}));


disp(lme_output.Interaction);

clearvars allAIC
allAIC(1,1)=lme.Intercept.ModelCriterion{1,1};
allAIC(1,2)=lme.Additive.ModelCriterion{1,1};
allAIC(1,3)=lme.Interaction.ModelCriterion{1,1};
allAIC(2,:)=allAIC(1,:)-allAIC(1,1);
disp(allAIC);
clearvars allLL
allLL(1,1)=lme.Intercept.ModelCriterion{1,3};
allLL(1,2)=lme.Additive.ModelCriterion{1,3};
allLL(1,3)=lme.Interaction.ModelCriterion{1,3};
allLL(2,:)=allLL(1,:)-allLL(1,1);
disp(allLL);




%% Study 2 analysis

fileidx=dir([homefolder 'TRFoutput/Study2/']);
fileidx = fileidx(~[fileidx.isdir]);
clearvars TRF_study
for f=1:size(fileidx,1)
    TRF_study.(strrep(fileidx(f).name,'.mat',''))=load([fileidx(f).folder '/' fileidx(f).name]);
end
filenames=fieldnames(TRF_study);
outputnames=fieldnames(TRF_study.(filenames{1}).TRF_master);
disp(outputnames);
outputnames(end)=[];

outputnames={'dempXYZspeed','dempXYZvel','dempXYZaccel'};
predtypes={'test','test_xflip','test_yflip'};
clearvars tbl
current_row=0;
for f=1:length(filenames)
    for t=1:length(outputnames)
        fns=fieldnames(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).test);
        for fn=1:length(fns)
            for p=1:length(predtypes)
                current_row=current_row+1;
                tbl{current_row,1}=filenames{f};
                tbl{current_row,2}=outputnames{t};
                tbl{current_row,3}=fns{fn};
                tbl{current_row,4}=TRF_study.(filenames{f}).TRF_master.(outputnames{t}).(predtypes{p}).(fns{fn}).r;
                [~,tbl{current_row,5}]=max(squeeze(sum(abs(TRF_study.(filenames{f}).TRF_master.(outputnames{t}).forwardmodel.(fns{fn}).w),2)));
                if contains(outputnames{t},'vel')
                    tbl{current_row,6}='velocity';
                end
                if contains(outputnames{t},'speed')
                    tbl{current_row,6}='speed';
                end
                if contains(outputnames{t},'accel')
                    tbl{current_row,6}='accel';
                end
                tbl{current_row,7}=predtypes{p};

            end
        end
    end
end
tbl=cell2table(tbl,'VariableNames',{'File','Bmottype','Fold','test_r','numbodyparts','Encoding','Orientation'});
tbl_accel=tbl(ismember(tbl.Encoding(:),'accel'),:);


clearvars lme lme_output modelnames modelequations results lmeresults
modelnames={'Intercept','Additive'};
modelequations={'test_r ~ 1  + (1|File)',...
    'test_r ~ 1  + Orientation + (1|File)'};% Loop to run lmers
for m=1:length(modelnames)
    lme.(modelnames{m}) = fitlme(tbl_accel,modelequations{m}, 'FitMethod','REML'); % or +(1|itemID_A)
    lme_output.(modelnames{m})=dataset2cell(lme.(modelnames{m}).Coefficients);
end
clearvars allAIC
allAIC(1,1)=lme.Intercept.ModelCriterion{1,1};
allAIC(1,2)=lme.Additive.ModelCriterion{1,1};
allAIC(2,:)=allAIC(1,:)-allAIC(1,1);
disp(allAIC);
clearvars allLL
allLL(1,1)=lme.Intercept.ModelCriterion{1,3};
allLL(1,2)=lme.Additive.ModelCriterion{1,3};
allLL(2,:)=allLL(1,:)-allLL(1,1);
disp(allLL);


clearvars lme lme_output modelnames modelequations results lmeresults
lme1 = fitlme(tbl_accel,'test_r ~ 1   + (1|File)', 'FitMethod','REML'); % or +(1|itemID_A)
lme2 = fitlme(tbl_accel,'test_r ~ 1  + Orientation + (1|File)', 'FitMethod','REML'); % or +(1|itemID_A)
lme1_output=dataset2cell(lme1.Coefficients);
lme2_output=dataset2cell(lme2.Coefficients);
lmeresults=compare(lme1,lme2);




tbl_accel_summary=groupsummary(tbl_accel,'Orientation','mean','test_r');
tbl_accel_error=groupsummary(tbl_accel,'Orientation','std','test_r');
tbl_accel_error.se=tbl_accel_error.std_test_r./sqrt(tbl_accel_error.GroupCount(1));
close all
figure; clc
    h=bar(1:3,tbl_accel_summary.mean_test_r);
    h.FaceColor=[0.5 0.6 0.6];
    hold on
    er = errorbar(1:3,tbl_accel_summary.mean_test_r,tbl_accel_error.se,tbl_accel_error.se);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth=1.5;
set(gca,'xtick',1:3,'xticklabel',{'Preserved','Horizontal Inversion','Vertical Inversion'},...
    'ytick',0:0.05:0.2), grid on, axis square
ylim([-0.02 0.2])
ylabel('Pearson Correlation')
xlabel('Body Oreintation')
title('Model Performance Means')

