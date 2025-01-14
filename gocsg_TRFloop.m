


function TRF_master=gocsg_TRFloop(readme)


skeld_dbI_folder=readme.Mfolder;
fileidxM=dir(skeld_dbI_folder);
fileidxM = fileidxM(~[fileidxM.isdir]);
fileidxM=struct2cell(fileidxM)';
fileidxA=readme.Afolder;
fileidxA=dir(fileidxA);
fileidxA = fileidxA(~[fileidxA.isdir]);
fileidxA=struct2cell(fileidxA)';
load([readme.setupfolder 'dbI_info.mat']);
load([readme.setupfolder 'jointsegmentkeys.mat']);



for s=1:23

    if s==1
        load([fileidxA{s,2} '/' fileidxA{s,1}]);
    else
        audiodta=load([fileidxA{s,2} '/' fileidxA{s,1}]);
    end
    mocapdta=load([fileidxM{s,2} '/' fileidxM{s,1}]);

    mocapdta.d.frameLength=mocapdta.d.freq;
    mocapdta.d.freq=1/mocapdta.d.freq;
    readme.oldfs=round(mocapdta.d.freq);

    readme.rawmocaplength=mocapdta.d.nFrames/mocapdta.d.freq;
    readme.rawaudiolength=length(audiodta.speech.wav_untrimmed)/audiodta.speech.wavSR;
    readme.trimaudiostart(1,1)=tbl_dbI{s,4};
    readme.trimaudiostart(1,2)=tbl_dbI{s,4}*audiodta.speech.wavSR;
    audiodta.speech.wavtrim=audiodta.speech.wav_untrimmed(readme.trimaudiostart(1,2):end);
    readme.trimaudiolength=length(audiodta.speech.wavtrim)/audiodta.speech.wavSR;
    audiodta.speech.wavtrim_env1=env1(audiodta.speech.wavtrim,audiodta.speech.wavSR,120,10);
    readme.newenv1length=length(audiodta.speech.wavtrim_env1)/120;
    readme.oldenv1length=length(audiodta.speech.wav_env1)/120;


    % trim ends to match sample lengths
    if length(audiodta.speech.wav_env1) > mocapdta.d.nFrames % if audio is longer
        readme.cuttomatchtime='speech';
        for i=1:abs(length(audiodta.speech.wav_env1)-mocapdta.d.nFrames)
            audiodta.speech.wav_env1(length(audiodta.speech.wav_env1),:)=[];
        end
    end
    if length(audiodta.speech.wav_env1) < mocapdta.d.nFrames % if mocap is longer
        readme.cuttomatchtime='mocap';
        mocapdta.d=mctrim(mocapdta.d,1,(mocapdta.d.nFrames-abs(length(audiodta.speech.wav_env1)-mocapdta.d.nFrames)),'frame');
    end

    % drop missing markers
    readme.markerdatasizecheck=length(mocapdta.d.markerName)*3==size(mocapdta.d.data,2);
    mocapdta.d.dim1horz =  mocapdta.d.data(:,1:3:end); % x
    mocapdta.d.dim2depth =  mocapdta.d.data(:,2:3:end); % y
    mocapdta.d.dim3vert =  mocapdta.d.data(:,3:3:end); % z
    findnan1 =find(isnan(mocapdta.d.dim1horz(1,:)));
    findnan2 = find(isnan(mocapdta.d.dim2depth(1,:)));
    findnan3 =find(isnan(mocapdta.d.dim3vert(1,:)));
    findnan =unique(horzcat(findnan1,findnan2,findnan3));
    mocapdta.d = mcrmmarker(mocapdta.d, findnan);
    readme.removemarkers=findnan;
    mocapdta.d.dim1horz =  mocapdta.d.data(:,1:3:end); % x
    mocapdta.d.dim2depth =  mocapdta.d.data(:,2:3:end); % y
    mocapdta.d.dim3vert =  mocapdta.d.data(:,3:3:end); % z
    clearvars findnan findnan1 findnan2 findnan3

    % identify if data needs to be flipped & save trace of flip decision
    xposlims=zeros(mocapdta.d.nFrames,2);
    yposlims=zeros(mocapdta.d.nFrames,2);
    zposlims=zeros(mocapdta.d.nFrames,2);
    for i=1:mocapdta.d.nFrames
        [xposlims(i,1),xposlims(i,2)]  =  bounds(mocapdta.d.dim1horz(i,:)); % x
        [yposlims(i,1),yposlims(i,2)]  =  bounds(mocapdta.d.dim2depth(i,:)); % y
        [zposlims(i,1),zposlims(i,2)]   =  bounds(mocapdta.d.dim3vert(i,:)); % z
    end

    readme.meanXdim1horzrange=mean(xposlims(:,2))-mean(xposlims(:,1));
    readme.meanYdim2depthrange=mean(yposlims(:,2))-mean(yposlims(:,1));
    readme.meanZdim3vertrange=mean(zposlims(:,2))-mean(zposlims(:,1));

    if readme.meanZdim3vertrange<readme.meanYdim2depthrange
        mocapdta.d = mcreorderdims(mocapdta.d, [1 3 2]);% flipping z and y axes - do for all of dbII except 3 outliers
        readme.flipdims_123='yes';
    else
        readme.flipdims_123='no';
    end

    mocapdta.d.dim1horz =  mocapdta.d.data(:,1:3:end); % x
    mocapdta.d.dim2depth =  mocapdta.d.data(:,2:3:end); % y
    mocapdta.d.dim3vert =  mocapdta.d.data(:,3:3:end); % z

    clearvars d2j
    d2j = mcm2j(mocapdta.d, m2jpar_dbI);
    d2s = mcj2s(d2j, j2spar_dbI);

    spar = mcgetsegmpar('Dempster',[0 0 8 7 6 0 8 7 6 13 12 10 11 3 2 1 11 3 2 1]);
    mocapdta.d.demppot = mcpotenergy(d2j, d2s, spar);
    [mocapdta.d.demptrans, mocapdta.d.demprot] = mckinenergy(d2j, d2s, spar);
    clearvars d2s spar

    mocapdta.d.dempkinetic=mocapdta.d.demprot+mocapdta.d.demptrans;
    mocapdta.d.demplagrangian=mocapdta.d.demppot-mocapdta.d.dempkinetic;

    % get kinetic derivatives
    % currently calculates gradient, then std both d and dv vectors separately
    for i=1:size(mocapdta.d.demppot,2)
        mocapdta.d.demppotdv(:,i)=gradient(mocapdta.d.demppot(:,i));
    end
    for i=1:size(mocapdta.d.demprot,2)
        mocapdta.d.dempkineticdv(:,i)=gradient(mocapdta.d.dempkinetic(:,i));
    end
    for i=1:size(mocapdta.d.demptrans,2)
        mocapdta.d.dempactiondv(:,i)=gradient(mocapdta.d.demplagrangian(:,i));
    end

    % document kinetics with missing bodypart values
    killpot=find(sum(mocapdta.d.demppot,1)==0);
    readme.killpot=killpot;
    readme.droppot=m2jpar_dbI.markerName(find(sum(mocapdta.d.demppot ,1)==0));
    readme.droprot=m2jpar_dbI.markerName(find(sum(mocapdta.d.demprot ,1)==0));
    readme.droptrans=m2jpar_dbI.markerName(find(sum(mocapdta.d.demptrans ,1)==0));

    % position data across dimensions
    mocapdta.d.demmpmarkerName=d2j.markerName;
    mocapdta.d.dempXpos1=d2j.data(:,1:3:end);
    mocapdta.d.dempYpos2=d2j.data(:,2:3:end);
    mocapdta.d.dempZpos3=d2j.data(:,3:3:end);
    % velocity/acceleration data
    mocapdta.d.d2jvel = mctimeder(d2j, 1,'acc'); % velocity
    mocapdta.d.d2jaccel = mctimeder(d2j, 2,'acc'); % acceleration
    mocapdta.d.d2jvel= mocapdta.d.d2jvel.data;
    mocapdta.d.d2jaccel= mocapdta.d.d2jaccel.data;

    % velocity/acceleration data across dimensions
    mocapdta.d.dempXvel1=mocapdta.d.d2jvel(:,1:3:end);
    mocapdta.d.dempYvel2=mocapdta.d.d2jvel(:,2:3:end);
    mocapdta.d.dempZvel3=mocapdta.d.d2jvel(:,3:3:end);
    mocapdta.d.dempXaccel1=mocapdta.d.d2jaccel(:,1:3:end);
    mocapdta.d.dempYaccel2=mocapdta.d.d2jaccel(:,2:3:end);
    mocapdta.d.dempZaccel3=mocapdta.d.d2jaccel(:,3:3:end);
    mocapdta.d.dempXspeed1=abs(mocapdta.d.dempXvel1);
    mocapdta.d.dempYspeed2=abs( mocapdta.d.dempYvel2);
    mocapdta.d.dempZspeed3=abs(  mocapdta.d.dempZvel3);
    for i=1:size(mocapdta.d.demppot,2)
        mocapdta.d.dempXspeed1dv(:,i)=gradient(mocapdta.d.dempXspeed1(:,i));
        mocapdta.d.dempYspeed2dv(:,i)=gradient(mocapdta.d.dempYspeed2(:,i));
        mocapdta.d.dempZspeed3dv(:,i)=gradient(mocapdta.d.dempZspeed3(:,i));
    end
    % match response columns across kinematic and kinetic variables

    mocapdta.d= rmfield(mocapdta.d,'demptrans');
    mocapdta.d= rmfield(mocapdta.d,'demprot');

    respnames=fieldnames(mocapdta.d);
    respnames(~contains(fieldnames(mocapdta.d),'demp'))=[];
    % normalize each biomotion parameter
    for i=1:length(respnames)
        mocapdta.d.(respnames{i})=mocapdta.d.(respnames{i})/std(mocapdta.d.(respnames{i})(:));
    end

    [readme.rho_posXpot,readme.pval_posXpot]=corr(mean(mocapdta.d.demppot,2),mean(mocapdta.d.dempXpos1,2));
    [readme.rho_posYpot,readme.pval_posYpot]=corr(mean(mocapdta.d.demppot,2),mean(mocapdta.d.dempYpos2,2));
    [readme.rho_posZpot,readme.pval_posZpot]=corr(mean(mocapdta.d.demppot,2),mean(mocapdta.d.dempZpos3,2));

    mocapdta.d.dempXYZvel=horzcat(mocapdta.d.dempXvel1,mocapdta.d.dempYvel2,mocapdta.d.dempZvel3);
    mocapdta.d.dempXYZspeed=horzcat(mocapdta.d.dempXspeed1,mocapdta.d.dempYspeed2,mocapdta.d.dempZspeed3);
    mocapdta.d.dempXYZaccel=horzcat(mocapdta.d.dempXaccel1,mocapdta.d.dempYaccel2,mocapdta.d.dempZaccel3);

    respnames=readme.bmotlist;

    % resample
    if readme.fs~=readme.oldfs
        for i=1:length(respnames)
            mocapdta.d.(respnames{i}) = resample(mocapdta.d.(respnames{i}), readme.fs, readme.oldfs);
        end
        audiodta.speech.wav_env1 = resample( audiodta.speech.wav_env1, readme.fs, readme.oldfs);
    end

    for i=1:length(respnames)
        respall.(respnames{i})=mocapdta.d.(respnames{i});
    end
    stimall.wav_env1=audiodta.speech.wav_env1;

    nfold=round((length(audiodta.speech.wav_env1)/readme.fs)/readme.foldlength);
    fs=readme.fs;
    Dir=readme.Dir;
    tmax=readme.tmax;
    tmin=readme.tmin;
    lambda=readme.lambda;
    zeropad=readme.zeropad;

    % clear any bmotion parameters that didn't compute (nans)
    clearvars killmiss
    for resptype=1:length(respnames)
        resp=respall.(respnames{resptype});
        killmiss(resptype)=  sum(find(isnan(resp)));
    end
    readme.containednansbmot=respnames(killmiss>0);
    readme.containednansrespnum=find(killmiss>0);
    for i=1:length(killmiss)
        if killmiss(i)>0
            respall=rmfield(respall,(respnames{i}));
        else
        end
    end

    respnames=fieldnames(respall);
    stimnames=fieldnames(stimall);

    for resptype=1:length(respnames)
        resp=respall.(respnames{resptype});
        stim=stimall.(stimnames{1});
        if size(resp,1)~=size(stim,1)
            stim=stim';
        end
        
        for testTrial=1:nfold
            [strain,rtrain,stest,rtest] = mTRFpartition(stim,resp,nfold,testTrial);
            cv.(['fold_' num2str(testTrial)]) = mTRFcrossval(strain,rtrain,fs,Dir,tmin,tmax,lambda,'zeropad',zeropad,'fast',1);
            [rmax.(['fold_' num2str(testTrial)]),idx.(['fold_' num2str(testTrial)])] = max(mean(cv.(['fold_' num2str(testTrial)]).r,1));
            model.(['fold_' num2str(testTrial)]) = mTRFtrain(strain,rtrain,fs,Dir,tmin,tmax,lambda(idx.(['fold_' num2str(testTrial)])),'zeropad',zeropad);
            [pred.(['fold_' num2str(testTrial)]),test.(['fold_' num2str(testTrial)])] = mTRFpredict(stest,rtest,model.(['fold_' num2str(testTrial)]),'zeropad',zeropad);
            forwardmodel.(['fold_' num2str(testTrial)]) = mTRFtransform(model.(['fold_' num2str(testTrial)]),rtrain);
            stest_all.(['fold_' num2str(testTrial)]) =stest;
            rtest_all.(['fold_' num2str(testTrial)]) =rtest;
            if isequal(readme.flippred,'yes')
                if size(rtest,2)>20
                    rtest_yflip= rtest;
                    rtest_xflip=rtest;
                    rtest_yflip(:,41:60)=rtest_yflip(:,41:60)*-1;
                    rtest_xflip(:,1:20)=rtest_xflip(:,1:20)*-1;
                    [pred_yflip.(['fold_' num2str(testTrial)]),test_yflip.(['fold_' num2str(testTrial)])] = mTRFpredict(stest,rtest_yflip,model.(['fold_' num2str(testTrial)]),'zeropad',zeropad);
                    [pred_xflip.(['fold_' num2str(testTrial)]),test_xflip.(['fold_' num2str(testTrial)])] = mTRFpredict(stest,rtest_xflip,model.(['fold_' num2str(testTrial)]),'zeropad',zeropad);
                end
            end


            clearvars rtest rtrain stest strain rtest_yflip rtest_xflip



        end %end testTrial loop

        TRF_master.(respnames{resptype}).stim=(stimnames{1});
        TRF_master.(respnames{resptype}).model=model;
        TRF_master.(respnames{resptype}).rmax=rmax;
        TRF_master.(respnames{resptype}).lambdaidx=idx;
        TRF_master.(respnames{resptype}).cv=cv;
        TRF_master.(respnames{resptype}).pred=pred;
        TRF_master.(respnames{resptype}).test=test;
        TRF_master.(respnames{resptype}).forwardmodel=forwardmodel;
        TRF_master.(respnames{resptype}).stest=stest_all;
        TRF_master.(respnames{resptype}).rtest=rtest_all;
        

        if exist('pred_xflip','var') == 1
            TRF_master.(respnames{resptype}).test_xflip=test_xflip;
            TRF_master.(respnames{resptype}).test_yflip=test_yflip;
            TRF_master.(respnames{resptype}).pred_xflip=pred_xflip;
            TRF_master.(respnames{resptype}).pred_yflip=pred_yflip;
        end

        clearvars model forwardmodel test pred cv rmax idx  stest_all rtest_all
        clearvars resp stim pred_xflip pred_yflip  test_yflip test_xflip
    end % end resp/stim type loop


    %  TRF_master.readme=readme;
    filename=[fileidxA{s,1}];
    filename=strrep(filename,'speech_','');
    filename=strrep(filename,'.mat','');
    filename=strrep(filename,'-','_');
    save([readme.outputfolder readme.saveTRFname 'rndperm_' filename ],'rndpermstats');
   % save([readme.outputfolder readme.saveTRFname 'rndperm_long_' filename ],'rndpermstats_long');
    % save([readme.outputfolder readme.saveTRFname filename ],'TRF_master');
    clearvars filename nfold

end % end file loop


