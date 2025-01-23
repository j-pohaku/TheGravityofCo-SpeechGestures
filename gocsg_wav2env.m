



function envfiles=gocsg_wav2env(readme)


fileidxA=readme.Afolder;
fileidxA=dir(fileidxA);
fileidxA = fileidxA(~[fileidxA.isdir]);
fileidxA=struct2cell(fileidxA)';
load([readme.setupfolder 'dbI_info.mat']);
load([readme.setupfolder 'jointsegmentkeys.mat']);
savefolder= readme.outputfolder;


              
for s=1:size(fileidxA,1)
    [mywav,wavSR]=readAudio([fileidxA{s,2} '/' fileidxA{s,1}]);

    wav_env1=env1(mywav,wavSR,120,10);
    [audiodta.speech.F0,audiodta.speech.strength,audiodta.speech.T_indF0,audiodta.speech.wflag] = getF0(mywav,wavSR,(1/120));
   
    audiodta.speech.F0=audiodta.speech.F0';

    audiodta.speech.wav_env1=wav_env1';

    save([savefolder '/' fileidxA{s,1}],'audiodta');
    clearvars mywav wavSR wav_env1

end

envfiles='Done';