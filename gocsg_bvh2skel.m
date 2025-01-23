

function skelfiles=gocsg_bvh2skel(readme)




fileidxM=dir(readme.Mfolder);
fileidxM = fileidxM(~[fileidxM.isdir]);
load([readme.setupfolder 'dbI_info.mat']);
load([readme.setupfolder 'jointsegmentkeys.mat']);
savefolder= readme.outputfolder;

for s=1:size(fileidxM,1)
    file = fileidxM(s).name;
    fn=[convertStringsToChars(readme.Mfolder) file];
    % read the data
    [skel, channels, frameLength] = bvhReadFile(fn);

    frameLength=1/59.94;
    fps=59.94;


    % create MoCap structure
    d.type = 'MoCap data';
    d.filename = skel.name;
    d.nFrames = size(channels,1);
    d.nCameras = 0;
    d.nMarkers = size(skel.tree,2);
    d.freq = fps;
    d.nAnalog = [];
    d.anaFreq= [];
    d.timederOrder = 0;
    d.markerName = [];
    d.data = [];
    d.analogdata = [];
    d.other = [];
    for i=1:size(skel.tree,2)
        d.markerName{i,1} = skel.tree(i).name;
    end
    for i = 1:size(channels,1) % for each frame
        vals=skel2xyz(skel, channels(i,:));
        tmp=vals(1,:);
        for j = 2:size(skel.tree,2)
            tmp=[tmp vals(j,:)];
        end
        d.data(i,:)=tmp;
    end
    d.other = [];
    d.framestep=fps;
    d.filelength_min =  (d.nFrames/d.framestep)/60;
    d.filelength_sec =  (d.nFrames/d.framestep);
    d.freq = 1/d.freq; %or 59.94

    disp(strcat(fn,' loaded'))
    save([savefolder 'skel_' file '.mat' ],'d','skel', 'channels', 'frameLength','fps');
    clear skel channels frameLength file fn d fps
    toc
end


skelfiles='Done';