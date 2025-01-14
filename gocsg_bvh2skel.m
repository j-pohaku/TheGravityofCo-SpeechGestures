

function skelfiles=gocsg_bvh2skel(readme)




% Load bvh and transform data to mc format

%parentfolder2= "D:\Fdrive_852022\JPM_dissertation\Trinitydb\TSGDB_II\allBVH";

 stimfolder='/home/AD/jmomsen/Documents/Diss_study23/stims/';
% DB1 training data
dblist=dir([ stimfolder 'Trinitydb/bvh/']);
for i=1:length(dblist)
kill(i)=dblist(i).isdir;
end
dblist(kill)=[];
clearvars kill i



parentfolder= ([ stimfolder 'Trinitydb/bvh/']);
savefolder= '/home/AD/jmomsen/Documents/Diss_study23/stims/TCDfiles/skel_files/';
savefolder2= '/home/AD/jmomsen/Documents/Diss_study23/stims/TCDfiles/d_files/';

%cd "E:\Fdrive_852022\JPM_dissertation\Trinitydb\TSGDB_I\AllBVHFromClap\";
%         opts = delimitedTextImportOptions("NumVariables", 1);
%         % Specify range and delimiter
%         opts.DataLines =[2, Inf] ;
%         opts.Delimiter = ",";
%         % Specify column names and types
%         opts.VariableNames = "HIERARCHY";
%         opts.VariableTypes = "string";
%         % Specify file level properties
%         opts.ExtraColumnsRule = "ignore";
%         opts.EmptyLineRule = "read";
%         opts.ConsecutiveDelimitersRule = "join";
%         % Specify variable properties
%         opts = setvaropts(opts, "HIERARCHY", "WhitespaceRule", "preserve");
%         opts = setvaropts(opts, "HIERARCHY", "EmptyFieldRule", "auto");

for s=1:length(dblist)
        % tic
                      file = dblist(s).name;
                         fn=[convertStringsToChars(parentfolder) file];

                    %[d japar] = mcread(fn);
                    
%                     uiimport(fn)
%                     load(fn)


%             % Import the data
%             myfile = table2cell(readtable(fn, opts));
% 
%             for i=1:size(myfile,1)
%                 if strcmp(myfile{i,1},"MOTION")==1;
%                motionrow=i;
%                 else
%                 end
%             end      
%             
%            framenum_bvh(s,1)=myfile{motionrow,1};
%            framenum_bvh(s,2)=myfile{motionrow+1,1};
%            framenum_bvh(s,3)=myfile{motionrow+2,1};
%            clearvars myfile motionrow i
% end
                    % read the data 
                    [skel, channels, frameLength] = bvhReadFile(fn);
                    bvhload=loadbvh(fn);
                    
                    % test second matrix computation 
                    %[skeleton,time] = loadbvh(fn);
                    frameLength=1/59.94;
                    fps=59.94;
                    %save([ 'skeleton' file '.mat' ],'skeleton', 'time');
                    save([savefolder 'skel_' file '.mat' ],'skel', 'channels', 'frameLength','fps');
                    save([savefolder 'bvhload_' file '.mat' ],'bvhload');
                    clearvars bvhload

       %end

%  skel=skel_jpm;
%  channels=channels_jpm;
%  frameLength=frameLength_jpm;
% load skel_jpmRecording_001.bvh
% load mcdemodata.mat


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

                               % OPTION 1 - mc_data maker
                               % Creates data structure exactly 3 times as
                               % big as channels per frame (end sites
                               % included)
                                        for i = 1:size(channels,1) % for each frame
                                            vals=skel2xyz(skel, channels(i,:));
                                            tmp=vals(1,:);
                                            for j = 2:size(skel.tree,2)    
                                                tmp=[tmp vals(j,:)];
                                            end
                                            d.data(i,:)=tmp;
                                        end
                               % OPTION 2 - copy "channels" from bvhReadFile
                                        % d.data=channels;
%{
                               % OPTION 3
                               % Fix 9/26 : try looping through channels
                               % and deleting end site markers
%                                skel_trim=skel;
%                                skelindex=zeros(1,1);
%                                %skel_trim = rmfield( skel_trim.tree , 'Hips' );
%                                for ch=1:length(skel.tree)
%                                             if strcmp(skel.tree(ch).name,'Site')==1
%                                                 skelindex(1,length(skelindex)+1)=ch;
%                                             else
%                                             end
%                                end 
%                                skelindex = nonzeros(skelindex');
%                                skel_trim.tree(skelindex)=[]; %skel_trim is skel without end site markers
% 
%                                for i = 1:size(channels,1) % for each frame
%                                             vals=skel2xyz(skel_trim, channels(i,:));
%                                             tmp=vals(1,:);
%                                             for j = 2:size(skel_trim.tree,2)    
%                                                 tmp=[tmp vals(j,:)];
%                                             end
%                                             d.data(i,:)=tmp;
%                                end

%}
                                        d.other = [];
                                        d.framestep=fps;
                                        d.filelength_min =  (d.nFrames/d.framestep)/60;
                                        d.filelength_sec =  (d.nFrames/d.framestep);
                                        d.freq = 1/d.freq; %or 59.94
                                        
                                        disp(strcat(fn,' loaded'))

                                save([savefolder2 'mcdata_' file '.mat' ],'d');
                    clear skel channels frameLength file fn d fps 
                    toc
end
%% save time by converting to mocap data

% DB1 training data
% DB1 training bvh files
dblistV=dir("D:\Fdrive_852022\JPM_dissertation\TCD_files\DB1_training\");
dblistV(1)=[];
dblistV(1)=[];
dblistV(1)=[];
dblistV(1)=[];
dblistV(1)=[];
parentfolderV= 'D:\Fdrive_852022\JPM_dissertation\TCD_files\DB1_training\';
savefolder= 'D:\Fdrive_852022\JPM_dissertation\TCD_files\DB1_training\Mocapversions\';



for s=1:length(dblistV)

          
        fileV = dblistV(s).name;
        fnV=[convertStringsToChars(parentfolderV) fileV];

       
         % preproc video
         load([fnV]); % load channels, skel, and frameLength
         
         
          % create MoCap structure and convert 
          % BVH format to 3d position coordinates
    
    d.type = 'MoCap data';        
    d.filename = skel.name;
    d.nFrames = size(channels,1);
    d.nCameras = 0;
    d.nMarkers = size(skel.tree,2);
    d.framestep = frameLength;
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
    d.filelength_min =  (d.framestep*d.nFrames)/60;
    d.filelength_sec =  (d.framestep*d.nFrames);
    d.freq = 1/d.framestep; 
    disp(strcat(fileV,' loaded'))
    clear skel channels frameLength

save([savefolder 'mocapdta_' fileV '.mat' ],'d');

end