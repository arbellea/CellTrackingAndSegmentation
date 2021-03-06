function main(txtFileName)
rng(1)
dbclear if caught error
myCluster = parcluster('local');
myCluster.NumWorkers
%poolobj = gcp('nocreate');
%if isempty(poolobj)
%poolobj = parpool(min(1,myCluster.NumWorkers),'IdleTimeout', 240);
%end
if ~isunix
dbstop if caught error
end
addpath(genpath(fullfile('.','SourceCode')));
if isunix
    addpath(genpath(fullfile('.','MitodixGit')));
else
    addpath(genpath(fullfile('..','Mitodix')));
end
warning off;

SendEMail;
emailAddress = [];

try
Params = read_parameters_txt([],txtFileName,false);
if isfield(Params.General,'emailAddress')
    emailAddress = Params.General.emailAddress;
else
    emailAddress = [];
end
if isfield(Params.Flags,'profile')&&Params.Flags.profile
    profile on;
else
    Params.Flags.profile = false;
end
fprintf('Loading Images... Please Wait...\n');
if isfield(Params.load_data_params,'expr');
   expr  = Params.load_data_params.expr ;
else
    expr  = '\w+(\d+).*';
end
if isfield(Params.load_data_params,'tagged_expr');
    tagged_expr = Params.load_data_params.tagged_expr;
else 
    tagged_expr  = '\w+(\d+).*';
end
 
Data = Load_Data(Params.load_data_params.data_path,Params.load_data_params.extention,expr);
fprintf('Loading Manual Segmentaion Images... Please Wait...\n');
TaggedData = Load_Data(Params.load_data_params.tagged_data_path,Params.load_data_params.tagged_extention,tagged_expr);
if isfield(Params.parameters,'stop_frame')
    Params.stop_frame = min(Data.Frame_Num,Params.parameters.stop_frame);
else
    Params.stop_frame = Data.Frame_Num;
end
Params.start_frame = 1;
disp('Initialize!!')
%which Initialize_Tracker
Tracking = Initialize_Tracker(Data,TaggedData,Params);
Tracking.stop_frame = Params.stop_frame;
disp('Lets Go!!')
saveDir = Track_and_Segment(Data,Tracking,Params,TaggedData);
msg = sprintf('Done Analyzing %s\n Results saved at: %s',txtFileName,saveDir);

if ~isempty(emailAddress)
sendmail(emailAddress,'BGU Cluster Report #1', msg);
end

if Params.Flags.WriteVideo
    for i = 1:TaggedData.Frame_Num
    [~,fname,fext]=fileparts(Data.Frame_name{i});
    segfile = fullfile(saveDir,'Results',sprintf('Seg_%s%s',fname,fext));
    copyfile(TaggedData.Frame_name{i},segfile,'f');
    end
    
    CreateOutputVideo(Data,saveDir);
    if ~isempty(emailAddress)
    msg = sprintf('Done Creating Video for %s\n Video saved at: %s',txtFileName,saveDir);
    
    sendmail(emailAddress,'BGU Cluster Report #2', msg);
    end
end

catch err
errmsg = getReport(err);
if ~isempty(emailAddress)
msg = sprintf('Error analyzing  %s!!!\n %s',txtFileName,errmsg);
sendmail(emailAddress,'BGU Cluster Report #0', msg);
end
%delete(poolobj);
if Params.Flags.profile
if isunix
    profile off
    profsave
else
    profile off
    profile viewer
end
end
rethrow(err);

end

if Params.Flags.profile
if isunix
    profile off
    profsave
else
    profile off
    profile viewer
end
end



%delete(poolobj);
%profile off;profile viewer;

