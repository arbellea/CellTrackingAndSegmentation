function [OutStruct] = read_parameters_txt(pathname,fname,GUI,varargin)
if GUI
h = waitbar(0,'Please Wait...');
else
    fprintf('Reading parameters file... Please wait...\n');
end

if length(varargin)==1
    OutStruct = varargin{1};
elseif length(varargin)>1
    error('Too many inputs')
elseif isempty(varargin)
    OutStruct = struct;
end
if ~isempty(pathname)
    fname = fullfile(pathname,fname);
end
fid = fopen(fname);  

while ~feof(fid)
    tline = fgetl(fid);
    tline = strtrim(tline);
    if isempty(tline)|| tline(1) == '%'
        continue;
    end
    splt_line = strsplit(tline);
    if splt_line{1} == '*'
        fldname = replace_wspace(splt_line,2);
    else
        splt_line = strsplit(tline,'=');
        splt_line = strtrim(splt_line);
        sub_splt = strsplit(splt_line{1});
        subfldname = replace_wspace(sub_splt,1);
        val = eval(splt_line{2});
        OutStruct.(fldname).(subfldname) = val;
    end    
end
if GUI
close(h);
else
   fprintf('Done Reading parameters file\n');
end

        
        
        
        
