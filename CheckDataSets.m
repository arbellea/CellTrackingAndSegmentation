%%
HOME = getenv('HOME');

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s101_tracking.txt') ;
data_path = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s101_t(\d+).TIF';


%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_CheckPoints','Results'); Name = 's101-Original'%ORIGINAL
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June6_CheckPoints','Results'); %GMMG
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_May31_Mitosis_CheckPoints','Results'); %MITOSIS
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June10_CheckPoints','Results'); Name = 's101';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June14_Mitosis_CheckPoints','Results'); Name = 's101-Mitodix'; 
resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June14_CheckPoints','Results'); Name = 's101'; 
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June17_CheckPoints','Results'); Name = 's101_Only Geodesic'; 
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June18_CheckPoints','Results'); Name = 's101_TPT';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June25_CheckPoints','Results'); Name = 's101 - Regularized'; 
resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June27_CheckPoints','Results'); Name = 's101 - Regularized'; 

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s101_t(\d+).TIF';

CheckJose(txtfilePath,data_path,extention,expr,resPath,resExt,resExp,Name);
%% s102

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s102_tracking.txt') ;

%D:\DataFromJose\M20150416_p21_20X
data_path = fullfile('D:','DataFromJose','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s102_t(\d+).TIF';


%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s102_June10_CheckPoints','Results'); 
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s102_June14_Mitosis_CheckPoints','Results'); Name = 's102-Mitodix';
resPath = fullfile('C:\Users\arbellea\Google Drive BGU\Cell_Tracking','Outputs','Results_LahavLab-Jose3_20X_s102_June14_CheckPoints','Results'); Name = 's102';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s102_June17_CheckPoints','Results'); Name = 's102-Geo Only';

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s102_t(\d+).TIF';

CheckJose(txtfilePath,data_path,extention,expr,resPath,resExt,resExp,Name);
%% s103

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s103_tracking.txt') ;
data_path = fullfile('D:','DataFromJose','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s103_t(\d+).TIF';


%resPath = fullfile('.','Outputs','Results_LahavLab-Jose3_20X_s103_June14_Mitosis_CheckPoints','Results'); Name = 's103 - Mitodix';
resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s103_June14_CheckPoints','Results'); Name = 's103';
resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s103_t(\d+).TIF';

CheckJose(txtfilePath,data_path,extention,expr,resPath,resExt,resExp,Name);
%% s104

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s104_tracking.txt') ;
data_path = fullfile('D:','DataFromJose','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s104_t(\d+).TIF';



%resPath = fullfile('.','Outputs','Results_LahavLab-Jose3_20X_s104_June8_CheckPoints','Results'); 
resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s104_June14_CheckPoints','Results'); Name = 's104';

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s104_t(\d+).TIF';

CheckJose(txtfilePath,data_path,extention,expr,resPath,resExt,resExp,Name);
%% s105

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s105_tracking.txt') ;
data_path = fullfile('D:','DataFromJose','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s105_t(\d+).TIF';


resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s105_CheckPoints','Results'); 

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s105_t(\d+).TIF';

CheckJose(txtfilePath,data_path,extention,expr,resPath,resExt,resExp);
