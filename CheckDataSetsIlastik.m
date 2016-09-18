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
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June27_2_CheckPoints','Results'); Name = 's101 - Regularized';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_July9_CheckPoints','Results'); Name = 's101-Regularized';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_July10_Mitosis_CheckPoints','Results'); Name = 's101-Reg-Mitosis';

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s101_t(\d+).TIF';

resPath = '/Users/assafarbelle/Documents/PhD/Ilastik/Raw_Data_ALL_Object-Identities.h5'; Name = 's102-Ilastik';

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s102_t(\d+).TIF';

CheckJoseIlastik(txtfilePath,data_path,extention,expr,resPath,Name);



%%

HOME = getenv('HOME');

txtfilePath = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_IR00GY_tracking','M20150416_IR00GY_s102_tracking.txt') ;
data_path = fullfile(HOME,'Input','Galit_Lahav','Jose3','M20150416_p21_20X');
extention = '*.TIF';
expr = 'M20150416_w2CFP_s102_t(\d+).TIF';


%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_CheckPoints','Results'); Name = 's101-Original'%ORIGINAL
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June6_CheckPoints','Results'); %GMMG
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_May31_Mitosis_CheckPoints','Results'); %MITOSIS
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June10_CheckPoints','Results'); Name = 's101';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June14_Mitosis_CheckPoints','Results'); Name = 's101-Mitodix';
resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June14_CheckPoints','Results'); Name = 's101';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June17_CheckPoints','Results'); Name = 's101_Only Geodesic';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June18_CheckPoints','Results'); Name = 's101_TPT';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June25_CheckPoints','Results'); Name = 's101 - Regularized';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_June27_2_CheckPoints','Results'); Name = 's101 - Regularized';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_July9_CheckPoints','Results'); Name = 's101-Regularized';
%resPath = fullfile(HOME,'Outputs','Results_LahavLab-Jose3_20X_s101_July10_Mitosis_CheckPoints','Results'); Name = 's101-Reg-Mitosis';

resExt = '*.TIF';
resExp = 'Seg_M20150416_w2CFP_s101_t(\d+).TIF';

resPath = 'D:\DataFromJose\Raw_Data_s102_Object-Identities.h5'; Name = 's102-Ilastik';


CheckJoseIlastik(txtfilePath,data_path,extention,expr,resPath,Name);
