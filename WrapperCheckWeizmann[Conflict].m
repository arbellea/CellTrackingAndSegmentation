HOME = getenv('HOME');

load(fullfile(HOME,'ManualTrackingAndSeg','Weizmann','ManualSeg','SaveDataWithSeg.mat'))


data_path = fullfile(HOME,'..','Cell_Tracking','Input','Weizmann','Input','calibrate2');
extention = '*.TIF';
expr = 'calibrate2-P01.(\d+).TIF';


resPath = fullfile(HOME,'..','Cell_Tracking','Outputs','Results_AlonLab_CheckPoints','Results'); Name = 'Weizmann-Ours'; %ORIGINAL
resPath = fullfile(HOME,'..','Cell_Tracking','Outputs','Results_AlonLab_07-Feb-2016_11-28-01','Results'); Name = 'Weizmann-Ours'; %ORIGINAL
%resPath = fullfile(HOME,'..','Cell_Tracking','Outputs','Results_AlonLab_29-Feb-2016_10-21-26','Results'); Name = 'Weizmann-Ours'; %ORIGINAL

resExt = '*.TIF';

resExp = 'Seg_calibrate2-P01.(\d+).TIF';





CheckWeizman(cellData,Link,data_path,extention,expr,resPath,resExt,resExp,Name)
%%
resPath = '/Users/assafarbelle/IlastikProject/Weizmann/Raw Data_Object-Identities.h5'; Name = 'H1299-Ilastik';
csvPath = '/Users/assafarbelle/IlastikProject/Weizmann/Raw Data-exported_data_table.csv';
%resPath = '/Users/assafarbelle/Downloads/Raw Data_Object-Identities.h5'; Name = 'Weizmann-Ilastik';
CheckWeizmanIlastik(cellData,Link,data_path,extention,expr,resPath,csvPath,Name)

