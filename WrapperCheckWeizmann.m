HOME = getenv('HOME');

load(fullfile(HOME,'ManualTrackingAndSeg','Weizmann','ManualSeg','SaveDataWithSeg.mat'))


data_path = fullfile(HOME,'..','Cell_Tracking','Input','Weizmann','Input','calibrate2');
extention = '*.TIF';
expr = 'calibrate2-P01.(\d+).TIF';


resPath = fullfile(HOME,'Outputs','Results_AlonLab_CheckPoints','Results'); Name = 'Weizmann-Ours'; %ORIGINAL

resExt = '*.TIF';

resExp = 'Seg_calibrate2-P01.(\d+).TIF';





CheckWeizman(cellData,Link,data_path,extention,expr,resPath,resExt,resExp,Name)

%resPath = '/Users/assafarbelle/Documents/PhD/Ilastik/Raw_Data_ALL_Object-Identities.h5'; Name = 's102-Ilastik';
