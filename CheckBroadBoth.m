HOME = getenv('HOME');

load(fullfile(HOME,'ManualTrackingAndSeg','Broad','SaveDataWithSeg.mat'))


data_path = fullfile(HOME,'Input','Broad','Input','SubPop1');
extention = '*.tif';
expr = 's17_sb_(\d+).tif';
%%

resPath = fullfile(HOME,'Outputs','Results_Broad_CheckPoints','Results'); Name = 'MCF-10A-Ours'; %ORIGINAL
resPath = fullfile(HOME,'Outputs','Results_Broad-KDE_CheckPoints','Results'); Name = 'MCF-10A-Ours'; %ORIGINAL


resExt = '*.tif';

resExp = 'Seg_s17_sb_(\d+).tif';





CheckWeizman(cellData,Link,segData,data_path,extention,expr,resPath,resExt,resExp,Name)

%resPath = '/Users/assafarbelle/Documents/PhD/Ilastik/Raw_Data_ALL_Object-Identities.h5'; Name = 's102-Ilastik';
%%
%csvPathMitosis = '/Users/assafarbelle/IlastikProject/Weizmann/Raw Data-exported_data_divisions.csv'
Name = 'MCF-10A-Ilastik'
csvPath = '/Users/assafarbelle/IlastikProject/Broad/RawData-exported_data_table.csv';
resPath = '/Users/assafarbelle/IlastikProject/Broad/RawData_Object-Identities.h5';
CheckWeizmanIlastik(cellData,Link,segData,data_path,extention,expr,resPath,csvPath,Name)