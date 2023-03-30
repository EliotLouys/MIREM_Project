function userInfo = UserSessionInfo_MIREM(userName)

% Collect user's info
%
% JB Eichenlaub, 2023 || jb.eichenlaub@gmail.com (MATLAB2022a)

% userInfo
switch userName
    case 'jb1'
        rootDirProject   = '';
        rootDirUtilities = '';
    case 'jb2'
        rootDirProject   = 'D:\JBE_local\MIREM\';
        rootDirUtilities = 'D:\JBE_local\utilities\';
    case 'el1'
        rootDirProject   ='C:\Users\eliot\OneDrive\Documents\Stage\MIREM_Project\MIREM\';
        rootDirUtilities ='C:\Users\eliot\OneDrive\Documents\Stage\MIREM_Project\utilities\';
    case 'el2'
        rootDirProject   ='D:\Documents\MIREM\MIREM_Project\MIREM\';
        rootDirUtilities ='D:\Documents\MIREM\MIREM_Project\utilities\';        
end
userInfo = extractUserInfo(rootDirProject, rootDirUtilities);

% add folders to search path
addpath([userInfo.matlabScripts 'subFunctions' filesep]);
addpath([userInfo.matlabScripts 'fromOthers'   filesep]);

% Fieldtrip
addpath(userInfo.fieldtripDir);
ft_defaults;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function userInfo = extractUserInfo(rootDirProject, rootDirUtilities)

userInfo.dataDir       = [rootDirProject   'data\raw_data\'];
userInfo.scorDir       = [rootDirProject   'data\scoring\'];
userInfo.analysisDir   = [rootDirProject   'analysis\'];
userInfo.matlabScripts = [rootDirProject   'matlabScripts_MIREM\'];
userInfo.notebook      = [rootDirProject   'notebook_MIREM\'];
userInfo.fieldtripDir  = [rootDirUtilities 'FieldTrip\fieldtrip-20230118\'];