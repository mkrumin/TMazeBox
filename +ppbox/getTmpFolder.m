function str = getTmpFolder()

% This location will be used to keep all the processed data of the ppbox
% package
% Make sure you have enough space for all the processed files there
% Do not use network drives, this will slow down the process. If you have
% SSD - use it, it will make things faster. The ppbox code tries to be
% memory-efficient (this is required for large datasets), which means it is
% writing to disk more than minimally needed.

str = 'C:\Temp2p\';

[~, hostname] = system('hostname');
hostname = hostname(1:end-1);

switch hostname
    case 'ZERO'
        str = 'G:\Processing\';
    case 'zpike'
        str = 'C:\Users\Federico\Documents\Data\2P';
    case 'zufolo'
        str = 'D:\Data\2P';
    case 'zigzag'
        str = 'F:\elad\Tmp';
    case 'zunder'
%         str = 'J:\Processing\';
%         str = 'C:\STORAGE\OneDriveData\OneDrive for Business\DATA\InfoStructs';
        str = 'C:\STORAGE\OneDrive - University College London\Lab\DATA\InfoStructs';
    case 'ZEALOT'
        str = 'C:\Users\Julie\OneDrive - University College London\04_Data';
    case 'DESKTOP-OH1FRI9'
        str = 'C:\DATA\InfoStructs';
end


