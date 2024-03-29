function [data] = ImportSmartPhoneData(filename,filename2,startRow, endRow,startRow2, endRow2)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [ACCX,ACCY,ACCZ,TIMESTAMP] = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   [ACCX,ACCY,ACCZ,TIMESTAMP] = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [AccX,AccY,AccZ,timestamp] = importfile('Acc_Wed Mar 11 140802 GMT+0000 2015.txt',1, 5614);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/04/08 10:48:06

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
AccX = dataArray{:, 1};
AccY = dataArray{:, 2};
AccZ = dataArray{:, 3};
timestamps = dataArray{:, 4};
%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow2 = 1;
    endRow2 = inf;
end

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID2 = fopen(filename2,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID2, formatSpec, endRow2(1)-startRow2(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow2(1)-1, 'ReturnOnError', false);
for block=2:length(startRow2)
    frewind(fileID2);
    dataArrayBlock = textscan(fileID2, formatSpec, endRow2(block)-startRow2(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines', startRow2(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID2);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
gyroX = dataArray{:, 1};
gyroY = dataArray{:, 2};
gyroZ = dataArray{:, 3};
timestamps2 = dataArray{:, 4};
samplingfreq=200;
data=struct('accX',AccX,'accY',AccY,'accZ',AccZ,'gyroX',gyroX,'gyroY',gyroY,'gyroZ',gyroZ,'acc_timestamp',timestamps,'gyr_timestamp',timestamps2,'fs',samplingfreq);



