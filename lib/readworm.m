function adjacency = readworm(filename, startRow, endRow)
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

formatSpec = '%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray));
for col=1:length(dataArray)
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

rawNumericColumns = str2double(raw(:, 4));
rawStringColumns = string(raw(:, [1,2,3]));

for catIdx = [1,2]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

adjacency = table;
adjacency.Neuron1 = categorical(rawStringColumns(:, 1));
adjacency.Neuron2 = categorical(rawStringColumns(:, 2));
adjacency.EMSection = rawStringColumns(:, 3);
adjacency.Weight = rawNumericColumns(:, 1);

