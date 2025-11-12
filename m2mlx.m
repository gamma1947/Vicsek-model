% convert .m file to .mlx file
% This code allow you converting multiple files
%ex-  fileName=["evenodd","ievenodd"];
% m2mlx(fileName)

function []=m2mlx(model.m)
tic

%preallocate the file size
mlx_file=fileNameInString;
m_file=fileNameInString;

% To make the .mlx and.m files
for i=1:length(fileNameInString)
    mlx_file(i) =strcat(fileNameInString(i),'.mlx');
    m_file(i) = strcat(fileNameInString(i),'.m');
end

%  Convert the .mlx file to .m file and dlete the .mlx file
for i=1:length(fileNameInString)
    matlab.internal.liveeditor.openAndSave(convertStringsToChars(m_file(i)),convertStringsToChars(mlx_file(i))); % convert .m 2 .mlx
    delete (convertStringsToChars(m_file(i))); %delete the .m file
end
disp('All .m files convereted into .mlx files and deleted all .m files')
toc
end