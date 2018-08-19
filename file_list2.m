function [file_array,num_file,date_array]=file_list2(folder)
a='2*data.ciemo';
files = dir(fullfile(folder,a));
files1={files.name};
out=regexp(files1,'\d+','match');
num_file=size(out,2);
file_array=cell(1,num_file);
date_array=zeros(1,num_file);
for ii=1:num_file
    time=out{1,ii}{1,2};
    file_array{1,ii}=time;
    date_array(1,ii)=str2double(out{1,ii}{1,1});
end

end
