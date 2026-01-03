%% Generate experimental data for videos
clear all;
close all;
%%

clc

data = load('AL0625_1_dam.mat');


%%

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_data = wave_data(:,:,1:150);

wave_data_un = wave_data;

Trans = ones(size(wave_data));
Trans(6:90,20:90,:) = 0;

wave_data(6:90,20:90,:) = 0;

trans = reshape(Trans,size(Trans,1)*size(Trans,2)*size(Trans,3),1);
%%
% 1. Please put in the name of the video accordingly
video=wave_data_un;

% 2. Run the lines
video=permute(video,[1,2,4,3]);
sz=size(video);
num=prod(sz(1:2));
ratio=0.9;
num=floor(num*ratio);   % number of lost entries
ind=randperm(prod(sz(1:2)));
mi=ind(1:num);
kn=ind(:,num+1:end);

Mi=mi;
Kn=kn;
for i=1:prod(sz(3:4))-1
    Mi=[Mi,i*prod(sz(1:2))+mi];   % the missing index for all three channel
    Kn=[Kn,i*prod(sz(1:2))+kn];
end

% 3. Choose your favorite name and name it!
save('favorite.mat','video','ratio','Mi','Kn','mi','kn')