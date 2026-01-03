
%%
%clear
clc

data = load('AL0625_1_dam.mat');


%%

wave_data_ = data.AL0625_1_dam;
wave_data = reshape(wave_data_,sqrt(size(wave_data_,1)),sqrt(size(wave_data_,1)),size(wave_data_,2));

wave_power = sum(wave_data.^2,3);

%%




mx = max(wave_power(:));
[ii,jj] = find(wave_power==mx);

wave_data = wave_data(:,:,1:100);

%%

lx = linspace(0,1,size(wave_data,1));
ly = linspace(0,1,size(wave_data,2));

axisFontsize = 15;

%%

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_power/max(wave_power(:)));
xlabel('\fontsize{16}Plate width [m]')
ylabel('\fontsize{16}Plate length [m]')
axis xy;
axis square;
ax = gca;
colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_energy.png','Resolution',600) 


%%

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data(:,:,5)/max(wave_data(:,:,5),[],'all'));
xlabel('\fontsize{16}Plate width [m]')
ylabel('\fontsize{16}Plate length [m]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_5.png','Resolution',600) 

%%

figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data(:,:,10)/max(wave_data(:,:,10),[],'all'));
xlabel('\fontsize{16}Plate width [m]')
ylabel('\fontsize{16}Plate length [m]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_10.png','Resolution',600) 



%%


figure,
set(gcf, 'Units', 'centimeters', 'Position', [10, 10, 13, 10])
imagesc(lx,ly,wave_data(:,:,13)/max(wave_data(:,:,13),[],'all'));
xlabel('\fontsize{16}Plate width [m]')
ylabel('\fontsize{16}Plate length [m]')
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize;
exportgraphics(gcf,'images/wave_at_13.png','Resolution',600) 

%%

load('positions.mat')
load('Positions_1.mat')
centers = [1;1;1;1]*[trueI,trueJ];
for i = 1:5
    
[a1,b1] = mids(1,1,centers(1,1),centers(1,2));
[a2,b2] = mids(1,100,centers(2,1),centers(2,2));
[a3,b3] = mids(100,100,centers(3,1),centers(3,2));
[a4,b4] = mids(100,1,centers(4,1),centers(4,2));

centers(1,1) = a1; centers(1,2) = b1;
centers(2,1) = a2; centers(2,2) = b2;
centers(3,1) = a3; centers(3,2) = b3;
centers(4,1) = a4; centers(4,2) = b4;
wave_data_ = wave_data;
wave_data_(a1:a4,b1:b2,:) = 0;

Trans = ones(size(wave_data));
Trans(a1:a4,b1:b2,:) = 0;
    
figure(i),
set(gcf, 'Units', 'centimeters', 'Position', [10, 5, 13*2, 10*2])
hold on;
imagesc(lx,ly,wave_data(:,:,5)/max(wave_data(:,:,5),[],'all'),'AlphaData',0.3);
imagesc(lx,ly,Trans(:,:,5)/max(Trans(:,:,5),[],'all'),'AlphaData',0.3);


H(1) = plot(trueJ/100,trueI/100,'-s','MarkerSize',8,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','red');
%text(trueJ/100+0.01,trueI/100,'GT')
%hold on;
H(2) = plot(Positions_HaLRTC{i}(2)/100,Positions_HaLRTC{i}(1)/100,'-p','MarkerSize',8,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','red');
%text(Positions_HaLRTC{i}(2)/100,Positions_HaLRTC{i}(1)/100,'HaLRTC')

H(3) = plot(Positions_lin_tubal{i}(2)/100,Positions_lin_tubal{i}(1)/100,'-d','MarkerSize',8,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','green');

%text(Positions_lin_tubal{i}(2)/100,Positions_lin_tubal{i}(1)/100,'TNN-DCT')


H(4) = plot(Positions_low_high_TIP{i}(2)/100,Positions_low_high_TIP{i}(1)/100,'-d','MarkerSize',8,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','blue');

%text(Positions_low_high_TIP{i}(2)/100,Positions_low_high_TIP{i}(1)/100+0.05,'HTNN-DCT')


H(5) = plot(Positions_PSTNN{i}(2)/100,Positions_PSTNN{i}(1)/100,'-d','MarkerSize',8,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor','black');



H(6) = plot(Positions_TC{i}(2)/100,Positions_TC{i}(1)/100,'-d','MarkerSize',8,...
    'MarkerEdgeColor','yellow',...
    'MarkerFaceColor','red');

P_wimf = Positions_WIMF{i}(6,:);
P_wimf_1 = Positions_WIMF{i}(1,:);
P_wimf_2 = Positions_WIMF{i}(2,:);
P_wimf_3 = Positions_WIMF{i}(3,:);
P_wimf_4 = Positions_WIMF{i}(4,:);
P_wimf_5 = Positions_WIMF{i}(5,:);

H(7) = plot(P_wimf(2)/100,P_wimf(1)/100,'-p','MarkerSize',8,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor','black');

H(8) = plot(P_wimf_1(2)/100,P_wimf_1(1)/100,'>','MarkerSize',8,'MarkerFaceColor','black');
H(9) = plot(P_wimf_2(2)/100,P_wimf_2(1)/100,'>','MarkerSize',8,'MarkerFaceColor','red');
H(10) = plot(P_wimf_3(2)/100,P_wimf_3(1)/100,'>','MarkerSize',8,'MarkerFaceColor','yellow');
H(11) = plot(P_wimf_4(2)/100,P_wimf_4(1)/100,'>','MarkerSize',8,'MarkerFaceColor','green');
H(12) = plot(P_wimf_5(2)/100,P_wimf_5(1)/100,'>','MarkerSize',8,'MarkerFaceColor','magenta');


%text(P_wimf_1(2)/100,P_wimf_1(1)/100,num2str(1))
if i==1
legend(H,{'GT','HaLRTC','TNN-DCT','HTNN-DCT','PSTNN','TNN','WIRL-full','WIRL-1',...
    'WIRL-2','WIRL-3','WIRL-4','WIRL-5'},'fontSize',12,'Location','NorthWest')
end
xlim([0,1])
ylim([0,1])
axis xy;
axis square;
ax = gca;
%colorbar;
ax.FontSize=axisFontsize*2;
xlabel('\fontsize{32}Plate width [m]')
ylabel('\fontsize{32}Plate length [m]')
grid on;
grid minor;

exportgraphics(gcf,strcat('images/recon',string(i),'.png'),'Resolution',600) 

end

