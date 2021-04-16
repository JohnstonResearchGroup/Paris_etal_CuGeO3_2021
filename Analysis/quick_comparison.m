clear;
clc;
close all;
A = dlmread('../Data/highf_timetrace.txt');
B = dlmread('../Code/zrs_intensity_vs_t_ppol.dat');

%normalize the data
I = find(A(:,1)<0);
norm = mean(A(I,2));
A(:,2) = A(:,2)/norm;
A(:,3) = A(:,3)/norm;

xx = [-1:0.1:0];
yy = ones(size(xx));

figure(1); hold on;
set(0,'defaulttextInterpreter','latex');
plot(B(:,1),B(:,4),'-k',xx,yy,'-k');
errorbar(A(:,1),A(:,2),A(:,3),'or','MarkerFaceColor','r')
axis([-1,10,0.5,1.1])
xlabel('Time (ps)','FontSize',30);
ylabel('$I_\mathrm{ZRS}(t)/I_\mathrm{ZRS}(t=0)$','FontSize',30)
set(gcf,'color','white')
set(gca,'FontSize',25,'XTick',[-2:2:10],'YTick',[0.5:0.1:1.2])
box on;
