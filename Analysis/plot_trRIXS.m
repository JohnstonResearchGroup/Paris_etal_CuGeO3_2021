clear;
clc;
close all;

Data = dlmread('../Data/zrs_intensity_vs_t_ppol.dat');
Data_Teffonly = dlmread('../Data/zrs_intensity_vs_t_ppol_Teffonly.dat');
Data_Phonly = dlmread('../Data/zrs_intensity_vs_t_ppol_phonly.dat');
RIXSpp = dlmread('../Data/rixs_Cu3O8_Ipp_Spump_Amp5.000_win0.500_geq0.000_gam0.200_000.dat');
RIXSps = dlmread('../Data/rixs_Cu3O8_Ips_Spump_Amp5.000_win0.500_geq0.000_gam0.200_000.dat');

fh = figure('Renderer', 'painters', 'Position', [10 10 700 800]);
set(fh,'color','white')

hbar = 0.658211951;
omega = 0.15d0*0.0136d0/hbar;

H = 0.38;
W = 0.38;
B = 0.12;
L = 0.1;
Offset = 0.12;

delays = RIXSpp(find(RIXSpp(:,2) == RIXSpp(1,2)),1);
omega = RIXSpp(find(RIXSpp(:,1)==RIXSpp(1,1)),2);

Intensity = reshape(RIXSpp(:,3),length(omega),length(delays)) + ...
            reshape(RIXSps(:,3),length(omega),length(delays));

subplot('position',[L,B+0.1+H,W,H]);
imagesc(delays/1000,omega,Intensity)
axis([0,4.5,3.2,4.2])
set(gca,'FontName','Arial','FontSize',14)
box on;
caxis([0,0.5])
xlabel('Time delay (ps)','FontSize',16,'FontName','Arial')
ylabel('Energy Loss (eV)','FontSize',16,'FontName','Arial')
text(-1.18,3.2,'(a)','FontSize',16','FontName','Arial')

subplot('position',[L+Offset+W,B+0.1+H,W,H]); hold on
plot(omega,Intensity(:,1),'-k','DisplayName',[num2str(delays(1),2) ' ps']);
plot(omega,Intensity(:,6),'-r','DisplayName',[num2str(delays(6)/1000,2) ' ps']);
plot(omega,Intensity(:,19),'-b','DisplayName',[num2str(delays(19)/1000,2) ' ps']);
plot(omega,Intensity(:,61),'-m','DisplayName',[num2str(delays(61)/1000,2) ' ps']);
plot(omega,Intensity(:,120),'-c','DisplayName',[num2str(delays(120)/1000,2) ' ps']);
axis([3.2,4.4,0,0.5])
set(gca,'FontName','Arial','FontSize',14)
box on;
xlabel('Energy Loss (eV)','FontSize',16,'FontName','Arial')
ylabel('Intensity (a.u.)','FontSize',16,'FontName','Arial')
text(2.88,0.5,'(b)','FontSize',16','FontName','Arial')
legend('location','northeast')
legend boxoff;

subplot('position',[L,B,W,H]); hold on;
x = [-0.5:0.01:0];
plot([x,Data(:,1)'],[zeros(size(x)),100*Data(:,2)'],'-k');
axis([-0.5,8,-2,4])
set(gca,'FontName','Arial','FontSize',14)
box on;
xlabel('Time delay (ps)','FontSize',16,'FontName','Arial')
ylabel(['Displacement (',char(215),'10^{-2} ',char(197),')'],'FontSize',16,'FontName','Arial')
text(-2.75,4,'(c)','FontSize',16','FontName','Arial')

subplot('position',[L+Offset+W,B,W,H]); hold on;
x = [-0.5:0.01:0];
plot([x,Data(:,1)'],[ones(size(x)),Data(:,4)'./Data(1,4)],'-k','DisplayName','Full Model');
plot([x,Data_Teffonly(:,1)'],[ones(size(x)),Data_Teffonly(:,4)'./Data_Teffonly(1,4)],'-.r','DisplayName','T_{eff} only');
plot([x,Data_Phonly(:,1)'],[ones(size(x)),Data_Phonly(:,4)'./Data_Phonly(1,4)],'--b','DisplayName','Ph only');
plot(Data( 1,1),Data( 1,4)./Data(1,4),'ok','MarkerFaceColor','k','HandleVisibility','off');
plot(Data( 4,1),Data( 4,4)./Data(1,4),'or','MarkerFaceColor','r','HandleVisibility','off')
plot(Data(13,1),Data(13,4)./Data(1,4),'ob','MarkerFaceColor','b','HandleVisibility','off')
plot(Data(40,1),Data(40,4)./Data(1,4),'om','MarkerFaceColor','m','HandleVisibility','off')
plot(Data(79,1),Data(79,4)./Data(1,4),'oc','MarkerFaceColor','c','HandleVisibility','off')
legend('location','southwest')
legend boxoff;

axis([-0.5,8,0.5,1.1])
set(gca,'FontName','Arial','FontSize',14)
box on;
xlabel('Time delay (ps)','FontSize',16,'FontName','Arial')
ylabel('ZRS Intensity (a.u.)','FontSize',16,'FontName','Arial')
text(-2.75,1.1,'(d)','FontSize',16','FontName','Arial')

saveas(gcf,'figure4.eps','epsc')
