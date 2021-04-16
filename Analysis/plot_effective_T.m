A = load('../Data/ZRS_temp_20201227.txt');

set(0,'defaulttextInterpreter','latex') %latex axis labels

I = find(A(:,1)>=0);
x = A(:,1);
y = 21.301302*ones(size(x));          %Numerical values were obtained by fits with cftool
y(I) = y(I) + 181.2*(1-exp(-x(I)/1.949));

fh = figure(17); 
hold on;
plot(A(:,1),A(:,2),'.')
plot(x,y,'--r')
xlabel('Time (ps)','FontSize',30);
ylabel('$T_\mathrm{eff}$ (K)','FontSize',30);
set(fh,'color','white');
set(gca,'XTick',[-1:1:10],'YTick',[0:50:225],'FontSize',25);
axis([-1,10,0,225])
box on;


