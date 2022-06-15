% chemotaxis OU model simulaiton

clear all
close all

%% parameters
sr = 1;                        
N_total = round(1/sr*10*10^4);   % 10~20*10^4 cell in the wave (exp data);

TB_mean = 0.325;
TB_std = 0.098;
TB = TB_mean + TB_std*randn(N_total,1);     TB(TB<0)=0;

g1=28;g2=g1;
Kd      =  2.9;
w       =  3.8;
Yp=g2*Kd./(g2-g1/2+log(1./TB-1))-Kd;
deltaG  =  g1/4-g2/2*(Yp./(Kd+Yp));
lambdaT   =  w*exp(deltaG);
lambdaR  =  w*exp(-deltaG);
Drot = 0.062;

v = 31;                         % um/s run speed on slide
theta=1-0.1564;
Diff=v.^2.*(1-TB)./3./(lambdaR*theta+Drot);
ki=22*Diff./(1+0.05./TB);
ki_std = std(ki);

epsilon = sqrt(2*Diff) ;  

dx = 10;                        % um
csa = 600*20;                   % section area of the channel  = 600um *20um
x = [0:dx:20000]';
Lx = length(x);
s = 200*ones(size(x));          % uM
s(x<-500)=0;

rho = zeros(size(x));
gamma0 =6*10^-12/60 ;           % maximum aspartate consumption rate 9.3*10^-12 umol/cell/s
Ks = 0.1;                       % uM, MM constant for consumtion

Koff = 1;                       % uM
Kon = 100;                      % uM 


xx = zeros(size(ki));           % initial posiiton 0
x_ind = floor((xx-x(1))/dx)+1;
for ix = 1:Lx
    rho(ix) = sum(x_ind==ix)*sr;
end
Ds = 1000;

%% main loop

dt = 0.05;       % s
ot = 60;         % s
TT = 60*60;      % s
ot2 = 50*60;     % time to start tracking aquisition
dt2 = 1;         % s

ox = nan(N_total,ceil(TT/ot)); 
os = nan(Lx,ceil(TT/ot)); 
orho = nan(Lx,ceil(TT/ot));
ox2 = nan(N_total,ceil((TT-ot2)/dt2)); 

io = 1;
io2 = 1;
sa=s;sb=s;sc=s;sd=s;

for it = 1:TT/dt
    
    s_gamma = gamma0*rho.*s./(s+Ks)/(csa*dx*10^-15);
    s1=s;s1(s1<=0)=0;
    sa(2:end-1)=Ds*(s1(1:end-2)+s1(3:end)-2*s1(2:end-1))/dx^2-s_gamma(2:end-1);      sa(1)=sa(2); sa(end)=sa(end-1);
    s1=s+sa*dt/2;s1(s1<=0)=0;
    sb(2:end-1)=Ds*(s1(1:end-2)+s1(3:end)-2*s1(2:end-1))/dx^2-s_gamma(2:end-1);      sb(1)=sb(2); sb(end)=sb(end-1);
    s1=s+sb*dt/2;s1(s1<=0)=0;
    sc(2:end-1)=Ds*(s1(1:end-2)+s1(3:end)-2*s1(2:end-1))/dx^2-s_gamma(2:end-1);      sc(1)=sc(2); sc(end)=sc(end-1);
    s1=s+sc*dt;s1(s1<=0)=0;
    sd(2:end-1)=Ds*(s1(1:end-2)+s1(3:end)-2*s1(2:end-1))/dx^2-s_gamma(2:end-1);      sd(1)=sd(2); sd(end)=sd(end-1);
    s=s+1/6*(sa+2*sb+2*sc+sd)*dt;   s(s<=0)=0;
    s(1)=s(2); s(end)=s(end-1); 

    
    f = log((1+s/Koff)./(1+s/Kon));
    dfdx = diff(f)/dx;
    
    xx(xx>x(end))=x(end);
    xx(xx<x(1))=x(1);
    x_ind = floor((xx-x(1))/dx)+1;
    vx = ki.*dfdx(x_ind)*dt ;
    dvx = epsilon.*randn(N_total,1)*sqrt(dt);
    xx = xx + vx + dvx;
    
    for ix = 1:Lx
        rho(ix) = sum(x_ind==ix)*sr;
    end
    
    if rem(it*dt,ot)==0
        disp(['time = ' num2str(it*dt/60) 'min']);
        yyaxis left
        plot(x,smooth(rho,3)/(csa*dx*10^-12),'linewidth',2);
        yyaxis right
        plot(x,s,'linewidth',2);
        ylim([0 300]);
        title(['Time = ' num2str(it*dt/60) 'min']);
        drawnow;
        ox(:,io) = xx;
        orho(:,io) = rho;
        os(:,io) = s;
        io = io+1;
    end
    
    if it*dt>ot2
        if rem(it*dt,dt2)==0
            ox2(:,io2) = xx;
            io2 = io2+1;
        end
    end
end


save('chemotaxis_OU_consumption_multi_type.mat','-v7.3');

%% movie of sub populations

ki_bin = [min(ki) ki_mean-2*ki_std ki_mean-ki_std ki_mean ki_mean+ki_std ki_mean+2*ki_std max(ki)];
ki_center = (ki_bin(1:end-1) + ki_bin(2:end))/2;
sect = length(ki_center);
Yki = discretize(ki,ki_bin);
cj = colormap('jet');
mc = cj(1:floor(64/sect):64,:);
orho_sub = nan(sect,Lx,ceil(TT/ot));

% v = VideoWriter([savingname '.avi']);
% open(v);
figure('position',[100 100 1100 400]);

for iot = 1:io-1
    yyaxis left
    for ik = 1:sect
        x_ind = floor((ox(Yki==ik,iot)-x(1))/dx)+1;
        for ix = 1:Lx
            orho_sub(ik,ix,iot) = sum(x_ind==ix);
        end
        plot(x,orho_sub(ik,:,iot)/(csa*dx*10^-12)*sr,'-','linewidth',3,'color',mc(ik,:)); hold all
    end
    ylabel('density cell/mL');
    hold off
    yyaxis right
    plot(x,os(:,iot),'k-','linewidth',3);
    xlabel('x [um]');
    ylabel('Asp [uM]');
    title(['time = ' num2str(iot*ot/60) 'min']);
    drawnow;
    frame = getframe();
    writeVideo(v,frame);
end
% close(v);
%% dynamic figure
h1=figure('position',[100 100 1100 400]);
dctime = floor(64/(io-1));
mctime = cj([1:io-1]*dctime,:);
for iot = 1:io-1
	plot(x,orho(:,iot)/(csa*dx*10^-12),'-','linewidth',2,'color',mctime(iot,:)); hold all
    legends_time{iot} = ['t=' num2str(iot)*ot/60 'min'];
end
% legend(legends_time);
xlabel('x [um]');
ylabel('density [cell/mL]');

%% peak pos
h2 = figure('position',[200 100 900 400]);
o_time = [1:io-1]*ot;
ot_stable = 50;
subplot(1,2,1)
for j = 1:io-1
    [~,peak_ind_all(j)] = max(squeeze(sum(orho_sub(:,:,j),1)));
    peak_pos_all(j) = x(peak_ind_all(j));

end
plot(o_time/60,peak_pos_all,'o','linewidth',2); hold all
pf = polyfit(o_time(ot_stable:end), peak_pos_all(ot_stable:end),1);
c = pf(1);
plot(o_time(ot_stable:end)/60,polyval(pf,o_time(ot_stable:end)),'k-','linewidth',2);
legend('peak pos all',['c = ' num2str(c) 'um/s'],'location','best');
xlabel('time [min]');
ylabel('peak position [um]');
set(gca,'fontsize',14);

subplot(1,2,2)
for i = 1:sect
    for j = 1:io-1
        [~,peak_ind(i,j)] = max(smooth(orho_sub(i,:,j),5));
        peak_pos(i,j) = x(peak_ind(i,j));
    end
    plot(o_time/60,peak_pos(i,:),'o','linewidth',2,'color',mc(i,:)); hold all
    legends{i} = sprintf('mean chi = %.0f',ki_center(i));
end
xlabel('time [min]');
ylabel('peak position [um]');
legend(legends{1:end-1},'location','best');
set(gca,'fontsize',14);
drawnow;


%% density all in mc compare
load('exp_plotv202003_compare.mat','xcenter_MC','density_profile_WT_track_mean_s');
cl = colormap('lines');
x_z = -4000:dx:4000;
rho_z_all = nan(length(x_z),ot-ot_stable);
for j = ot_stable:io-1
    rho_z_all(:,j-ot_stable+1) = orho(peak_ind_all(j)-400:peak_ind_all(j)+400,j);
end
rho_z_all_mean = nanmean(rho_z_all,2);
rho_z_all_std = nanstd(rho_z_all,[],2);


h3 = figure();
shadedErrorBar(x_z, rho_z_all_mean/100*2*10^9, rho_z_all_std/100*2*10^9,'lineprops',{'-','color',cl(1,:)}); hold all
plot(x_z, smooth(rho_z_all_mean,5)/100*2*10^9, 'k-', 'linewidth', 2,'color',cl(1,:)); hold all

shadedErrorBar(xcenter_MC*1000, mean(density_profile_WT_track_mean_s), std(density_profile_WT_track_mean_s),'lineprops',{'-','color',cl(2,:)}); hold all
plot(xcenter_MC*1000,mean(density_profile_WT_track_mean_s),'linewidth',2,'color',cl(2,:)); hold all
xlabel('z (\mum)');
ylabel('cell counts');


%% density sub in mc
rho_z = nan(sect,length(x_z),io-1-ot_stable);
for i = 1:sect
    for j = ot_stable:io-1
        rho_z(i,:,j-ot_stable+1) = orho_sub(i,peak_ind_all(j)-400:peak_ind_all(j)+400,j);
        
    end
end

s_z = nan(length(x_z),io-1-ot_stable);
for j = ot_stable:io-1
    s_z(:,j-ot_stable+1) = os(peak_ind_all(j)-400:peak_ind_all(j)+400,j);
end
rho_z_mean = nanmean(rho_z,3);
rho_z_std = nanstd(rho_z,[],3);
s_z_mean = nanmean(s_z,2);
s_z_std = nanstd(s_z,[],2);
rho_z_mean_smooth = nan(size(rho_z_mean));
for i = 1:sect
    rho_z_mean_smooth(i,:) = smooth(rho_z_mean(i,:),5);
end

h4 = figure();
yyaxis left
for i = 1:sect
    shadedErrorBar(x_z, rho_z_mean(i,:), rho_z_std(i,:),'lineprops',{'-','color',mc(i,:)}); hold all
    plot(x_z, rho_z_mean_smooth(i,:), '-', 'linewidth', 2 ,'color',mc(i,:)); hold all
end
xlabel('z (\mum)');
ylabel('cell counts');
ax = gca();
ylim([0 ax.YLim(2)]);
% legend(legends,'location','best');
yyaxis right
shadedErrorBar(x_z',s_z_mean,s_z_std,'lineprops',{'-','color',[0.5 0.5 0.5]}); hold all
plot(x_z',s_z_mean,'k-','linewidth',2); hold all
ylim([0 200])
xlabel('z (\mum)');
ylabel('s');
set(gca,'fontsize',14);

axes('Position',[.2 .6 .25 .25],'box','on');
[ki_y, ki_x] = hist(ki,30);
plot(ki_x, ki_y/sum(ki_y)/(ki_x(2)-ki_x(1)),'linewidth',2); hold all
ax = gca();
for i = 1:length(ki_bin)
    plot([ki_bin(i) ki_bin(i)], [ax.YLim],'k--','linewidth',2); hold all
end
xlabel('\chi');
ylabel('PDF');
drawnow;


[~,peak_ind_z] = max(rho_z_mean_smooth,[],2);
peak_pos_z = x_z(peak_ind_z);

%% 
% xbin = 100;
% xcenter_MC = -3000 : xbin : 3000 ;
% xedge_MC = -3000-xbin/2 : xbin : 3000+xbin/2;
% rho_MC = nan(length(OUki),length(xcenter_MC));
% x_MC = nan(length(OUki),N_type,length(o_time)-ot_stable+1);
% 
% for j = 1:length(OUki)
%     x_MC(j,:,:) = squeeze(ox(j,:,ot_stable:end)) - repmat(peak_pos_all(ot_stable:end),N_type,1);
%     ytemp = squeeze(discretize(x_MC(j,:,:),xedge_MC));
%     for ip = 1:length(xcenter_MC)
%         rho_MC(j,ip) = nansum(nansum(ytemp==ip));
%     end
% end
% 
% h6 = figure;  
% for j = 1:length(OUki)
%     plot(xcenter_MC,rho_MC(j,:),'-','linewidth',3,'color',mc(j,:)); hold all
% end
% xlabel('z (\mum)');
% ylabel('cell counts');
% legend(legends{1:end-1},'location','nw');
% % xlim([-2000 2000]);
% set(gca,'fontsize',14);
% drawnow;

% [~,peak_ind_MC] = max(rho_MC,[],2);
% peak_pos_MC = xcenter_MC(peak_ind_MC);

%% chi*g(z)

f_z_mean = log((1+s_z_mean/Koff)./(1+s_z_mean/Kon));
g_z_mean = diff(f_z_mean)/dx;
chig_cross = nan(sect,1);
for i = 1:sect
    try
    temp_ind = find(ki_center(i)*g_z_mean>c);
    chig_cross(i) = x_z(temp_ind(end));
    catch error
        continue;
    end
end


h5 = figure();
yyaxis left
for i = 1:sect
    plot(x_z, rho_z_mean_smooth(i,:), '-', 'linewidth', 2 ,'color',mc(i,:)); hold all
end
ylabel('density');
yyaxis right
for i = 1:sect
    plot(x_z(2:end), ki_center(i)*g_z_mean', '-', 'linewidth', 2 ,'color',mc(i,:)); hold all
end
ylabel('\chi*g(z)');
xlabel('z (um)');
% legend(legends{1:end-1},'location','best');
%% compare peaks to cross point
for i = 1:sect
    x_z_average(i) = weighted(x_z,rho_z_mean(i,:));
end

h6 = figure();
plot(ki_center,peak_pos_z,'o','linewidth',2); hold all
plot(ki_center,x_z_average,'x','linewidth',2); hold all
plot(ki_center,chig_cross,'-d','linewidth',2); hold all
xlabel('\chi');
ylabel('relative peak pos z [um]');
legend('peak pos','average pos','chi*g cross pos','location','best');
% title('peak posistion in moving coordinate');
% xlim([2000 10000]);
set(gca,'fontsize',14);

%% get MSD

OU_track_time = ot2+dt:dt2:TT;
ox2_MC = ox2 - repmat(c*OU_track_time,N_total,1);
OU_MSD_xMC = (ox2_MC - repmat(ox2_MC(:,1),1,length(OU_track_time))).^2 ;
OU_MSD_mean_xMC = nanmean(OU_MSD_xMC);

%% plot MSD
h7 = figure('position',[100 400 900 400]);
subplot(1,2,1)
plot(OU_track_time-OU_track_time(1),OU_MSD_mean_xMC,'-','linewidth',3); hold all
xlabel('time (s)');
ylabel('msd (\mum^2)');
title('MSD linear plot');
set(gca,'fontsize',14);
subplot(1,2,2)
loglog(OU_track_time-OU_track_time(1),OU_MSD_mean_xMC,'-','linewidth',3); hold all
xlabel('time (s)');
ylabel('msd (\mum^2)');
title('MSD log plot');
set(gca,'fontsize',14);
drawnow;

%% chi(z) & g(z)
f_z = log((1+s_z/Koff)./(1+s_z/Kon));
g_z = diff(f_z)/dx;
g_z_mean = mean(g_z,2);
g_z_std = std(g_z,[],2);
x_z_center = (x_z(1:end-1) + x_z(2:end))/2;

x_z_bin = x_z(1:10:end);
chi_z_mean = nan(size(x_z_bin));
chi_z_std = nan(size(x_z_bin));
% Yx = discretize(ox2_MC,x_z_bin);
Yx = discretize(ox2_MC,x_z);
% kiki = repmat(ki,1,length(OU_track_time));
kiki = repmat(ki,1,length(OU_track_time));
for ix = 1:length(x_z_bin)
    chi_z_mean(ix) = nanmean(kiki(Yx==ix));
    chi_z_std(ix) = std(kiki(Yx==ix));
end


h8 = figure();
yyaxis left
shadedErrorBar(x_z_center',g_z_mean,g_z_std,'lineprops',{'-','color',cl(1,:)}); hold all
xlabel('rel coord (um)');
ylabel('g(z)');
yyaxis right
shadedErrorBar(x_z_bin,chi_z_mean,chi_z_std/sqrt(length(OU_track_time)),'lineprops',{'-','color',cl(2,:)}); hold all
ylabel('\chi(z)');
xlim([-1500 1500])

%% save
