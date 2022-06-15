% chemotaxis OU model simulaiton

clear all
close all

%% parameters

OUki = 3000;
N_type = 80000;

epsilon = 25 ;

dx = 10;
x = -2000:dx:20000;
Lx = length(x);
s = 200*ones(size(x));
s(1:151) = 0; 
rho_all = zeros(size(x));
gamma = 0.7*10^-4;

Koff = 3.5; 
Kon = 1000; 


xx = zeros(length(OUki),N_type); 
x_ind = floor((xx-x(1))/dx)+1;
for ix = 1:Lx
    rho_all(ix) = sum(sum(x_ind==ix));
end
s = s - gamma*rho_all;
s = smooth(s,3)';
f = log((1+s/Koff)./(1+s/Kon));
dfdx = diff(f)/dx;

%% main loop

dt = 0.1;           % s
ot = 60;            % s
TT = 60*60;         % s
ot2 = 50*60;        % time to start tracking aquisition


ox = nan(length(OUki),N_type,ceil(TT/ot)); 
os = nan(length(x),ceil(TT/ot)); 
o_rho = nan(length(OUki),Lx,ceil(TT/ot));
ox2 = nan(length(OUki),N_type,ceil((TT-ot2)/dt)); 

kiki = repmat(OUki',1,N_type);
io = 1;
io2 = 1;

for it = 1:TT/dt
    
    xx(xx>x(end))=x(end);
    xx(xx<x(1))=x(1);
    x_ind = floor((xx-x(1))/dx)+1;
    vx = kiki.*dfdx(x_ind)*dt ;
    dvx = epsilon*randn(length(OUki),N_type)*sqrt(dt);
    xx = xx + vx + dvx;
    
    for ix = 1:Lx
        rho_all(ix) = sum(sum(x_ind==ix));
    end
    s = s - gamma*rho_all;
    s = smooth(s,3)';
    s(s<0)=0;
    f = log((1+s/Koff)./(1+s/Kon));
    dfdx = diff(f)/dx;
    
    if rem(it*dt,ot)==0
        disp(['time = ' num2str(it*dt/60) 'min']);
        yyaxis left
        plot(x,smooth(rho_all,3),'linewidth',2);
        yyaxis right
        plot(x,s,'linewidth',2);
        ylim([0 300]);
        title(['Time = ' num2str(it*dt/60) 'min']);
        drawnow;
        ox(:,:,io) = xx;
        o_rho(:,:,io) = rho_all;
        os(:,io) = s';
        io = io+1;
    end
    
    if it*dt>ot2
        ox2(:,:,io2) = xx;
        io2 = io2+1;
    end
end


for iot = 1:io-1
    x_ind = floor((ox(:,:,iot)-x(1))/dx)+1;
    for ix = 1:Lx
        o_rho(:,ix,iot) = sum(x_ind==ix,2);
    end
end

save('chemotaxis_OU_consumption_single_type.mat','-v7.3');



