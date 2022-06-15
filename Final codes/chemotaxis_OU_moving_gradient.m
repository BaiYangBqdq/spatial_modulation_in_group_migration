% chemotaxis OU model simulaiton

clear all
close all


load moving_gradient.mat
c = Vd; clear Vd;

%% get ki and s

OUki = [2547.09901242214,3831.74099920417,5215.00382964141,6422.80979213398,7747.80755133456,8934.49543167301,10319.2304624969,11293.5444486250];
Koff = 3.5;                         % uM 
Kon = 1000;                         % uM 
f = log((1+s/Koff)./(1+s/Kon));
dfdz = diff(f)/dz;
dsdz = diff(s)/dz;


%% chemotaxis OU model simulaiton
N_type = 1000;
dt = 0.1;                           % s
ot = 60;                            % s
TT = 60*60;                         % s
ot2 = 50*60;                        % time to start tracking aquisition
xx = zeros(length(OUki),N_type);    % initial posiiton 0
ox = nan(length(OUki),N_type,ceil(TT/ot)); 
ox2 = nan(length(OUki),N_type,ceil((TT-ot2)/dt)); 
x_ind = nan(length(OUki),N_type);
kiki = repmat(OUki',1,N_type);
dfdx = dfdz;
io = 1;
io2 = 1;

for it = 1:TT/dt

    x = z + c*it*dt;
    
    xx(xx>x(end))=x(end);
    xx(xx<x(1))=x(1);
    x_ind = floor(xx-x(1))/dz+1;
    vx = kiki.*dfdx(x_ind) + 100*randn(length(OUki),N_type);
    xx = xx + vx*dt;
    
    if rem(it*dt,ot)==0
        disp(['time = ' num2str(it*dt) 's']);
        ox(:,:,io) = xx;
        io = io+1;
    end
    
    if it*dt>ot2
        ox2(:,:,io2) = xx;
        io2 = io2+1;
    end
end

save('chemotaxis_OU_moving_gradient.mat','-v7.3');
