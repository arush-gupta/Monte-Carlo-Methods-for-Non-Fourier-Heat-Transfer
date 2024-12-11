% Arush Gupta Thursday 4th and 11th July, 2019
% Jean-Philippe Peraud - December 3, 2015
close all
clear all
%%% this script does the following:
% -loads silicon data from dataSi.txt
% -solves the LINEARIZED phonon transport problem between four walls (2D)
% with given prescribed temperatures T_l, T_r, T_d, and T_u
 
 
 
 
%%%%% load the data
load 'dataSi.txt'; % Note: instead of having two TA branches in those data,
                   % we accounted for the double effect by multiplying the
                   % density of states by two for the TA mode.
% the data is organized as follows:
% - each line represents a phonon frequency
% - column 1 is the frequency
% - column 2 is the density of states
% - column 3 is the group velocity
% - column 4 is the width of the frequency cell (often a constant but not
% necessarily
% - column 5 is the relaxation time for this particular mode
% - column 6 is the polarization (1 if LA, 2 if TA). Not strictly needed
% for this code
 
% define reduced Planck constant and Bolzmann constant
hbar=1.054517e-34;
boltz=1.38065e-23;
 
% define linearization temperature (also referred to as "equilibrium"
% temperature)
Teq = 300;
 
Tinit = 300; % initial temperature
Tl = 310; % left wall
Tr = 300; % right wall
Tu = 300; % top wall
Td = 300; % bottom wall
dORs='d'; % 'd' for diffuse walls, 's' for specular walls (degree = 1)
spec = 0; %degree of specularity

LX = 1e-7; % System length
LY = 1e-7;
N = 20000000; % number of particles
 
NX = 25; % number of spatial cells
NY = 25;
Ntt = 100; % number of measurement times
tmax = 9.9e-10; % maximum time
 
tt = 0:tmax/(Ntt-1):tmax; %times of measurement
%Ntt = length(tt); % in some cases one may want to define tt differently than above. 
                  % Important for Ntt to be exactly the length of it
                   
 
xx = 1e9*(LX/2/NX:LX/NX:LX-LX/2/NX); % centroids of cells for measuring temperature
                    % not really needed but useful for plotting results
yy = 1e9*(LY/2/NY:LY/NY:LY-LY/2/NY);                    
                     
 
Dx = LX/NX; % length of a cell
Dy = LY/NY;
 
% dataSi(:,5) = 4e-11*ones(Nmodes,1);  %% <-- uncomment this for single relaxation time
SD = dataSi(:,2);  % Density of states
V = dataSi(:,3); % velocities
Dom = dataSi(:,4); % Delta frequencies
tau_inv = 1./dataSi(:,5); % relaxation rates
tau = dataSi(:,5); % relaxation times
F = dataSi(:,1); % frequencies
de_dT = (hbar*F/Teq).^2/boltz.*exp(hbar*F/boltz/Teq)./(exp(hbar*F/boltz/Teq)-1).^2; %derivative of Bose-Einstein
 
% plot(F,SD,'sb')
% hold on
% plot(F(1:1000),SD(1:1000),'ok')
% hold on
% plot(F(1001:end),SD(1001:end),'or')
 
 
Nmodes = length(F); 
 
T  = zeros(Ntt,NX,NY); % array storing the temperature solution
Qx = zeros(Ntt,NX,NY); % array storing the heat flux solution
Qy = zeros(Ntt,NX,NY); 

% cumulative distribution functions
cumul_base = zeros(Nmodes,1);
cumul_coll = zeros(Nmodes,1);
cumul_vel  = zeros(Nmodes,1);
 
cumul_base(1) = SD(1)*de_dT(1)*Dom(1);
cumul_coll(1) = SD(1)*de_dT(1)*Dom(1)*tau_inv(1);
cumul_vel(1) = SD(1)*de_dT(1)*Dom(1)*V(1);
 
for i=2:Nmodes
    cumul_base(i) = cumul_base(i-1)+SD(i)*de_dT(i)*Dom(i);
    cumul_coll(i) = cumul_coll(i-1)+SD(i)*de_dT(i)*tau_inv(i)*Dom(i);
    cumul_vel(i) = cumul_vel(i-1) + SD(i)*de_dT(i)*V(i)*Dom(i);
end
C = cumul_base(Nmodes); % Heat capacity at Teq
 
% plot(F,cumul_base)
 
%Deviational energy from the different sources
enrgInit = LX*LY*C*abs(Tinit-Teq);
enrgLeft = LY*cumul_vel(Nmodes)*tt(Ntt)*abs(Tl-Teq)/4;
enrgRight = LY*cumul_vel(Nmodes)*tt(Ntt)*abs(Tr-Teq)/4;
enrgUp = LX*cumul_vel(Nmodes)*tt(Ntt)*abs(Tu-Teq)/4;
enrgDown = LX*cumul_vel(Nmodes)*tt(Ntt)*abs(Td-Teq)/4;
 
 
% enrgInit = LX*LY*C*abs(Tinit-Teq);
% enrgLeft = pi*cumul_velsqr(Nmodes)*tt(Ntt)*tt(Ntt)*abs(Tl-Teq)/4;   %%%%%%%%%denominator to be reviewed
% enrgRight = pi*cumul_velsqr(Nmodes)*tt(Ntt)*tt(Ntt)*abs(Tr-Teq)/4;  %%%%%%%%%denominator to be reviewed
% enrgUp = pi*cumul_velsqr(Nmodes)*tt(Ntt)*tt(Ntt)*abs(Tu-Teq)/4;    %%%%%%%%%denominator to be reviewed
% enrgDown = pi*cumul_velsqr(Nmodes)*tt(Ntt)*tt(Ntt)*abs(Td-Teq)/4;  %%%%%%%%%denominator to be reviewed
% 
% Total deviational energy
enrgTot = enrgInit + enrgLeft + enrgRight + enrgDown + enrgUp;
Eeff = enrgTot/N;
 
% enrgInitX = LX*C*abs(Tinit-Teq);
% enrgLeftX = cumul_vel(Nmodes)*tt(Ntt)*abs(Tl-Teq)/4;
% enrgRightX = cumul_vel(Nmodes)*tt(Ntt)*abs(Tr-Teq)/4;
% 
% enrgInitY = LY*C*abs(Tinit-Teq);
% enrgDownY = cumul_vel(Nmodes)*tt(Ntt)*abs(Td-Teq)/4;
% enrgUpY = cumul_vel(Nmodes)*tt(Ntt)*abs(Tu-Teq)/4;
 
% effective energy
% EeffX = (enrgInitX+enrgLeftX+enrgRightX)/N;
% EeffY = (enrgInitY+enrgDownY+enrgUpY)/N;
 
 
% calculate thermal conductivity (optional. just to make sure parameters are realistic)
ktest = sum(SD.*Dom.*V.*V.*tau.*de_dT)/3;
 
 
driftvelx=0;
driftvely=0;
meanphononlifetime=0;
tic
for i=1:N %loop over the N particles   
    flag = 0;
%%%% decide the origin of the particle   
    while ~flag
        Ri = rand(); %draw random number to decide is particle is emitted from left or right wall
                % or from initial condition 
                         
        if Ri < enrgInit/enrgTot % this case: emitted from initial source
            x0 = LX*rand();
            y0 = LY*rand();
            ind_mod = select_mode(cumul_base,Nmodes);
             
            Rx = 2*rand()-1;
            Ry = 2*rand()-1;
            Rz = 2*rand()-1;
             
            cosa = Rx/sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
            cosb = Ry/sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
             
                     
            vx = V(ind_mod)*cosa; %! not be confused between V and vx
            vy = V(ind_mod)*cosb;
         
            t0 = 0; % initial source => initial time of the particle is 0
            psign = sign(Tinit-Teq); % particle sign
             
            flag = 1;
         
        elseif Ri > enrgInit/enrgTot && Ri < (enrgInit+enrgLeft)/enrgTot  
                                      % emitted by left wall
            x0 = 0;
            y0 = LY*rand();
            ind_mod = select_mode(cumul_vel,Nmodes);
                         
            costhetax = sqrt(rand());
            vx = V(ind_mod)*costhetax;
             
            azy = 2*pi*rand();
            vy = V(ind_mod)*sqrt(1-costhetax*costhetax)*cos(azy);
             
            Vxy = sqrt(vx*vx+vy*vy);
             
            psign = sign(Tl-Teq);
            
            t0 = rand()*tt(Ntt); % emission by a wall => initial time of the particle 
                             % is anything between 0 and tmax
            flag = 2;
                              
        elseif Ri > (enrgInit+enrgLeft)/enrgTot && Ri < (enrgInit+enrgLeft+enrgRight)/enrgTot  % emitted by right wall
            x0 = LX;
            y0 = LY*rand();
            ind_mod = select_mode(cumul_vel,Nmodes);
            
            costhetax = -sqrt(rand());
            vx = V(ind_mod)*costhetax;
             
            azy = 2*pi*rand();
            vy = V(ind_mod)*sqrt(1-costhetax*costhetax)*cos(azy);
             
            Vxy = sqrt(vx*vx+vy*vy);
             
            psign = sign(Tr-Teq);
            
            t0 = rand()*tt(Ntt); % emission by a wall => initial time of the particle 
                             % is anything between 0 and tmax
            flag = 3;     
              
        elseif Ri > (enrgInit+enrgLeft+enrgRight)/enrgTot && Ri < (enrgInit+enrgLeft+enrgRight+enrgDown)/enrgTot % emitted by bottom wall
            x0 = LX*rand();
            y0 = 0;
            ind_mod = select_mode(cumul_vel,Nmodes);
             
            costhetay = sqrt(rand());
            vy = V(ind_mod)*costhetay;
             
            azx = 2*pi*rand();
            vx = V(ind_mod)*sqrt(1-costhetay*costhetay)*cos(azx);
             
            Vxy = sqrt(vx*vx+vy*vy);
             
            
            psign = sign(Td-Teq);
            t0 = rand()*tt(Ntt); % emission by a wall => initial time of the particle 
                             % is anything between 0 and tmax
            flag = 4;    
             
        else  % emitted by top wall
            x0 = LX*rand();
            y0 = LY;
            ind_mod = select_mode(cumul_vel,Nmodes);
             
            costhetay = -sqrt(rand());
            vy = V(ind_mod)*costhetay;
             
            azx = 2*pi*rand();
            vx = V(ind_mod)*sqrt(1-costhetay*costhetay)*cos(azx);
             
            Vxy = sqrt(vx*vx+vy*vy);
             
             
            psign = sign(Tu-Teq);  
            t0 = rand()*tt(Ntt); % emission by a wall => initial time of the particle 
                             % is anything between 0 and tmax
            flag = 5;                    
         
        end
    end
%%%%    
         
         
    phlegm = -1;
   
    finished = false;  % as long as "false", the current particle is active
    im = 1; % index for tracking measurement times
    while tt(im)<t0
        im = im+1;
    end
    
    phlegm = 0;
    %The while loop below is for propagating the particle till its end
    
    while ~finished
        Delta_t = -tau(ind_mod)*log(rand()); % time to next scattering event
        t1 = t0 + Delta_t; % time of next scattering event
        x1 = x0 + Delta_t*vx; % position of next scattering event
        y1 = y0 + Delta_t*vy;
        % --- this part handles the contribution of the current particle to
        % the final estimates ------------------- %
        while (im<Ntt+1 && t0<=tt(im) && t1>tt(im))
            x_ = x0 + (tt(im)-t0)*vx; %% position at time tt(im)
            y_ = y0 + (tt(im)-t0)*vy; 
%             if (x_<0 || x_>LX || y_<0 || y_>LY)
%                 %in=specRef(spec,vx,vy,x0,y0,x_,y_,LX,LY); %for chosen degree of specularity
%                 in=adiaRef(V(ind_mod),dORs,vx,vy,x0,y0,x_,y_,LX,LY);
%                 x_ = in(1);
%                 y_ = in(2);
%                 vx = in(3);
%                 vy = in(4);
%             endt
            indxx = floor((x_/LX)*NX)+1; %% index of the current spatial cell
            indxy = floor((y_/LY)*NY)+1; %
            if (indxx>0 && indxx<NX+1 && indxy>0 && indxy<NY+1)
                T(im,indxx,indxy) = T(im,indxx,indxy) + psign*Eeff/(C*Dx*Dy); % temperature
                driftvelx=driftvelx+vx;
                driftvely=driftvely+vy;
                %                 Qx(im,indxx,indxy) = Qx(im,indxx,indxy) + psign*Eeff*vx/(Dx*Dy); % heat flux
                %                 Qy(im,indxx,indxy) = Qy(im,indxx,indxy) + psign*Eeff*vy/(Dx*Dy);
            end
            im = im + 1;
        end
        % --------------------------------------- %
         
        % select post-collision mode
        ind_mod = select_mode(cumul_coll,Nmodes);
         
        % update particle parameters
        Rx = 2*rand()-1;
        Ry = 2*rand()-1;
        Rz = 2*rand()-1;
         
        ran = rand();
             
        cosa = Rx/sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
        cosb = Ry/sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
                             
        vx = V(ind_mod)*cosa; %! not be confused between V and vx
        vy = V(ind_mod)*cosb;
         
         
        x0 = x1;
        t0 = t1;
        y0 = y1;
                 
        if (t0>tt(Ntt) || x0<0 || x0>LX || y0<0 || y0>LY) % particle is terminated if it exits 
                                        % the system or if it overshoots the 
                                        % largest time of measurement
            finished = true;
        end
        meanphononlifetime= meanphononlifetime + Delta_t; %keep adding all the scattering times till the particle terminates
    end
    phlegm = 1;
     
end  
toc
driftvelx=driftvelx/N;
driftvely=driftvely/N;
meanphononlifetime=meanphononlifetime/N;
 
[xx,yy]=meshgrid(xx,yy);
 
 
 
% un-comment this for visualizing the solution
movieLength = 5; % duration of movie in seconds
figure();
 
tl = num2str(Tl);
tl = strcat(tl,'K');
td = num2str(Td);
td = strcat(td,'K');
tu = num2str(Tu);
tu = strcat(tu,'K');
tr = num2str(Tr);
tr = strcat(tr,'K');
tinit = strcat(num2str(Tinit),'K');
tinit = strcat('Initial temperature =  ',tinit);
 
unit_time = tmax/Ntt/(1e-9);
 
vid = VideoWriter('Isothermal__25by25_20mil.avi');
vid.FrameRate = 15;
open(vid);
xx = [0 100];
yy = [100 0];
clims = [300 310];
y=1e9*linspace(LY/2/NY,LY-LY/2/NY,NY);
for i=1:Ntt
   %fourierT[i]=fourierT(i);
   Tfp = reshape(T(i,:,:)+Teq,[NX,NY]);
   %contour(yy,xx,rot90(Tfp(:,:) + Teq));
   clf
   %subplot(2,2,1);
   imagesc(xx,yy,Tfp,clims);
   shading interp;
   colormap(jet);
   colorbar();
   set(gca,'XTick',[0:20:100]);
   set(gca,'Xticklabel',{'0','20','40','60','80','100'},'FontSize',10);
   set(gca,'YTick',[0:20:100],'FontSize',10);
   set(gca,'Yticklabel',{'100','80','60','40','20','0'},'FontSize',10);
   set(gcf,'Color',[1,1,1]);
   set(gca,'FontSize',10);
   xlabel('y (nm)','FontSize',13);
   ylabel('x (nm)','FontSize',13);
   set(gca,'XMinorTick','on','YMinorTick','on');
   
   %Create line
annotation('line',[0.792982456140351 0.792982456140351],...
    [0.10852380952381 0.926190476190476]);% % Create line
annotation('line',[0.130326355554293 0.130326355554293],...
    [0.10852380952381 0.926190476190476]);
% Create line
annotation('line',[0.293705797142614 0.293705797142614],...
    [0.10852380952381 0.926190476190476]);
% Create line
annotation('line',[0.458164569545829 0.458164569545829],...
    [0.10852380952381 0.926190476190476]);
% Create line
annotation('line',[0.628426165251233 0.628426165251233],...
    [0.10852380952381 0.926190476190476]);
   
   
   
   %title('Solving BTE in 2D using Deviational Monte Carlo Simulation')
%    annotation('textbox',[0.492187500000001 0.944017563117453 0.0312499999999994 0.0363216245883645],'String',tr,'EdgeColor','none')
%    annotation('textbox',[0.906547619047619 0.508730158730159 0.0904761904761904 0.0873015873015883],'String',tu,'EdgeColor','none')
%    annotation('textbox',[0.00873809523809524 0.538888888888889 0.0870952380952381 0.0603174603174608],'String',td,'EdgeColor','none')
%    annotation('textbox',[0.485895833333333 0.0197585071350165 0.0333750000000002 0.0384193194291987],'String',tl,'EdgeColor','none')
%    annotation('textbox',[0.00992857142857141 0.00714285714285714 0.347809523809524 0.0555555555555557],'String',tinit,'EdgeColor','none')
   %    axis([0 L Teq+min(min(T)) Teq+max(max(T))]);
   %    annotation('textbox',[0.199214285714286 0.945238095238095 0.235904761904762 0.0460317460317463],'String','Arush, Prof Dipanshu, Peraud','Edg','off');
   %    colorbar;
   time = i*unit_time;
   time = strcat(num2str(time),'ns');
   time = strcat('Time =',time);
   annotation('textbox',[0.667071428571429 0.93888888888889 0.341857142857143 0.0507936507936545],'String',time,'EdgeColor','none');
   %    annotation('textbox',[0.80397619047619 0.0579365079365081 0.139476190476191 0.0476190476190478],'String',{'(in nm)'},'EdgeColor','none');
   %    annotation('textbox',[0.0444523809523806 0.942063492063494 0.13947619047619 0.0476190476190479],'String',{'(in nm)'},'EdgeColor','none');
   %    annotation('textbox',[0.906547619047619 0.492857142857143 0.0904761904761902 0.0873015873015883],'String',{'(in K)'},'EdgeColor','none');
   %    annotation('textbox',[0.345833333333333 0.942063492063494 0.251190476190476 0.0460317460317461],'String','Top wall  = 300 K','EdgeColor','none');
   %    ylabel({'Left wall = 300 K'});
   %    xlabel({'Bottom wall = 310 K'});
    
    
    
%    clf
%    plot(y,Tfp(:,13),'Linewidth',2);
%    ylim([300 310]);
%    xlabel('x (nm)','FontSize',13);
%    ylabel('Temperature','FontSize',13);
%    time = i*unit_time;
%    time = strcat(num2str(time),'ns');
%    time = strcat('Time =',time)
%    annotation('textbox',[0.667071428571429 0.93888888888889 0.341857142857143 0.0507936507936545],'String',time,'EdgeColor','none');
 
         
 
   frame = getframe(gcf) ;
   writeVideo(vid, frame);
   pause(0.0001);
    
end
 



% close the writer object
close(vid);
a=0;
for i=1:Ntt
    for j=1:NX
        for k=1:NY
            if(T(i,j,k)+Teq>=310)
                a=a+T(i,j,k)+Teq-310;
            end
        end
    end
end
 
%Tfourier=fourierT();
%plot(xx,T(Ntt,:)+Teq-Tfourier,'Linewidth',2);