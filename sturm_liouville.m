                                                                     
                                             
function precession_orbits
%
% orbits - Kelperian except that L is not conserved - driving harmonic
% force
%
global pert m AngL E x z
syms r(t) th(t) w yy L 
%
% first see if you have a symbolic solution
%

%
% no - so setup for numeric solution
%
irun = 1;
iloop = 0;
%
while irun > 0
    kk = menu('Pick Another Strength of Perturbation?','Yes','No');
    if kk == 2
        irun = -1;
        break
    end
    if kk == 1
        %
        a = input('Enter Initial radius (~ 4): ');
        m = input('Enter angular frequency of force: ');
        LL = input('Enter Initial Angular Velocity 1/a^3/2 Circular: ');
       % fprintf('Angular velocity < circular is not stable and hits r = 0 \n')
        %
        % circular orbits
        %
        Lc = sqrt(a);
        vc = 1.0 ./(a .^1.5);  % angular velocity 1/a^2/3
        tauc = (0.5 .*pi) ./0.01 ;        
        %
        fprintf('Circular Orbits for Initial Radius a \n')
        
        %
        tspan = linspace(0, 50.0 .*tauc ,10000); %linspace(0,10.0.*tauc,100000); % time frame - 4 circular periods
        %
        [t,y] = ode45(@mech,tspan,[0 a LL 0 0]);   % vr, r vth, th - initial conditions 
        %for unperturbed ellipse
        %
   %     iloop = iloop + 1;
   %     figure(iloop)
    %    plot(t,y(:,2))
  %      title('r vs t')
   %     xlabel('t')
    %    ylabel('r')
        %
 %       iloop = iloop + 1;
%        figure(iloop)
 %       plot(t,y(:,4))
  %      title('\theta vs t')
%        xlabel('t')
 %       ylabel('\theta')
        %
        xx = y(:,2) .*cos(y(:,4));
        yy = y(:,2) .*sin(y(:,4));
        iloop = iloop + 1;
        figure(iloop)
        plot(xx,yy,'b',0,0,'r*')
        title('x vs y')
        xlabel('x')
        ylabel('y')
        %
        AngL = y(:,2) .*y(:,2) .*y(:,3);
        iloop = iloop +1;
%        figure(iloop)
%        plot(t,y(:,2) .*y(:,2) .*y(:,3))
        title('angular momentum vs t')
        xlabel('t')
        ylabel('angular momentum')
        %
        V = y(:,2) .*y(:,3);
        iloop = iloop+1;
%        figure(iloop)
%        plot(t,y(:,2) .*y(:,3))
        title('velocity vs t')
        xlabel('t')
        ylabel('velocity')
        %
   
         iloop = iloop +1;
     %    figure(iloop)
      %  plot(y(:,4),y(:,2))
        xlabel('t')
        ylabel('r vs theta')
        %
        iloop=iloop +1;
%        figure(iloop)
%        plot(t,y(:,3))
%        xlabel('t')
%        ylabel('angular velocity')
        %
 
        
        
        
        % iloop = iloop +1;
       % figure(iloop)
       % AngL = y(:,2) .*y(:,2) .*y(:,3);
       % D2A = diff(AngL,t,2);
       % x = D2A + y(:,3) .*y(:,3) .*AngL;
       % plot(t,x);   
        % 
      %  iloop = iloop +1;
       % figure(iloop)
        %z = y(:,3).*y(:,2) .*y(:,3) .*y(:,2).*(1.0 ./y(:,2) .^2) +y(:,2) .*y(:,2) .*y(:,3) .*0.5 .*m .*m .*cos(m .*t);
        %x =   y(:,3) .*y(:,3) .*AngL ;
        %plot(t,z,'r',t,x,'b')
    %
    %   iloop = iloop +1;
       %figure(iloop)
       %plot(t,y(:,4))
    end
end
%---------------------------------------------------------------------------------
function dy = mech(t,y)
%
global  m AngL V E
dy = zeros(5,1);
% y1 = drdt, 2 = r, 4 = dthetadt, 5 = theta 3 = domegadt
%
dy(1) = (y(3) .^2) .*y(2) - 1.0 ./y(2) .^2;
%dy(1) =  -(0.05+m.*sin(m.*t)).*y(3).* y(2)   + m .*m .*(0 .*cos(m .*t) - 0 .*cos(2 .*m .*t)); 
%dy(1) = (y(3) .^2) .*y(2) - 1.0 ./y(2) .^2 + 1. * m .*m .*(2.5 .*cos(m .*t) - 0.5.*cos(2 .*m .*t)); 
dy(3) = ((-2 .*y(1) .*y(3) .*y(2) + 1 .*m .*0.5 .*sin( m .*t)) ./(y(2) .^2));
%dy(3) = ((-2 .*y(1) .*y(3) .*y(2)+(0.5+m.* sin(m.*t)).*y(1).*y(2)) ./(y(2) .^2)); 
%dy(4) = (1 .*m .*m .*sin (m .*t)-y(4)-2 .*y(1) .*y(1) .*y(4)-2 .*y(2) .*((y(4) .^2) .*y(2) - 1.0 ./y(2) .^2) .*y(4)-4 .*y(2) .*y(1) .*y(4)) ./(y(2) .^2);
dy(2) = y(1);
dy(4) = y(3);
%y(2) = y(3) .^(-1.5); 



