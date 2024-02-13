% PLEASE UNCOMMENT (ctrl + T) ONLY THE SCRIPT SECTION YOU WOULD LIKE TO RUN 
% AND LET THE OTHER SCRIPTS COMMENTED (ctrl + R)

%% script with constant arc voltage assumption

clear all;
clc;

global f omega R L Tau dt I1 I2 AngI;

%PLEASE CHOOSE N OF STEPS OF S.C. PROSPECTIVE CURRENT BETWEEN I1 AND I2
stepIrms=2;

E=400; %voltage in rms
%to get i max should be AngE=pi+AngZ so i have max transient
AngE=-pi/2; %condition for assimetrical prospective current e(0)=0

f=50;
omega=2*pi*f;
I1=1e3; %min short circuit prospective current
I2=100e3; %max short circuit prospective current
dt=2e-6;
t0=0.5e-3; %contact separation
te=1.8e-3;

%constant arc voltage
EA0=350; 


I_temp=linspace(I1,I2,stepIrms); %sequence of possible short circuit rms current

R_temp=[];
Z_temp=[];
L_temp=[];
I2t_temp=[];

for l=1:stepIrms
    I=I_temp(l);
    Z=E/I; %module
    cosFi=0.2;
    R=Z*cosFi;
    X=sqrt(Z.^2-R.^2);
    L=X/omega;
    Tau=L/R;
    AngZ=atan(X/R);
    AngZ_degrees=AngZ/pi*180;

    k=1.02+0.98*exp(-3*R/X);
    AngI=AngE-AngZ;
    AngI=0;
    AngI_degrees=AngI/pi*180;
    I0=EA0/R-I*sqrt(2)*cos(AngI)*exp(-t0/Tau);
   
    
 %if arc voltage is big enough, itot will go to zero immediately therefore
 %to find tf we have to search in 1/4 of period (b=0.25/f) otherwise 1/f
 
    a=t0;
    
if (I*sqrt(2).*cos(omega.*(t0*3)+AngI)-EA0./R+I0*exp(-(t0*3-t0)/Tau))<0;
    b=0.25/f;
else
    b=1/f;
end
    
while ((b-a)>dt)
    c=(a+b)/2;
    if ( i_total(I,EA0,I0,t0,a).*i_total(I,EA0,I0,t0,c) ) <0
            b=c;
    else
            a=c;
    end
end

tf=(a+b)/2;

%possible time ranges for plots
t_A=time_range(0,t0);
t_B=time_range(t0,tf);
t_B0=time_range(0,1/f);
t_C0=time_range(tf,1/f);
t_F=time_range(t0,1/f);
t_int=time_range(0,tf);


    figure (I_temp(l))
    plot(t_A, i1Gen(I,t_A), 'r', ...
            t_F, i1Gen(I,t_F), 'r--', ...
            t_B, i1(I, EA0, t_B), 'b', ...
            t_C0, i1(I, EA0, t_C0), 'b--', ...
            t_B, RL_transient(I0, t0, t_B), 'g', ...
            t_C0, RL_transient(I0, t0, t_C0), 'g--',...
            [t0 t0], [0 I0], 'y--', ...
            [t0 t0], [i1Gen(I, t0) i1(I, EA0, t0)], 'm--', ...
            t_B, -i1Gen(I,t_B)+i_total(I,EA0,I0,t0,t_B), 'c', ...
            t_C0, -i1Gen(I,t_C0)+i_total(I,EA0,I0,t0,t_C0), 'c--', ...
            [tf tf], [i1Gen(I,tf) 0], 'k');
        grid on

        hold on     
        plot(t_B, i_total(I,EA0,I0,t0,t_B), 'k','linewidth',2)
        hold on
        plot(t_C0, i_total(I,EA0,I0,t0,t_C0), 'k--','linewidth',2)
        hold off
            legend('i1Gen','i1Gen(t>t0)',...
            'i''(i1Gen-I1arc)','i''(t>tf)',...
            'i''''(RL transient)','i''''(t>tf)',...
            'I''''step','I''step',...
            'arc distortion(i_total-i1Gen)',...
            'arc distortion(t>tf)','tf setpoint','i total','i total(t>tf)')
 
% figure
% plot(t_F,i_total(I,EA0,I0,t0,t_F),'k')
% grid on
% legend(strcat ('itotal with Irms=',num2str(I_temp(l))))

 R_temp(l)=R;
 Z_temp(l)=Z;
 L_temp(l)=L;
 

 %find I2t
I2t=trapz(t_B,i_total(I,EA0,I0,t0,t_B).^2)+trapz(t_A,i1Gen(I, t_A).^2);
%alternative method
%  itot2=@(t) i_total(I,EA0,I0,t0,t).^2;
%  I2t=integral(itot2,0,tf)

%find I2t in function of Irms --> really time computing (approx. 5min)
 %syms Irms
 %itot2_Irms=i_total(Irms,EA0,I0,t0,t_int).^2;
 %I2t_Irms=int(itot2_Irms)

I2t_temp(l)=I2t;  
 
%arc energy 
Earc=abs(EA0*trapz(t_B,i_total(I,EA0,I0,t0,t_B)));

%alternative method
% itot=@(t) i_total(I,EA0,I0,t0,t);
% Earc=EA0*integral(itot,t0,tf)

%find Earc in function of Irms --> really time computing (approx. 2min)
% syms Irms
% itot_Irms=i_total(Irms,EA0,I0,t0,t_int);
% Earc=EA0*int(itot_Irms)

Earc_temp(l)=Earc; 


%find peak value

I_peak=max(i_total(I,EA0,I0,t0,t_B));

%alternative methods
% t_peak=t0;
% I_peak=i_total(I,EA0,I0,t0,t_peak)

%peak of not limited current
%I_peak1=-findpeaks(-i_total(I,EA0,I0,t0,t_B0)) 

I_peak_mod=abs(I_peak);
I_peak_temp(l)=I_peak_mod; 

I_peak_prospective_temp(l)=k*sqrt(2)*I;

end


%collect all data: vector of length stepIrms
R_values=R_temp;
Z_values=Z_temp;
L_values=L_temp;
I2t_values=I2t_temp;
Earc_values=Earc_temp;
I_peak_values=I_peak_temp;
I_peak_prospective=I_peak_prospective_temp; 

figure
loglog(I_temp,I_peak_values,'m')
grid on
xlabel('Irms values [A]') 
ylabel('I peak values [A]') 
hold on
loglog(I_temp,I_peak_prospective,'r')
hold off
legend('I peak(Irms)','I peak prospective(Irms)')




figure
loglog(I_temp,I2t_values,'y')
grid on
xlabel('Irms values [A]') 
ylabel('I2t values [J/ohm]')
legend('I2t values(Irms)')

figure
loglog(I_temp,Earc_values,'r')
grid on
xlabel('Irms values [A]') 
ylabel('Earc values [V]')
legend('Earc values(Irms)')




function y=time_range(t1,t2)
% this function is uset to generate time samples 
% which are equally spaced, with an approximate
% sample interval equal to dt, which span the 
% range from t1 to t2 in such a way that values 
% t1 and t2 are both included in the results set
global dt
y=t1:((t2-t1)/ceil((t2-t1)/dt)):t2;
end


function y=i1Gen(I, t)
global omega Tau AngI
y=I*sqrt(2).*cos(omega.*t+AngI)-I*sqrt(2)*cos(AngI)*exp(-t/Tau);
end

function y=i1Gen1(I, t)
global omega Tau AngI
y=I*sqrt(2).*cos(omega.*t+AngI);
end

function y=i1Arc(EA, t)
global R
y=ones(size(t))*EA./R;
end

function y=i1(I, EA, t)
y=i1Gen1(I,t)-i1Arc(EA,t);
end

function y=i_total(I,EA,I0,t0,t)
y=i1(I, EA, t)+RL_transient(I0, t0, t);
end


function y=RL_transient(I0, t0, t)
global Tau
y=I0*exp(-(t-t0)/Tau);
end


 
%% script to see the comparison with a HV non limiting circuit breaker

% select the script and type ctrl + t to uncomment
% remember to comment the scripts you don't need (ctrl + r)

% clear all;
% clc;
% 
% global f omega R L Tau dt I1 I2 AngI;
% 
% %PLEASE CHOOSE N OF STEPS OF S.C. PROSPECTIVE CURRENT BETWEEN I1 AND I2
% stepIrms=2;
% 
% E=400; %voltage in rms
% AngE=-pi/2; %condition for assimetrical prospective current e(0)=0
% f=50;
% omega=2*pi*f;
% I1=1e3; %min short circuit prospective current
% I2=100e3; %max short circuit prospective current
% dt=2e-6;
% t0=0.5e-3; %contact separation
% te=1.8e-3;
% 
% %constant arc voltage
% EA0=10; 
% 
% 
% I_temp=linspace(I1,I2,stepIrms); %sequence of possible short circuit rms current
% 
% R_temp=[];
% Z_temp=[];
% L_temp=[];
% I2t_temp=[];
% 
% for l=1:stepIrms
%     I=I_temp(l);
%     Z=E/I; %module
%     cosFi=0.2;
%     R=Z*cosFi;
%     X=sqrt(Z.^2-R.^2);
%     L=X/omega;
%     Tau=L/R;
%     AngZ=atan(X/R);
%     AngZ_degrees=AngZ/pi*180;
%     k=1.02+0.98*exp(-3*R/X);
%     AngI=AngE-AngZ;
%     I0=EA0/R-I*sqrt(2)*cos(AngI)*exp(-t0/Tau);
%    
%     
%  %if arc voltage is big enough, itot will go to zero immediately therefore
%  %to find tf we have to search in 1/4 of period (b=0.25/f) otherwise 1/f
%  
%     a=t0;
%     
% if (I*sqrt(2).*cos(omega.*(t0*4)+AngI)-EA0./R+I0*exp(-(t0*4-t0)/Tau))<0;
%     b=0.25/f;
% else
%     b=1/f;
% end
%     
% while ((b-a)>dt)
%     c=(a+b)/2;
%     if ( i_total(I,EA0,I0,t0,a).*i_total(I,EA0,I0,t0,c) ) <0
%             b=c;
%     else
%             a=c;
%     end
% end
% 
% tf=(a+b)/2;
% 
% %possible time ranges for plots
% t_A=time_range(0,t0);
% t_B=time_range(t0,tf);
% t_B0=time_range(0,1/f);
% t_C0=time_range(tf,1/f);
% t_F=time_range(t0,1/f);
% t_int=time_range(0,tf);
% 
% 
%     figure (I_temp(l))
%     plot(t_A, i1Gen(I,t_A), 'r', ...
%             t_F, i1Gen(I,t_F), 'r--', ...
%             t_B, i1(I, EA0, t_B), 'b', ...
%             t_C0, i1(I, EA0, t_C0), 'b--', ...
%             t_B, RL_transient(I0, t0, t_B), 'g', ...
%             t_C0, RL_transient(I0, t0, t_C0), 'g--',...
%             [t0 t0], [0 I0], 'y--', ...
%             [t0 t0], [i1Gen(I, t0) i1(I, EA0, t0)], 'm--', ...
%             t_B, -i1Gen(I,t_B)+i_total(I,EA0,I0,t0,t_B), 'c', ...
%             t_C0, -i1Gen(I,t_C0)+i_total(I,EA0,I0,t0,t_C0), 'c--', ...
%             [tf tf], [i1Gen(I,tf) 0], 'k');
%         grid on
% 
%         hold on     
%         plot(t_B, i_total(I,EA0,I0,t0,t_B), 'k','linewidth',2)
%         hold on
%         plot(t_C0, i_total(I,EA0,I0,t0,t_C0), 'k--','linewidth',2)
%         hold off
%             legend('i1Gen','i1Gen(t>t0)',...
%             'i''(i1Gen-I1arc)','i''(t>tf)',...
%             'i''''(RL transient)','i''''(t>tf)',...
%             'I''''step','I''step',...
%             'arc distortion(i_total-i1Gen)',...
%             'arc distortion(t>tf)','tf setpoint','i total','i total(t>tf)')
%  
% % figure
% % plot(t_F,i_total(I,EA0,I0,t0,t_F),'k')
% % grid on
% % legend(strcat ('itotal with Irms=',num2str(I_temp(l))))
% 
%  R_temp(l)=R;
%  Z_temp(l)=Z;
%  L_temp(l)=L;
%  
% 
%  %find I2t
% I2t=trapz(t_int,i_total(I,EA0,I0,t0,t_int).^2);
% %alternative method
% %  itot2=@(t) i_total(I,EA0,I0,t0,t).^2;
% %  I2t=integral(itot2,0,tf)
% 
% %find I2t in function of Irms --> really time computing (approx. 5min)
%  %syms Irms
%  %itot2_Irms=i_total(Irms,EA0,I0,t0,t_int).^2;
%  %I2t_Irms=int(itot2_Irms)
% 
% I2t_temp(l)=I2t;  
%  
% %arc energy 
% Earc=abs(EA0*trapz(t_B,i_total(I,EA0,I0,t0,t_B)));
% 
% %alternative method
% % itot=@(t) i_total(I,EA0,I0,t0,t);
% % Earc=EA0*integral(itot,t0,tf)
% 
% %find Earc in function of Irms --> really time computing (approx. 2min)
% % syms Irms
% % itot_Irms=i_total(Irms,EA0,I0,t0,t_int);
% % Earc=EA0*int(itot_Irms)
% 
% Earc_temp(l)=Earc; 
% 
% 
% %find peak value
% 
% I_peak=max(i_total(I,EA0,I0,t0,t_B));
% 
% %alternative methods
% % t_peak=t0;
% % I_peak=i_total(I,EA0,I0,t0,t_peak)
% 
% %peak of not limited current
% %I_peak1=-findpeaks(-i_total(I,EA0,I0,t0,t_B0)) 
% 
% I_peak_mod=abs(I_peak);
% I_peak_temp(l)=I_peak_mod; 
% 
% I_peak_prospective_temp(l)=k*sqrt(2)*I;
% 
% end
% 
% 
% %collect all data: vector of length stepIrms
% R_values=R_temp;
% Z_values=Z_temp;
% L_values=L_temp;
% I2t_values=I2t_temp;
% Earc_values=Earc_temp;
% I_peak_values=I_peak_temp;
% I_peak_prospective=I_peak_prospective_temp; 
% 
% figure
% loglog(I_temp,I_peak_values,'m')
% grid on
% xlabel('Irms values [A]') 
% ylabel('I peak values [A]') 
% hold on
% loglog(I_temp,I_peak_prospective,'r')
% hold off
% legend('I peak(Irms)','I peak prospective(Irms)')
% 
% 
% 
% 
% figure
% loglog(I_temp,I2t_values,'y')
% grid on
% xlabel('Irms values [A]') 
% ylabel('I2t values [J/ohm]')
% legend('I2t values(Irms)')
% 
% figure
% loglog(I_temp,Earc_values,'r')
% grid on
% xlabel('Irms values [A]') 
% ylabel('Earc values [V]')
% legend('Earc values(Irms)')
% 
% 
% 
% 
% function y=time_range(t1,t2)
% % this function is uset to generate time samples 
% % which are equally spaced, with an approximate
% % sample interval equal to dt, which span the 
% % range from t1 to t2 in such a way that values 
% % t1 and t2 are both included in the results set
% global dt
% y=t1:((t2-t1)/ceil((t2-t1)/dt)):t2;
% end
% 
% 
% function y=i1Gen(I, t)
% global omega Tau AngI
% y=I*sqrt(2).*cos(omega.*t+AngI)-I*sqrt(2)*cos(AngI)*exp(-t/Tau);
% end
% 
% function y=i1Gen1(I, t)
% global omega Tau AngI
% y=I*sqrt(2).*cos(omega.*t+AngI);
% end
% 
% function y=i1Arc(EA, t)
% global R
% y=ones(size(t))*EA./R;
% end
% 
% function y=i1(I, EA, t)
% y=i1Gen1(I,t)-i1Arc(EA,t);
% end
% 
% function y=i_total(I,EA,I0,t0,t)
% y=i1(I, EA, t)+RL_transient(I0, t0, t);
% end
% 
% 
% function y=RL_transient(I0, t0, t)
% global Tau
% y=I0*exp(-(t-t0)/Tau);
% end


