% THE VOLTAGE RECORDINGS OF THE BIOLOGICAL NEURONS
%MelibeData1=xlsread('MelibeData1.xlsx');
time1 = MelibeData1(:,1)-MelibeData1(1,1);
si1L = MelibeData1(:,2);
Iapp = MelibeData1(:,3);
si1R = MelibeData1(:,4);
si2L = MelibeData1(:,5);
si2R = MelibeData1(:,6);

g13=0.35;
alphae=.0061;
betae=.002;

noise=0.4;

Ca_shift1=  -85.8;
Ca_shift2 = -100.6;
x_shift =  -3.;

alpha = 0.03;
beta =  0.002;
%THE CONDUCTANCES OF THE BLENDED SYNAPSES
g12=0.015;
g21=0.015;
gelec=0.0000;
taum=600;

t1=10*1000;
t2=230*1000;
      
% H-current        
    gh    = 0.0001;
    Vhh   = -53;

% Intergration    
t_final=time1(end)-24;   
step=1;
step=1000*time1(end)/length(time1);

si3L = si1R;
si3R = si1L;
s3R_Array = zeros(length(time1),1);
s3L_Array = zeros(length(time1),1);
peaks_si3L = findpeaks(si1L,'MinPeakHeight',10, 'MinPeakDistance',5);
peaks_si3R = findpeaks(si1R,'MinPeakHeight',10, 'MinPeakDistance',5);    

% Initial values the V1, vv1,... represent the variables of model-neuron 1
% with V2, vv2,... represent the varibels of model-neuron 2
V1= 0; V2= 0; Ca1=.5; Ca2=.5; h1 =0; h2 =0; n1 =0; n2 =0; x1 =0.8; x2 =0.85;
y1 =0; y2 =0; s1 =0; s2 =0; m1 =0; m2 =0;i=1; tt=1; seL=0; seR=0; s3L=0;s3R=0;

clear ss1 mm1  clear time vv1 vv2 seeL seeR ss2 mm2 clear Caa1 Caa2 xx1 xx2

time=zeros(length(time1),1);vv1=time;vv2=time;ss1=time;ss2=time;mm1=time;mm2=time;Caa1=time;
Caa2=time;xx1=time;xx2=time;seeL=time; seeR=time;siiL=time; siiR=time;

tic
for i=1:length(time1) 
   
%Blended synapse model 
seL=seL+step*(alphae*(1-seL)/(1+exp(-10*(si1L(i)+20)))-betae*seL);
seeL(i)=seL;
seR=seR+step*(alphae*(1-seR)/(1+exp(-10*(si1R(i)+20)))-betae*seR);
seeR(i)=seR;

%Si3 neuron synapse
s3L=s3L+step*(alpha*(1-s3L)/(1+exp(-10*(si3L(i)+20)))-beta*s3L);
s3R=s3R+step*(alpha*(1-s3R)/(1+exp(-10*(si3R(i)+20)))-beta*s3R);
s3R_Array(i)=s3R;
s3L_Array(i)=s3L;

% THE MODEL NEURONS WITH BLENDED SYNAPTIC CURRENTS (can be seen as the terms with g12 and g21)
V1 =V1 +step*(4*((0.1*(50-(127*V1/105+8265/105))/(exp((50 - ...
    (127*V1/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V1/105 + 8265/105))/(exp((50 - (127*V1/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V1/105 + 8265/105))/18))))^3*h1*(30 - V1) + 0.3*n1^4*(-75 - V1)+0.01*x1*(30-V1) +0.03*Ca1/(.5 + Ca1)*(-75 - V1)...
    +0.003*(-40 - V1)   +gh*((1/(1+exp(-(V1+63)/7.8)))^3)*y1*(+120-V1)...
    -g13*(V1+65)*seL-g21*(V1+80)*s2*m2+gelec*(si1L(i)-V1) +noise*(rand-0.5) );
V2= V2 +step*(4*((0.1*(50-(127*V2/105+8265/105))/(exp((50 - ...
    (127*V2/105 + 8265/105))/10) - 1))/((0.1*(50 - (127*V2/105 + 8265/105))/(exp((50 - (127*V2/105 + 8265/105))/10) - 1))+...
    (4*exp((25 - (127*V2/105 + 8265/105))/18))))^3*h2*(30 - V2) + 0.3*n2^4*(-75 - V2)+0.01*x2*(30-V2) +0.03*Ca2/(.5 + Ca2)*(-75 - V2)...
    +0.003*(-40 - V2) +gh*((1/(1+exp(-(V2+63)/7.8)))^3)*y2*(-V2+120)...
    -g13*(V2+60)*seR-g12*(V2+80)*s1^2*m1+gelec*(si1R(i)-V2)+noise*(rand-0.5));

Ca1=Ca1+step*(0.0001*(0.0085*x1*(140-V1+Ca_shift1)-Ca1) );
Ca2=Ca2+step*(0.0001*(0.0085*x2*(140-V2+Ca_shift2)-Ca2));
x1 =x1+step*(((1/(exp(0.15*(-V1-50+x_shift))+1))-x1)/235);
x2 =x2+step*(((1/(exp(0.15*(-V2-50+x_shift))+1))-x2)/235);
h1 =h1+step*(((1-h1)*(0.07*exp((25 - (127*V1/105 + 8265/105))/20))-h1*(1.0/(1 + exp((55 - (127*V1/105 + 8265/105))/10))))/12.5);
h2 =h2+step*(((1-h2)*(0.07*exp((25 - (127*V2/105 + 8265/105))/20))-h2*(1.0/(1 + exp((55 - (127*V2/105 + 8265/105))/10))))/12.5);
n1 =n1+step*(((1-n1)*(0.01*(55 - (127*V1/105 + 8265/105))/(exp((55 - (127*V1/105 + 8265/105))/10) - 1))-n1*(0.125*exp((45 - (127*V1/105 + 8265/105))/80)))/12.5 );
n2 =n2+step*(((1-n2)*(0.01*(55 - (127*V2/105 + 8265/105))/(exp((55 - (127*V2/105 + 8265/105))/10) - 1))-n2*(0.125*exp((45 - (127*V2/105 + 8265/105))/80)))/12.5);
y1 =y1+step*(.5*((1/(1+exp(10*(V1-Vhh))))-y1)/(7.1+10.4/(1+exp((V1+68)/2.2))));
y2 =y2+step*(.5*((1/(1+exp(10*(V2-Vhh))))-y2)/(7.1+10.4/(1+exp((V2+68)/2.2))));

s1 =s1+step*(alpha*(1-s1)/(1+exp(-10*(V1+30)))-beta*s1 +0.01*noise*(rand-0.5));
s2 =s2+step*(alpha*(1-s2)/(1+exp(-10*(V2+30)))-beta*s2 +0.012*noise*(rand-0.5));
m1 =m1+step*((1/(1+exp(-(V1+30)))-m1)/taum);
m2 =m2+step*((1/(1+exp(-(V2+30)))-m2)/taum);

time(i)=(i-1)*step;
vv1(i)=V1;
vv2(i)=V2;
ss1(i)=s1;
ss2(i)=s2;
mm1(i)=m1;
mm2(i)=m2;
Caa1(i)=Ca1;
Caa2(i)=Ca2;
xx1(i)=x1;
xx2(i)=x2;
end  
toc    