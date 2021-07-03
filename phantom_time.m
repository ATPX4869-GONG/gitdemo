
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 9);

% 指定范围和分隔符
opts.DataLines = [7, Inf];
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = ["ChannelTitle1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 导入数据
tbl = readtable("F:\研究\20210610\result\20210610_rE03_sE20_2.txt", opts);
[extention,contraction]=setsudotime('re0.3se2.0_2.csv');

ChannelTitle = tbl.ChannelTitle1;
Fxv = tbl.VarName2;
Fyv = tbl.VarName3;
Fzv = tbl.VarName4;
Mxv = tbl.VarName5;
Myv = tbl.VarName6;
anglev = tbl.VarName7;
Mzv = tbl.VarName8;

clear opts tbl



%% 電圧データから角度（rad）／力(N)にする
start_time = 1;
over_time = 99000;

Data_M = [Fxv';Fyv';Fzv';Mxv';Myv';Mzv'];

Hoseix=[1.190113405 -1.02086084 -0.632202425 43.85204566  -43.50735015 0.121485036];
Hoseiy=[0.510630059 0.346199264 -1.388111543 25.13755688 24.76078099 -47.62953003];

FHx = (Hoseix * Data_M)';
FHy = (Hoseiy * Data_M)' ; 
Fn=sqrt(FHx.^2+FHy.^2);

% F4(1:200000)= 4;
% F = F4(start_time:over_time)';   %使F=4N  必要时需要改回
time=ChannelTitle(start_time:over_time);

% d_0=2.031;
% d_90=2.712; 
% angleradinx= (anglev-d_0)/((d_90-d_0)/90)/180*pi;
% angleradin = angleradinx(start_time:over_time);

% 
d_180=2.031;% calibration angle 135 degree
d_90=2.712; %90degree 
zerov=d_180-180*(d_180-d_90)/90;
angleradinx=(anglev-zerov)/((d_180-d_90)/90)/180*pi;
angleradin = angleradinx(start_time:over_time);



%% filter/data election

Fn = Fn(start_time:over_time);
[Fn,time_F] = resample(Fn,time,100);
[angleradin,time_A] = resample(angleradin,time,100);
[b,a] = butter(1,50/500,'low');
Fbutter=filtfilt(b,a,Fn);
[x,y] = butter(1,50/500,'low');
anglebutter=filtfilt(x,y,angleradin);
anglezuihou=filtfilt(x,y,angleradin);


%% get the data from phantom
Mix=[extention,contraction];
angle_ex(1:length(extention))=anglezuihou(extention);
angle_flx(1:length(contraction))=anglezuihou(contraction);




%%  calculation
str='Mix';
str1=eval(str);
% Stiffness=zeros(length(str1),1);
% Viscosity=zeros(length(str1),1)

Stiffness_100ms=zeros(length(str1),1);
Viscosity_100ms=zeros(length(str1),1);
Stiffness_150ms=zeros(length(str1),1);
Viscosity_150ms=zeros(length(str1),1);
Stiffness_200ms=zeros(length(str1),1);
Viscosity_200ms=zeros(length(str1),1);

plus1=0;
for plus2=(10+plus1):5:(20+plus1)
   for counterCal=1:length(str1)
   %データの開始時刻を定義
   forcec=Fbutter(str1(counterCal)+plus1:str1(counterCal)+plus2)-Fbutter(str1(counterCal)+plus1);
   if counterCal > length(extention)
       forcec= -forcec;
   end
   timec=0:0.01:(plus2-plus1)*0.01;
   timec=timec';
   anglec=anglezuihou(str1(counterCal)+plus1:str1(counterCal)+plus2)-anglezuihou(str1(counterCal)+plus1);
%--------------------------------------------------------------------


    f=fittype('poly4');                             %(fitting)
    fit1=fit(timec,anglec,f);                        %time -angle fitting curve
    fit2=fit(timec,forcec,f);                       %time-force fitting curve
    [a1,a2]=differentiate(fit1,timec);              %time-angle 1st 2st differentiation
    a0=polyfit(timec,anglec,4);                      %coefficient of time- angle fitting 
    f0=polyfit(timec,forcec,4);                     %coefficient of time- force fitting
    anglekin=polyval(a0,timec);                            %fitted value of time-angle 
    Forcekin=polyval(f0,timec);                           %fitted value of time-force  

    weight=67;                                           %the weight of body
    aw=weight*0.022+0.109;                               %the weight of arm
    arml=0.32;                                        %the length of hand
    moment=aw*arml*arml/3;                                   %Moment of Inertia
    Torque=Forcekin*arml;                                         %torque(fiitting)
    P1=Torque-(moment*a2);                                    %torque-inertia*accreclation
    T1=Torque;                                           %torque
    A1=a2;                                          %angular acceleration Matrix
    M=[anglekin a1];                                       %angle and angular velocity
    N=pinv(M);                                      % pseudo-inverse matrix
    O=P1;                                           %torque-inertia*accreclation 
    A01=N*O;                                         %calculate the stiffness and viscocity Matrix

        if plus2==(10+plus1)
        Stiffness_100ms(counterCal)=abs(A01(1));
        Viscosity_100ms(counterCal)=abs(A01(2));
        elseif plus2==(15+plus1)
            Stiffness_150ms(counterCal)=abs(A01(1));
            Viscosity_150ms(counterCal)=abs(A01(2));
        elseif plus2==(20+plus1)
            Stiffness_200ms(counterCal)=abs(A01(1));
            Viscosity_200ms(counterCal)=abs(A01(2));
        end
    end
end


%% Fig output  
% for mode=100:50:200
% mode=150;
% name1=strcat('re0.7se0.5_',num2str(mode),'_ms(2)','.fig');%change the name to store the different filename
% name2=strcat('Stiffness_',num2str(mode),'ms');
% name3=strcat('Viscosity_',num2str(mode),'ms');
% Figuresti=eval(name2);
% Figurevis=eval(name3);

figure(1)
subplot(2,1,1)
plot(time_F,anglebutter)
hold on
plot(extention/100,anglebutter(extention),'rx')
hold on
plot(contraction/100,anglebutter(contraction),'gx')
hold on 
%xlabel('time/s') 
ylabel('Angle/rad')
%title(name1,'Interpreter','none')

subplot(2,1,2)
plot(time_F,Fbutter)
hold on
plot(extention/100,Fbutter(extention),'rx') 
hold on 
plot(contraction/100,Fbutter(contraction),'gx') 

 hold on 
 ylabel('Force/N')

 
% subplot(4,1,3)
% plot(extention/100,Stiffness_150ms(1:length(extention)),'ro') 
% hold on 
% plot(contraction/100,Stiffness_150ms(length(extention)+1:length(Mix)),'go') 
% hold on
% ylabel('Stiffness/Nm/Rad ') 
% legend('extensor','flexsor')
% 
% %  
% subplot(4,1,4)
% plot(extention/100,Viscosity_150ms(1:length(extention)),'ro') 
% hold on 
% plot(contraction/100,Viscosity_150ms(length(extention)+1:length(Mix)),'go') 
% hold on 
% ylabel('Visc/Nms/rad')


figure(2)
subplot(4,1,1)
plot(fit1,timec,anglec) % cfit plot method
xlabel('Time/s') 
ylabel('Angle/rad')

subplot(4,1,2)
plot(timec,a1,'m') % double plot method
grid on
legend('1st derivative of time-angle')

subplot(4,1,3)
plot(timec,a2,'c') % double plot method
grid on
legend('2nd derivative of time-angle')

subplot(4,1,4)
plot(fit2,timec,forcec)
xlabel('Time/s') 
ylabel('Force/N')


%% save program
% zhuangtai=str;
% qianzhui1='Derivative_';
% qianzhui2='status_';
% bianhao=num2str(counterCal);
% name1=[zhuangtai,qianzhui1,bianhao];
% name2=[zhuangtai,qianzhui2,bianhao];
% name3='all';
% 
% saveas(1,name1,'jpeg')
% saveas(2,name2,'jpeg')
% saveas(3,name3,'jpeg')
% close all


%% Outlier
% for i = 1:length(angle_ex)
%     if abs(mean(angle_ex)-angle_ex(i))>0.2
%         angle_ex(i)= 0;
%          Stiffness_150ms(i) = 0;
%          Viscosity_150ms(i) = 0;
%     end
% end
% 
% for i = 1:length(angle_flx)
%     if abs(mean(angle_flx)-angle_flx(i))>0.2
%         angle_flx(i)= 0;
%          Stiffness_150ms(length(angle_ex)+i) = 0;
%          Viscosity_150ms(length(angle_ex)+i) = 0;
%     end
% end
%         

Stie_150 = Stiffness_150ms(1:length(angle_ex))';
Stic_150 = Stiffness_150ms(length(angle_ex)+1:end)';
Ve_150 = Viscosity_150ms(1:length(angle_ex))';
Vc_150 = Viscosity_150ms(length(angle_ex)+1:end)';

 
% %根据角度平稳性去除离群值 
% checke = isoutlier(angle_ex,'median','ThresholdFactor',0.8);
% checkf = isoutlier(angle_flx,'median','ThresholdFactor',0.8);
% 
%  Stie_150(checke==1)=[];    
%  Stic_150(checkf==1)=[];   
% Ve_150(checke==1)=[];  
% Vc_150(checkf==1)=[];  
% angle_ex(checke==1)=[]; 
% angle_flx(checkf==1)=[];

%根据刚性值去除离群值
% checkse = isoutlier(Stie_150,'median','ThresholdFactor',0.7);
% checksc = isoutlier(Stic_150,'median','ThresholdFactor',0.7);
%  
% Stie_150(checkse==1)=[];      
% Stic_150(checksc==1)=[]; 
% Ve_150(checkse==1)=[];
% Vc_150(checksc==1)=[]; 
% angle_ex(checkse==1)=[]; 
% angle_flx(checksc==1)=[]; 
%wdnmd dafen