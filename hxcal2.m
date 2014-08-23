%%hxcal2.m 2007-08-20: result as log(k_ex) instead of deltaG
%%hxcal.m 2007-07-10: to calculate(normalize)/fit/plot HX rate from gCOSY NMR data (vol measured in Felix):  
%%current version updated at 2007-07-25 

disp('check! first import corresponding "data", "timeLog", "peakName", and "kc"[129X1]')
N=input('input number of time points (.fid files): '); %number of time points (.fid files) (21)
sizer=size(timeLog);
if (sizer(1)~=N*2) error('wrong input of time points or timeLog'); end
M=input('input number of peaks: '); %number of all selected peaks (90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculate intensity; subtract blank; and reference peak correction:
vol=data(:,2:N+1);
area=data(:,1);

for i=1:N
    inten(:,i)=vol(:,i)./area;
end

disp('check! currently peak(row) #63, 64, 65, 68, 69 as blank')
for i=1:N
    blank(i)=(inten(63,i)+inten(64,i)+inten(65,i)+inten(68,i)+inten(69,i))/5;
end

for i=1:N
    intenSub(:,i)=inten(:,i)-blank(i);
end

disp('check! currently peak(row) #1, 2, 62 as non-exchange reference')
for i=1:N
    ref(i)=(intenSub(1,i)/intenSub(1,1)+intenSub(2,i)/intenSub(2,1)+intenSub(62,i)/intenSub(62,1))/3;
end

% intenCorr=zeros(M,N);
for i=1:N
    intenCorr(:,i)=intenSub(:,i)/ref(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cal time points from timeLog:
t1=int2str(timeLog); 
sizer=size(timeLog);
for i=1:sizer(1)
    dd(i)=str2double(t1(i,1:2)); if(t1(i,1:2)==' ') dd(i)=0; end
    hh(i)=str2double(t1(i,5:6)); if(t1(i,5:6)==' ') hh(i)=0; end
    mm(i)=str2double(t1(i,7:8)); if(t1(i,7:8)==' ') mm(i)=0; end
    ss(i)=str2double(t1(i,9:10)); if(t1(i,9:10)==' ') ss(i)=0; end
end

time0=input('input the HX time before NMR starting: [mm ss] '); %[1X2]

%time of each .fid file: (nt=2 => intv=0.3819 hr)
intv= 24*(dd(2)-dd(1)) + (hh(2)-hh(1)) + (mm(2)-mm(1))/60 + (ss(2)-ss(1))/3600;

time(1)= time0(1)/60 + time0(2)/3600 + intv/2;
for i=2:N
    time(i)=time(1) - intv + 24*(dd(2*i)-dd(1)) + (hh(2*i)-hh(1)) + (mm(2*i)-mm(1))/60 + (ss(2*i)-ss(1))/3600;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%single-exponential fitting and free energy calculation:
for i=1:M
    if(intenCorr(i,1)/intenCorr(i,N)<1) 
        iniK=0;
    else
        iniK=log(intenCorr(i,1)/intenCorr(i,N))/((time(N)-time(1))*3600);
    end
    iniPara=[intenCorr(i,1); iniK];  
    options = optimset('TolX', 1e-9, 'TolFun', 1e-15);
    [fitPara,r1,r2,exitFlag,output]=lsqnonlin(@hxfit, iniPara, [0;0],[], options, time, i, intenCorr);
    disp(i);
    fitA(i)=fitPara(1);
    fitK(i)=fitPara(2);
end

% deltaG=zeros(1,M);
% temp=input('input current experimental temperature: ("C)');
% for i=1:M
%     str=char(peakName(i));
%     ResidueNum=str2double(str(2:4));
%     j=0; if(ResidueNum>0&&ResidueNum<130) j=ResidueNum; end  %lysozyme has 129 residues 
%     if(j>0&&kc(j)~=0) deltaG(i)=-8.314*(273.15+temp)*log(fitK(i)/kc(j))/4200; end %kcal/mol
% end
% format short
% deltaG
k=-log(fitK)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%plot results. adjust scale and label when necessary:
spi=input('input subplot row number i: ');
spj=input('input subplot column number j: ');
pk=input('plot peak# (input 0) or peakName (input 1)? ');
if(pk~=0&&pk~=1) error('input 1 or 0 for peak plot option'); end
for i=1:M
    subplot(spi,spj,i)
    semilogx(time,intenCorr(i,:),'.')
    hold on
    timeFit=[0:0.1:500];
    semilogx(timeFit,fitA(i)*exp(-fitK(i)*timeFit*3600),'r','LineWidth',1.1)
    if(pk==1) text(0.3, 18000, peakName(i), 'EdgeColor','green'); end
    if(pk==0) text(0.3, 18000, num2str(i), 'EdgeColor','green'); end
    text(10, 10000, num2str(-log(fitK(i)),'%2.2f'))
    axis([0.2 500 -1000 24000])
    set(gca,'XTick',[0.1 1 10 100])
    if i==M
        set(gca,'XTickLabel',{'0.1','1','10','100'})
    else
        set(gca,'XTickLabel',{'','','',''})
    end
    set(gca,'YTick',0:5000:20000)
    if i==1
        set(gca,'YTickLabel',{'0','1','2','3','4'}) % *5000
    else
        set(gca,'YTickLabel',{'','',''})
    end
    hold on
end
