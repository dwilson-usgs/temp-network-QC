% TempNetQC.m
% 
% Matlab script to QC temporary networks using event based metrics.
%
% This script uses the irisFetch matlab function. The latest version is
% available here:
% http://ds.iris.edu/ds/nodes/dmc/software/downloads/irisfetch.m/
% The IRIS WS jar file is also required, it is available here:
% http://ds.iris.edu/ds/nodes/dmc/software/downloads/IRIS-WS/ 
% make sure the java jar is in the path
% The version used here is IRIS-WS-2.0.15.jar
%
%**Disclaimer:**
%
%This software is preliminary or provisional and is subject to revision. It is 
%being provided to meet the need for timely best science. The software has not 
%received final approval by the U.S. Geological Survey (USGS). No warranty, 
%expressed or implied, is made by the USGS or the U.S. Government as to the 
%functionality of the software and related material nor shall the fact of release 
%constitute any such warranty. The software is provided on the condition that 
%neither the USGS nor the U.S. Government shall be held liable for any damages 
%resulting from the authorized or unauthorized use of the software.

clear
javaaddpath('include\IRIS-WS-2.0.15.jar');
path(path,'include')

% Load network information  ----------------------------------
stas = [{'GS','OK044','00','HHZ','HH1','HH2'};
    {'GS','OK045','00','HHZ','HH1','HH2'};
    {'GS','OK046','00','HHZ','HH1','HH2'};
    {'GS','OK047','00','HHZ','HH1','HH2'};
    {'GS','OK048','00','HHZ','HH1','HH2'};
    {'GS','OK049','00','HHZ','HH1','HH2'};
    {'GS','OK050','00','HHZ','HH1','HH2'};
    {'GS','OK051','00','HHZ','HH1','HH2'}];

stcoords = [ 36.372700 -96.796330 290 ;
 36.447670 -96.924950 284  ;
 36.397250 -96.907410 294  ;
 36.468870 -96.723620 253  ;
 36.416220 -96.943740 297  ;
 36.406880 -97.304410 327  ;
 36.394500 -96.981800 315  ;
 36.505190 -96.836040 263];

% coordinate of monitoring target (mainshock, etc.)
xcoord=[36.42591 -96.92931 -5.6];

% title for plots
titl='Pawnee';

% radius for grabbing events, default is max station dist from xcoord.
Radkm = [];

% time to compute QC
tstart='2016-09-15 00:00:00';
tend='2016-09-16 00:00:00';

% frequency bin to analyze for s/n (Hz)
Fbin = [ 5 20];

% velocity for time shift
Vp=4.5;

debug=0;
pubplot=1;

% --------------------------------------------------------------


lnst=length(stas(:,1));
taxis = datenum(tstart):1/86400:datenum(tend);
if length(taxis) > 3600,
    taxis2 = datenum(tstart):3600/86400:(datenum(tend)-3600/86400);
else
    taxis2 = datenum(tstart);
end
% compute station distances
sdist=[];
for n = 1:lnst,
    sdist(n,1)=edist(xcoord(1),xcoord(2),stcoords(n,1),stcoords(n,2))*111.1949;
end
if length(Radkm)<1,
    Radkm=max(sdist);
end
sdist=((sdist).^2 + (stcoords(:,3)./1000 - xcoord(3)).^2).^.5;

DataMatrixz=zeros(lnst,length(taxis));
DataMatrixh1=zeros(lnst,length(taxis));
DataMatrixh2=zeros(lnst,length(taxis));
 for m=1:lnst,
     for n = 1:length(taxis2),
           tdelay= sdist(m)/Vp;
           Mytstart=taxis2(n)+tdelay/86400; Mytend=Mytstart+3600/86400;
           for cn=4:6,
               tic;
                mytracez=irisFetch.Traces(stas(m,1),stas(m,2),stas(m,3),stas(m,cn),datestr(Mytstart,31),datestr(Mytend,31),'includePZ');
                disp(sprintf('%2.2f seconds fetching %s-%s-%s-%s',toc,char(stas(m,1)),char(stas(m,2)),char(stas(m,3)),char(stas(m,cn))));
                tic;
                if ~isempty(mytracez),
                  for x=1:length(mytracez),
                    [B,A]=butter(2,Fbin/(mytracez(x).sampleRate));
                    datz=abs(hilbert(filtfilt(B,A,RemoveResp(mytracez(x),1,.001))));
                    if mytracez(x).sampleRate>1,
                      datz = meanfilter(datz,mytracez(x).sampleRate);
                    end
                    
                    datz_dec = datz(1:mytracez(x).sampleRate:length(datz));
                    [ yy, Istart] = min(abs(taxis -(datenum(mytracez(x).startTime)-tdelay/86400)));
                    Iend = Istart+length(datz_dec)-1;
                    switch cn,
                      case 4,
                        DataMatrixz(m,Istart:Iend)=datz_dec';
                      case 5,
                        DataMatrixh1(m,Istart:Iend)=datz_dec';
                      case 6,
                        DataMatrixh2(m,Istart:Iend)=datz_dec';
                    end
                  end
                end
                disp(sprintf('%2.2f seconds processing',toc));
            end
            if debug==1,
              figure(1); clf
              subplot(3,1,1)
              plot(DataMatrixz(m,:),'k')
              title(sprintf('%s-%s %i of %i',char(stas(m,1)),char(stas(m,2)),m,lnst))
              hold on
              subplot(3,1,2)
              plot(DataMatrixh1(m,:),'k')
              hold on
              subplot(3,1,3)
              plot(DataMatrixh2(m,:),'k')
              hold on
              pause(.1)
            end
        end
     end



%compute network median stacks
medz=median(DataMatrixz,1);
DataMatrixH=(DataMatrixh1.^2 + DataMatrixh2.^2).^.5;
medh=median(DataMatrixH,1);


% identify network events
[zmedfilt15, zmedstd15]=meanfilter(medz,15,'B');
[zmedfilt5, zmedstd5]=meanfilter(medz,5,'F');
Zsrcfn=(zmedfilt5-zmedfilt15)./zmedstd15;
Zsrcfn=max(0,Zsrcfn);
[ Imax, maxval ] = localmax(Zsrcfn);
II=find(maxval>=3);
Imax2=Imax(II);
maxval2=maxval(II);
zvels=zmedfilt5(Imax2);

[hmedfilt15, hmedstd15]=meanfilter(medh,15,'B');
[hmedfilt5, hmedstd5]=meanfilter(medh,5,'F');
Hsrcfn=(hmedfilt5-hmedfilt15)./hmedstd15;
Hsrcfn=max(0,Hsrcfn);
[ Imax, maxval ] = localmax(Hsrcfn);
II=find(maxval>=3);
Imax2h=Imax(II);
maxval2h=maxval(II);
hvels=hmedfilt5(Imax2h);



% compute individual station signal to noise
dathSN=[]; datzSN=[]; 
for n = 1:lnst,
    [hmedfilt15, hmedstd15]=meanfilter(DataMatrixH(n,:),15,'B');
    [hmedfilt5, hmedstd5]=meanfilter(DataMatrixH(n,:),5,'F');
    Hsrcfn2=(hmedfilt5-hmedfilt15)./hmedstd15;
    
        dathSN(:,n)=Hsrcfn2(Imax2h);
         
        [hmedfilt15, hmedstd15]=meanfilter(DataMatrixz(n,:),15,'B');
    [hmedfilt5, hmedstd5]=meanfilter(DataMatrixz(n,:),5,'F');
    Hsrcfn2=(hmedfilt5-hmedfilt15)./hmedstd15;
    
        datzSN(:,n)=Hsrcfn2(Imax2);
       
end
      


    figure(3); clf
    plot(Zsrcfn,'k')
    hold on
    plot(Hsrcfn,'r')
    
 % compute estimated magnitudes based on coefficitents from best fit of 
 % 1 week of events from the Fairview, Feb 2016 deployment.  
zmag_predicted=log10(zvels).*.9308 + 7.0077;
hmag_predicted=log10(hvels).*.9052 + 6.4078;

    figure(4); clf
   
     semilogy(zmag_predicted,maxval2,'kx','linewidth',1.5,'markersize',9)
    hold on
    semilogy(hmag_predicted,maxval2h,'rx','linewidth',1.5,'markersize',9)
  
% fetch events for time period
MyEvents = irisFetch.Events('startTime',tstart,'endTime',tend,'radialcoordinates',[ xcoord(1) xcoord(2) Radkm/111.1949]);
nh=0; nz=0; magh=[]; magz=[]; maghvels=[]; magzvels=[];
if ~isempty(MyEvents),
    % find matching times, and plot on figure 4
    for x = 1:length(MyEvents),
         [YYh,IIh]=min(abs(datenum(MyEvents(x).PreferredTime)-taxis(Imax2h)));
         [YY,II]=min(abs(datenum(MyEvents(x).PreferredTime)-taxis(Imax2)));
         
         figure(4);
         if YYh*86400 < 5,
            %loglog(maxval2h(IIh),hvels(IIh),'ro')
            semilogy(hmag_predicted(IIh),maxval2h(IIh),'ro','linewidth',1.4,'markersize',9)
            nh=nh+1;
            magh(nh)=MyEvents(x).PreferredMagnitudeValue;
            %semilogy(magh(nh),maxval2h(IIh),'rd')
            maghvels(nh)=hvels(IIh);
         end
         if YY*86400 < 5,
            %loglog(maxval2(II),zvels(II),'ko')
            semilogy(zmag_predicted(II),maxval2(II),'ko','linewidth',1.4,'markersize',9)
            nz=nz+1;
            magz(nz)=MyEvents(x).PreferredMagnitudeValue;
            %semilogy(magz(nz),maxval2(II),'kd')
            magzvels(nz)=zvels(II);
         end
      
    end
else
    disp('No catalog events found for this time period')
end
ylabel('signal/noise (# of std from noise mean)','fontsize',14)
xlabel('estimated magnitude','fontsize',14)
%title(sprintf('%s %s-%s',titl,tstart,tend),'fontsize',14)
legend(gca,sprintf('Vertical (%i events)',length(maxval2)),sprintf('Horizontal (%i events)',length(maxval2h)),'Location','NorthWest')
set(gca,'fontsize',14)

% create histogram of estimated magnitudes.
minX = min([min(zmag_predicted) min(hmag_predicted)])-.1;
maxX = max([max(zmag_predicted) max(hmag_predicted)])+.1;
edges=floor(minX*5)/5:.2:ceil(maxX*5)/5;
zhist = histc(zmag_predicted,edges);
hhist = histc(hmag_predicted,edges);


figure(6); clf
h=[];
h(1)=semilogy(edges,zhist,'kx','linewidth',1.5,'markersize',9);
hold on
h(2)=semilogy(edges,hhist,'rx','linewidth',1.5,'markersize',9);
xlabel('magnitude','fontsize',14)
ylabel('number of events','fontsize',14)
set(gca,'fontsize',14)
legend('1 day vertical','1 day horizontals')
% title(sprintf('%s %s-%s',titl,tstart,tend))
 axis([minX maxX 1 max(max([zhist hhist]))*1.5])
grid on


dathSN=max(dathSN ,.1);
figure(7); clf
colo=(jet(lnst)+.75)./2;
colo1=jet(lnst);
h=[];

for n = 1:lnst,
    semilogy(hmag_predicted,dathSN(:,n),'.','MarkerSize',9,'LineWidth',1.5,'Color',colo(n,:));
    hold on
    
end
for n = 1:lnst,
    mcnt=0;
    aveM=[];
    for m=edges,
        [I]=find(abs(m-hmag_predicted)<=.2);
        if ~isempty(I),
            mcnt=mcnt+1;
            aveM(mcnt,:)=[m mean(dathSN(I,n))];
        end
    end
    semilogy(aveM(:,1),aveM(:,2),'-','LineWidth',2.2,'Color',colo1(n,:));
    h(n)=semilogy(aveM(:,1),aveM(:,2),'-','LineWidth',1.5,'Color',colo(n,:));
end
ylabel('Horizontal Signal to noise (# std)','fontsize',14)
xlabel('estimated magnitude','fontsize',14)
semilogy([min(hmag_predicted) max(hmag_predicted)],[3 3],'k--','LineWidth',2)
legend(h,stas(:,2),'Location','EastOutside')
axis([min(zmag_predicted)*.9 max(zmag_predicted)*1.1 .5 5000])
set(gca,'fontsize',14)
orient Landscape
%title(sprintf('%s %s-%s',titl,tstart,tend),'fontsize',14)

datzSN=max(datzSN ,.1);
figure(8); clf
h=[];

for n = 1:lnst,
    semilogy(zmag_predicted,datzSN(:,n),'.','MarkerSize',9,'LineWidth',1.5,'Color',colo(n,:));
    hold on
    
end
for n = 1:lnst,
    mcnt=0;
    aveM=[];
    for m=edges,
        [I]=find(abs(m-zmag_predicted)<=.2);
        if ~isempty(I),
            mcnt=mcnt+1;
            aveM(mcnt,:)=[m mean(datzSN(I,n))];
        end
    end
    semilogy(aveM(:,1),aveM(:,2),'-','LineWidth',2.2,'Color',colo1(n,:));
    h(n)=semilogy(aveM(:,1),aveM(:,2),'-','LineWidth',1.5,'Color',colo(n,:));
end
ylabel('Vertical Signal to noise (# std)','fontsize',14)
xlabel('estimated magnitude','fontsize',14)
semilogy([min(zmag_predicted) max(zmag_predicted)],[3 3],'k--','LineWidth',2)
legend(h,stas(:,2),'Location','EastOutside')
axis([min(zmag_predicted)*.9 max(zmag_predicted)*1.1 .2 5000])
set(gca,'fontsize',12)
orient Landscape
%title(sprintf('%s %s-%s',titl,tstart,tend),'fontsize',14)

colo=jet(lnst);
figure(9); clf
for n=1:lnst,
          subplot(3,1,1)
          semilogy(DataMatrixz(n,:),'linewidth',1.5,'Color',colo(n,:))
          hold on
          subplot(3,1,2)
          semilogy(DataMatrixH(n,:),'linewidth',1.5,'Color',colo(n,:))
          hold on
end    
          subplot(3,1,1)
          title('Vertical Station Energy Envelopes','fontsize',14)
          ylabel('m/s','fontsize',12)
          set(gca,'fontsize',12)
          axis([0 2*3600 min(min(DataMatrixz)) 1e-5])
          
          subplot(3,1,2)
          title('Horizontal Station Energy Envelopes','fontsize',14)
          ylabel('m/s','fontsize',12)
          set(gca,'fontsize',12)
          axis([0 2*3600 min(min(DataMatrixH)) 1e-5])
          subplot(3,1,3)
          semilogy(medz,'k','linewidth',1.5)
          hold on
          semilogy(medh,'r','linewidth',1.5)
          legend('Z','H','Location','NorthEast')
          title('Network Median Energy Envelopes','fontsize',14)
          axis([0 2*3600 min(medz) 1e-5])
          ylabel('m/s','fontsize',12)
          set(gca,'fontsize',12)
          xlabel(sprintf('Seconds after %s',tstart),'fontsize',12)
       
          
          


% plot envelope example
figure(19); clf;
  subplot(2,1,1)
          semilogy(medz,'k','linewidth',1.5)
          hold on
          semilogy(medh,'-','Color',[.5 .5 .5],'linewidth',1.5)
          legend('Z','H','Location','NorthEast')
          %title('Network Median Energy Envelopes','fontsize',14)
          axis([5200 2*3600 min(medz) 1e-5])
          %ax=axis; axis([-20 20 ax(3:4)])
          ylabel('m/s','fontsize',14)
          set(gca,'fontsize',14)
          yticks([1e-7 1e-6 1e-5])
          
          %xlabel(sprintf('Seconds after %s',tstart),'fontsize',14)
          subplot(2,1,2)
          
          semilogy(max(.5,Zsrcfn),'k','linewidth',1.5)
    hold on
    semilogy(max(.5,Hsrcfn),':','Color',[.5 .5 .5],'linewidth',1.5)
    semilogy([5200 7200],[3 3],'k--','linewidth',1.2)
    axis([5200 2*3600 .9 100])
    set(gca,'fontsize',14)
    ylabel('Signal/Noise')
          xlabel(sprintf('Seconds after %s',tstart),'fontsize',14)
          legend('Z','H','S/N cutoff','Location','NorthWest')
          yticks([1 10 100])
          
          