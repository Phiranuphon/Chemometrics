function [normdata,dpeak_sort,dstrstp_avg]=intenorm2(data,start,peak,stop);


[n m]=size(data);
dstr=data(:,start);
dpeak=data(:,peak);
[dpn dpm]=size(data(:,peak));
dstp=data(:,stop);
normdata=zeros(n,m);
dpeak_sort=zeros(n,dpm);
dstrstp_avg=zeros(n,3);
dmaxind=zeros(n,3);
dpeak2=zeros(n,3);

for x=1:n
    %dstrstp_avg(x,1)=mean(dstr(x,:));
    %dstrstp_avg(x,2)=mean(dstp(x,:));
    dstrstp_avg(x,3)=mean(data(x,[start,stop]));
    
    dpeak_sort(x,:)=sort(dpeak(x,:),'descend');
    %dmaxind(x,1)=find(dpeak(x,:)==dpeak_sort(x,1));
    %dmaxind(x,2)=find(dpeak(x,:)==dpeak_sort(x,2));
    %dmaxind(x,3)=find(dpeak(x,:)==dpeak_sort(x,3));
    
    %dpeak2(x,1)=dpeak(x,dmaxind(x,1));
    %dpeak2(x,2)=dpeak(x,dmaxind(x,2));
    %dpeak2(x,3)=dpeak(x,dmaxind(x,3));
    
    dpeak_avg(x,1)=mean(dpeak_sort(x,1:3));
    peaksize(x,1)=dpeak_avg(x,1)-dstrstp_avg(x,1);
    %peaksize(x,1)=dpeak_sort(x,1)-dstrstp_avg(x,3);
    normdata(x,:)=data(x,:)/peaksize(x,1);
end

end

%[ICAsub.IC1_5subintenorm2]
%=intenorm2(ICAsub.IC1_5sub,1172:1180,1181:1188,1189:1190);