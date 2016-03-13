function [normdata]=intenorm2max(data,start,peak,stop);

%data=data';
[n m]=size(data);
dstr=data(:,start);
dpeak=data(:,peak);
dstp=data(:,stop);



normdata=zeros(n,m);
    for x=1:n
        
        dstrstp_avg(x,1)=mean(dstr(x,:));
        dstrstp_avg(x,2)=mean(dstp(x,:));
        dstrstp_avg(x,3)=mean(dstrstp_avg(x,1:2));
        
      
        
        peaksize(x,1)=max(dpeak(x,:))-dstrstp_avg(x,3);
        
        
        normdata(x,:)=data(x,:)/peaksize(x,1);
    end


end