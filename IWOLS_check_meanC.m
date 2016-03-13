function [IWRMSEall,OLS_percent,RMSEGC]= IWOLS_check_meanC(GC,sample,lipidlibrary,compnum,range,step);
WNstr=60;
WNstp=WNstr+step;
WNstp1=WNstr+step;
GC2=GC;
count=1;

[n m]=size(sample)
j=compnum;
Cmix=eye(j);

    if compnum==3
        GCdata(1,1:j)=GC2(:,2:4);
    elseif compnum==4
        GCdata(1,1:j)=GC2(:,1:4);
    end

while(WNstp<range)
    if compnum==4
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[1:4],WNstr:WNstp);
    elseif compnum==3
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[2:4],WNstr:WNstp);%2:4 only for POL, C10 must use 1:4 and j must be 4  
    end

OLS(1:n,:)=C(j+1:n+j,:);%%create OLS conc. matrix

    for y=1:j
        for x=1:n
          OLS_percent(x,y)=OLS(x,y)/sum(OLS(x,1:j));%re-calculate conc. into %
        end
    end


   OLS_percent(n+2,1:j)=mean(OLS_percent(1:n,1:j));
      
   for y=1:j
   RMSEGC(1,y)=GCdata(1,y)-OLS_percent(n+2,y);
   end
  
  

    for x=1
        for y=1:j%for y=1:4
            RMSEGC1(x,y)= sqrt(RMSEGC(x,y)^2); %  RMSEGC1(x,1)= RMSEGC(x,1)^2;    
        end    
    end

    for x=1
        RMSEGC2(x,1)=sum(RMSEGC1(x,1:j))/j;% RMSEGC2(x,1)=sqrt(RMSEGC1(x,1));
    end

RMSEall(1,count)=RMSEGC2(1,1);
RMSEall(1+1,count)=WNstp;

count=count+1;
WNstp=WNstp1+(step*count);


end

    for x=1
        testmin(1,1)=min(RMSEall(x,:));
        [testa1,testb1]=find(RMSEall(x,:)==testmin(1,1));
        [tb1n tb1m]=size(testb1);
        if tb1m==1
            RMSEall(x+2,1)=RMSEall(x,testb1);
            RMSEall(x+2,2)=RMSEall(1+1,testb1);
        elseif tb1m>1
            for i=1:tb1m
                RMSEall(x+2,i+i-1)=RMSEall(x,testb1(1,i));
                RMSEall(x+2,i*2)=RMSEall(1+1,testb1(1,i));  
            end 
        end
    end


IWRMSEall=RMSEall;

end