function [MWRMSEall,RMSEall,OLS_percent,RMSEGC]= MWOLS_GCcheck(GC,sample,lipidlibrary,compnum,range,step);
WNstr=80;
%WNstp=273;
count=1;
WNstr_step=WNstr+step;
GC2=GC;

[n m]=size(sample);
j=compnum;
Cmix=eye(j);

if compnum==3
        GCdata(1,1:j)=GC2(:,2:4);
    elseif compnum==4
        GCdata(1,1:j)=GC2(:,1:4);
    end

while WNstr_step<range%range is WN limit

 if compnum==4
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[1:4],WNstr:WNstr_step);
    elseif compnum==3
        [C,A,R,ssq]=ALS(sample,lipidlibrary,j,[2:4],WNstr:WNstr_step);%2:4 only for POL, C10 must use 1:4 and j must be 4  
    end

OLS(1:n,:)=C(j+1:n+j,:);%%create OLS conc. matrix

    for y=1:j
        for x=1:n
            OLS_percent(x,y)=OLS(x,y)/sum(OLS(x,1:j));%re-calculate conc. into %
        end
    end

    %for x=1:j
     %   OLS_percent(48,x)=std(OLS_percent([1:15],x));

     %  OLS_percent(51,x)=std(OLS_percent(16:30,x));

     % OLS_percent(54,x)=std(OLS_percent(31:45,x));

     %   OLS_percent(47,x)=mean(OLS_percent([1:15],x));

     %   OLS_percent(50,x)=mean(OLS_percent(16:30,x));

     %   OLS_percent(53,x)=mean(OLS_percent(31:45,x));
 
    %end
    
    %ต้องทำnormalizeใหม่ตรงนี้โดย 1.ใช้1440 2.local max(บริเวณที่ทำ)
    %**เฉพาะการหาerrorจาก spectra เท่านั้นถ้าหาจาก C ไม่ต้อง

    for x=1:n
        for y=1:j
            RMSEGC(x,y)=GCdata(1,y)-OLS_percent(x,y);
        end
    end

    for x=1:n
        for y=1:j%for y=1:4
            RMSEGC(x,y+5)= sqrt(RMSEGC(x,y)^2);   
        end    
    end

    for x=1:n
        RMSEGC(x,11)=sum(RMSEGC(x,1+5:j+5))/j;
    end

RMSEall(1:n,count)=RMSEGC(1:n,11);
RMSEall(n+1,count)=WNstr_step;




count=count+1;
WNstr=WNstr_step;
WNstr_step=WNstr_step+step;
%WNstp=120+(step*count);
if WNstr_step<range
    fin_count=count;

end    

end

 nancheck(1:n,:)=isnan(RMSEall(1:n,:));
 [nc1 nc2]=size(nancheck);
    for x=1:nc1
        for y=1:nc2
            if nancheck(x,y)==1
                RMSEall(x,y)=0;
            end   
        end
     end   
 




 for x=1:n
        testmean=[];
        testmean(1,1)=mean(RMSEall(x,1:fin_count));
        [testa1,testb1]=find(RMSEall(x,1:fin_count)<testmean(1,1)/2);%selection criteria
        [tb1n tb1m]=size(testb1);
        if tb1m==1
        RMSEall(x+n+2,1)=RMSEall(x,testb1);
        RMSEall(x+n+2,2)=RMSEall(n+1,testb1);
        elseif tb1m>1
            for i=1:tb1m
                RMSEall(x+x+n+2-1,i)=RMSEall(x,testb1(1,i));
                RMSEall(x+x+n+2,i)=RMSEall(n+1,testb1(1,i));  
            end 
        end
    end

    
   MWRMSEall=RMSEall;
end

