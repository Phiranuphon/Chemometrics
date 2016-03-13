function [MWRMSEallspec,OLS_percent,RMSEGC]= MWOLS_check_spec(sample,lipidlibrary,compnum,range,step);
WNstr=80;   
%WNstp=273;
count=1;
WNstr_step=WNstr+step;

[n m]=size(sample);
j=compnum;
Cmix=eye(j);

while WNstr_step<range

[C,A,R,ssq]=ALS(sample,n,lipidlibrary,j,[2:4],WNstr:WNstr_step);%2:4 only for POL, C10 must use 1:4 and j must be 4

OLS(1:n,:)=C(j+1:n+j,:);%%create OLS conc. matrix

    for y=1:j
        for x=1:n
            OLS_percent(x,y)=OLS(x,y)/sum(OLS(x,1:j));%re-calculate conc. into %
        end
    end


    
    %ต้องทำnormalizeใหม่ตรงนี้โดย 1.ใช้1440 2.local max(บริเวณที่ทำ)
    %**เฉพาะการหาerrorจาก spectra เท่านั้นถ้าหาจาก C ไม่ต้อง
RMSEGC =zeros(n,m);
   
    for x=1:n
    for y=1:j
    simspec(x,:)=OLS_percent(x,y)*lipidlibrary(y,:);
    
    end
    end



RMSEall(1:n,count)=RMSEGC(1:n,11);
RMSEall(n+1,count)=WNstr_step;




count=count+1;
WNstr=WNstr_step;
WNstr_step=WNstr_step+step;
%WNstp=120+(step*count);


end
%check NaN entry
 nancheck(1:3,:)=isnan(RMSEall(1:3,:));
 [nc1 nc2]=size(nancheck);
    for x=1:nc1
        for y=1:nc2
            if nancheck(x,y)==1
                RMSEall(x,y)=0;
            end   
        end
     end   
 for y=1:j
     if RMSEall(y,:)
     
     
 end


%find min of each row

   for y=1:j
        testmin=find(RMSEall(y,:)<mean(RMSEall(y,:)));
        [tn tm]=size(testmin(1,:));
            for x=1:tm
                RMSEall(6+y*2,x)=RMSEall(y,testmin(1,x));
                RMSEall(7+y*2,x)=RMSEall(4,testmin(1,x));
            end
   end
    
   MWRMSEallspec=RMSEall;
end

