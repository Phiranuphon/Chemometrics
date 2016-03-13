function [MWRMSEall,OLS_percent,RMSEGC,reglist,allC]=MWrun(GCdata,sample,lipidlibrary,compnum,range,step,type);%need correction
%this function has a problem
%function MWOLS_GCcheck and MWOLS_GCcheckmean are correct but the
%reconstruction of ALS calculation in line 14 has a problem. I suggest to
%use IWrun instead of this one
if type==1
[MWRMSEallcheck,OLS_percent,RMSEGC]= MWOLS_GCcheck(GCdata,sample,lipidlibrary,compnum,range,step);
elseif type==2
[MWRMSEallcheckmean,OLS_percent,RMSEGC]=MWOLS_GCcheckmean(GCdata,sample,lipidlibrary,compnum,range,step);
end

[n m]=size(sample);
allC=zeros(n,compnum);
    if type==1
        reglist(:,1)=MWRMSEallcheck(n+3:2*n+2,2);
    
         if compnum==4
             for x=1:n
                 [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,1:4,[81:reglist(x,1)]);
                    allC(x,:)=C(7,1:4);
             end   
         elseif compnum==3
                for x=1:n
                    [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,2:4,[81:reglist(x,1)]);    
                allC(x,:)=C(8,1:3);
                end
         end
 elseif type==2
         reglist(1,1)=MWRMSEallcheckmean(3,2);
          if compnum==4
             for x=1:n
                 [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,1:4,[81:reglist(1,1)]);
                    allC(x,:)=C(7,1:4);
             end   
         elseif compnum==3
                for x=1:n
                    [C,A,R,ssq]=ALS(sample(x,:),lipidlibrary,compnum,2:4,[81:reglist(1,1)]);    
                allC(x,:)=C(8,1:3);
                end
          end     
         
    end
       
    if compnum==4        
        figure
        plot(GCdata(1,:)','b')
        hold on
        [Cn Cm]=size(allC);
        allC(Cn+1,:)=mean(allC);
        plot(allC(Cn+1,:)','r')
  
        
elseif compnum==3
        figure
        plot(GCdata(1,2:4)','b')
        hold on
        [Cn Cm]=size(allC);
        allC(Cn+1,:)=mean(allC);
        plot(allC(Cn+1,:)','r')% plot(allC(Cn+1,:)','r')

    end

end


%check 0 entry
%for x=1:15
%[n m]=find(allC(x,:)==0)
%if isempty(n)
%test(x,:)=allC(x,:)
%end
%end





  