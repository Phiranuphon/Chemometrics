function [C,residual,sim_spec]=MWLS(lipidlibrary,sample,windowsize,region,component);
lib=lipidlibrary;
WS=windowsize;
lipsam=sample;
NC=component;

[m n]=size(sample);

for z=2:NC
for x=1 
     rep=1;
     y_old=1;
 for y=1:WS:region
  
   C(rep,z)=lipsam(x,y_old:y)*pinv(lib(z,y_old:y));
   sim_spec(rep,y_old:y)=C(rep,z)*lib(z,y_old:y);
  R=lipsam(x,y_old:y)-sim_spec(rep,y_old:y);
  ssq=sum(sum(R.*R));
  residual(rep,z)=ssq;
  rep=rep+1;
  y_old=y;
 end 

end

end

%RMSE


%for x=1:m

%sim_spec(x,:)=C(x,:)*lib(:,region);
%RMSE(x,:)=sqrt(sum((lipsam(x,regino)-sim_spec(x,:))^2)/m)
%result(x,m+2)=RMSE(x,:);
%end
%result(1:m,1:NC)=C;
figure
plot(residual)
xlabel(WS)
end



