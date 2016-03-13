function eqinterval_FTIR(a1,a2,a3,a4,a5,a6);

curfol=cd;
%for j=str:stp%change the number of folder
for j={a1,a2,a3,a4,a5,a6}
%foldernum=num2str(j{:});
dishfol= [curfol,'\',j{:}];
cd(dishfol);    

mkdir ('eqinterval')
imported=dir('*.txt');   
[m n]=size(imported);
for i=1:m
raw=imported(i,1).name;
fid = fopen(raw);
tscan= textscan(fid,'%n%n','headerlines',4,'delimiter','\t');
fclose(fid);
c(1,:)=tscan{1,1}';
c(2,:)=tscan{1,2}';

%%wavenumber
newWN(1,1)=c(1,1);
flnum=floor(c(1,2075));
decinum=c(1,1)-floor(c(1,1));
flnum=flnum+decinum;
for x=1:2075
if newWN<flnum
newWN(1,x+1)=newWN(1,x)+2;
end
end



%interpol by Lagrange

yi(1,1) = c(2,1);
yi(1,2) = LagrangeInter(c(1,1:3),c(2,1:3),newWN(1,2));
yi(1,3) = LagrangeInter(c(1,1:3),c(2,1:3),newWN(1,3));
for x=4:2001
for y=1:2075
flsel(1,y)=abs(newWN(1,x)-c(1,y));
end
[idx idx] = min(flsel(1,:));
flselidx(1,x)=idx;
idx2=idx-2;
yi(1,x) = LagrangeInter(c(1,idx2:idx),c(2,idx2:idx),newWN(1,x));
end
newWN(2,:)=yi;


prtfile(:,1)=num2cell(newWN(1,:));
prtfile(:,2)=num2cell(newWN(2,:));
dlmcell('C:\Users\Kin\Desktop\test interval\test.txt',prtfile,'\t','w');%% change back to text
newname=[dishfol,'\eqinterval\',raw,'.txt'];%%new destination
oldname= ('C:\Users\Kin\Desktop\test interval\test.txt');%%old destination
copyfile(oldname,newname);


end
end

cd(curfol);    
end

