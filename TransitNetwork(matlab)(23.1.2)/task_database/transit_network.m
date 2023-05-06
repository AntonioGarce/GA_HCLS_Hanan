clc,close all,clear all
%% read excel file
[num1,txt1,input] = xlsread('KRAKOW - FULL DATA.xlsx','TRAFFIC-DEMAND');

%% get sort matrix
y1=diff(num1(:,1));
[ind,y2]=find(y1>0);
y3=num1(ind);
y3(length(y3)+1)=num1(ind(end)+1);
j=length(y3);
%% make output 
x1=num1(:,1);x2=num1(:,2);x3=num1(:,3);
z=zeros(j+1);
z(1,2:end)=y3;
z(2:end,1)=y3';
output=cell(j+1);
output(1,2:end)=num2cell(y3);
output(2:end,1)=num2cell(y3');
%  for i=1:200%length(num1)
%       [indx,xx]=find(z(1,:)==x1(i));
%       [indy,yy]=find(z(:,1)==x2(i));
%       if x3(i)==0
%          x3(i)=inf; 
%       end
%       z(indx+1,indy)=x3(i);      
%  end
%     
 ind(end+1) = length(x1);
 inds = 1;
 for i = 1:length(y3)
    inde = ind(i);
    for j = inds:inde
        [indy,yy]=find(z(:,1)==x2(j));
        z(i+1,indy)=x3(j);
        if(x3(j)==0)
            z(i+1,indy)=inf;
        end
        
    end
    inds = inde;
 end
 xlswrite('output1',z)