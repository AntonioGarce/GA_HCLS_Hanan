function [Parent]=FFTT(Parentge,No_Routs,P1,P2)
    %load('Mandl_data')
    load('KRAKOW_data_final.mat')        
    MandlTravelTimes=KRAKOWTravelTimes;
    MandlDemand=KRAKOWDemand;
    Rout=Parentge;ssdo=0;sstode=0;
    for rou=1:No_Routs
         Ro=Rout{rou};
         roulen = length(Ro);
         pto=zeros(roulen,1); Total_trans=zeros(roulen,1);
        
         for i=1:roulen-1
%              pto(i)= MandlTravelTimes(Ro(i),Ro(i+1))*MandlDemand(Ro(i),Ro(i+1));  
             for ir=i+1:roulen
                 ssdo=ssdo+ MandlDemand(Ro(i),Ro(ir));   
                 sstode=sstode+ MandlDemand(Ro(i),Ro(ir)) ; 
             end
             Total_trans(i)=MandlTravelTimes(Ro(i),Ro(i+1))*MandlDemand(Ro(i),Ro(i+1));    
         end
         ddo=  ssdo; 
         total_demand = sstode;
         
%          sum_pto(rou)=sum(pto);   
         sum_ddo(rou)=(ddo); 
         sum_totra(rou)=sum(Total_trans);
    end
    
    Parent.ddo=sum(sum_ddo);
    Parent.PTo=sum(sum_ddo);
    Parent.TIVT=sum(sum_totra);
    Parent.Total_trans=sum(total_demand);
    
    %% one transfer
    Rout=[];sum_pt1=[];dd1=[];sum_t1=[];sum_dd1=[];
    co=0;
    Rout=Parentge;
    One_trans=[];
    for rou=1:No_Routs-1
        R= [];R=        {rou};
        for rou1=rou+1:No_Routs
            R1= [];R1=Rout{rou1}; 
            [ chek_intr,locR]=intersect(R,R1,'stable');
            locR1=find(chek_intr==R1(:));
        
             %[ ~,locR1]=intersect(R1,R,'stable');
            T1=[];pt1=[];dd1=[];
            if isempty(chek_intr)==0
                new_rout_end=[];new_rout_end= [R(1:locR(1)) R1(locR1(1)+1:end)] ;  
                uniq_new_end=[];uniq_new_end=new_rout_end;%unique( new_rout_end,'stable');   
                new_rout_start=[];new_rout_start= [(R(1:locR(1))) flip(R1(1:locR1(1)))]  ;
                uniq_new_start=[];uniq_new_start=new_rout_start;%unique( new_rout_start,'stable');   
                 
                for ii=1:length(uniq_new_end)-1
                    T1(ii)=MandlDemand(uniq_new_end(ii),uniq_new_end(ii+1)); 
                    pt1(ii)= MandlTravelTimes(uniq_new_end(ii),uniq_new_end(ii+1))*MandlDemand(uniq_new_end(ii),uniq_new_end(ii+1));   
                    dd1(ii)=MandlDemand(uniq_new_end(ii),uniq_new_end(ii+1)); 
                end
          
                for ii=1:length(uniq_new_start)-1
                   T1(end+ii)=MandlDemand(uniq_new_start(ii),uniq_new_start(ii+1)); 
                   pt1(end+ii)= MandlTravelTimes(uniq_new_start(ii),uniq_new_start(ii+1))*MandlDemand(uniq_new_start(ii),uniq_new_start(ii+1));    
                   dd1(end+ii)=MandlDemand(uniq_new_start(ii),uniq_new_start(ii+1)); 
                end
             else
                 T1=[];
                 pt1=[];
                 dd1=[];
                 new_rout_end=[];
                 uniq_new_end=[]; 
                 new_rout_start=[];  
                 uniq_new_start=[];  
             end
             co=co+1;
             sum_t1(co)=sum(T1);
             sum_pt1(co)=sum(pt1);
             sum_dd1(co)=sum(dd1);
             One_trans(co).inter2end=uniq_new_end;      
             One_trans(co).inter2start=flip(uniq_new_start);
        end 
    end
    Parent.ontrans=One_trans;
    %Parent(ge).T1=sum(sum_t1);
    Parent.P_T1=sum(sum_dd1)*P1;
    if sum(sum_dd1)<=0
        Parent.PT1=0;
    else
        Parent.PT1=sum(sum_pt1);
    end
        Parent.dd1=sum(sum_dd1);
    
    %% two transfer
    two_trans=[];
    co=0;sum_t2=[];sum_dd2=[];
    for rou=1:No_Routs-2
        R= [];R=Rout{rou}; 
        for rou1=rou+1:No_Routs-1
            R1= [];R1=Rout{rou1}; 
            [ chek_intr,locR]=intersect(R,R1,'stable');
            locR1=find(chek_intr==R1(:));
    
            %[ ~,locR1]=intersect(R1,R,'stable');
            if isempty(chek_intr)==0
                R1= [];R1=Rout{rou1}; 
                new_rout_end=[];new_rout_end=  R1(locR1(1):end);   
         
                for rou2=rou1+1:No_Routs
                    R2= [];R2=Rout{rou2}; 
                    [ chek_intr2,locR2]=intersect(new_rout_end,R2,'stable');
                    locR3=find(chek_intr2==R2(:));
    
                    %[ ~,locR3]=intersect(R2,new_rout_end,'stable');
                    T2=[];pt2=[];dd2=[];
    
                    if isempty(chek_intr2)==0
                        new_rout_end2=[];new_rout_end2= [R(1:locR(1)-1) new_rout_end(1:locR2(1)) R2(locR3(1)+1:end)] ;  
                        uniq_new_end2=[];uniq_new_end2=new_rout_end2;%unique( new_rout_end2,'stable');   
                        new_rout_start2=[];new_rout_start2= [R(1:locR(1)-1) new_rout_end(1:locR2(1)) flip(R2(1:locR3(1)))];
                        uniq_new_start2=[];uniq_new_start2=new_rout_start2;%unique( new_rout_start2,'stable');   
    
                        for ii=1:length(uniq_new_end2)-1
                            T2(ii)= MandlDemand(uniq_new_end2(ii),uniq_new_end2(ii+1)); 
                            pt2(ii)= MandlTravelTimes(uniq_new_end2(ii),uniq_new_end2(ii+1))*MandlDemand(uniq_new_end2(ii),uniq_new_end2(ii+1));   
                            dd2(ii)=MandlDemand(uniq_new_end2(ii),uniq_new_end2(ii+1)); 
                        end
                        for ii=1:length(new_rout_start2)-1
                            T2(end+ii)= MandlDemand(new_rout_start2(ii),new_rout_start2(ii+1)); 
                            pt2(end+ii)= MandlTravelTimes(new_rout_start2(ii),new_rout_start2(ii+1))*MandlDemand(new_rout_start2(ii),new_rout_start2(ii+1));   
                            dd2(end+ii)=MandlDemand(new_rout_start2(ii),new_rout_start2(ii+1)); 
                        end
                    else
                        T2=[];
                        pt2=[];
                        dd2=[];
                        new_rout_end2=[];
                        uniq_new_end2=[]; 
                        new_rout_start2=[];  
                        uniq_new_start2=[];  
     
                    end
                    co=co+1;
                    sum_t2(co)=sum(T2);
                    sum_pt2(co)=sum(pt2);
                    sum_dd2(co)=sum(dd2);
                    two_trans(co).inter2end=uniq_new_end2;
                    two_trans(co).inter2start=flip(uniq_new_start2);
                end
            else
                T2=[];
                pt2=[];
                dd2=[];
                new_rout_end2=[];
                uniq_new_end2=[]; 
                new_rout_start2=[];  
                uniq_new_start2=[];  
            end
        end
    end

    Parent.tworans=two_trans;
    Parent.T2=sum(sum_t2);
    Parent.P_T2=sum(sum_dd2)*P2;
    Parent.dd2=sum(sum_dd2);

    if sum(sum_dd2)<=0
        Parent.PT2=0;
    else
        Parent.PT2=sum(sum_pt2);
    end
end



