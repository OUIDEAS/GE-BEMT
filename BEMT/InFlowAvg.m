for i = 1:20
   FlowAvgs(i) = sum(Inflow2(:,i))/60;
    
end

flowAvgFin = sum(FlowAvgs)/20;

