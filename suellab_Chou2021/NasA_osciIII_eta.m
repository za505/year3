function dy = NasA_osciIII_eta(t,y,p)

    aq = p(1);
    bq = p(2);
    aF = p(3);
    bF = p(4);
    TT = p(5);
    bT = p(6);
    aA = p(7);
    bA = p(8);

    n1 = p(9);
    n2 = p(10);
    
    cq = p(11);
    aS = p(12);
    bS = p(13);
    aE = p(14);
    aN = p(15);
    ao = p(16);
    co = p(17);
    KT = p(18);
    KF = p(19);
    aT = p(20);
    cE = p(21);
    aR = p(22);
    bR = p(23);
    KR = p(24);
    e = p(25);
    k = p(26);
    cR = p(27);
    
    

    dy  =  zeros(10,1); 

    % glutamine production = supply from GS - gltAB using up glutamine - 
    % glutamine used to make biomass.
    dy(1) = e/(1+y(10))*(aq*y(5)*y(6)*y(7) - bq*y(4)*y(1)*y(8) - cq*y(1));
    
    % several downstream products of glutamine binding to GS (S) forming FBI-GS (y(2)) 
    dy(2) = e/(1+y(10))*(aF*y(5)*y(1)^n1/(KF^n1 + y(1)^n1) - bF*y(2));
    
    % Total TnrA = FBI-GS-bound inactive complex (FT) + free TnrA (T).
    dy(3) = e/(1+y(10))*(aT*TT  - bT*y(2)*y(3) - aT*y(3));
    %tnrA(M96A)
    %dy(3) = 0;

    % expression of gltAB = a constant synthesis that can be repressed by
    % the presence of free TnrA (T) - degradation of gltAB.
    dy(4) = e/(1+y(10))*(aA*(KT^n2)/(KT^n2+y(3)^n2) -  bA*y(4));

    % GS is being constantly produced and degraded, but also inhibited by
    % glutamine
    dy(5) = e/(1+y(10))*(aS - aF*y(5)*y(1)^n1/(KF^n1 + y(1)^n1) + bF*y(2) - bS*y(5));
    
    % glutamate has a constant production while also being consumed by
    % GS and GDH and resupplied by GltAB 
    dy(6) = e/(1+y(10))*(aE - aq*y(5)*y(6)*y(7) + bq*y(4)*y(1)*y(8)- aN*y(6) - cE*y(6));

    % ammonium is being made by rocG and consumed by GS
    dy(7) = e/(1+y(10))*(aN*y(6) - aq*y(5)*y(6)*y(7));
    
    % alpha-ketoglutarate is being supplied from the TCA cycle, being made
    % by rocG, and being consumed by gltAB
    dy(8) = e/(1+y(10))*(ao + aN*y(6) - bq*y(4)*y(1)*y(8) - co*y(8));  
    
    % Fluorescent reporter PnasA-yfp
    dy(9) = e/(1+y(10))*(aR*(y(3)^n2)/(KR^n2+y(3)^n2) -  bR*y(9)) - cR*y(9); 

    % The freeze factor (in h-1)
    dy(10) = k;
end