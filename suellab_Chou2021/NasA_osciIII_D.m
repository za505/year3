function dy = NasA_osciIII_nD(t,y,p)

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
    
    % parameters for the alternative negative feedback loop
    %dq = p(26);
    %aX = p(27);
    %bX = p(28);
    
    

    dy  =  zeros(9,1);% For alternative negative feedback there are 10 valuables

    % glutamine production = supply from GS - gltAB using up glutamine - 
    % glutamine used to make biomass.
    dy(1) = aq*y(5)*y(6)*y(7) - bq*y(4)*y(1)*y(8) - cq*y(1); % - dq*y(10)*y(1) additional term for the alternative negative feedback
    
    % several downstream products of glutamine binding to GS (S) forming FBI-GS (y(2)) 
    dy(2) = aF*y(5)*y(1)^n1/(KF^n1 + y(1)^n1) - bF*y(2);
    
    % Total TnrA = FBI-GS-bound inactive complex (FT) + free TnrA (T).
    dy(3) = aT*(TT - y(3)) - bT*y(2)*y(3);
    %dy(3) = 0;

    % expression of gltAB = a constant synthesis that can be repressed by
    % the presence of free TnrA (T) - degradation of gltAB.
    
    dy(4) = aA*(KT^n2)/(KT^n2+y(3)^n2) -  bA*y(4);
    %dy(4) = aA - bA*y(4); alternative negative feedback

    % GS is being constantly produced and degraded, but also inhibited by
    % glutamine
    dy(5) = aS - bS*y(5) - aF*y(5)*y(1)^n1/(KF^n1 + y(1)^n1) + bF*y(2);
    
    % glutamate has a constant supply while also being consumed by
    % GS and GDH and resupplied by GltAB 
    dy(6) = aE - aN*y(6) - aq*y(5)*y(6)*y(7) + bq*y(4)*y(1)*y(8) - cE*y(6);
    
    % ammonium is being made by rocG and consumed by GS
    dy(7) = aN*y(6) - aq*y(5)*y(6)*y(7);
    
    % alpha-ketoglutarate is being supplied from the TCA cycle, being made
    % by rocG, and being consumed by gltAB
    dy(8) = ao + aN*y(6) - co*y(8) - bq*y(4)*y(1)*y(8); 
    
    % Fluorescent reporter PnasA-yfp
    dy(9) = aR*(y(3)^n2)/(KR^n2+y(3)^n2) -  bR*y(9);
    
    % Gene "X", consumes glutamine and inhibited by TnrA in the alternative
    % negative feedback loop
    %dy(10) = aX*(KT^n2)/(KT^n2+y(3)^n2) -  bX*y(10);
    
    dy = dy*e;
end