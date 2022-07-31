function [lb,ub,dim] = func_bound(Func_num)
%func_bound 输入CEC2019函数编号，返回相应函数的上下界和维度
switch Func_num
    case 1
        lb=-8192;ub=8192;dim=9;
    
    case 2
        lb=-16384;ub=16384;dim=16;
    
    case 3
        lb=-4;ub=4;dim=18;
    
    case 4
        lb=-100;ub=100;dim=10;
    
    case 5
        lb=-100;ub=100;dim=10;
    
    case 6
        lb=-100;ub=100;dim=10;
    
    case 7
        lb=-100;ub=100;dim=10;
    
    case 8
        lb=-100;ub=100;dim=10;
    
    case 9
        lb=-100;ub=100;dim=10;
    
    case 10
        lb=-100;ub=100;dim=10;
end
end

