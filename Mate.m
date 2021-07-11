function [contractiona,contractionf] = Mate(contractionax,contractionfx)
% Because of the standard of the lable algorithm , the purpose of the 
% function is to get the same length matirx of angle and force.

% The input data is the start point of the angle and force
% The output data is the respective
%此函数用于函数对齐/校准  使得峰检出不一致的矩阵元素为同矩阵
counter4 = 1;
counter5 = 1;
if length(contractionax)>length(contractionfx)
    contractionf = contractionfx;
    for i = 1:length(contractionax)
        for counter1 =1:length(contractionfx)
             if abs(contractionax(i)-contractionfx(counter1))<100% 100 is the based on time frame
                contractiona(counter4) = contractionax(i);
                counter4 = counter4+1;
             end
        end
    end
    
else
    contractiona = contractionax;
    for i = 1:length(contractionfx)
        for counter1 = 1:length(contractionax)
            if abs(contractionfx(i) - contractionax(counter1)) <100
               contractionf(counter5) = contractionfx(i);
               counter5 = counter5+1;
            end
        end
    end
end
    
            
       

end

