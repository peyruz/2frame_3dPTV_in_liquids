function [contSort]=arrangeContours(contours, ixCol,CalCamMaskAll)

imNum=max(size(contours));

contSort=cell(size(contours));

for ii=1:imNum
    lag=0;
    for jj=1:size(CalCamMaskAll,1)
        if CalCamMaskAll(jj)==1
                contSort{ii}{jj-lag}=contours{ii}{ixCol{ii}(jj)};
            else
            lag=lag+1;
        end
    end
end

end


