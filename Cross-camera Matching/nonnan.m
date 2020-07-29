function [NonNaNoutput]=nonnan(NaNinput)
NonNaNoutput=NaNinput(~isnan(NaNinput));
end