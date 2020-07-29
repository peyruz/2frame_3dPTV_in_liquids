function out1=InvRefract(i,c,s1,s2,s3,n)
out1=(s2-c(2))/(c(3)/sqrt(1-i^2-((s1-c(1))/(s2-c(2)))^2*i^2)-s3/sqrt(n^2-i^2-((s1-c(1))/(s2-c(2)))^2*i^2))-i;
end