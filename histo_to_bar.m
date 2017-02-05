function [ newbins,newhisto ] = histo_to_bar( bins,histovals )
%makes continuous bar plot from normal histogram
binw=abs((bins(2)-bins(1)))/2;
newbins=zeros(1,2*length(bins)+1);
newhisto=zeros(1,2*length(bins)+1);
counter=1;
newbins(1,1)=(bins(1)-binw);
for j=2:2:length(newbins)-1
    newbins(1,j)=(bins(1)-binw)+(j-2)*(binw);
    newbins(1,j+1)=(bins(1)+binw)+(j-2)*(binw);
    newhisto(1,j)=histovals(counter);
    newhisto(1,j+1)=histovals(counter);
    counter=counter+1;
end
newbins(1,j+2)=(bins(1)+binw)+(j-2)*(binw);
newhisto(1,j+2)=0;

